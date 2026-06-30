//! Recursive consensus gene prediction — trellis + recursion on tail/intergenic regions.

use crate::algo::coding_scores::CodingScores;
use crate::algo::filter::filter_predictions_low_support;
use crate::algo::intergenic::{get_intergenic_regions, IntergenicScores};
use crate::algo::introns::{make_intron_key, IntronEvidenceMap, IntronScoreMap, IntronVec};
use crate::algo::trellis::{build_trellis, traverse_path};
use crate::types::evidence::EvWeightMap;
use crate::types::exon::{Exon, ExonType, LinkageTables};
use crate::types::genome::MaskVec;
use crate::types::prediction::{EvmPrediction, PredMode};
use anyhow::Result;

/// Parameters controlling the recursive search.
pub struct ConsensusParams<'a> {
    pub exons: &'a mut Vec<Exon>,
    pub introns_to_score: &'a IntronScoreMap,
    pub introns_to_evidence: &'a IntronEvidenceMap,
    pub ev_weights: &'a EvWeightMap,
    pub mask: &'a MaskVec,
    /// Augmented intergenic scores (start/stop peak augmentation applied) — used
    /// by the trellis on the first build.
    pub intergenic_scores: &'a IntergenicScores,
    /// Base (non-augmented) intergenic scores — used by the low-support filter.
    pub base_intergenic_scores: &'a IntergenicScores,
    pub linkages: &'a LinkageTables,
    pub stop_codons: &'a [[u8; 3]],
    pub coding_scores: &'a CodingScores,
    pub fwd_intron_vec: &'a IntronVec,
    pub rev_intron_vec: &'a IntronVec,
    pub max_prev_exons_compare: usize,
    pub min_intergenic_size_on_re_search: u32,
    pub min_gene_length_size_on_re_search: u32,
    pub report_elm: bool,
    pub recursion_limit: usize,
}

/// Generate consensus gene predictions for [range_lend, range_rend].
///
/// Outputs formatted prediction text (matching the Perl output format) to `output`.
pub fn generate_consensus_gene_predictions(
    range_lend: u32,
    range_rend: u32,
    mode: PredMode,
    params: &mut ConsensusParams,
    recursion_count: &mut usize,
    output: &mut Vec<String>,
) -> Result<()> {
    *recursion_count += 1;

    let region_len = range_rend.saturating_sub(range_lend) + 1;
    if region_len < params.min_gene_length_size_on_re_search && mode != PredMode::Standard {
        *recursion_count -= 1;
        return Ok(());
    }

    if *recursion_count > params.recursion_limit {
        log::warn!("Recursion limit ({}) exceeded", params.recursion_limit);
        *recursion_count -= 1;
        return Ok(());
    }

    // Collect exons within range
    let exon_indices_in_range: Vec<usize> = {
        let exons_ref = &*params.exons;
        (0..exons_ref.len())
            .filter(|&i| {
                let (l, r) = exons_ref[i].coords_sorted();
                l >= range_lend && r <= range_rend
            })
            .collect()
    };

    if exon_indices_in_range.is_empty() {
        log::debug!("No exons in range {}-{}", range_lend, range_rend);
        *recursion_count -= 1;
        return Ok(());
    }

    // Build a local copy of exons for the trellis (subset in range)
    let mut local_exons: Vec<Exon> = exon_indices_in_range
        .iter()
        .map(|&i| params.exons[i].clone())
        .collect();
    local_exons.sort_by_key(|e| e.end5);

    let top_idx = build_trellis(
        &mut local_exons,
        range_lend,
        range_rend,
        params.linkages,
        params.introns_to_score,
        params.intergenic_scores,
        params.stop_codons,
        params.max_prev_exons_compare,
    );

    let top_idx = match top_idx {
        Some(i) => i,
        None => {
            *recursion_count -= 1;
            return Ok(());
        }
    };

    let mut predictions = traverse_path(&local_exons, top_idx);
    if predictions.is_empty() {
        *recursion_count -= 1;
        return Ok(());
    }

    // Compute prediction_score / orientation / intron coords (Perl _init parity).
    for pred in predictions.iter_mut() {
        pred.finalize(&local_exons, params.introns_to_score);
    }

    // Filter low-support predictions
    filter_predictions_low_support(
        &mut predictions,
        &local_exons,
        params.base_intergenic_scores,
        params.fwd_intron_vec,
        params.rev_intron_vec,
        params.introns_to_evidence,
        params.ev_weights,
        params.mask,
        mode.as_str(),
    );

    // Convert 5'-partial genes to complete genes where an overlapping initial
    // exon exists (Perl convert_5prime_partials_to_complete_genes_where_possible,
    // runs after the filter, before reporting / tail recursion).
    convert_5prime_partials_to_complete_genes(
        &mut predictions,
        &mut local_exons,
        params.exons,
        params.introns_to_score,
    );

    let preds_remain: Vec<&EvmPrediction> =
        predictions.iter().filter(|p| !p.is_eliminated).collect();

    if preds_remain.is_empty() && !params.report_elm {
        *recursion_count -= 1;
        return Ok(());
    }

    // Determine span of predictions. Perl uses the full set of predictions
    // (including eliminated models when --report_ELM is on) for both the
    // "!! Predictions spanning range" line and tail-recursion boundaries.
    let (pred_span_lend, pred_span_rend) = predictions
        .iter()
        .fold((range_rend, range_lend), |(l, r), p| {
            (l.min(p.lend), r.max(p.rend))
        });

    // Emit predictions: one "!!" range line for the call, then each prediction's
    // block followed by a blank line (Perl prints `toString() . "\n"`).
    output.push(format!(
        "!! Predictions spanning range {} - {} [R{}]\n",
        pred_span_lend, pred_span_rend, *recursion_count
    ));
    for pred in &predictions {
        if pred.is_eliminated && !params.report_elm {
            continue;
        }
        let text = format_prediction(
            pred,
            &local_exons,
            params.introns_to_evidence,
            mode.as_str(),
        );
        output.push(text);
        output.push("\n".to_string());
    }

    // Recursion: tail regions
    if mode != PredMode::Intron {
        let left_len = pred_span_lend.saturating_sub(range_lend);
        if left_len >= params.min_intergenic_size_on_re_search {
            generate_consensus_gene_predictions(
                range_lend,
                pred_span_lend - 1,
                mode.clone(),
                params,
                recursion_count,
                output,
            )?;
        }
        let right_len = range_rend.saturating_sub(pred_span_rend);
        if right_len >= params.min_intergenic_size_on_re_search {
            generate_consensus_gene_predictions(
                pred_span_rend + 1,
                range_rend,
                mode.clone(),
                params,
                recursion_count,
                output,
            )?;
        }

        // Recursion: intergenic regions between predictions
        if params.min_gene_length_size_on_re_search > 0 {
            // Perl get_intergenic_regions operates on all predictions, including
            // eliminated models when --report_ELM is enabled.
            let spans: Vec<(u32, u32)> = predictions.iter().map(|p| (p.lend, p.rend)).collect();
            for (ig_l, ig_r) in get_intergenic_regions(&spans) {
                let ig_len = ig_r.saturating_sub(ig_l) + 1;
                if ig_len >= params.min_gene_length_size_on_re_search {
                    generate_consensus_gene_predictions(
                        ig_l,
                        ig_r,
                        mode.clone(),
                        params,
                        recursion_count,
                        output,
                    )?;
                }
            }
        }
    }

    *recursion_count -= 1;
    Ok(())
}

/// A prediction is 5'-partial if none of its exons is an initial or single exon
/// (Perl `EVM_prediction::is_5prime_partial`).
fn is_5prime_partial(pred: &EvmPrediction, exons: &[Exon]) -> bool {
    !pred
        .exon_indices
        .iter()
        .any(|&i| matches!(exons[i].exon_type, ExonType::Initial | ExonType::Single))
}

/// Faithful port of Perl `convert_5prime_partials_to_complete_genes_where_possible`.
///
/// For each multi-exon 5'-partial prediction whose gene-start exon is `internal`,
/// search the global exon pool for an overlapping `initial` exon with the same
/// orientation, identical 3' coordinate (`end3`) and end frame, and replace the
/// internal gene-start exon with the best-scoring such initial exon. The
/// prediction is then re-initialised (span / score / introns recomputed).
fn convert_5prime_partials_to_complete_genes(
    predictions: &mut [EvmPrediction],
    local_exons: &mut Vec<Exon>,
    search_pool: &[Exon],
    introns_to_score: &IntronScoreMap,
) {
    for pred in predictions.iter_mut() {
        if !is_5prime_partial(pred, local_exons) {
            continue;
        }
        // Only multi-exon genes (Perl skips single-exon).
        if pred.exon_indices.len() < 2 {
            continue;
        }

        let orient = pred.orient;
        // exon_indices are sorted by end5 ascending (finalize). The gene-start
        // exon is the first for '+' and the last for '-' (Perl reverses for '-').
        let gene_start_pos = if orient == '-' {
            pred.exon_indices.len() - 1
        } else {
            0
        };
        let gs_idx = pred.exon_indices[gene_start_pos];
        let gs = &local_exons[gs_idx];

        // Perl only converts when the gene-start exon is internal.
        if gs.exon_type != ExonType::Internal {
            continue;
        }

        let (gs_l, gs_r) = gs.coords_sorted();
        let gs_end3 = gs.end3;
        let gs_end_frame = gs.end_frame;

        // Find the best overlapping initial exon (Perl find_overlapping_exons:
        // end5 < rendRange && end3 > lendRange — strict overlap on sorted coords).
        let mut best: Option<&Exon> = None;
        for e in search_pool {
            let (el, er) = e.coords_sorted();
            if !(el < gs_r && er > gs_l) {
                continue;
            }
            if e.exon_type != ExonType::Initial {
                continue;
            }
            if e.orientation.as_char() != orient {
                continue;
            }
            if e.end3 != gs_end3 {
                continue;
            }
            if e.end_frame != gs_end_frame {
                continue;
            }
            match best {
                Some(b) if b.base_score >= e.base_score => {}
                _ => best = Some(e),
            }
        }

        if let Some(b) = best {
            // Splice the replacement exon into the local pool and the prediction,
            // then re-init (Perl replace_exons -> _init).
            let new_idx = local_exons.len();
            local_exons.push(b.clone());
            pred.exon_indices[gene_start_pos] = new_idx;
            pred.finalize(local_exons, introns_to_score);
        }
    }
}

/// Format a prediction as EVM output text — faithful port of Perl
/// `EVM_prediction::toString` (header + coordinate-sorted exon/intron rows).
fn format_prediction(
    pred: &EvmPrediction,
    exons: &[Exon],
    introns_to_evidence: &IntronEvidenceMap,
    mode: &str,
) -> String {
    let orient = pred.orient;

    // Header line. score_ratio is stored pre-rounded (Perl sprintf "%.2f");
    // the remaining numeric fields are formatted with two decimals here.
    let mut s = format!(
        "# EVM prediction: Mode:{} S-ratio: {:.2} {}-{} orient({}) score({:.2}) \
noncoding_equivalent({:.2}) raw_noncoding({:.2}) offset({:.2}) ",
        mode,
        pred.score_ratio,
        pred.lend,
        pred.rend,
        orient,
        pred.total_score,
        pred.noncoding_equivalent,
        pred.raw_noncoding,
        pred.offset_noncoding,
    );
    if pred.is_eliminated {
        s.push_str(" *** ELIMINATED *** ");
    }
    s.push('\n');

    // Build the interleaved, coordinate-sorted component list (Perl orders by
    // the first stored coordinate of each exon/intron).
    enum Comp<'a> {
        Exon(&'a Exon),
        Intron(u32, u32, String),
    }
    let mut components: Vec<(u32, Comp)> = Vec::new();

    for &(intron_lend, intron_rend) in &pred.intron_coords {
        // Display coords: 5'→3' for the strand.
        let (intron_end5, intron_end3) = if orient == '+' {
            (intron_lend, intron_rend)
        } else {
            (intron_rend, intron_lend)
        };
        // Evidence key: donor/acceptor-adjusted (Perl: '+' end3--, '-' end3++).
        let (key5, key3) = if orient == '+' {
            (intron_end5, intron_end3 - 1)
        } else {
            (intron_end5, intron_end3 + 1)
        };
        let key = make_intron_key(key5, key3);
        let ev = introns_to_evidence.get(&key).cloned().unwrap_or_default();
        // Perl emits evidence in hash-iteration order, which is non-deterministic
        // across runs; sort canonically so Rust output is reproducible.
        let mut toks: Vec<String> = ev
            .iter()
            .map(|(acc, et)| format!("{{{};{}}}", acc, et))
            .collect();
        toks.sort();
        let ev_str = toks.join(",");
        components.push((intron_end5, Comp::Intron(intron_end5, intron_end3, ev_str)));
    }

    for &idx in &pred.exon_indices {
        let exon = &exons[idx];
        components.push((exon.end5, Comp::Exon(exon)));
    }

    components.sort_by_key(|(c, _)| *c);

    for (_, comp) in &components {
        match comp {
            Comp::Exon(exon) => {
                let mut row = format!(
                    "{}\t{}\t{}{}\t{}\t{}\t",
                    exon.end5,
                    exon.end3,
                    exon.exon_type.as_str(),
                    exon.orientation.as_char(),
                    exon.start_frame,
                    exon.end_frame,
                );
                // Sort evidence canonically (see intron note above).
                let mut toks: Vec<String> = exon
                    .evidence
                    .iter()
                    .map(|(acc, et)| format!("{{{};{}}}", acc, et))
                    .collect();
                toks.sort();
                row.push_str(&toks.join(","));
                s.push_str(&row);
                s.push('\n');
            }
            Comp::Intron(e5, e3, ev_str) => {
                s.push_str(&format!("{}\t{}\tINTRON\t\t\t{}\n", e5, e3, ev_str));
            }
        }
    }

    s
}
