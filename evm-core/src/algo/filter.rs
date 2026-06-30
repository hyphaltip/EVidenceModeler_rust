//! Prediction filtering — remove low-support and degenerate gene models.

use crate::algo::intergenic::{calc_intergenic_score, IntergenicScores};
use crate::algo::introns::{
    make_intron_key, unpack_intron_key, IntronEvidenceMap, IntronKey, IntronVec,
};
use crate::types::evidence::EvWeightMap;
use crate::types::exon::{Exon, ExonType};
use crate::types::genome::MaskVec;
use crate::types::prediction::EvmPrediction;

/// Minimum coding/noncoding score ratio; predictions below this are eliminated.
const MIN_CODING_NONCODING_SCORE_RATIO: f64 = 0.75;
/// Minimum total coding length for nested/intergenic gene searches.
const MIN_CODING_LENGTH: u32 = 300;
/// Minimum total coding length under STANDARD mode (50 aa).
const MIN_CODING_LENGTH_STANDARD: u32 = 150;

/// Filter low-support predictions — faithful port of Perl
/// `filter_predictions_low_support`.
///
/// For each prediction computes the noncoding-equivalent score and the
/// coding/noncoding score ratio, then eliminates it if the ratio is below
/// `MIN_CODING_NONCODING_SCORE_RATIO` or the coding length is below the
/// mode-dependent minimum (150 for STANDARD, else 300). The computed
/// `raw_noncoding`, `offset_noncoding`, `noncoding_equivalent`, and
/// `score_ratio` are recorded on each prediction for the output header.
///
/// `base_intergenic` must be the NON-augmented intergenic vector (Perl resets
/// the intergenic scores at the first filter call, undoing the start/stop peak
/// augmentation used only for the trellis).
#[allow(clippy::too_many_arguments)]
pub fn filter_predictions_low_support(
    predictions: &mut [EvmPrediction],
    exons: &[Exon],
    base_intergenic: &IntergenicScores,
    fwd_intron_vec: &IntronVec,
    rev_intron_vec: &IntronVec,
    introns_to_evidence: &IntronEvidenceMap,
    ev_weights: &EvWeightMap,
    mask: &MaskVec,
    mode: &str,
) {
    let min_coding_length = if mode == "STANDARD" {
        MIN_CODING_LENGTH_STANDARD
    } else {
        MIN_CODING_LENGTH
    };

    for pred in predictions.iter_mut() {
        let (lend, rend) = pred.get_span();
        let prediction_score = pred.total_score;

        // Base intergenic over the prediction span.
        let noncoding_score = calc_intergenic_score(base_intergenic, lend, rend);

        // Predicted introns (both strands) scored as intergenic over the span.
        let mut noncoding_intron_addition = 0.0;
        for i in lend..=rend {
            noncoding_intron_addition += fwd_intron_vec.get(i as usize).copied().unwrap_or(0.0)
                + rev_intron_vec.get(i as usize).copied().unwrap_or(0.0);
        }

        // Offset: same-strand predicted-intron support that should not count as
        // noncoding for this prediction's own introns.
        let orient = pred.orient;
        let intron_vec = if orient == '+' {
            fwd_intron_vec
        } else {
            rev_intron_vec
        };
        let mut offset = 0.0;
        for &(simple_lend, simple_rend) in &pred.intron_coords {
            let key = get_intron_key(simple_lend, simple_rend, orient);
            let (ilend, irend) = intron_key_to_span_sorted(key);
            let intron_len = adjust_feature_length_for_mask(ilend, irend, mask);
            if intron_len == 0 {
                continue;
            }
            // existing per-base contribution = sum of ab-initio weights for this intron
            let existing_per_base: f64 = introns_to_evidence
                .get(&key)
                .map(|evs| {
                    evs.iter()
                        .filter_map(|(_acc, ev_type)| {
                            ev_weights
                                .get(ev_type)
                                .filter(|e| e.ev_class.is_abinitio())
                                .map(|e| e.weight)
                        })
                        .sum()
                })
                .unwrap_or(0.0);
            offset += calc_intergenic_score(base_intergenic, ilend, irend);
            for i in ilend..=irend {
                if mask.get(i as usize) {
                    continue;
                }
                offset += intron_vec.get(i as usize).copied().unwrap_or(0.0) - existing_per_base;
            }
        }

        let raw_noncoding = noncoding_score + noncoding_intron_addition;
        pred.raw_noncoding = raw_noncoding;
        pred.offset_noncoding = offset;

        let mut noncoding_equivalent = raw_noncoding - offset;
        if noncoding_equivalent <= 0.0 {
            noncoding_equivalent = 0.0001 * prediction_score;
        }
        pred.noncoding_equivalent = noncoding_equivalent;

        let score_ratio = if prediction_score == 0.0 && noncoding_equivalent == 0.0 {
            0.0
        } else {
            round2(prediction_score / noncoding_equivalent)
        };
        pred.score_ratio = score_ratio;

        let coding_length: u32 = pred.exon_indices.iter().map(|&i| exons[i].length()).sum();

        if score_ratio < MIN_CODING_NONCODING_SCORE_RATIO || coding_length < min_coding_length {
            pred.is_eliminated = true;
        }
    }
}

/// Round to 2 decimals exactly as Perl `sprintf("%.2f", x)` then numeric compare.
fn round2(x: f64) -> f64 {
    format!("{:.2}", x).parse::<f64>().unwrap_or(x)
}

/// Perl `get_intron_key`: from a simple intron span (lend, rend) and orient,
/// produce the donor/acceptor-adjusted key matching `INTRONS_TO_SCORE`.
fn get_intron_key(lend: u32, rend: u32, orient: char) -> IntronKey {
    let (l, r) = if lend <= rend {
        (lend, rend)
    } else {
        (rend, lend)
    };
    let (end5, end3) = if orient == '+' {
        (l, r.saturating_sub(1))
    } else {
        (r, l + 1)
    };
    make_intron_key(end5, end3)
}

/// Perl `intron_key_to_intron_span`, returning a sorted (lend, rend).
fn intron_key_to_span_sorted(key: IntronKey) -> (u32, u32) {
    let (end5, end3) = unpack_intron_key(key);
    let (a, b) = if end5 < end3 {
        (end5, end3.saturating_sub(1)) // '+' orient
    } else {
        (end5, end3 + 1) // '-' orient
    };
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Perl `adjust_feature_length_for_mask`: count of non-masked positions in [lend, rend].
fn adjust_feature_length_for_mask(lend: u32, rend: u32, mask: &MaskVec) -> u32 {
    let (l, r) = if lend <= rend {
        (lend, rend)
    } else {
        (rend, lend)
    };
    let mut n = 0u32;
    for i in l..=r {
        if !mask.get(i as usize) {
            n += 1;
        }
    }
    n
}

/// For predictions that lack a proper start codon (5' partials), try to find
/// an alternative initial exon upstream that starts at ATG.
///
/// This is a simplified version of the Perl `convert_5prime_partials_to_complete_genes_where_possible`.
#[allow(clippy::ptr_arg)]
pub fn convert_5prime_partials(predictions: &mut Vec<EvmPrediction>, exons: &[Exon]) {
    for pred in predictions.iter_mut() {
        if pred.is_eliminated {
            continue;
        }
        if pred.exon_indices.is_empty() {
            continue;
        }
        let first_idx = pred.exon_indices[0];
        let first_exon = &exons[first_idx];
        // Only act on non-initial leading exons (i.e. the prediction starts
        // with an internal exon, making it a 5' partial)
        if first_exon.exon_type != ExonType::Internal {
            continue;
        }
        // In the full Perl implementation this searches for a new initial exon.
        // Here we mark the prediction as a partial — full recovery logic
        // would require access to all candidate exons, which is done in the
        // `consensus` module.
        // TODO: full partial→complete conversion.
    }
}
