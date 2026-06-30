//! Load gene prediction GFF3 data and create exon candidates.

use crate::algo::coding_scores::{add_match_coverage, CodingScores};
use crate::algo::introns::{
    add_introns, IntronEvidenceMap, IntronScoreBuilder, PredictedIntronMap,
};
use crate::algo::phases::determine_good_phases;
use crate::io::gff3::Gff3Record;
use crate::types::evidence::{EvClass, EvWeightMap};
use crate::types::exon::{end_frame, Exon, ExonPhase, ExonType, Orientation};
use crate::types::genome::{FeatureVec, MaskVec, FEAT_ACCEPTOR, FEAT_DONOR, FEAT_START, FEAT_STOP};
use anyhow::Result;
use std::collections::HashMap;

/// Load CDS exon records from gene-prediction GFF3 and populate the exon pool.
///
/// Returns the updated exon pool. Intron maps are updated in-place.
#[allow(clippy::too_many_arguments)]
pub fn load_prediction_data(
    genomic_strand: char,
    records: &[Gff3Record],
    genome_features: &FeatureVec,
    genome_seq: &[u8],
    ev_weights: &EvWeightMap,
    min_intron_length: u32,
    mask: &MaskVec,
    coding_scores: &mut CodingScores,
    introns_to_score: &mut IntronScoreBuilder,
    introns_to_evidence: &mut IntronEvidenceMap,
    predicted_introns: &mut PredictedIntronMap,
    exons: &mut Vec<Exon>,
    exons_via_coords: &mut HashMap<String, usize>,
    genomic_seq_len: usize,
) -> Result<()> {
    // Group CDS records by model id (ev_type + Parent id)
    let mut model_to_coords: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut model_to_entry: HashMap<String, (String, EvClass, f64)> = HashMap::new(); // model_id → (ev_type, class, weight)

    for rec in records {
        if rec.feature != "CDS" {
            continue;
        }
        if rec.strand != genomic_strand {
            continue;
        }

        let ev_type = &rec.source;
        let entry = match ev_weights.get(ev_type) {
            Some(e) => e,
            None => {
                log::warn!("Skipping ev_type {} not in weights file", ev_type);
                continue;
            }
        };

        // Only process prediction types
        if !entry.ev_class.is_prediction() {
            continue;
        }

        let parent = rec
            .raw_attributes
            .split("Parent=")
            .nth(1)
            .and_then(|s| s.split([';', ' ']).next())
            .unwrap_or("")
            .to_string();
        if parent.is_empty() {
            continue;
        }

        let model_id = format!("{}_{}", ev_type, parent);

        // Convert to forward-strand end5/end3
        let (end5, end3) = if genomic_strand == '+' {
            (rec.start, rec.end)
        } else {
            (rec.end, rec.start)
        };

        let (fwd_end5, fwd_end3) = if genomic_strand == '-' {
            let rc5 = genomic_seq_len as u32 - end5 + 1;
            let rc3 = genomic_seq_len as u32 - end3 + 1;
            (rc3.min(rc5), rc3.max(rc5))
        } else {
            (end5.min(end3), end5.max(end3))
        };

        model_to_coords
            .entry(model_id.clone())
            .or_default()
            .push((fwd_end5, fwd_end3));
        model_to_entry
            .entry(model_id)
            .or_insert_with(|| (ev_type.clone(), entry.ev_class.clone(), entry.weight));
    }

    // Process each model
    for (model_id, mut coordsets) in model_to_coords {
        coordsets.sort_by_key(|&(a, _)| a);
        let (ev_type, ev_class, weight) = model_to_entry[&model_id].clone();

        // Add coding coverage
        for &(e5, e3) in &coordsets {
            add_match_coverage(coding_scores, mask, e5, e3, weight, &ev_class);
        }

        // Try to classify exons
        if coordsets.len() == 1 {
            let (end5, end3) = coordsets[0];
            // Single exon: needs start at end5 and stop at end3-2
            if genome_features.get(end5 as usize) == FEAT_START
                && genome_features.get((end3 as usize).saturating_sub(2)) == FEAT_STOP
                && (end3 - end5 + 1) >= 3
            {
                add_or_update_exon(
                    &model_id,
                    end5,
                    end3,
                    ExonType::Single,
                    1,
                    &ev_type,
                    &ev_class,
                    weight,
                    genome_seq,
                    genome_features,
                    mask,
                    coding_scores,
                    exons,
                    exons_via_coords,
                    genomic_seq_len,
                );
            }
        } else {
            // Multi-exon: classify each
            let mut valid = true;
            let mut cds_len = 0u32;
            let mut classified: Vec<(u32, u32, ExonType, ExonPhase)> = Vec::new();

            for (i, &(end5, end3)) in coordsets.iter().enumerate() {
                let exon_len = end3 - end5 + 1;
                let start_frame = (cds_len % 3 + 1) as ExonPhase;
                cds_len += exon_len;

                let has_start = genome_features.get(end5 as usize) == FEAT_START;
                let has_stop = genome_features.get((end3 as usize).saturating_sub(2)) == FEAT_STOP;
                let has_acceptor =
                    genome_features.get((end5 as usize).saturating_sub(2)) == FEAT_ACCEPTOR;
                let has_donor = genome_features.get((end3 + 1) as usize) == FEAT_DONOR;

                let exon_type = if i == 0 && has_start && has_donor && exon_len >= 3 {
                    Some(ExonType::Initial)
                } else if i == coordsets.len() - 1 && has_acceptor && has_stop && exon_len >= 3 {
                    Some(ExonType::Terminal)
                } else if i != 0 && i != coordsets.len() - 1 && has_acceptor && has_donor {
                    Some(ExonType::Internal)
                } else {
                    None
                };

                match exon_type {
                    Some(t) => classified.push((end5, end3, t, start_frame)),
                    None => {
                        valid = false;
                        break;
                    }
                }
            }

            if valid {
                for (end5, end3, exon_type, start_frame) in &classified {
                    add_or_update_exon(
                        &model_id,
                        *end5,
                        *end3,
                        *exon_type,
                        *start_frame,
                        &ev_type,
                        &ev_class,
                        weight,
                        genome_seq,
                        genome_features,
                        mask,
                        coding_scores,
                        exons,
                        exons_via_coords,
                        genomic_seq_len,
                    );
                }
                add_introns(
                    &model_id,
                    &coordsets,
                    genomic_strand,
                    weight,
                    &ev_type,
                    &ev_class,
                    min_intron_length,
                    genome_features,
                    mask,
                    introns_to_score,
                    introns_to_evidence,
                    predicted_introns,
                    genomic_seq_len,
                );
            } else {
                // Try to recover partials
                recover_partial_prediction(
                    &model_id,
                    &coordsets,
                    genomic_strand,
                    weight,
                    &ev_type,
                    &ev_class,
                    genome_seq,
                    genome_features,
                    mask,
                    coding_scores,
                    introns_to_score,
                    introns_to_evidence,
                    predicted_introns,
                    exons,
                    exons_via_coords,
                    min_intron_length,
                    genomic_seq_len,
                );
            }
        }
    }

    Ok(())
}

/// Add a new exon to the pool, or update evidence on an existing matching exon.
#[allow(clippy::too_many_arguments)]
fn add_or_update_exon(
    accession: &str,
    end5: u32,
    end3: u32,
    exon_type: ExonType,
    start_frame: ExonPhase,
    ev_type: &str,
    _ev_class: &EvClass,
    _weight: f64,
    genome_seq: &[u8],
    genome_features: &FeatureVec,
    _mask: &MaskVec,
    _coding_scores: &mut CodingScores,
    exons: &mut Vec<Exon>,
    exons_via_coords: &mut HashMap<String, usize>,
    _genomic_seq_len: usize,
) {
    // Validate phase
    let check_end3 = match exon_type {
        ExonType::Terminal | ExonType::Single => end3.saturating_sub(3),
        _ => end3,
    };
    let good_phases = determine_good_phases(genome_features, end5, check_end3);
    if !good_phases.contains(&start_frame) {
        log::warn!(
            "add_exon: {} {}-{} invalid phase {}",
            accession,
            end5,
            end3,
            start_frame
        );
        return;
    }

    let ef = end_frame(start_frame, end3 - end5 + 1);
    let coord_key = format!("{}_{}_{:?}_{}", end5, end3, exon_type, start_frame);

    if let Some(&idx) = exons_via_coords.get(&coord_key) {
        // Existing exon — append evidence
        exons[idx].append_evidence(accession, ev_type);
    } else {
        let mut exon = Exon::new(end5, end3);
        exon.exon_type = exon_type;
        exon.start_frame = start_frame;
        exon.end_frame = ef;
        exon.orientation = if end5 <= end3 {
            Orientation::Fwd
        } else {
            Orientation::Rev
        };
        exon.refresh_type_orient();
        // Store 2-char boundary sequences for stop-codon junction check
        if end5 >= 2 && (end5 - 2) as usize + 2 <= genome_seq.len() {
            let lb = &genome_seq[(end5 - 2) as usize..(end5) as usize];
            exon.left_seq_boundary = [lb[0], lb[1]];
        }
        if (end3 as usize) < genome_seq.len() {
            let rb = &genome_seq[(end3 - 1) as usize..(end3 + 1) as usize];
            exon.right_seq_boundary = [rb[0], if rb.len() > 1 { rb[1] } else { 0 }];
        }
        exon.append_evidence(accession, ev_type);
        let idx = exons.len();
        exons.push(exon);
        exons_via_coords.insert(coord_key, idx);
    }
}

#[allow(clippy::too_many_arguments)]
fn recover_partial_prediction(
    model_id: &str,
    coordsets: &[(u32, u32)],
    genomic_strand: char,
    weight: f64,
    ev_type: &str,
    ev_class: &EvClass,
    genome_seq: &[u8],
    genome_features: &FeatureVec,
    mask: &MaskVec,
    coding_scores: &mut CodingScores,
    introns_to_score: &mut IntronScoreBuilder,
    introns_to_evidence: &mut IntronEvidenceMap,
    predicted_introns: &mut PredictedIntronMap,
    exons: &mut Vec<Exon>,
    exons_via_coords: &mut HashMap<String, usize>,
    min_intron_length: u32,
    genomic_seq_len: usize,
) {
    for &(end5, end3) in coordsets {
        let exon_len = end3 - end5 + 1;
        let has_start = genome_features.get(end5 as usize) == FEAT_START;
        let has_stop = genome_features.get((end3 as usize).saturating_sub(2)) == FEAT_STOP;
        let has_acceptor = genome_features.get((end5 as usize).saturating_sub(2)) == FEAT_ACCEPTOR;
        let has_donor = genome_features.get((end3 + 1) as usize) == FEAT_DONOR;

        if has_start && has_donor && exon_len >= 3 {
            add_or_update_exon(
                model_id,
                end5,
                end3,
                ExonType::Initial,
                1,
                ev_type,
                ev_class,
                weight,
                genome_seq,
                genome_features,
                mask,
                coding_scores,
                exons,
                exons_via_coords,
                genomic_seq_len,
            );
        }
        if has_acceptor && has_stop && exon_len >= 3 {
            let num_prev = exon_len % 3;
            let phase: ExonPhase = match num_prev {
                0 => 1,
                1 => 3,
                2 => 2,
                _ => 1,
            };
            add_or_update_exon(
                model_id,
                end5,
                end3,
                ExonType::Terminal,
                phase,
                ev_type,
                ev_class,
                weight,
                genome_seq,
                genome_features,
                mask,
                coding_scores,
                exons,
                exons_via_coords,
                genomic_seq_len,
            );
        }
        // Candidate internal exons (Perl try_recover_partial_prediction line 1076):
        // acceptor at end5-2 and donor at end3+1. Internal exons are added for
        // every valid reading phase (the Perl add_exon validates phase).
        if has_acceptor && has_donor {
            for phase in determine_good_phases(genome_features, end5, end3) {
                add_or_update_exon(
                    model_id,
                    end5,
                    end3,
                    ExonType::Internal,
                    phase,
                    ev_type,
                    ev_class,
                    weight,
                    genome_seq,
                    genome_features,
                    mask,
                    coding_scores,
                    exons,
                    exons_via_coords,
                    genomic_seq_len,
                );
            }
        }
    }
    add_introns(
        model_id,
        coordsets,
        genomic_strand,
        weight,
        ev_type,
        ev_class,
        min_intron_length,
        genome_features,
        mask,
        introns_to_score,
        introns_to_evidence,
        predicted_introns,
        genomic_seq_len,
    );
}
