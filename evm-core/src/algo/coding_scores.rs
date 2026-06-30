//! Per-base coding score accumulation.

use crate::types::evidence::EvClass;
use crate::types::exon::Exon;
use crate::types::genome::MaskVec;

/// Coding score vector: one f64 per genome position (1-indexed, index 0 unused).
pub type CodingScores = Vec<f64>;

/// Create a zeroed coding-score vector of length `seq_len + 2`.
pub fn new_coding_scores(seq_len: usize) -> CodingScores {
    vec![0.0f64; seq_len + 2]
}

/// Increment the coding score vector for a range [end5..=end3] by `weight`,
/// but only for positions that are not masked.
/// Only PROTEIN and ABINITIO_PREDICTION classes contribute to the coding vector.
/// Negative weights are floored at 0 (prevent negative coding score).
pub fn add_match_coverage(
    scores: &mut CodingScores,
    mask: &MaskVec,
    end5: u32,
    end3: u32,
    weight: f64,
    ev_class: &EvClass,
) {
    if !ev_class.contributes_to_coding_vec() {
        return;
    }
    if end5 > end3 {
        return; // only called with forward-strand coords
    }
    for i in end5..=end3 {
        let i = i as usize;
        if !mask.get(i) {
            if weight < 0.0 {
                let cur = scores[i];
                scores[i] = f64::max(cur + weight, 0.0);
            } else {
                scores[i] += weight;
            }
        }
    }
}

/// Assign base_score to each exon by summing the coding scores over its span
/// plus evidence-specific contributions for TRANSCRIPT and OTHER_PREDICTION exons.
#[allow(clippy::ptr_arg)]
pub fn score_exons(
    exons: &mut Vec<Exon>,
    coding_scores: &CodingScores,
    mask: &MaskVec,
    ev_weight_fn: &dyn Fn(&str) -> Option<f64>, // ev_type → weight
    ev_class_fn: &dyn Fn(&str) -> Option<EvClass>, // ev_type → class
) {
    for exon in exons.iter_mut() {
        let (end5, end3) = exon.coords_sorted();
        let mut coding_score = 0.0f64;

        // Sum coding vector contributions (PROTEIN + ABINITIO already in vector)
        for i in end5..=end3 {
            let v = coding_scores[i as usize];
            if v > 0.0 {
                coding_score += v;
            }
        }

        // Add exon-specific contributions for TRANSCRIPT and OTHER_PREDICTION
        for (accession, ev_type) in &exon.evidence {
            let cls = ev_class_fn(ev_type);
            let matches_class = cls
                .as_ref()
                .is_some_and(|c| matches!(c, EvClass::Transcript | EvClass::OtherPrediction));
            if matches_class {
                if let Some(w) = ev_weight_fn(ev_type) {
                    let mut contrib = 0.0;
                    for i in end5..=end3 {
                        if !mask.get(i as usize) {
                            contrib += w;
                        }
                    }
                    coding_score += contrib;
                    let _ = accession; // used implicitly
                }
            }
        }

        exon.base_score = coding_score;
        exon.sum_score = coding_score;
    }
}
