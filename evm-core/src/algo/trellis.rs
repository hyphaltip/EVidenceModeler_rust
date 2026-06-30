//! Dynamic-programming trellis: finding the highest-scoring path through exons.

use crate::algo::intergenic::{calc_intergenic_score, IntergenicScores};
use crate::algo::introns::{make_intron_key, IntronScoreMap};
use crate::algo::phases::is_stop_codon;
use crate::types::exon::{
    Exon, ExonType, LinkageTables, TYPE_ORIENT_BOUND, TYPE_ORIENT_INITIAL_REV,
    TYPE_ORIENT_TERMINAL_REV,
};
use crate::types::prediction::EvmPrediction;

/// Result of compatibility check between two exons.
pub enum CompatResult {
    /// Incompatible — cannot link A → B.
    Incompatible,
    /// Compatible; the score bonus (intron score or intergenic score) to add.
    Compatible(f64),
}

/// Test whether exon A (left / upstream) can be followed by exon B (right / downstream)
/// in the trellis, and compute the join score.
///
/// Returns `Compatible(score)` on success, `Incompatible` otherwise.
#[allow(clippy::too_many_arguments)]
pub fn are_compatible_exons(
    exon_a: &Exon,
    exon_b: &Exon,
    linkages: &LinkageTables,
    introns_to_score: &IntronScoreMap,
    intergenic_scores: &IntergenicScores,
    stop_codons: &[[u8; 3]],
) -> CompatResult {
    let key_a = exon_a.type_orient;
    let key_b = exon_b.type_orient;

    // Check linkage is allowed
    if !linkages.contains_acceptable(key_a, key_b) {
        return CompatResult::Incompatible;
    }

    let (a_lend, a_rend) = exon_a.coords_sorted();
    let (b_lend, b_rend) = exon_b.coords_sorted();

    // No overlap allowed
    if a_lend <= b_rend && a_rend >= b_lend {
        return CompatResult::Incompatible;
    }

    let reverse_strand = (TYPE_ORIENT_INITIAL_REV..=TYPE_ORIENT_TERMINAL_REV).contains(&key_a);

    if linkages.contains_phased(key_a, key_b) {
        // Check intron validity
        let (intron_end5, intron_end3) = if reverse_strand {
            // Reverse strand: intron is upstream of A in genomic terms
            (b_lend - 1, a_rend + 2)
        } else {
            (a_rend + 1, b_lend - 2)
        };

        let intron_key = make_intron_key(intron_end5, intron_end3);
        let intron_score = match introns_to_score.get(&intron_key) {
            Some(&s) => s,
            None => return CompatResult::Incompatible,
        };

        // Check frame compatibility
        let (before, after) = if reverse_strand {
            (exon_b, exon_a)
        } else {
            (exon_a, exon_b)
        };
        if !linkages.contains_frame_pair(before.end_frame, after.start_frame) {
            return CompatResult::Incompatible;
        }

        // Check no stop codon created across the junction
        let end_frame = before.end_frame % 3;
        let seq_junction: [u8; 4] = [
            before.right_seq_boundary[0],
            before.right_seq_boundary[1],
            after.left_seq_boundary[0],
            after.left_seq_boundary[1],
        ];
        let potential_stop = match end_frame {
            1 => Some(&seq_junction[1..4]),
            2 => Some(&seq_junction[0..3]),
            _ => None, // frame 0 / 3: no partial codon at junction
        };
        if let Some(codon) = potential_stop {
            if is_stop_codon(codon, stop_codons) {
                return CompatResult::Incompatible;
            }
        }

        CompatResult::Compatible(intron_score)
    } else if linkages.contains_intergenic(key_a, key_b) {
        let score = calc_intergenic_score(intergenic_scores, a_rend + 1, b_lend - 1);
        CompatResult::Compatible(score)
    } else {
        CompatResult::Incompatible
    }
}

/// Score the connection between a boundary node and an exon (or vice-versa).
///
/// Mirrors the Perl `score_boundary_condition`: a feature may extend to the
/// sequence terminus, scoring the flanking region as intergenic.
/// - left bound → exon B: intergenic over (bound_end3, exonB_lend - 1)
/// - exon A → right bound: intergenic over (exonA_rend + 1, bound_end5)
/// - both bound: not allowed (incompatible).
///
/// `exon_a` is the upstream (left) node, `exon_b` the downstream (right) node.
pub fn score_boundary_condition(
    exon_a: &Exon,
    exon_b: &Exon,
    intergenic_scores: &IntergenicScores,
) -> CompatResult {
    let a_bound = exon_a.exon_type == ExonType::Bound;
    let b_bound = exon_b.exon_type == ExonType::Bound;

    if a_bound && b_bound {
        return CompatResult::Incompatible;
    }

    if a_bound {
        // left boundary → exon B; bound has end5 == end3 == range_lend.
        let (b_lend, _) = exon_b.coords_sorted();
        let s = calc_intergenic_score(intergenic_scores, exon_a.end3, b_lend.saturating_sub(1));
        CompatResult::Compatible(s)
    } else {
        // exon A → right boundary; bound has end5 == end3 == range_rend.
        let (_, a_rend) = exon_a.coords_sorted();
        let s = calc_intergenic_score(intergenic_scores, a_rend + 1, exon_b.end5);
        CompatResult::Compatible(s)
    }
}

/// Build the trellis over `exons` restricted to [range_lend, range_rend].
///
/// Exons must be sorted by end5 ascending before calling this function.
/// Returns the index of the highest-scoring exon (or None if empty).
#[allow(clippy::too_many_arguments)]
pub fn build_trellis(
    exons: &mut Vec<Exon>,
    range_lend: u32,
    range_rend: u32,
    linkages: &LinkageTables,
    introns_to_score: &IntronScoreMap,
    intergenic_scores: &IntergenicScores,
    stop_codons: &[[u8; 3]],
    max_prev_exons_compare: usize,
) -> Option<usize> {
    if exons.is_empty() {
        return None;
    }

    // Sort real exons by 5' end (Perl re-sorts inside build_trellis).
    exons.sort_by_key(|e| e.end5);

    // Add boundary sentinel nodes in Perl order: left bound at the FRONT,
    // right bound at the END. Layout: [left_bound, exons…, right_bound].
    let mut left_bound = Exon::new(range_lend, range_lend);
    left_bound.exon_type = ExonType::Bound;
    left_bound.type_orient = TYPE_ORIENT_BOUND;
    left_bound.start_frame = 1;
    left_bound.end_frame = 1;
    exons.insert(0, left_bound);

    let mut right_bound = Exon::new(range_rend, range_rend);
    right_bound.exon_type = ExonType::Bound;
    right_bound.type_orient = TYPE_ORIENT_BOUND;
    right_bound.start_frame = 1;
    right_bound.end_frame = 1;
    exons.push(right_bound);

    // Reset link and sum_score for every node (bounds included).
    for exon in exons.iter_mut() {
        exon.sum_score = exon.base_score;
        exon.link = None;
    }

    let num_exons = exons.len();
    let mut highest_score = 0.0f64;
    let mut highest_idx = 0usize;

    for i in 1..num_exons {
        let base_score_i = exons[i].base_score;
        let mut best_sum = exons[i].sum_score;

        let mut compare_count = 0usize;
        let mut found_compatible = false;

        let i_type_bound = exons[i].exon_type == ExonType::Bound;

        let mut j = (i - 1) as isize;

        while j >= 0 && (compare_count < max_prev_exons_compare || !found_compatible) {
            let ji = j as usize;
            compare_count += 1;
            j -= 1;

            let j_type_bound = exons[ji].exon_type == ExonType::Bound;

            let join_score = if i_type_bound || j_type_bound {
                match score_boundary_condition(&exons[ji], &exons[i], intergenic_scores) {
                    CompatResult::Compatible(s) => Some(s),
                    CompatResult::Incompatible => None,
                }
            } else {
                match are_compatible_exons(
                    &exons[ji],
                    &exons[i],
                    linkages,
                    introns_to_score,
                    intergenic_scores,
                    stop_codons,
                ) {
                    CompatResult::Compatible(s) => Some(s),
                    CompatResult::Incompatible => None,
                }
            };

            if let Some(js) = join_score {
                found_compatible = true;
                let candidate = base_score_i + exons[ji].sum_score + js;
                if candidate > best_sum {
                    best_sum = candidate;
                    exons[i].link = Some(ji);
                    exons[i].sum_score = best_sum;
                }
            }
        }

        if best_sum >= highest_score {
            highest_score = best_sum;
            highest_idx = i;
        }
    }

    Some(highest_idx)
}

/// Traverse the trellis from the highest-scoring exon and assemble
/// gene predictions.  Returns a Vec of `EvmPrediction` objects.
pub fn traverse_path(exons: &[Exon], top_idx: usize) -> Vec<EvmPrediction> {
    let mut chain: Vec<usize> = Vec::new();
    let mut cur = Some(top_idx);
    while let Some(idx) = cur {
        chain.push(idx);
        cur = exons[idx].link;
    }
    chain.reverse(); // left to right

    // Split chain into individual gene predictions at terminal/single boundaries
    let mut predictions: Vec<EvmPrediction> = Vec::new();
    let mut current: Vec<usize> = Vec::new();

    for &idx in &chain {
        let exon = &exons[idx];
        let type_orient = exon.type_orient_key();
        if type_orient.contains("bound") {
            if !current.is_empty() {
                let pred = EvmPrediction::new(current.clone(), exons);
                predictions.push(pred);
                current.clear();
            }
            continue;
        }
        current.push(idx);
        // Terminate prediction at terminal+ or single, or initial- (reverse strand terminal)
        if type_orient == "terminal+" || type_orient.contains("single") || type_orient == "initial-"
        {
            let pred = EvmPrediction::new(current.clone(), exons);
            predictions.push(pred);
            current.clear();
        }
    }
    if !current.is_empty() {
        let pred = EvmPrediction::new(current, exons);
        predictions.push(pred);
    }

    predictions
}
