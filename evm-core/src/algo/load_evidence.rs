//! Load protein/transcript alignment evidence and create evidence-based exons.

use crate::algo::coding_scores::{add_match_coverage, CodingScores};
use crate::algo::introns::{
    add_introns, IntronEvidenceMap, IntronScoreBuilder, PredictedIntronMap,
};
use crate::algo::phases::determine_good_phases;
use crate::io::gff3::Gff3Record;
use crate::types::evidence::{EvClass, EvWeightMap, EvidenceChain};
use crate::types::exon::{end_frame, Exon, ExonType, Orientation};
use crate::types::genome::{FeatureVec, MaskVec, FEAT_ACCEPTOR, FEAT_DONOR};
use std::collections::HashMap;

/// Parse evidence chains from a GFF3 file containing protein or transcript alignments.
pub fn parse_evidence_chains(
    genomic_strand: char,
    records: &[Gff3Record],
    ev_weights: &EvWeightMap,
    genomic_seq_len: usize,
) -> Vec<EvidenceChain> {
    let mut acc_to_chain: HashMap<String, EvidenceChain> = HashMap::new();

    for rec in records {
        let ev_type = &rec.source;
        if ev_weights.get(ev_type).is_none() {
            log::warn!("Skipping ev_type {} not in weights file", ev_type);
            continue;
        }

        let ev_class = ev_weights[ev_type].ev_class.clone();
        let orient = rec.strand;

        // Filter by strand
        if genomic_strand != '?' && orient != genomic_strand {
            continue;
        }

        let lend = rec.start.min(rec.end);
        let rend = rec.start.max(rec.end);

        // Parse chain ID and parent
        let chain_id = rec.attr("ID").unwrap_or("").to_string();
        let parent_id = rec.attr("Parent").unwrap_or("").to_string();
        let target = rec
            .attr("Target")
            .or_else(|| rec.attr("Query"))
            .map(|s| s.to_string());

        let is_child = !parent_id.is_empty();
        let key_chain_id = if is_child {
            parent_id.clone()
        } else {
            chain_id.clone()
        };
        let key = format!("ev_type:{}/ID={}", ev_type, key_chain_id);

        let (mut end5, mut end3) = if orient == '+' {
            (lend, rend)
        } else {
            (rend, lend)
        };

        if genomic_strand == '-' {
            // Transpose coordinates to forward-strand reference
            end5 = genomic_seq_len as u32 - end5 + 1;
            end3 = genomic_seq_len as u32 - end3 + 1;
        }

        let applied_orient = if genomic_strand == '-' { '+' } else { orient };

        let chain = acc_to_chain
            .entry(key.clone())
            .or_insert_with(|| EvidenceChain {
                accession: key.clone(),
                target: None,
                ev_type: ev_type.clone(),
                ev_class: ev_class.clone(),
                lend: u32::MAX,
                rend: 0,
                links: Vec::new(),
                gaps: Vec::new(),
                applied_orient,
            });

        if let Some(t) = target {
            chain.target = Some(t);
        }

        chain.links.push((end5, end3));
    }

    // Finalise each chain
    let mut chains: Vec<EvidenceChain> = Vec::new();
    for (_, mut chain) in acc_to_chain {
        chain.links.sort_by_key(|&(a, _)| a);
        let all_coords: Vec<u32> = chain.links.iter().flat_map(|&(a, b)| [a, b]).collect();
        chain.lend = *all_coords.iter().min().unwrap_or(&0);
        chain.rend = *all_coords.iter().max().unwrap_or(&0);

        // Build gap list. Perl re-sorts each link's coordinates (min,max) before
        // computing the gap, so this is correct for reverse-orient chains in the
        // '?' parse where links are stored end5>end3.
        for w in chain.links.windows(2) {
            let prev_max = w[0].0.max(w[0].1);
            let curr_min = w[1].0.min(w[1].1);
            let gap_l = prev_max + 1;
            let gap_r = curr_min.saturating_sub(1);
            if gap_r > gap_l {
                chain.gaps.push((gap_l, gap_r));
            }
        }
        chains.push(chain);
    }
    chains
}

/// Protein-alignment gaps at most `INTRON_MEDIAN_FACTOR × median gap` are treated
/// as inferred (imperfect) introns and decrement the coding coverage.
const INTRON_MEDIAN_FACTOR: f64 = 2.0;
/// Alignment gaps shorter than this are simple alignment gaps, not introns, and
/// are excluded from the median-gap estimate.
const MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON: u32 = 30;

/// Median of a slice of gap lengths (Perl `median`: mean of the two central
/// values for an even count, the central value for an odd count; 0 if empty).
fn median(nums: &[u32]) -> f64 {
    if nums.is_empty() {
        return 0.0;
    }
    let mut v: Vec<u32> = nums.to_vec();
    v.sort_unstable();
    let n = v.len();
    let mid = n / 2;
    if n.is_multiple_of(2) {
        (v[mid - 1] as f64 + v[mid] as f64) / 2.0
    } else {
        v[mid] as f64
    }
}

/// Faithful port of Perl `decrement_coding_using_protein_alignment_introns`.
///
/// Protein alignments whose internal gaps look like introns (length at most
/// `INTRON_MEDIAN_FACTOR × median gap`, where the median is taken over *all*
/// strand gaps ≥ `MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON`) subtract the alignment's
/// weight from the coding-score vector over the gap span — preventing protein
/// alignment coverage from scoring across introns as if it were coding. Runs
/// per strand after protein evidence is loaded, before peak analysis.
pub fn decrement_coding_using_protein_alignment_introns(
    records: &[Gff3Record],
    ev_weights: &EvWeightMap,
    mask: &MaskVec,
    genomic_seq_len: usize,
    genomic_strand: char,
    coding_scores: &mut CodingScores,
) {
    // Median gap length is computed over all chains (strand '?'), matching Perl.
    let all_chains = parse_evidence_chains('?', records, ev_weights, genomic_seq_len);
    let mut gap_lengths: Vec<u32> = Vec::new();
    for chain in &all_chains {
        for &(l, r) in &chain.gaps {
            let len = r - l + 1;
            if len >= MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON {
                gap_lengths.push(len);
            }
        }
    }
    let max_gap_length = INTRON_MEDIAN_FACTOR * median(&gap_lengths);

    // Apply decrements using strand-specific chains.
    let chains = parse_evidence_chains(genomic_strand, records, ev_weights, genomic_seq_len);
    for chain in &chains {
        let weight = ev_weights
            .get(&chain.ev_type)
            .map(|e| e.weight)
            .unwrap_or(0.0);
        for &(end5, end3) in &chain.gaps {
            let gap_length = end3.abs_diff(end5) + 1;
            if (gap_length as f64) <= max_gap_length {
                add_match_coverage(coding_scores, mask, end5, end3, -weight, &chain.ev_class);
            }
        }
    }
}

/// Instantiate evidence-based exons from the parsed chains.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::ptr_arg)]
pub fn instantiate_evidence_based_exons(
    chains: &[EvidenceChain],
    begins: &mut Vec<f64>,
    ends: &mut Vec<f64>,
    genome_features: &FeatureVec,
    _genome_seq: &[u8],
    mask: &MaskVec,
    ev_weights: &EvWeightMap,
    genomic_strand: char,
    coding_scores: &mut CodingScores,
    introns_to_score: &mut IntronScoreBuilder,
    introns_to_evidence: &mut IntronEvidenceMap,
    predicted_introns: &mut PredictedIntronMap,
    exons: &mut Vec<Exon>,
    exons_via_coords: &mut HashMap<String, usize>,
    min_intron_length: u32,
    genomic_seq_len: usize,
) {
    for chain in chains {
        let weight = ev_weights
            .get(&chain.ev_type)
            .map(|e| e.weight)
            .unwrap_or(1.0);
        let accession = if let Some(t) = &chain.target {
            format!("{}/Target={}", chain.accession, t)
        } else {
            chain.accession.clone()
        };

        // PROTEIN chains increment begin/end peaks
        if chain.ev_class == EvClass::Protein {
            if (chain.lend as usize) < begins.len() {
                begins[chain.lend as usize] += weight;
            }
            if (chain.rend as usize) < ends.len() {
                ends[chain.rend as usize] += weight;
            }
        }

        let links = &chain.links;
        let num_links = links.len();

        for (link_idx, &(end5, end3)) in links.iter().enumerate() {
            // Internal alignment segments (not first/last) can contribute exons
            // if they have proper splice boundaries
            let got_acceptor = genome_features.get((end5 as usize).saturating_sub(2))
                == FEAT_ACCEPTOR
                && link_idx != 0
                && link_idx != num_links - 1;
            let got_donor = genome_features.get((end3 + 1) as usize) == FEAT_DONOR
                && link_idx != 0
                && link_idx != num_links - 1;

            if got_donor && got_acceptor {
                let good_phases = determine_good_phases(genome_features, end5, end3);
                for phase in &good_phases {
                    let ef = end_frame(*phase, end3 - end5 + 1);
                    let coord_key = format!("{}_{}_{:?}_{}", end5, end3, ExonType::Internal, phase);
                    if let Some(&idx) = exons_via_coords.get(&coord_key) {
                        exons[idx].append_evidence(&accession, &chain.ev_type);
                    } else {
                        let mut exon = Exon::new(end5, end3);
                        exon.exon_type = ExonType::Internal;
                        exon.orientation = Orientation::Fwd;
                        exon.refresh_type_orient();
                        exon.start_frame = *phase;
                        exon.end_frame = ef;
                        exon.append_evidence(&accession, &chain.ev_type);
                        let i = exons.len();
                        exons.push(exon);
                        exons_via_coords.insert(coord_key, i);
                    }

                    // TRANSCRIPT internal exons with ORF also contribute to coding
                    if chain.ev_class == EvClass::Transcript {
                        add_match_coverage(
                            coding_scores,
                            mask,
                            end5,
                            end3,
                            weight,
                            &chain.ev_class,
                        );
                    }
                }
            }

            // PROTEIN chains always contribute to coding coverage
            if chain.ev_class == EvClass::Protein {
                add_match_coverage(coding_scores, mask, end5, end3, weight, &chain.ev_class);
            }
        }

        // Add introns from chain
        add_introns(
            &accession,
            links,
            genomic_strand,
            weight,
            &chain.ev_type,
            &chain.ev_class,
            min_intron_length,
            genome_features,
            mask,
            introns_to_score,
            introns_to_evidence,
            predicted_introns,
            genomic_seq_len,
        );
    }
}
