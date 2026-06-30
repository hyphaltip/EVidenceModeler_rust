//! Coordinates the per-strand evidence loading pipeline.

use crate::algo::coding_scores::score_exons;
use crate::algo::coding_scores::{new_coding_scores, CodingScores};
use crate::algo::introns::{
    populate_intron_vectors, IntronEvidenceMap, IntronScoreBuilder, IntronVec, PredictedIntronMap,
};
use crate::algo::load_evidence::{
    decrement_coding_using_protein_alignment_introns, instantiate_evidence_based_exons,
    parse_evidence_chains,
};
use crate::algo::load_predictions::load_prediction_data;
use crate::algo::peaks::analyze_peaks;
use crate::algo::splice_sites::populate_genome_features;
use crate::io::gff3::Gff3Record;
use crate::types::evidence::EvWeightMap;
use crate::types::exon::Exon;
use crate::types::genome::{GenomeSequence, MaskVec};
use anyhow::Result;
use std::collections::HashMap;

/// State accumulated across one strand's processing pass.
pub struct StrandState {
    pub exons: Vec<Exon>,
    pub exons_via_coords: HashMap<String, usize>,
    pub coding_scores: CodingScores,
    pub introns_to_score: IntronScoreBuilder,
    pub introns_to_evidence: IntronEvidenceMap,
    pub predicted_introns: PredictedIntronMap,
    pub begins: Vec<f64>,
    pub ends: Vec<f64>,
    /// Peaks as (forward-strand position, score, strand) — strand recorded per
    /// Perl analyze_gene_boundaries so augmentation can split by strand.
    pub start_peaks: Vec<(u32, f64, char)>,
    pub end_peaks: Vec<(u32, f64, char)>,
    pub fwd_intron_vec: IntronVec,
    pub rev_intron_vec: IntronVec,
}

impl StrandState {
    pub fn new(seq_len: usize) -> Self {
        StrandState {
            exons: Vec::new(),
            exons_via_coords: HashMap::new(),
            coding_scores: new_coding_scores(seq_len),
            introns_to_score: IntronScoreBuilder::default(),
            introns_to_evidence: IntronEvidenceMap::default(),
            predicted_introns: PredictedIntronMap::default(),
            begins: vec![0.0f64; seq_len + 2],
            ends: vec![0.0f64; seq_len + 2],
            start_peaks: Vec::new(),
            end_peaks: Vec::new(),
            fwd_intron_vec: Vec::new(),
            rev_intron_vec: Vec::new(),
        }
    }
}

/// Configuration for one EVM run on a single partition.
pub struct ProcessConfig<'a> {
    pub genome_seq: &'a GenomeSequence,
    pub stop_codons: &'a [[u8; 3]],
    pub ev_weights: &'a EvWeightMap,
    pub gene_pred_records: &'a [Gff3Record],
    pub protein_records: Option<&'a [Gff3Record]>,
    pub transcript_records: Option<&'a [Gff3Record]>,
    pub min_intron_length: u32,
    pub mask: &'a MaskVec,
    pub sum_genepred_weights: f64,
    /// Window size for peak analysis (CHAIN_TERMINI_WINDOW = 250 in Perl).
    pub chain_termini_window: usize,
}

/// Run the complete per-strand feature-processing pipeline.
///
/// Returns a `StrandState` containing all exons, introns, and peaks for that strand.
pub fn process_features(genomic_strand: char, cfg: &ProcessConfig) -> Result<StrandState> {
    let seq_len = cfg.genome_seq.len();
    let seq = cfg.genome_seq.as_bytes();

    // Build feature vector for this strand
    let genome_features = populate_genome_features(seq, cfg.stop_codons);

    let mut state = StrandState::new(seq_len);

    // Load gene predictions
    load_prediction_data(
        genomic_strand,
        cfg.gene_pred_records,
        &genome_features,
        seq,
        cfg.ev_weights,
        cfg.min_intron_length,
        cfg.mask,
        &mut state.coding_scores,
        &mut state.introns_to_score,
        &mut state.introns_to_evidence,
        &mut state.predicted_introns,
        &mut state.exons,
        &mut state.exons_via_coords,
        seq_len,
    )?;

    // Load transcript evidence (before protein, as per Perl ordering)
    if let Some(recs) = cfg.transcript_records {
        let chains = parse_evidence_chains(genomic_strand, recs, cfg.ev_weights, seq_len);
        instantiate_evidence_based_exons(
            &chains,
            &mut state.begins,
            &mut state.ends,
            &genome_features,
            seq,
            cfg.mask,
            cfg.ev_weights,
            genomic_strand,
            &mut state.coding_scores,
            &mut state.introns_to_score,
            &mut state.introns_to_evidence,
            &mut state.predicted_introns,
            &mut state.exons,
            &mut state.exons_via_coords,
            cfg.min_intron_length,
            seq_len,
        );
    }

    // Load protein evidence
    if let Some(recs) = cfg.protein_records {
        let chains = parse_evidence_chains(genomic_strand, recs, cfg.ev_weights, seq_len);
        instantiate_evidence_based_exons(
            &chains,
            &mut state.begins,
            &mut state.ends,
            &genome_features,
            seq,
            cfg.mask,
            cfg.ev_weights,
            genomic_strand,
            &mut state.coding_scores,
            &mut state.introns_to_score,
            &mut state.introns_to_evidence,
            &mut state.predicted_introns,
            &mut state.exons,
            &mut state.exons_via_coords,
            cfg.min_intron_length,
            seq_len,
        );

        // Decrement coding coverage over inferred protein-alignment introns
        // (Perl decrement_coding_using_protein_alignment_introns, runs after the
        // protein evidence load, before gene-boundary peak analysis).
        decrement_coding_using_protein_alignment_introns(
            recs,
            cfg.ev_weights,
            cfg.mask,
            seq_len,
            genomic_strand,
            &mut state.coding_scores,
        );
    }

    // Analyse gene boundary peaks
    let start_peaks_raw = analyze_peaks(
        &state.begins,
        seq_len,
        cfg.chain_termini_window,
        cfg.sum_genepred_weights,
    );
    let end_peaks_raw = analyze_peaks(
        &state.ends,
        seq_len,
        cfg.chain_termini_window,
        cfg.sum_genepred_weights,
    );

    // Transpose peaks back to forward-strand coordinates if processing reverse
    for (pos, score) in &start_peaks_raw {
        let fwd_pos = if genomic_strand == '-' {
            seq_len as u32 - pos + 1
        } else {
            *pos
        };
        state.start_peaks.push((fwd_pos, *score, genomic_strand));
    }
    for (pos, score) in &end_peaks_raw {
        let fwd_pos = if genomic_strand == '-' {
            seq_len as u32 - pos + 1
        } else {
            *pos
        };
        state.end_peaks.push((fwd_pos, *score, genomic_strand));
    }

    // Score exons
    let ev_weight_fn =
        |ev_type: &str| -> Option<f64> { cfg.ev_weights.get(ev_type).map(|e| e.weight) };
    let ev_class_fn = |ev_type: &str| -> Option<crate::types::evidence::EvClass> {
        cfg.ev_weights.get(ev_type).map(|e| e.ev_class.clone())
    };
    score_exons(
        &mut state.exons,
        &state.coding_scores,
        cfg.mask,
        &ev_weight_fn,
        &ev_class_fn,
    );

    // Build intron vectors
    let (fwd, rev) = populate_intron_vectors(&state.predicted_introns, cfg.mask, seq_len);
    state.fwd_intron_vec = fwd;
    state.rev_intron_vec = rev;

    Ok(state)
}
