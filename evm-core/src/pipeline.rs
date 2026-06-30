//! Single-partition EVM driver shared by the `evidence_modeler` shim and the
//! `EVidenceModeler` orchestrator's per-partition step.
//!
//! This is the one place the forward/reverse `process_features` runs are merged,
//! the intergenic vectors are built, and the trellis/consensus is invoked — so
//! the two binaries can never drift apart (a risk flagged while building Phase C).

use crate::algo::consensus::{generate_consensus_gene_predictions, ConsensusParams};
use crate::algo::intergenic::{
    augment_intergenic_from_start_stop_peaks, populate_intergenic_scores,
};
use crate::algo::introns::{
    IntronEvidenceMap, IntronScoreBuilder, IntronScoreMap, PredictedIntronMap,
};
use crate::algo::process::{process_features, ProcessConfig, StrandState};
use crate::algo::splice_sites::parse_stop_codons;
use crate::io::fasta::read_fasta_file;
use crate::io::gff3::read_gff3_file;
use crate::io::weights::read_weights_file;
use crate::types::exon::{build_acceptable_linkages, Exon, Orientation};
use crate::types::genome::{GenomeSequence, MaskVec};
use crate::types::prediction::PredMode;
use anyhow::Result;
use std::io::Write;

/// Tunables for one single-partition EVM run.
pub struct SinglePartitionParams {
    pub stop_codons: String,
    pub min_intron_length: u32,
    pub forward_only: bool,
    pub reverse_only: bool,
    pub report_elm: bool,
    pub terminal_intergenic_re_search: u32,
    pub intergenic_adjust: f64,
    pub max_prev_exons_compare: usize,
    /// Optional repeats GFF3 file; positions covered by repeats are masked.
    pub repeats: Option<String>,
}

impl Default for SinglePartitionParams {
    fn default() -> Self {
        SinglePartitionParams {
            stop_codons: "TAA,TGA,TAG".to_string(),
            min_intron_length: 20,
            forward_only: false,
            reverse_only: false,
            report_elm: false,
            terminal_intergenic_re_search: 10000,
            intergenic_adjust: 1.0,
            max_prev_exons_compare: 500,
            repeats: None,
        }
    }
}

/// Map a forward reading frame (1,2,3) to its reverse equivalent (4,5,6) when
/// transposing reverse-strand exons back to forward coordinates.
fn fwd_frame_to_rev(frame: crate::types::exon::ExonPhase) -> crate::types::exon::ExonPhase {
    match frame {
        1 => 4,
        2 => 5,
        3 => 6,
        other => other,
    }
}

/// Run EVM on a single partition's input files and return the `evm.out` text
/// blocks (one `String` per consensus prediction block, already formatted).
pub fn run_single_partition(
    genome_path: &str,
    weights_path: &str,
    gene_pred_path: &str,
    protein_path: Option<&str>,
    transcript_path: Option<&str>,
    params: &SinglePartitionParams,
) -> Result<Vec<String>> {
    let fasta_records = read_fasta_file(genome_path)?;
    if fasta_records.is_empty() {
        anyhow::bail!("No sequences found in genome file: {}", genome_path);
    }
    let genome_seq = GenomeSequence::new(&fasta_records[0].sequence);
    let seq_len = genome_seq.len();

    let ev_weights = read_weights_file(weights_path)?;
    let stop_codons = parse_stop_codons(&params.stop_codons)?;

    let gene_pred_records = read_gff3_file(gene_pred_path)?;
    let protein_records = protein_path.map(|p| read_gff3_file(p).unwrap_or_default());
    let transcript_records = transcript_path.map(|p| read_gff3_file(p).unwrap_or_default());

    let mut mask = MaskVec::new(seq_len);
    // Mask N positions (always) and repeat regions (if provided).
    apply_n_mask(&genome_seq, &mut mask);
    if let Some(repeats_path) = &params.repeats {
        apply_repeats_mask(repeats_path, &mut mask, seq_len)?;
    }

    let sum_pred_weights: f64 = ev_weights
        .values()
        .filter(|e| e.ev_class.is_prediction())
        .map(|e| e.weight)
        .sum();

    let cfg = ProcessConfig {
        genome_seq: &genome_seq,
        stop_codons: &stop_codons,
        ev_weights: &ev_weights,
        gene_pred_records: &gene_pred_records,
        protein_records: protein_records.as_deref(),
        transcript_records: transcript_records.as_deref(),
        min_intron_length: params.min_intron_length,
        mask: &mask,
        sum_genepred_weights: sum_pred_weights,
        chain_termini_window: 250,
    };

    let fwd_state = if !params.reverse_only {
        Some(process_features('+', &cfg)?)
    } else {
        None
    };

    let rev_genome = genome_seq.to_reverse_complement();
    let rev_cfg = ProcessConfig {
        genome_seq: &rev_genome,
        ..cfg
    };
    let rev_state = if !params.forward_only {
        let mut s = process_features('-', &rev_cfg)?;
        // transpose_exons_back_to_forward_strand: revcomp coords, remap reading
        // frames 1→4/2→5/3→6, flip orientation.
        for exon in s.exons.iter_mut() {
            let new_e5 = seq_len as u32 - exon.end5 + 1;
            let new_e3 = seq_len as u32 - exon.end3 + 1;
            exon.end5 = new_e5;
            exon.end3 = new_e3;
            exon.start_frame = fwd_frame_to_rev(exon.start_frame);
            exon.end_frame = fwd_frame_to_rev(exon.end_frame);
            exon.orientation = Orientation::Rev;
            exon.refresh_type_orient();
            exon.refresh_coords();
        }
        Some(s)
    } else {
        None
    };

    // Debug dumps — per-strand (enabled via EVM_DUMP_DIR env var)
    let dump_dir_opt = std::env::var("EVM_DUMP_DIR").ok();
    if let Some(ref dump_dir) = dump_dir_opt {
        let dd = std::path::Path::new(dump_dir);
        std::fs::create_dir_all(dd).ok();
        if let Some(ref s) = fwd_state {
            dump_coding_vec(&dd.join("coding_vector.+.dat"), s);
            dump_vec(
                &dd.join("introns_decomposed_to_vec.+.dat"),
                s.fwd_intron_vec.as_slice(),
            );
            dump_vec(&dd.join("begins.+.coords"), &s.begins);
            dump_vec(&dd.join("ends.+.coords"), &s.ends);
        }
        if let Some(ref s) = rev_state {
            dump_coding_vec(&dd.join("coding_vector.-.dat"), s);
            dump_vec(
                &dd.join("introns_decomposed_to_vec.-.dat"),
                s.fwd_intron_vec.as_slice(),
            );
            dump_vec(&dd.join("begins.-.coords"), &s.begins);
            dump_vec(&dd.join("ends.-.coords"), &s.ends);
        }
    }

    // Merge forward + (transposed) reverse states.
    let mut all_exons: Vec<Exon> = Vec::new();
    let mut all_introns_to_score = IntronScoreBuilder::default();
    let mut all_introns_to_evidence: IntronEvidenceMap = IntronEvidenceMap::default();
    let mut all_predicted_introns = PredictedIntronMap::default();
    let mut all_coding_scores = crate::algo::coding_scores::new_coding_scores(seq_len);
    let mut all_start_peaks = Vec::new();
    let mut all_end_peaks = Vec::new();
    let mut all_fwd_intron_vec = vec![0.0f64; seq_len + 2];
    let mut all_rev_intron_vec = vec![0.0f64; seq_len + 2];

    for state in [fwd_state, rev_state].into_iter().flatten() {
        all_exons.extend(state.exons);
        for (k, v) in state.introns_to_score {
            *all_introns_to_score.entry(k).or_insert(0.0) += v;
        }
        for (k, v) in state.introns_to_evidence {
            all_introns_to_evidence.entry(k).or_default().extend(v);
        }
        for (k, v) in state.predicted_introns {
            *all_predicted_introns.entry(k).or_insert(0.0) += v;
        }
        for (i, &v) in state.coding_scores.iter().enumerate() {
            all_coding_scores[i] += v;
        }
        all_start_peaks.extend(state.start_peaks);
        all_end_peaks.extend(state.end_peaks);
        for (i, &v) in state.fwd_intron_vec.iter().enumerate() {
            if i < all_fwd_intron_vec.len() {
                all_fwd_intron_vec[i] += v;
            }
        }
        for (i, &v) in state.rev_intron_vec.iter().enumerate() {
            if i < all_rev_intron_vec.len() {
                all_rev_intron_vec[i] += v;
            }
        }
    }

    // Base intergenic (for the low-support filter) + a start/stop-peak-augmented
    // copy (for the trellis).
    let ig_base = populate_intergenic_scores(
        seq_len,
        &gene_pred_records,
        &ev_weights,
        &mask,
        params.intergenic_adjust,
    );
    let mut ig_scores = ig_base.clone();
    augment_intergenic_from_start_stop_peaks(
        &mut ig_scores,
        &all_start_peaks,
        &all_end_peaks,
        &all_exons,
        &mask,
        seq_len as u32,
        sum_pred_weights,
        500,
    );

    // Debug dumps — merged data
    if let Some(ref dump_dir) = dump_dir_opt {
        let dd = std::path::Path::new(dump_dir);
        dump_vec(&dd.join("intergenic.bps"), ig_base.per_base());
        dump_peaks(&dd.join("start_peaks"), &all_start_peaks);
        dump_peaks(&dd.join("end_peaks"), &all_end_peaks);
        dump_exon_list(&dd.join("exon_list.out"), &all_exons);
        dump_vec(&dd.join("pred_fwd_intron_vec.dat"), &all_fwd_intron_vec);
        dump_vec(&dd.join("pred_rev_intron_vec.dat"), &all_rev_intron_vec);
        log::info!("Debug dumps written to {}", dump_dir);
    }

    let mut output: Vec<String> = Vec::new();
    let mut recursion_count = 0usize;

    let linkage_tables = build_acceptable_linkages();
    let all_introns_to_score = IntronScoreMap::from_hashmap(all_introns_to_score);
    let mut consensus_params = ConsensusParams {
        exons: &mut all_exons,
        introns_to_score: &all_introns_to_score,
        introns_to_evidence: &all_introns_to_evidence,
        ev_weights: &ev_weights,
        mask: &mask,
        intergenic_scores: &ig_scores,
        base_intergenic_scores: &ig_base,
        linkages: &linkage_tables,
        stop_codons: &stop_codons,
        coding_scores: &all_coding_scores,
        fwd_intron_vec: &all_fwd_intron_vec,
        rev_intron_vec: &all_rev_intron_vec,
        max_prev_exons_compare: params.max_prev_exons_compare,
        min_intergenic_size_on_re_search: params.terminal_intergenic_re_search,
        min_gene_length_size_on_re_search: 0,
        report_elm: params.report_elm,
        recursion_limit: 10000,
    };

    generate_consensus_gene_predictions(
        1,
        seq_len as u32,
        PredMode::Standard,
        &mut consensus_params,
        &mut recursion_count,
        &mut output,
    )?;

    Ok(output)
}

/// Mask positions that are 'N' in the genome sequence.
fn apply_n_mask(genome_seq: &GenomeSequence, mask: &mut MaskVec) {
    let bytes = genome_seq.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'N' {
            mask.set(i + 1, true);
        }
    }
}

/// Mask positions covered by repeats in a GFF3 file.
///
/// Mirrors Perl `repeatMask`: reads all non-comment, non-blank lines, takes
/// columns 4 (start) and 5 (end), and marks every position in [start, end]
/// in the mask. Coordinates are 1-based and clipped to the sequence length.
fn apply_repeats_mask(repeats_path: &str, mask: &mut MaskVec, seq_len: usize) -> Result<()> {
    let records = read_gff3_file(repeats_path)?;
    for rec in records {
        let start = rec.start.max(1) as usize;
        let end = (rec.end as usize).min(seq_len);
        for pos in start..=end {
            mask.set(pos, true);
        }
    }
    Ok(())
}

// ── Debug dump helpers ──────────────────────────────────────────────────

fn dump_vec(path: &std::path::Path, data: &[f64]) {
    if let Ok(mut f) = std::fs::File::create(path) {
        for (i, &v) in data.iter().enumerate() {
            writeln!(f, "{}\t{}", i, v).ok();
        }
    }
}

fn dump_coding_vec(path: &std::path::Path, state: &StrandState) {
    if let Ok(mut f) = std::fs::File::create(path) {
        let cs = &state.coding_scores;
        for (i, val) in cs.iter().enumerate() {
            writeln!(f, "{}\t{}", i, val).ok();
        }
    }
}

fn dump_peaks(path: &std::path::Path, peaks: &[(u32, f64, char)]) {
    if let Ok(mut f) = std::fs::File::create(path) {
        for &(pos, score, strand) in peaks {
            writeln!(f, "{}\t{}\t{}", pos, score, strand).ok();
        }
    }
}

fn dump_exon_list(path: &std::path::Path, exons: &[Exon]) {
    if let Ok(mut f) = std::fs::File::create(path) {
        for e in exons {
            let (l, r) = e.coords_sorted();
            let mut evs: Vec<String> = e
                .evidence
                .iter()
                .map(|(a, t)| format!("{{{};{}}}", a, t))
                .collect();
            evs.sort();
            writeln!(
                f,
                "{}\t{}\t{}{}\t{}\t{}\t{}\tbase score: {}, score_per_base: {:.2}",
                l,
                r,
                e.exon_type.as_str(),
                e.orientation.as_char(),
                e.start_frame,
                e.end_frame,
                evs.join(","),
                e.base_score,
                if e.length() > 0 {
                    e.base_score / e.length() as f64
                } else {
                    0.0
                },
            )
            .ok();
        }
    }
}
