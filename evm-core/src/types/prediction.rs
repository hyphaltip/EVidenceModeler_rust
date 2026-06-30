//! EVM prediction objects.

use crate::algo::introns::{make_intron_key, IntronScoreMap};
use crate::types::exon::Exon;

/// Run mode for a prediction search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredMode {
    Standard,
    Intron,
}

impl PredMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            PredMode::Standard => "STANDARD",
            PredMode::Intron => "INTRON",
        }
    }
}

/// A complete or partial consensus gene prediction produced by EVM.
#[derive(Debug, Clone)]
pub struct EvmPrediction {
    /// Indices into the exon pool for the exons making up this prediction,
    /// ordered 5′ → 3′ along the genomic sequence.
    pub exon_indices: Vec<usize>,
    /// True if the prediction was eliminated by the low-support filter.
    pub is_eliminated: bool,
    /// Mode under which the prediction was generated.
    pub mode: PredMode,
    /// Leftmost coordinate of the prediction span (1-based).
    pub lend: u32,
    /// Rightmost coordinate of the prediction span (1-based).
    pub rend: u32,
    /// Total evidence-weighted score of this prediction
    /// (Σ exon base scores + Σ intron scores), matching Perl prediction_score.
    pub total_score: f64,
    /// Strand orientation of the prediction ('+'/'-'), from its first exon.
    pub orient: char,
    /// Simple intron gap coordinates (lend, rend) between consecutive exons,
    /// ordered 5'→3'. Used by the low-support filter and the output writer.
    pub intron_coords: Vec<(u32, u32)>,
    /// Noncoding scores computed by the low-support filter (Perl parity).
    pub raw_noncoding: f64,
    pub offset_noncoding: f64,
    pub noncoding_equivalent: f64,
    pub score_ratio: f64,
}

impl EvmPrediction {
    pub fn new(exon_indices: Vec<usize>, exons: &[Exon]) -> Self {
        let (lend, rend) = span_of_indices(&exon_indices, exons);
        EvmPrediction {
            exon_indices,
            is_eliminated: false,
            mode: PredMode::Standard,
            lend,
            rend,
            total_score: 0.0,
            orient: '+',
            intron_coords: Vec::new(),
            raw_noncoding: 0.0,
            offset_noncoding: 0.0,
            noncoding_equivalent: 0.0,
            score_ratio: 0.0,
        }
    }

    /// Compute `total_score`, `orient`, and `intron_coords` from the member
    /// exons and the intron score map — mirrors the Perl `EVM_prediction::_init`
    /// (prediction_score = Σ exon coding scores + Σ intron scores; introns keyed
    /// by donor/acceptor-adjusted coordinates).
    pub fn finalize(&mut self, exons: &[Exon], introns_to_score: &IntronScoreMap) {
        // Perl sorts the prediction's exons by end5 ascending.
        self.exon_indices.sort_by_key(|&i| exons[i].end5);
        if self.exon_indices.is_empty() {
            return;
        }

        self.orient = exons[self.exon_indices[0]].orientation.as_char();

        let mut score = 0.0;
        let mut intron_coords = Vec::new();
        for w in self.exon_indices.windows(2) {
            let a = &exons[w[0]];
            let b = &exons[w[1]];
            // Donor/acceptor-adjusted key coordinates (Perl _get_exon_pair_intron_score).
            let (i5, i3) = if self.orient == '+' {
                (a.end3 + 1, b.end5.saturating_sub(2))
            } else {
                (b.end3.saturating_sub(1), a.end5 + 2)
            };
            let key = make_intron_key(i5, i3);
            if let Some(&s) = introns_to_score.get(&key) {
                score += s;
            }
            // Simple intron gap coordinates (exonA_rend+1, exonB_lend-1).
            let a_rend = a.end5.max(a.end3);
            let b_lend = b.end5.min(b.end3);
            intron_coords.push((a_rend + 1, b_lend.saturating_sub(1)));
        }
        for &i in &self.exon_indices {
            score += exons[i].base_score;
        }
        self.total_score = score;
        self.intron_coords = intron_coords;

        // Perl _init also recomputes the prediction span from member exons.
        let (lend, rend) = span_of_indices(&self.exon_indices, exons);
        self.lend = lend;
        self.rend = rend;
    }

    pub fn is_eliminated(&self) -> bool {
        self.is_eliminated
    }

    pub fn get_span(&self) -> (u32, u32) {
        (self.lend, self.rend)
    }
}

fn span_of_indices(indices: &[usize], exons: &[Exon]) -> (u32, u32) {
    let mut lend = u32::MAX;
    let mut rend = 0u32;
    for &idx in indices {
        let (el, er) = exons[idx].coords_sorted();
        if el < lend {
            lend = el;
        }
        if er > rend {
            rend = er;
        }
    }
    (lend, rend)
}

/// Intermediate prediction structure used during partition recombination.
#[derive(Debug, Clone)]
pub struct PartitionPred {
    pub lend: u32,
    pub rend: u32,
    pub class: PredClass,
    /// Serialised EVM output text for this prediction.
    pub text: String,
    /// Span length used for DP scoring.
    pub length: u32,
    /// Cumulative path score.
    pub path_score: u32,
    /// Index of predecessor prediction in the sorted array (None = no link).
    pub prev_link: Option<usize>,
    /// Predictions nested within introns of this one.
    pub intronic_preds: Vec<PartitionPred>,
    /// Set to true if this pred is encapsulated inside another pred's intron.
    pub encaps: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredClass {
    Complete,
    Partial,
}
