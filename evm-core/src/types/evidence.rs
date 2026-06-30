//! Evidence types, weights, and chains.

use std::collections::HashMap;

/// Evidence class categories recognised by EVM.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum EvClass {
    Protein,
    Transcript,
    AbinitioPrediction,
    OtherPrediction,
}

impl EvClass {
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<EvClass> {
        match s {
            "PROTEIN" => Some(EvClass::Protein),
            "TRANSCRIPT" => Some(EvClass::Transcript),
            "ABINITIO_PREDICTION" => Some(EvClass::AbinitioPrediction),
            "OTHER_PREDICTION" => Some(EvClass::OtherPrediction),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            EvClass::Protein => "PROTEIN",
            EvClass::Transcript => "TRANSCRIPT",
            EvClass::AbinitioPrediction => "ABINITIO_PREDICTION",
            EvClass::OtherPrediction => "OTHER_PREDICTION",
        }
    }

    /// True if this class contributes to the per-base coding score vector
    /// (PROTEIN and ABINITIO_PREDICTION only).
    pub fn contributes_to_coding_vec(&self) -> bool {
        matches!(self, EvClass::Protein | EvClass::AbinitioPrediction)
    }

    /// True if this class is a gene prediction type.
    pub fn is_prediction(&self) -> bool {
        matches!(self, EvClass::AbinitioPrediction | EvClass::OtherPrediction)
    }

    /// True if ab initio (scores intergenic regions).
    pub fn is_abinitio(&self) -> bool {
        matches!(self, EvClass::AbinitioPrediction)
    }
}

/// A single evidence type entry from the weights file.
#[derive(Debug, Clone)]
pub struct EvEntry {
    pub ev_class: EvClass,
    pub weight: f64,
}

/// Map from ev_type string → (EvClass, weight).
pub type EvWeightMap = HashMap<String, EvEntry>;

/// Alignment chain parsed from a GFF3 evidence file.
#[derive(Debug, Clone)]
pub struct EvidenceChain {
    /// Combined key used internally (ev_type/ID=chainID).
    pub accession: String,
    /// Target/Query hit name (if present).
    pub target: Option<String>,
    pub ev_type: String,
    pub ev_class: EvClass,
    /// Leftmost coordinate (1-based, forward strand).
    pub lend: u32,
    /// Rightmost coordinate (1-based, forward strand).
    pub rend: u32,
    /// Individual alignment blocks: each is (end5, end3) already transposed
    /// to the forward strand if the chain was on the minus strand.
    pub links: Vec<(u32, u32)>,
    /// Gap spans between consecutive links.
    pub gaps: Vec<(u32, u32)>,
    /// Orientation after any transposition.
    pub applied_orient: char,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ev_class_round_trip() {
        for s in &[
            "PROTEIN",
            "TRANSCRIPT",
            "ABINITIO_PREDICTION",
            "OTHER_PREDICTION",
        ] {
            let cls = EvClass::from_str(s).unwrap();
            assert_eq!(cls.as_str(), *s);
        }
    }
}
