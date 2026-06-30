//! Exon data type and exon pool.

/// Classification of an exon within a gene model.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum ExonType {
    Initial,
    Internal,
    Terminal,
    Single,
    /// Boundary sentinel node used only inside the trellis.
    Bound,
}

impl ExonType {
    pub fn as_str(&self) -> &'static str {
        match self {
            ExonType::Initial => "initial",
            ExonType::Internal => "internal",
            ExonType::Terminal => "terminal",
            ExonType::Single => "single",
            ExonType::Bound => "bound",
        }
    }
}

/// Strand orientation.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum Orientation {
    Fwd,
    Rev,
}

impl Orientation {
    pub fn as_char(&self) -> char {
        match self {
            Orientation::Fwd => '+',
            Orientation::Rev => '-',
        }
    }
}

/// Phase: 1-3 for forward strand, 4-6 for reverse strand.
/// Encodes reading frame continuity at exon boundaries.
pub type ExonPhase = u8;

/// Compact encoding of exon type + strand orientation used for linkage lookups.
/// Values are chosen to fit in a u8 and are stable across the codebase.
pub type TypeOrient = u8;

pub const TYPE_ORIENT_INITIAL_FWD: TypeOrient = 0;
pub const TYPE_ORIENT_INTERNAL_FWD: TypeOrient = 1;
pub const TYPE_ORIENT_TERMINAL_FWD: TypeOrient = 2;
pub const TYPE_ORIENT_SINGLE_FWD: TypeOrient = 3;
pub const TYPE_ORIENT_INITIAL_REV: TypeOrient = 4;
pub const TYPE_ORIENT_INTERNAL_REV: TypeOrient = 5;
pub const TYPE_ORIENT_TERMINAL_REV: TypeOrient = 6;
pub const TYPE_ORIENT_SINGLE_REV: TypeOrient = 7;
pub const TYPE_ORIENT_BOUND: TypeOrient = 8;

pub fn type_orient_from(exon_type: ExonType, orientation: Orientation) -> TypeOrient {
    match (exon_type, orientation) {
        (ExonType::Initial, Orientation::Fwd) => TYPE_ORIENT_INITIAL_FWD,
        (ExonType::Internal, Orientation::Fwd) => TYPE_ORIENT_INTERNAL_FWD,
        (ExonType::Terminal, Orientation::Fwd) => TYPE_ORIENT_TERMINAL_FWD,
        (ExonType::Single, Orientation::Fwd) => TYPE_ORIENT_SINGLE_FWD,
        (ExonType::Initial, Orientation::Rev) => TYPE_ORIENT_INITIAL_REV,
        (ExonType::Internal, Orientation::Rev) => TYPE_ORIENT_INTERNAL_REV,
        (ExonType::Terminal, Orientation::Rev) => TYPE_ORIENT_TERMINAL_REV,
        (ExonType::Single, Orientation::Rev) => TYPE_ORIENT_SINGLE_REV,
        (ExonType::Bound, _) => TYPE_ORIENT_BOUND,
    }
}

pub fn type_orient_to_str(t: TypeOrient) -> &'static str {
    match t {
        TYPE_ORIENT_INITIAL_FWD => "initial+",
        TYPE_ORIENT_INTERNAL_FWD => "internal+",
        TYPE_ORIENT_TERMINAL_FWD => "terminal+",
        TYPE_ORIENT_SINGLE_FWD => "single+",
        TYPE_ORIENT_INITIAL_REV => "initial-",
        TYPE_ORIENT_INTERNAL_REV => "internal-",
        TYPE_ORIENT_TERMINAL_REV => "terminal-",
        TYPE_ORIENT_SINGLE_REV => "single-",
        _ => "bound",
    }
}

/// Compute the end frame given a start frame and exon length.
/// Phase cycles 1→2→3→1 (fwd) or 4→5→6→4 (rev).
pub fn end_frame(start_frame: ExonPhase, exon_len: u32) -> ExonPhase {
    let base = if start_frame <= 3 { 0u8 } else { 3u8 };
    let zero_frame = (start_frame - 1 - base) as u32;
    let end_zero = (zero_frame + exon_len - 1) % 3;
    base + end_zero as u8 + 1
}

/// GFF3 phase conversion: (gff_phase_digit, strand_char) → ExonPhase [1..=6].
pub fn gff_phase_to_exon_phase(gff_phase: u8, orient: char) -> ExonPhase {
    match (gff_phase, orient) {
        (0, '+') => 1,
        (1, '+') => 2,
        (2, '+') => 3,
        (0, '-') => 4,
        (1, '-') => 5,
        (2, '-') => 6,
        _ => 1,
    }
}

/// EVM phase back to GFF3 phase (0-2).
pub fn exon_phase_to_gff_phase(phase: ExonPhase) -> u8 {
    match phase {
        1 | 4 => 0,
        2 | 5 => 1,
        3 | 6 => 2,
        _ => 0,
    }
}

/// A candidate exon stored in the exon pool.
#[derive(Debug, Clone)]
pub struct Exon {
    /// 5-prime coordinate (1-based, in the **forward** strand coordinate
    /// system after transposition).  For forward exons end5 < end3;
    /// for reverse exons end5 > end3 (stored that way for consistency
    /// with the Perl original).
    pub end5: u32,
    pub end3: u32,
    pub exon_type: ExonType,
    pub orientation: Orientation,
    pub start_frame: ExonPhase,
    pub end_frame: ExonPhase,
    /// Sorted genomic coordinates (lend <= rend), cached for fast overlap checks.
    pub lend: u32,
    pub rend: u32,
    /// Cumulative evidence-weighted score for the exon (set by `score_exons`).
    pub base_score: f64,
    /// Running best path score through this exon (set during trellis build).
    pub sum_score: f64,
    /// Two-character dinucleotide at the left (5') boundary of the exon
    /// (used for checking stop-codon creation across splice junctions).
    pub left_seq_boundary: [u8; 2],
    /// Two-character dinucleotide at the right (3') boundary of the exon.
    pub right_seq_boundary: [u8; 2],
    /// Evidence supporting this exon: Vec of (accession, ev_type).
    pub evidence: Vec<(String, String)>,
    /// Index into the exon pool of the best predecessor in the trellis.
    /// None = no predecessor / boundary.
    pub link: Option<usize>,
    /// Compact (exon_type, orientation) code for fast linkage lookups.
    pub type_orient: TypeOrient,
}

impl Exon {
    pub fn new(end5: u32, end3: u32) -> Self {
        let (lend, rend) = if end5 <= end3 {
            (end5, end3)
        } else {
            (end3, end5)
        };
        Exon {
            end5,
            end3,
            exon_type: ExonType::Single,
            orientation: Orientation::Fwd,
            start_frame: 1,
            end_frame: 1,
            lend,
            rend,
            base_score: 0.0,
            sum_score: 0.0,
            left_seq_boundary: [0, 0],
            right_seq_boundary: [0, 0],
            evidence: Vec::new(),
            link: None,
            type_orient: TYPE_ORIENT_SINGLE_FWD,
        }
    }

    /// Recompute the compact type+orientation code from fields.
    pub fn refresh_type_orient(&mut self) {
        self.type_orient = type_orient_from(self.exon_type, self.orientation);
    }

    /// Recompute cached sorted coordinates after a coordinate transposition.
    pub fn refresh_coords(&mut self) {
        if self.end5 <= self.end3 {
            self.lend = self.end5;
            self.rend = self.end3;
        } else {
            self.lend = self.end3;
            self.rend = self.end5;
        }
    }

    /// Sorted (lend, rend) regardless of strand orientation.
    pub fn coords_sorted(&self) -> (u32, u32) {
        (self.lend, self.rend)
    }

    pub fn length(&self) -> u32 {
        let (l, r) = self.coords_sorted();
        r - l + 1
    }

    pub fn type_orient_key(&self) -> String {
        type_orient_to_str(self.type_orient).to_string()
    }

    pub fn append_evidence(&mut self, accession: impl Into<String>, ev_type: impl Into<String>) {
        self.evidence.push((accession.into(), ev_type.into()));
    }
}

const NUM_TYPE_ORIENTS: usize = 9;
const NUM_PHASES: usize = 7; // phases 1..6

/// Fast O(1) linkage lookup tables indexed by compact type/orient or phase codes.
pub struct LinkageTables {
    pub acceptable: [[bool; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS],
    pub phased: [[bool; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS],
    pub intergenic: [[bool; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS],
    pub frame_pairs: [[bool; NUM_PHASES]; NUM_PHASES],
}

impl LinkageTables {
    pub fn contains_acceptable(&self, a: TypeOrient, b: TypeOrient) -> bool {
        self.acceptable[a as usize][b as usize]
    }
    pub fn contains_phased(&self, a: TypeOrient, b: TypeOrient) -> bool {
        self.phased[a as usize][b as usize]
    }
    pub fn contains_intergenic(&self, a: TypeOrient, b: TypeOrient) -> bool {
        self.intergenic[a as usize][b as usize]
    }
    pub fn contains_frame_pair(&self, a: ExonPhase, b: ExonPhase) -> bool {
        self.frame_pairs[a as usize][b as usize]
    }
}

/// Acceptable exon linkage table: (typeA_orient, typeB_orient) → phased?.
/// This mirrors the `@acceptableExonLinkages` table from the Perl.
pub fn build_acceptable_linkages() -> LinkageTables {
    fn tok(s: &str) -> TypeOrient {
        match s {
            "initial+" => TYPE_ORIENT_INITIAL_FWD,
            "internal+" => TYPE_ORIENT_INTERNAL_FWD,
            "terminal+" => TYPE_ORIENT_TERMINAL_FWD,
            "single+" => TYPE_ORIENT_SINGLE_FWD,
            "initial-" => TYPE_ORIENT_INITIAL_REV,
            "internal-" => TYPE_ORIENT_INTERNAL_REV,
            "terminal-" => TYPE_ORIENT_TERMINAL_REV,
            "single-" => TYPE_ORIENT_SINGLE_REV,
            _ => TYPE_ORIENT_BOUND,
        }
    }

    let mut all = [[false; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS];
    let mut phased = [[false; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS];

    let raw: &[(&str, &str, bool)] = &[
        // Forward strand
        ("initial+", "terminal+", true),
        ("initial+", "internal+", true),
        ("internal+", "internal+", true),
        ("internal+", "terminal+", true),
        ("terminal+", "initial+", false),
        ("terminal+", "single+", false),
        ("single+", "single+", false),
        ("single+", "initial+", false),
        // Reverse strand
        ("terminal-", "initial-", true),
        ("internal-", "initial-", true),
        ("internal-", "internal-", true),
        ("terminal-", "internal-", true),
        ("initial-", "terminal-", false),
        ("single-", "terminal-", false),
        ("single-", "single-", false),
        ("initial-", "single-", false),
        // Fwd → Rev transitions
        ("single+", "single-", false),
        ("single+", "terminal-", false),
        ("terminal+", "terminal-", false),
        ("terminal+", "single-", false),
        // Rev → Fwd transitions
        ("single-", "single+", false),
        ("single-", "initial+", false),
        ("initial-", "initial+", false),
        ("initial-", "single+", false),
    ];
    for (a, b, p) in raw {
        let ai = tok(a) as usize;
        let bi = tok(b) as usize;
        all[ai][bi] = true;
        if *p {
            phased[ai][bi] = true;
        }
    }

    let mut intergenic = [[false; NUM_TYPE_ORIENTS]; NUM_TYPE_ORIENTS];
    let intergenic_raw: &[(&str, &str)] = &[
        ("terminal+", "initial+"),
        ("terminal+", "single+"),
        ("single+", "single+"),
        ("single+", "initial+"),
        ("initial-", "terminal-"),
        ("single-", "terminal-"),
        ("single-", "single-"),
        ("initial-", "single-"),
        ("single+", "single-"),
        ("single+", "terminal-"),
        ("terminal+", "terminal-"),
        ("terminal+", "single-"),
        ("single-", "single+"),
        ("single-", "initial+"),
        ("initial-", "initial+"),
        ("initial-", "single+"),
    ];
    for (a, b) in intergenic_raw {
        intergenic[tok(a) as usize][tok(b) as usize] = true;
    }

    let mut frame_pairs = [[false; NUM_PHASES]; NUM_PHASES];
    for (a, b) in &[
        (1u8, 2u8),
        (2, 3),
        (3, 1), // fwd
        (4, 5),
        (5, 6),
        (6, 4), // rev
    ] {
        frame_pairs[*a as usize][*b as usize] = true;
    }

    LinkageTables {
        acceptable: all,
        phased,
        intergenic,
        frame_pairs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_frame_cycle() {
        // phase1, len=3: zero_frame=0, end_zero=(0+3-1)%3=2 → result=3
        assert_eq!(end_frame(1, 3), 3);
        // phase1, len=4: end_zero=(0+4-1)%3=0 → result=1
        assert_eq!(end_frame(1, 4), 1);
        // phase2, len=3: zero_frame=1, end_zero=(1+3-1)%3=3%3=0 → result=1
        assert_eq!(end_frame(2, 3), 1);
        // phase3, len=3: zero_frame=2, end_zero=(2+3-1)%3=4%3=1 → result=2
        assert_eq!(end_frame(3, 3), 2);
    }

    #[test]
    fn exon_coords_sorted() {
        let mut e = Exon::new(100, 50);
        let (l, r) = e.coords_sorted();
        assert_eq!(l, 50);
        assert_eq!(r, 100);
        e.end5 = 50;
        e.end3 = 100;
        e.refresh_coords();
        let (l, r) = e.coords_sorted();
        assert_eq!(l, 50);
        assert_eq!(r, 100);
    }
}
