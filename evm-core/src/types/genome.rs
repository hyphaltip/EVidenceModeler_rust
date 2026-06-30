//! Genome sequence and feature vector types.

/// 1-indexed feature codes stored in the feature vector parallel to the genome.
pub const FEAT_START: u8 = 1;
pub const FEAT_DONOR: u8 = 2;
pub const FEAT_ACCEPTOR: u8 = 3;
pub const FEAT_STOP: u8 = 4;

/// Uppercase DNA sequence stored as bytes.
/// All coordinate access uses 1-based indexing matching the Perl original.
#[derive(Clone)]
pub struct GenomeSequence {
    seq: Vec<u8>,
}

impl GenomeSequence {
    /// Create from a raw sequence string; uppercases the sequence.
    pub fn new(raw: &str) -> Self {
        let seq: Vec<u8> = raw.bytes().map(|b| b.to_ascii_uppercase()).collect();
        GenomeSequence { seq }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Return bytes slice (0-indexed internally).
    pub fn as_bytes(&self) -> &[u8] {
        &self.seq
    }

    /// Get a single base at 1-based coordinate.
    pub fn get(&self, pos: usize) -> u8 {
        self.seq[pos - 1]
    }

    /// Return a substring at 1-based [start, start+len).
    pub fn substr(&self, start: usize, length: usize) -> &[u8] {
        &self.seq[start - 1..start - 1 + length]
    }

    /// Reverse-complement coordinate: given a 1-based coord, return its
    /// revcomp coordinate in a sequence of this length.
    pub fn rev_coord(&self, pos: usize) -> usize {
        self.seq.len() - pos + 1
    }

    /// In-place reverse complement.
    pub fn reverse_complement(&mut self) {
        self.seq = reverse_complement_bytes(&self.seq);
    }

    /// Return a reversed-complement copy.
    pub fn to_reverse_complement(&self) -> Self {
        GenomeSequence {
            seq: reverse_complement_bytes(&self.seq),
        }
    }
}

/// Reverse-complement a byte slice of DNA.
pub fn reverse_complement_bytes(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement_base(b)).collect()
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'G' => b'C',
        b'g' => b'c',
        b'C' => b'G',
        b'c' => b'g',
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'M' => b'K',
        b'm' => b'k',
        b'K' => b'M',
        b'k' => b'm',
        b'S' => b'S',
        b's' => b's',
        b'W' => b'W',
        b'w' => b'w',
        b'H' => b'D',
        b'h' => b'd',
        b'D' => b'H',
        b'd' => b'h',
        b'B' => b'V',
        b'b' => b'v',
        b'V' => b'B',
        b'v' => b'b',
        b'N' => b'N',
        b'n' => b'n',
        b'X' => b'X',
        b'x' => b'x',
        other => other,
    }
}

/// Feature vector: 1-indexed parallel array to genome holding
/// `FEAT_*` constants at each position (0 = no feature).
#[derive(Clone)]
pub struct FeatureVec {
    data: Vec<u8>,
}

impl FeatureVec {
    pub fn new(len: usize) -> Self {
        FeatureVec {
            data: vec![0u8; len + 2], // extra padding for boundary checks
        }
    }

    /// Set feature at 1-based position.
    pub fn set(&mut self, pos: usize, feat: u8) {
        if pos < self.data.len() {
            self.data[pos] = feat;
        }
    }

    /// Get feature at 1-based position. Returns 0 if out of range.
    pub fn get(&self, pos: usize) -> u8 {
        self.data.get(pos).copied().unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// Boolean mask vector (1-indexed). True = masked (repeat / N region).
#[derive(Clone)]
pub struct MaskVec {
    data: Vec<bool>,
}

impl MaskVec {
    pub fn new(len: usize) -> Self {
        MaskVec {
            data: vec![false; len + 2],
        }
    }

    pub fn set(&mut self, pos: usize, val: bool) {
        if pos < self.data.len() {
            self.data[pos] = val;
        }
    }

    pub fn get(&self, pos: usize) -> bool {
        self.data.get(pos).copied().unwrap_or(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rev_complement() {
        let mut g = GenomeSequence::new("ATGCGT");
        let rc = g.to_reverse_complement();
        assert_eq!(rc.as_bytes(), b"ACGCAT");
        g.reverse_complement();
        assert_eq!(g.as_bytes(), b"ACGCAT");
    }

    #[test]
    fn rev_coord() {
        let g = GenomeSequence::new("ATGCGT"); // len=6
        assert_eq!(g.rev_coord(1), 6);
        assert_eq!(g.rev_coord(6), 1);
    }

    #[test]
    fn feature_vec_bounds() {
        let mut f = FeatureVec::new(10);
        f.set(1, FEAT_START);
        assert_eq!(f.get(1), FEAT_START);
        assert_eq!(f.get(0), 0);
        assert_eq!(f.get(100), 0);
    }
}
