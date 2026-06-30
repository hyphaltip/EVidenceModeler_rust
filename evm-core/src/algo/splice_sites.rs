//! Splice-site and codon position scanning.

use crate::types::genome::{FeatureVec, FEAT_ACCEPTOR, FEAT_DONOR, FEAT_START, FEAT_STOP};
use memchr::memmem;

/// Find all 0-based start positions of `pattern` in `seq`.
pub fn find_all_positions(seq: &[u8], pattern: &[u8]) -> Vec<usize> {
    let finder = memmem::Finder::new(pattern);
    finder.find_iter(seq).collect()
}

/// Populate a `FeatureVec` with splice-site and codon features.
///
/// * Donor sites:    `GT` and `GC` (stored at the first nucleotide of the dinucleotide + 1)
/// * Acceptor sites: `AG` (stored at the second nucleotide + 1, i.e. end of dinucleotide)
/// * Start codons:   `ATG`
/// * Stop codons:    configurable (default TAA, TGA, TAG)
///
/// Coordinate convention mirrors the Perl:
///   donor stored at position of second nt of the exon end + 1
///   acceptor stored at position of first nt of the next exon - 2
///
/// In the Perl code:
///   GT at 0-based pos p → DONOR stored at p+2 (1-based = p+2)
///   This is `GENOME_FEATURES[$pos+1] = $DONOR` where $pos is the 0-based
///   index from `index()`.  The same convention is applied here.
pub fn populate_genome_features(seq: &[u8], stop_codons: &[[u8; 3]]) -> FeatureVec {
    let mut fv = FeatureVec::new(seq.len());

    // Donors (GT and GC): position stored as 0-based pos + 1 (= 1-based pos of the G)
    for pattern in [b"GT".as_ref(), b"GC".as_ref()] {
        for p in find_all_positions(seq, pattern) {
            // In the Perl: `GENOME_FEATURES[$pos+1] = DONOR` where $pos is 0-based.
            // We store at 1-based position = p+1.
            fv.set(p + 1, FEAT_DONOR);
        }
    }

    // Acceptors (AG)
    for p in find_all_positions(seq, b"AG") {
        fv.set(p + 1, FEAT_ACCEPTOR);
    }

    // Start codons (ATG)
    for p in find_all_positions(seq, b"ATG") {
        fv.set(p + 1, FEAT_START);
    }

    // Stop codons
    for codon in stop_codons {
        for p in find_all_positions(seq, codon.as_ref()) {
            fv.set(p + 1, FEAT_STOP);
        }
    }

    fv
}

/// Default stop codons (universal code).
pub const DEFAULT_STOP_CODONS: [[u8; 3]; 3] = [*b"TAA", *b"TGA", *b"TAG"];

/// Parse a comma-separated stop-codon list (e.g. "TAA,TGA,TAG").
pub fn parse_stop_codons(arg: &str) -> anyhow::Result<Vec<[u8; 3]>> {
    let mut out = Vec::new();
    for s in arg.split(',') {
        let s = s.trim();
        if s.len() != 3 {
            anyhow::bail!("Stop codon must be 3 nucleotides: {}", s);
        }
        let bytes = s.as_bytes();
        out.push([
            bytes[0].to_ascii_uppercase(),
            bytes[1].to_ascii_uppercase(),
            bytes[2].to_ascii_uppercase(),
        ]);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::genome::{FEAT_ACCEPTOR, FEAT_DONOR, FEAT_START, FEAT_STOP};

    #[test]
    fn find_positions() {
        let seq = b"ATGATG";
        let pos = find_all_positions(seq, b"ATG");
        assert_eq!(pos, vec![0, 3]);
    }

    #[test]
    fn feature_vec_populated() {
        // seq: ATGGTAAGCCTAA  (0-based indices)
        //      0123456789...
        // ATG at 0 → START at 1-based pos 1
        // GT  at 3 → DONOR at 1-based pos 4
        // AG  at 6 (A=6,G=7) → ACCEPTOR stored at p+1 = 7
        // TAA at 10 → STOP at 1-based pos 11
        let seq = b"ATGGTAAGCCTAA";
        let stops = vec![*b"TAA", *b"TGA", *b"TAG"];
        let fv = populate_genome_features(seq, &stops);
        assert_eq!(fv.get(1), FEAT_START);
        assert_eq!(fv.get(4), FEAT_DONOR);
        assert_eq!(fv.get(7), FEAT_ACCEPTOR);
        assert_eq!(fv.get(11), FEAT_STOP);
    }
}
