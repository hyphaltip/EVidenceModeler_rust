//! Phase / reading-frame utilities.

use crate::types::exon::ExonPhase;
use crate::types::genome::{FeatureVec, FEAT_STOP};

/// Determine which reading phases (1, 2, 3) are valid for an exon spanning
/// [end5, end3] by checking whether any in-frame stop codon would be created.
///
/// Phase 1 means the first codon begins at end5.
/// Phase 2 means one extra base precedes the first codon (i.e. end5+1 starts a codon).
/// Phase 3 means two extra bases precede the first codon.
///
/// A stop codon at position i (1-based, pointing to the first nt of the stop)
/// eliminates the phase for which position i is in-frame.
pub fn determine_good_phases(genome_features: &FeatureVec, end5: u32, end3: u32) -> Vec<ExonPhase> {
    let mut phase_ok = [true; 3]; // index 0 = phase1, 1 = phase2, 2 = phase3

    // Walk through potential stop-codon positions within [end5, end3-2]
    // (stop codon occupies 3 bases so last possible start is end3-2)
    if end3 < end5 + 2 {
        // Exon too short to contain a stop
        return vec![1, 2, 3];
    }

    for i in end5..=(end3 - 2) {
        if genome_features.get(i as usize) == FEAT_STOP {
            let delta = (i - end5) as usize;
            let phase = match delta % 3 {
                0 => 0, // phase 1 affected
                1 => 2, // phase 3 affected (delta%3==1 ↔ phase 3)
                2 => 1, // phase 2 affected
                _ => unreachable!(),
            };
            phase_ok[phase] = false;
        }
    }

    let mut good = Vec::new();
    if phase_ok[0] {
        good.push(1u8);
    }
    if phase_ok[1] {
        good.push(2u8);
    }
    if phase_ok[2] {
        good.push(3u8);
    }
    good
}

/// Return true if the 3-byte slice represents a stop codon.
pub fn is_stop_codon(triplet: &[u8], stop_codons: &[[u8; 3]]) -> bool {
    if triplet.len() < 3 {
        return false;
    }
    let t = [
        triplet[0].to_ascii_uppercase(),
        triplet[1].to_ascii_uppercase(),
        triplet[2].to_ascii_uppercase(),
    ];
    stop_codons.contains(&t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algo::splice_sites::{populate_genome_features, DEFAULT_STOP_CODONS};

    #[test]
    fn no_stops_all_phases_good() {
        // ATGCCC — no stop codon
        let seq = b"ATGCCC";
        let fv = populate_genome_features(seq, &DEFAULT_STOP_CODONS);
        let phases = determine_good_phases(&fv, 1, 6);
        assert_eq!(phases, vec![1, 2, 3]);
    }

    #[test]
    fn stop_eliminates_phase() {
        // ATG TAA CCC
        // stop (TAA) at pos 4 → delta from end5=1 is 3 → 3%3=0 → phase1 eliminated
        let seq = b"ATGTAACCC";
        let fv = populate_genome_features(seq, &DEFAULT_STOP_CODONS);
        // end5=1, end3=9, stop at 4 (1-based)
        let phases = determine_good_phases(&fv, 1, 9);
        assert!(!phases.contains(&1), "phase 1 should be eliminated");
        assert!(phases.contains(&2));
        assert!(phases.contains(&3));
    }

    #[test]
    fn is_stop() {
        let stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
        assert!(is_stop_codon(b"TAA", &stops));
        assert!(is_stop_codon(b"tga", &stops));
        assert!(!is_stop_codon(b"ATG", &stops));
    }
}
