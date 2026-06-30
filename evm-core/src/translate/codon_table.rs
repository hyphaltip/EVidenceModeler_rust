//! Codon translation table and sequence translation utilities.

use std::collections::HashMap;

/// Build the standard genetic code codon table.
fn build_standard_table() -> HashMap<[u8; 3], char> {
    let raw: &[(&str, char)] = &[
        // Phe
        ("TTT", 'F'),
        ("TTC", 'F'),
        // Leu
        ("TTA", 'L'),
        ("TTG", 'L'),
        ("CTT", 'L'),
        ("CTC", 'L'),
        ("CTA", 'L'),
        ("CTG", 'L'),
        // Ile
        ("ATT", 'I'),
        ("ATC", 'I'),
        ("ATA", 'I'),
        // Met
        ("ATG", 'M'),
        // Val
        ("GTT", 'V'),
        ("GTC", 'V'),
        ("GTA", 'V'),
        ("GTG", 'V'),
        // Ser
        ("TCT", 'S'),
        ("TCC", 'S'),
        ("TCA", 'S'),
        ("TCG", 'S'),
        ("AGT", 'S'),
        ("AGC", 'S'),
        // Pro
        ("CCT", 'P'),
        ("CCC", 'P'),
        ("CCA", 'P'),
        ("CCG", 'P'),
        // Thr
        ("ACT", 'T'),
        ("ACC", 'T'),
        ("ACA", 'T'),
        ("ACG", 'T'),
        // Ala
        ("GCT", 'A'),
        ("GCC", 'A'),
        ("GCA", 'A'),
        ("GCG", 'A'),
        // Tyr
        ("TAT", 'Y'),
        ("TAC", 'Y'),
        // Stop
        ("TAA", '*'),
        ("TAG", '*'),
        ("TGA", '*'),
        // His
        ("CAT", 'H'),
        ("CAC", 'H'),
        // Gln
        ("CAA", 'Q'),
        ("CAG", 'Q'),
        // Asn
        ("AAT", 'N'),
        ("AAC", 'N'),
        // Lys
        ("AAA", 'K'),
        ("AAG", 'K'),
        // Asp
        ("GAT", 'D'),
        ("GAC", 'D'),
        // Glu
        ("GAA", 'E'),
        ("GAG", 'E'),
        // Cys
        ("TGT", 'C'),
        ("TGC", 'C'),
        // Trp
        ("TGG", 'W'),
        // Arg
        ("CGT", 'R'),
        ("CGC", 'R'),
        ("CGA", 'R'),
        ("CGG", 'R'),
        ("AGA", 'R'),
        ("AGG", 'R'),
        // Gly
        ("GGT", 'G'),
        ("GGC", 'G'),
        ("GGA", 'G'),
        ("GGG", 'G'),
    ];

    let mut map = HashMap::new();
    for &(codon, aa) in raw {
        let bytes = codon.as_bytes();
        map.insert([bytes[0], bytes[1], bytes[2]], aa);
    }
    map
}

/// Translate a CDS byte sequence to a protein string.
///
/// Translation starts at position 0 and proceeds in triplets.
/// Unknown codons are rendered as 'X'; stop codons as '*'.
/// `stop_codons` optionally overrides which triplets terminate translation.
/// If `stop_codons` is empty, the standard table stop codons are used.
pub fn translate(cds: &[u8], stop_codons: &[[u8; 3]]) -> String {
    let table = build_standard_table();
    let mut protein = String::new();
    let custom_stops = !stop_codons.is_empty();

    for chunk in cds.chunks(3) {
        if chunk.len() < 3 {
            break;
        }
        let codon = [
            chunk[0].to_ascii_uppercase(),
            chunk[1].to_ascii_uppercase(),
            chunk[2].to_ascii_uppercase(),
        ];

        if custom_stops {
            // Only the specified codons terminate translation
            if stop_codons.contains(&codon) {
                protein.push('*');
                break;
            }
            // Standard-table stops are reassigned to unknown when custom stops override
            let aa = table.get(&codon).copied().unwrap_or('X');
            if aa == '*' {
                // Reassigned stop: treat as unknown amino acid
                protein.push('X');
            } else {
                protein.push(aa);
            }
        } else {
            // Standard genetic code: stop at any '*' in the table
            let aa = table.get(&codon).copied().unwrap_or('X');
            if aa == '*' {
                protein.push('*');
                break;
            }
            protein.push(aa);
        }
    }

    protein
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translate_met() {
        let cds = b"ATG";
        let prot = translate(cds, &[]);
        assert_eq!(prot, "M");
    }

    #[test]
    fn translate_stop() {
        let cds = b"ATGTAA";
        let prot = translate(cds, &[]);
        assert_eq!(prot, "M*");
    }

    #[test]
    fn translate_custom_stop() {
        // With TGA as the only stop, TAA is reassigned to 'X' (not a terminator)
        let cds = b"ATGTAACCC";
        let stops = vec![*b"TGA"];
        let prot = translate(cds, &stops);
        assert_eq!(prot, "MXP"); // ATG=M, TAA=X (reassigned stop), CCC=P
                                 // TGA still stops translation
        let cds2 = b"ATGTGA";
        let prot2 = translate(cds2, &stops);
        assert_eq!(prot2, "M*");
    }

    #[test]
    fn translate_full_orf() {
        // ATG GGC TTT TAA
        let cds = b"ATGGGCTTTAA";
        let prot = translate(cds, &[]);
        // ATG=M, GGC=G, TTT=F; 11 chars => last codon AA is incomplete -> ignored
        assert_eq!(prot, "MGF");
        let cds2 = b"ATGGGCTTTTAA";
        let prot2 = translate(cds2, &[]);
        // ATG=M, GGC=G, TTT=F, TAA=* (stop, not included when in stop_codons list)
        assert_eq!(prot2, "MGF*");
    }
}
