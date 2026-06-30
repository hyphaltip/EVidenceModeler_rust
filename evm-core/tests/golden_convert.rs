//! Golden-file parity tests for the Phase C output converters
//! (`evm.out` → GFF3 → BED, and GFF3 → protein/CDS), against fixtures captured
//! from the Perl `EVidenceModeler` pipeline on the `testing/` Contig1 dataset.
//!
//! GFF3 and BED are asserted byte-for-byte. Protein/CDS FASTA are compared as a
//! set of (header, sequence) records because the Perl `gff3_file_to_proteins.pl`
//! emits records in hash order (non-deterministic); the per-record content is
//! exact.
//!
//! `multicontig_gff3_ordering_matches_perl` verifies that the multi-contig
//! concatenation (Contig1 then Contig2 in listing order) matches the Perl
//! `find … -regex .*evm.out.gff3 -exec cat {}` golden output.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use evm_core::gff3_convert::evm_to_gff3::evm_output_to_gff3;
use evm_core::gff3_convert::gff3_to_bed::gff3_to_bed;
use evm_core::gff3_convert::gff3_to_proteins::{extract_sequences, SeqType};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../tests/fixtures")
        .join(name)
}

fn tmp(prefix: &str) -> PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let id = std::process::id();
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!("evm_golden_{}_{}_{}", id, n, prefix))
}

/// Parse a (possibly line-wrapped) FASTA into header → concatenated sequence.
fn parse_fasta(text: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let mut header = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !header.is_empty() {
                map.insert(header.clone(), std::mem::take(&mut seq));
            }
            header = h.to_string();
        } else {
            seq.push_str(line);
        }
    }
    if !header.is_empty() {
        map.insert(header, seq);
    }
    map
}

#[test]
fn evm_out_to_gff3_matches_perl() {
    let out = tmp("evm.out.gff3");
    evm_output_to_gff3(
        fixture("Contig1.perl.evm.out").to_str().unwrap(),
        "Contig1",
        out.to_str().unwrap(),
    )
    .unwrap();
    let got = fs::read_to_string(&out).unwrap();
    let want = fs::read_to_string(fixture("Contig1.perl.EVM.gff3")).unwrap();
    assert_eq!(got, want, "GFF3 output differs from Perl golden");
    let _ = fs::remove_file(&out);
}

#[test]
fn gff3_to_bed_matches_perl() {
    // Build GFF3 first (the BED converter consumes it).
    let gff3 = tmp("for_bed.gff3");
    evm_output_to_gff3(
        fixture("Contig1.perl.evm.out").to_str().unwrap(),
        "Contig1",
        gff3.to_str().unwrap(),
    )
    .unwrap();

    let lines = gff3_to_bed(gff3.to_str().unwrap()).unwrap();
    let got = format!("{}\n", lines.join("\n"));
    let want = fs::read_to_string(fixture("Contig1.perl.EVM.bed")).unwrap();
    assert_eq!(got, want, "BED output differs from Perl golden");
    let _ = fs::remove_file(&gff3);
}

fn check_sequences(seq_type: SeqType, fixture_name: &str) {
    let gff3 = tmp("for_seq.gff3");
    evm_output_to_gff3(
        fixture("Contig1.perl.evm.out").to_str().unwrap(),
        "Contig1",
        gff3.to_str().unwrap(),
    )
    .unwrap();

    let stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
    let recs = extract_sequences(
        gff3.to_str().unwrap(),
        fixture("Contig1.genome.fasta").to_str().unwrap(),
        &seq_type,
        &stops,
    )
    .unwrap();

    let mut got: BTreeMap<String, String> = BTreeMap::new();
    for (hdr, seq) in recs {
        got.insert(hdr, seq);
    }

    let want = parse_fasta(&fs::read_to_string(fixture(fixture_name)).unwrap());
    assert_eq!(
        got, want,
        "{} records differ from Perl golden",
        fixture_name
    );
    let _ = fs::remove_file(&gff3);
}

#[test]
fn gff3_to_proteins_matches_perl() {
    check_sequences(SeqType::Prot, "Contig1.perl.EVM.pep");
}

/// Verify that eliminated predictions are emitted with source `EVM_elm`.
///
/// The EVM output format marks eliminated models with
/// ` *** ELIMINATED *** ` in the header; `EVM_to_GFF3.pl` sets the GFF3
/// source column to `EVM_elm` for those models. This fixture exercises a
/// file containing one regular and one eliminated prediction.
#[test]
fn evm_out_to_gff3_emits_eliminated_source() {
    let out = tmp("evm.out.elm.gff3");
    evm_output_to_gff3(
        fixture("Contig1.with_elm.perl.evm.out").to_str().unwrap(),
        "Contig1",
        out.to_str().unwrap(),
    )
    .unwrap();
    let got = fs::read_to_string(&out).unwrap();
    let want = fs::read_to_string(fixture("Contig1.with_elm.perl.EVM.gff3")).unwrap();
    assert_eq!(
        got, want,
        "GFF3 output differs for eliminated-model fixture"
    );
    let _ = fs::remove_file(&out);
}

#[test]
fn gff3_to_cds_matches_perl() {
    check_sequences(SeqType::Cds, "Contig1.perl.EVM.cds");
}

/// Verify multi-contig GFF3 concatenation order.
///
/// Simulates two contigs (Contig1 + Contig2, same underlying evm.out) and
/// checks that Rust produces the same GFF3 as Perl's
/// `find … -regex .*evm.out.gff3 -exec cat {}` (which returns them in
/// filesystem/alphabetical order — matching partition listing order).
#[test]
fn multicontig_gff3_ordering_matches_perl() {
    // Convert the same evm.out for both contigs into temp files.
    let gff3_c1 = tmp("multicontig_c1.gff3");
    let gff3_c2 = tmp("multicontig_c2.gff3");
    evm_output_to_gff3(
        fixture("Contig1.perl.evm.out").to_str().unwrap(),
        "Contig1",
        gff3_c1.to_str().unwrap(),
    )
    .unwrap();
    evm_output_to_gff3(
        fixture("Contig1.perl.evm.out").to_str().unwrap(),
        "Contig2",
        gff3_c2.to_str().unwrap(),
    )
    .unwrap();

    // Concatenate in listing order (Contig1 then Contig2), mirroring the
    // orchestrator's concatenate_gff3_outputs and Perl's find order.
    let got = format!(
        "{}{}",
        fs::read_to_string(&gff3_c1).unwrap(),
        fs::read_to_string(&gff3_c2).unwrap(),
    );
    let want = fs::read_to_string(fixture("multicontig.perl.EVM.gff3")).unwrap();
    assert_eq!(got, want, "multi-contig GFF3 differs from Perl golden");

    let _ = fs::remove_file(&gff3_c1);
    let _ = fs::remove_file(&gff3_c2);
}
