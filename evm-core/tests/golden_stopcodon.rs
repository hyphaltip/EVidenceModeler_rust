//! End-to-end golden test for alternate stop-codons parity.
//!
//! The `stopcodon` fixture contains a tiny genome whose only in-frame stop is a
//! TGA. Running with `--stop_codons TGA` (Perl-style underscores) should produce
//! the single-exon gene; running with the default TAA/TGA/TAG set would not
//! because the intended stop is ignored and the ORF runs off the contig end.

use std::fs;
use std::path::PathBuf;

use evm_core::gff3_convert::evm_to_gff3::evm_output_to_gff3;
use evm_core::pipeline::{run_single_partition, SinglePartitionParams};

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
    std::env::temp_dir().join(format!("evm_stop_{}_{}_{}", id, n, prefix))
}

#[test]
fn alternate_stop_codon_tga_produces_gene() {
    let params = SinglePartitionParams {
        stop_codons: "TGA".to_string(),
        min_intron_length: 20,
        forward_only: false,
        reverse_only: false,
        report_elm: false,
        terminal_intergenic_re_search: 10000,
        intergenic_adjust: 1.0,
        max_prev_exons_compare: 500,
        repeats: None,
    };

    let output = run_single_partition(
        fixture("stopcodon.genome.fasta").to_str().unwrap(),
        fixture("stopcodon.weights.txt").to_str().unwrap(),
        fixture("stopcodon.gene_predictions.gff3").to_str().unwrap(),
        None,
        None,
        &params,
    )
    .expect("run_single_partition should succeed on stopcodon fixture");

    let got_evm_out = tmp("stop.evm.out");
    fs::write(&got_evm_out, output.concat()).unwrap();

    let want_evm_out = fs::read_to_string(fixture("Contig_stop.perl.evm.out")).unwrap();
    assert_eq!(
        fs::read_to_string(&got_evm_out).unwrap(),
        want_evm_out,
        "Rust evm.out differs from Perl golden on --stop_codons TGA fixture"
    );

    let got_gff3 = tmp("stop.evm.out.gff3");
    evm_output_to_gff3(
        got_evm_out.to_str().unwrap(),
        "Contig_stop",
        got_gff3.to_str().unwrap(),
    )
    .unwrap();

    let want_gff3 = fs::read_to_string(fixture("Contig_stop.perl.EVM.gff3")).unwrap();
    assert_eq!(
        fs::read_to_string(&got_gff3).unwrap(),
        want_gff3,
        "Rust GFF3 differs from Perl golden on --stop_codons TGA fixture"
    );

    let _ = fs::remove_file(&got_evm_out);
    let _ = fs::remove_file(&got_gff3);
}
