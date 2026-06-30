//! End-to-end golden test for eliminated-model production.
//!
//! The `elm` fixture contains a tiny genome with one strong single-exon gene
//! (survives the low-support filter) and one short single-exon gene
//! (coding length < 150, eliminated under STANDARD mode). Running with
//! `--report_ELM` produces both predictions; the eliminated model is then
//! converted to GFF3 with source `EVM_elm`.

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
    std::env::temp_dir().join(format!("evm_elm_{}_{}_{}", id, n, prefix))
}

#[test]
fn single_partition_produces_eliminated_model() {
    let params = SinglePartitionParams {
        stop_codons: "TAA,TGA,TAG".to_string(),
        min_intron_length: 20,
        forward_only: false,
        reverse_only: false,
        report_elm: true,
        terminal_intergenic_re_search: 10000,
        intergenic_adjust: 1.0,
        max_prev_exons_compare: 500,
        repeats: None,
    };

    let output = run_single_partition(
        fixture("elm.genome.fasta").to_str().unwrap(),
        fixture("elm.weights.txt").to_str().unwrap(),
        fixture("elm.gene_predictions.gff3").to_str().unwrap(),
        None,
        None,
        &params,
    )
    .expect("run_single_partition should succeed on elm fixture");

    let got_evm_out = tmp("elm.evm.out");
    fs::write(&got_evm_out, output.concat()).unwrap();

    let want_evm_out = fs::read_to_string(fixture("Contig_elm.perl.evm.out")).unwrap();
    assert_eq!(
        fs::read_to_string(&got_evm_out).unwrap(),
        want_evm_out,
        "Rust evm.out differs from Perl golden on eliminated-model fixture"
    );

    // Convert the Rust evm.out to GFF3 and compare to the Perl-derived golden.
    let got_gff3 = tmp("elm.evm.out.gff3");
    evm_output_to_gff3(
        got_evm_out.to_str().unwrap(),
        "Contig_elm",
        got_gff3.to_str().unwrap(),
    )
    .unwrap();

    let want_gff3 = fs::read_to_string(fixture("Contig_elm.perl.EVM.gff3")).unwrap();
    assert_eq!(
        fs::read_to_string(&got_gff3).unwrap(),
        want_gff3,
        "Rust GFF3 differs from Perl golden on eliminated-model fixture"
    );

    let _ = fs::remove_file(&got_evm_out);
    let _ = fs::remove_file(&got_gff3);
}
