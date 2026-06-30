//! EVidenceModeler — Rust re-implementation of the main EVM orchestrator.
//!
//! Drop-in replacement for the `EVidenceModeler` Perl script.

use anyhow::Result;
use clap::Parser;
use log::{info, warn};
use rayon::prelude::*;
use std::fs;
use std::io::Write;
use std::path::Path;

use evm_core::gff3_convert::evm_to_gff3::convert_all_to_gff3;
use evm_core::gff3_convert::gff3_to_bed::gff3_to_bed;
use evm_core::gff3_convert::gff3_to_proteins::{extract_sequences, SeqType};
use evm_core::io::partitions::read_partitions_file;
use evm_core::io::weights::read_weights_file;
use evm_core::partition::partition::{run_partition, InputFile};
use evm_core::pipeline::{run_single_partition, SinglePartitionParams};
use evm_core::recombine::recombine::recombine_outputs;

const VERSION: &str = "EVidenceModeler-v2.1.0-rust";

/// CLI flags mirror the Perl `EVidenceModeler` driver `GetOptions` exactly,
/// including its underscore/camelCase long names (`--sample_id`,
/// `--gene_predictions`, `--segmentSize`, `--CPU`, ...) so this binary is a
/// drop-in replacement. clap's automatic kebab-case renaming is disabled per
/// field via explicit `long = "..."`.
#[derive(Parser, Debug)]
#[command(name = "EVidenceModeler", version = VERSION, about = "Evidence Modeler — Rust implementation")]
struct Cli {
    /// Sample ID (used for naming outputs)
    #[arg(long = "sample_id")]
    sample_id: String,

    /// Genome sequence in FASTA format
    #[arg(long = "genome", short = 'G')]
    genome: String,

    /// Weights file for evidence types
    #[arg(long = "weights", short = 'w')]
    weights: String,

    /// Gene predictions GFF3 file
    #[arg(long = "gene_predictions", short = 'g')]
    gene_predictions: String,

    /// Segment size for genome partitioning
    #[arg(long = "segmentSize")]
    segment_size: u32,

    /// Overlap size between partitions
    #[arg(long = "overlapSize")]
    overlap_size: u32,

    /// Protein alignments GFF3 file
    #[arg(long = "protein_alignments", short = 'p')]
    protein_alignments: Option<String>,

    /// Transcript alignments GFF3 file
    #[arg(long = "transcript_alignments", short = 'e')]
    transcript_alignments: Option<String>,

    /// Repeats GFF3 file
    #[arg(long = "repeats", short = 'r')]
    repeats: Option<String>,

    /// Terminal exons file from PASA
    #[arg(long = "terminalExonsFile", short = 't')]
    terminal_exons: Option<String>,

    /// Stop codons (default: TAA,TGA,TAG)
    #[arg(long = "stop_codons", default_value = "TAA,TGA,TAG")]
    stop_codons: String,

    /// Minimum intron length (default: 20)
    #[arg(long = "min_intron_length", default_value_t = 20)]
    min_intron_length: u32,

    /// Number of CPUs for parallel execution
    #[arg(long = "CPU", default_value_t = 4)]
    cpu: usize,

    /// Run only on the forward strand
    #[arg(long = "forwardStrandOnly")]
    forward_strand_only: bool,

    /// Run only on the reverse strand
    #[arg(long = "reverseStrandOnly")]
    reverse_strand_only: bool,

    /// Report eliminated EVM predictions
    #[arg(long = "report_ELM")]
    report_elm: bool,

    /// Verbose output
    #[arg(short = 'S')]
    verbose: bool,

    /// Debug mode
    #[arg(long = "debug")]
    debug: bool,

    /// Stitch ends (passed through to per-partition EVM; accepted for parity)
    #[arg(long = "stitch_ends")]
    stitch_ends: Option<String>,

    /// Extend to terminal (accepted for parity)
    #[arg(long = "extend_to_terminal")]
    extend_to_terminal: Option<String>,

    /// Per-partition execution dir (accepted for parity; unused by orchestrator)
    #[arg(long = "exec_dir")]
    exec_dir: Option<String>,

    /// Limit analysis to coordinates >= this (accepted for parity)
    #[arg(long = "limit_range_lend")]
    limit_range_lend: Option<u32>,

    /// Limit analysis to coordinates <= this (accepted for parity)
    #[arg(long = "limit_range_rend")]
    limit_range_rend: Option<u32>,

    /// Trellis search limit (max previous exons compared)
    #[arg(long = "trellis_search_limit")]
    trellis_search_limit: Option<usize>,

    /// Search for nested genes in long introns (0 = off)
    #[arg(long = "search_long_introns", default_value_t = 0)]
    search_long_introns: u32,

    /// Re-examine intergenic regions of minimum length (0 = off)
    #[arg(long = "re_search_intergenic", default_value_t = 0)]
    re_search_intergenic: u32,

    /// Re-examine terminal intergenic regions of minimum length
    #[arg(long = "terminal_intergenic_re_search", default_value_t = 10000)]
    terminal_intergenic_re_search: u32,

    /// Intergenic score adjustment factor (default: 1.0)
    #[arg(long = "INTERGENIC_SCORE_ADJUST_FACTOR", default_value_t = 1.0)]
    intergenic_adjust: f64,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    let segment_size = cli.segment_size;
    let overlap_size = cli.overlap_size;

    let forward_only = cli.forward_strand_only;
    let reverse_only = cli.reverse_strand_only;
    if forward_only && reverse_only {
        anyhow::bail!("--forwardStrandOnly and --reverseStrandOnly are mutually exclusive");
    }

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.cpu)
        .build_global()
        .ok();

    let sample_id = &cli.sample_id;
    let checkpts_dir = format!("__{}-EVM_chckpts", sample_id);
    fs::create_dir_all(&checkpts_dir).ok();

    let genome_path = &cli.genome;
    let genome_basename = Path::new(genome_path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("genome.fasta")
        .to_string();

    let partition_dir = format!("{}.partitions", sample_id);
    let partition_listing = format!("{}.partitions.listing", sample_id);

    // ─── Step 1: Partition inputs ────────────────────────────────────────────
    let partition_ckpt = format!("{}/partition_inputs.ok", checkpts_dir);
    if !Path::new(&partition_ckpt).exists() {
        info!("Partitioning inputs...");
        let mut input_files = vec![InputFile::new("gene_predictions", &cli.gene_predictions)];
        if let Some(p) = &cli.protein_alignments {
            input_files.push(InputFile::new("protein_alignments", p));
        }
        if let Some(t) = &cli.transcript_alignments {
            input_files.push(InputFile::new("transcript_alignments", t));
        }
        if let Some(r) = &cli.repeats {
            input_files.push(InputFile::new("repeats", r));
        }

        run_partition(
            &partition_dir,
            genome_path,
            &input_files,
            &genome_basename,
            segment_size,
            overlap_size,
            &partition_listing,
        )?;
        fs::write(&partition_ckpt, "")?;
        info!("Partitioning complete.");
    } else {
        info!("Partitioning already done — skipping.");
    }

    let entries = read_partitions_file(&partition_listing)?;

    // ─── Step 2: Run EVM on each partition in parallel ───────────────────────
    let evm_ckpt = format!("{}/run_evm_cmds.ok", checkpts_dir);
    if !Path::new(&evm_ckpt).exists() {
        info!("Running EVM on {} partitions...", entries.len());
        let _weights = read_weights_file(&cli.weights)?;
        let _stop_codons_parsed =
            evm_core::algo::splice_sites::parse_stop_codons(&cli.stop_codons)?;

        // Build list of per-partition work items
        let work_items: Vec<_> = entries
            .iter()
            .map(|e| {
                let data_dir = if e.is_partitioned {
                    e.partition_dir.clone().unwrap_or(e.base_dir.clone())
                } else {
                    e.base_dir.clone()
                };
                (e.accession.clone(), data_dir)
            })
            .collect();

        let results: Vec<Result<()>> = work_items
            .par_iter()
            .map(|(acc, data_dir)| {
                run_evm_on_partition(
                    acc,
                    data_dir,
                    &genome_basename,
                    &cli.gene_predictions,
                    cli.protein_alignments.as_deref(),
                    cli.transcript_alignments.as_deref(),
                    &cli.weights,
                    &cli.stop_codons,
                    cli.min_intron_length,
                    forward_only,
                    reverse_only,
                    cli.report_elm,
                    cli.search_long_introns,
                    cli.re_search_intergenic,
                    cli.terminal_intergenic_re_search,
                    cli.intergenic_adjust,
                    cli.trellis_search_limit.unwrap_or(500),
                    cli.repeats.as_deref(),
                )
            })
            .collect();

        let mut had_error = false;
        for r in results {
            if let Err(e) = r {
                warn!("EVM partition error: {:#}", e);
                had_error = true;
            }
        }
        if had_error {
            anyhow::bail!("One or more EVM partitions failed");
        }
        fs::write(&evm_ckpt, "")?;
        info!("All EVM partitions complete.");
    }

    // ─── Step 3: Recombine partial outputs ───────────────────────────────────
    let recombine_ckpt = format!("{}/recombined_EVM_partial_outputs.ok", checkpts_dir);
    if !Path::new(&recombine_ckpt).exists() {
        info!("Recombining partial outputs...");
        recombine_outputs(&entries, "evm.out")?;
        fs::write(&recombine_ckpt, "")?;
    }

    // ─── Step 4: Convert to GFF3 ─────────────────────────────────────────────
    let gff3_ckpt = format!("{}/evm_out_to_gff3_format.ok", checkpts_dir);
    if !Path::new(&gff3_ckpt).exists() {
        info!("Converting EVM output to GFF3...");
        convert_all_to_gff3(&entries, "evm.out")?;
        fs::write(&gff3_ckpt, "")?;
    }

    // ─── Step 5: Concatenate GFF3 outputs ────────────────────────────────────
    let concat_ckpt = format!("{}/concatenate_final_output_gff3.ok", checkpts_dir);
    let final_gff3 = format!("{}.EVM.gff3", sample_id);
    if !Path::new(&concat_ckpt).exists() {
        info!("Concatenating GFF3 outputs...");
        concatenate_gff3_outputs(&entries, &partition_dir, &final_gff3)?;
        fs::write(&concat_ckpt, "")?;
    }

    // ─── Step 6: Extract protein sequences ───────────────────────────────────
    let pep_ckpt = format!("{}/make_evm_pep.ok", checkpts_dir);
    let pep_file = format!("{}.EVM.pep", sample_id);
    if !Path::new(&pep_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Extracting protein sequences...");
        let default_stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
        let seqs = extract_sequences(&final_gff3, genome_path, &SeqType::Prot, &default_stops)?;
        let mut out = fs::File::create(&pep_file)?;
        for (hdr, seq) in &seqs {
            writeln!(out, ">{}", hdr)?;
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
            }
        }
        fs::write(&pep_ckpt, "")?;
    }

    // ─── Step 7: Extract CDS sequences ───────────────────────────────────────
    let cds_ckpt = format!("{}/make_evm_cds.ok", checkpts_dir);
    let cds_file = format!("{}.EVM.cds", sample_id);
    if !Path::new(&cds_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Extracting CDS sequences...");
        let default_stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
        let seqs = extract_sequences(&final_gff3, genome_path, &SeqType::Cds, &default_stops)?;
        let mut out = fs::File::create(&cds_file)?;
        for (hdr, seq) in &seqs {
            writeln!(out, ">{}", hdr)?;
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
            }
        }
        fs::write(&cds_ckpt, "")?;
    }

    // ─── Step 8: Make BED file ────────────────────────────────────────────────
    let bed_ckpt = format!("{}/make_bed.ok", checkpts_dir);
    let bed_file = format!("{}.EVM.bed", sample_id);
    if !Path::new(&bed_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Generating BED file...");
        let bed_lines = gff3_to_bed(&final_gff3)?;
        let mut out = fs::File::create(&bed_file)?;
        for line in &bed_lines {
            writeln!(out, "{}", line)?;
        }
        fs::write(&bed_ckpt, "")?;
    }

    info!("Done. See {}.EVM.* outputs", sample_id);
    Ok(())
}

/// Run the EVM algorithm on a single partition directory.
#[allow(clippy::too_many_arguments)]
fn run_evm_on_partition(
    _accession: &str,
    data_dir: &str,
    genome_basename: &str,
    gene_pred_global: &str,
    protein_global: Option<&str>,
    transcript_global: Option<&str>,
    weights_path: &str,
    stop_codons_str: &str,
    min_intron_length: u32,
    forward_only: bool,
    reverse_only: bool,
    report_elm: bool,
    _search_long_introns: u32,
    _re_search_intergenic: u32,
    terminal_intergenic_re_search: u32,
    intergenic_adjust: f64,
    max_prev_exons_compare: usize,
    repeats_global: Option<&str>,
) -> Result<()> {
    let output_path = format!("{}/evm.out", data_dir);

    // Skip if already done
    let ckpt = format!("{}/evm.done.ok", data_dir);
    if Path::new(&ckpt).exists() {
        return Ok(());
    }

    // Resolve the per-partition copies of each input file.
    let genome_path = format!("{}/{}", data_dir, genome_basename);
    let gene_pred_path = format!(
        "{}/{}",
        data_dir,
        Path::new(gene_pred_global)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
    );
    let protein_path = protein_global.map(|p| {
        format!(
            "{}/{}",
            data_dir,
            Path::new(p)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
        )
    });
    let transcript_path = transcript_global.map(|p| {
        format!(
            "{}/{}",
            data_dir,
            Path::new(p)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
        )
    });
    let repeats_path = repeats_global.map(|p| {
        format!(
            "{}/{}",
            data_dir,
            Path::new(p)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
        )
    });

    let params = SinglePartitionParams {
        stop_codons: stop_codons_str.to_string(),
        min_intron_length,
        forward_only,
        reverse_only,
        report_elm,
        terminal_intergenic_re_search,
        intergenic_adjust,
        max_prev_exons_compare,
        repeats: repeats_path,
    };

    let output = run_single_partition(
        &genome_path,
        weights_path,
        &gene_pred_path,
        protein_path.as_deref(),
        transcript_path.as_deref(),
        &params,
    )?;

    let mut out_file = fs::File::create(&output_path)?;
    for block in &output {
        write!(out_file, "{}", block)?;
    }

    fs::write(&ckpt, "")?;
    Ok(())
}

/// Concatenate all per-contig GFF3 outputs into the final output file.
fn concatenate_gff3_outputs(
    entries: &[evm_core::io::partitions::PartitionEntry],
    _partition_dir: &str,
    output_path: &str,
) -> Result<()> {
    use std::collections::HashSet;
    let mut seen_dirs: HashSet<String> = HashSet::new();
    let mut out = fs::File::create(output_path)?;

    for entry in entries {
        if seen_dirs.insert(entry.base_dir.clone()) {
            let gff3 = format!("{}/evm.out.gff3", entry.base_dir);
            if Path::new(&gff3).exists() {
                let content = fs::read_to_string(&gff3)?;
                write!(out, "{}", content)?;
            }
        }
    }
    Ok(())
}
