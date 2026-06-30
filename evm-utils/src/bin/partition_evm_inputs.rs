//! partition_evm_inputs — partition EVM inputs for parallel processing.
//!
//! Replaces EvmUtils/partition_EVM_inputs.pl

use anyhow::Result;
use clap::Parser;

use evm_core::partition::partition::{run_partition, InputFile};

#[derive(Parser, Debug)]
#[command(
    name = "partition_evm_inputs",
    about = "Partition EVM inputs for parallel execution"
)]
struct Cli {
    /// Output directory for partitions
    #[arg(long = "partition_dir")]
    partition_dir: String,

    /// Genome FASTA
    #[arg(long = "genome", short = 'g')]
    genome: String,

    /// Gene predictions GFF3
    #[arg(long = "gene_predictions")]
    gene_predictions: Option<String>,

    /// Protein alignments GFF3
    #[arg(long = "protein_alignments")]
    protein_alignments: Option<String>,

    /// Transcript alignments GFF3
    #[arg(long = "transcript_alignments")]
    transcript_alignments: Option<String>,

    /// PASA terminal exons file (partitioned for downstream compatibility)
    #[arg(long = "pasaTerminalExons")]
    pasa_terminal_exons: Option<String>,

    /// Repeats GFF3
    #[arg(long = "repeats")]
    repeats: Option<String>,

    /// Segment size (nt)
    #[arg(long = "segmentSize", default_value_t = 100000)]
    segment_size: u32,

    /// Overlap between segments (nt)
    #[arg(long = "overlapSize", default_value_t = 10000)]
    overlap_size: u32,

    /// Output partitions listing file
    #[arg(long = "partition_listing")]
    partition_listing: Option<String>,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    let genome_basename = std::path::Path::new(&cli.genome)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("genome.fasta")
        .to_string();

    let mut input_files: Vec<InputFile> = Vec::new();
    if let Some(p) = &cli.gene_predictions {
        input_files.push(InputFile::new("gene_predictions", p));
    }
    if let Some(p) = &cli.protein_alignments {
        input_files.push(InputFile::new("protein_alignments", p));
    }
    if let Some(p) = &cli.transcript_alignments {
        input_files.push(InputFile::new("transcript_alignments", p));
    }
    if let Some(p) = &cli.pasa_terminal_exons {
        input_files.push(InputFile::new("pasaTerminalExons", p));
    }
    if let Some(p) = &cli.repeats {
        input_files.push(InputFile::new("repeats", p));
    }

    let listing_path = cli
        .partition_listing
        .unwrap_or_else(|| format!("{}.listing", cli.partition_dir));

    run_partition(
        &cli.partition_dir,
        &cli.genome,
        &input_files,
        &genome_basename,
        cli.segment_size,
        cli.overlap_size,
        &listing_path,
    )?;

    eprintln!("Partitioning complete. See {}", listing_path);
    Ok(())
}
