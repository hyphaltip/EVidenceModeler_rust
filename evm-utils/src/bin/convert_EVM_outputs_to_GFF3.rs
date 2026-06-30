//! convert_EVM_outputs_to_GFF3 — convert per-contig EVM outputs to GFF3.
//!
//! Drop-in replacement for EvmUtils/convert_EVM_outputs_to_GFF3.pl, which just
//! drives EvmUtils/EVM_to_GFF3.pl for each entry in the partitions listing.

use anyhow::Result;
use clap::Parser;

use evm_core::gff3_convert::evm_to_gff3::convert_all_to_gff3;
use evm_core::io::partitions::read_partitions_file;

#[derive(Parser, Debug)]
#[command(
    name = "convert_EVM_outputs_to_GFF3",
    about = "Convert EVM outputs to GFF3 format",
    version
)]
struct Cli {
    /// Partitions listing file
    #[arg(long = "partitions")]
    partitions: String,

    /// EVM output filename within each base/partition dir (default: evm.out)
    #[arg(long = "output_file_name", short = 'O', default_value = "evm.out")]
    output_file_name: String,

    /// Genome FASTA (accepted for Perl CLI parity; not used by the Rust converter)
    #[arg(long = "genome")]
    _genome: Option<String>,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    let entries = read_partitions_file(&cli.partitions)?;
    convert_all_to_gff3(&entries, &cli.output_file_name)?;

    eprintln!("GFF3 conversion complete.");
    Ok(())
}
