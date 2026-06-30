//! recombine_evm_outputs — merge per-partition EVM outputs.
//!
//! Replaces EvmUtils/recombine_EVM_partial_outputs.pl

use anyhow::Result;
use clap::Parser;

use evm_core::io::partitions::read_partitions_file;
use evm_core::recombine::recombine::recombine_outputs;

#[derive(Parser, Debug)]
#[command(
    name = "recombine_evm_outputs",
    about = "Recombine partial EVM outputs"
)]
struct Cli {
    /// Partitions listing file
    #[arg(long = "partitions")]
    partitions_list: String,

    /// EVM output filename within each partition dir (default: evm.out)
    #[arg(long = "output_file_name", short = 'O', default_value = "evm.out")]
    evm_output_file: String,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    let entries = read_partitions_file(&cli.partitions_list)?;
    recombine_outputs(&entries, &cli.evm_output_file)?;

    eprintln!("Recombination complete.");
    Ok(())
}
