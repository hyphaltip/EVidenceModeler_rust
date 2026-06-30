//! gff3_file_to_proteins — extract protein/CDS/cDNA sequences from GFF3 + genome.
//!
//! Replaces EvmUtils/gff3_file_to_proteins.pl

use anyhow::Result;
use clap::Parser;
use std::fs;
use std::io::Write;

use evm_core::gff3_convert::gff3_to_proteins::{extract_sequences, SeqType};

#[derive(Parser, Debug)]
#[command(
    name = "gff3_file_to_proteins",
    about = "Extract protein/CDS/cDNA from GFF3",
    version
)]
struct Cli {
    /// GFF3 file with gene models
    gff3: String,

    /// Genome FASTA
    fasta: String,

    /// Output type: prot, cds, cdna, or gene (default: prot)
    #[arg(default_value = "prot")]
    seqtype: String,

    /// Upstream/downstream flank (accepted for Perl parity; currently ignored)
    flank: Option<String>,

    /// Stop codon list (default: TAA,TGA,TAG)
    #[arg(long, default_value = "TAA,TGA,TAG")]
    stop_codons: String,

    /// Output file (default: stdout)
    #[arg(long, short = 'o')]
    output: Option<String>,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let cli = Cli::parse();

    if cli.flank.is_some() {
        log::warn!("flank argument is accepted for Perl parity but not yet implemented");
    }

    let seq_type = match cli.seqtype.as_str() {
        "prot" | "protein" => SeqType::Prot,
        "cds" => SeqType::Cds,
        "cdna" | "gene" => SeqType::Cdna,
        other => anyhow::bail!("Unknown sequence type: {}. Use prot, cds, or cdna.", other),
    };

    let stop_codons = parse_stop_codons(&cli.stop_codons)?;
    let sequences = extract_sequences(&cli.gff3, &cli.fasta, &seq_type, &stop_codons)?;

    let output: Box<dyn Write> = match &cli.output {
        Some(path) => Box::new(fs::File::create(path)?),
        None => Box::new(std::io::stdout()),
    };
    let mut out = std::io::BufWriter::new(output);

    for (header, seq) in &sequences {
        writeln!(out, ">{}", header)?;
        for chunk in seq.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
        }
    }

    Ok(())
}

fn parse_stop_codons(s: &str) -> Result<Vec<[u8; 3]>> {
    let mut result = Vec::new();
    for codon in s.split(',') {
        let codon = codon.trim();
        if codon.len() != 3 {
            anyhow::bail!("Invalid codon: {} (must be 3 nt)", codon);
        }
        let b = codon.as_bytes();
        result.push([b[0], b[1], b[2]]);
    }
    Ok(result)
}
