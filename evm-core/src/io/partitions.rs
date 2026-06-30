//! Partition listing file parser/writer.
//!
//! Line format (tab-separated):
//!   accession  base_dir  Y|N  [partition_dir]
//!
//! If the third field is 'N' there is no fourth field (sequence not further partitioned).

use anyhow::Result;
use std::io::{BufRead, Write};

#[derive(Debug, Clone)]
pub struct PartitionEntry {
    pub accession: String,
    pub base_dir: String,
    pub is_partitioned: bool,
    /// Present only when `is_partitioned` is true.
    pub partition_dir: Option<String>,
}

/// Parse a partitions listing file.
pub fn read_partitions_file(path: &str) -> Result<Vec<PartitionEntry>> {
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("Cannot open partitions file {}: {}", path, e))?;
    parse_partitions(std::io::BufReader::new(f))
}

pub fn parse_partitions<R: BufRead>(reader: R) -> Result<Vec<PartitionEntry>> {
    let mut entries = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let is_partitioned = parts[2] == "Y";
        let partition_dir = if is_partitioned && parts.len() >= 4 {
            Some(parts[3].to_string())
        } else {
            None
        };
        entries.push(PartitionEntry {
            accession: parts[0].to_string(),
            base_dir: parts[1].to_string(),
            is_partitioned,
            partition_dir,
        });
    }
    Ok(entries)
}

/// Write a partitions listing.
pub fn write_partitions<W: Write>(writer: &mut W, entries: &[PartitionEntry]) -> Result<()> {
    for e in entries {
        if e.is_partitioned {
            let pd = e.partition_dir.as_deref().unwrap_or("");
            writeln!(writer, "{}\t{}\tY\t{}", e.accession, e.base_dir, pd)?;
        } else {
            writeln!(writer, "{}\t{}\tN", e.accession, e.base_dir)?;
        }
    }
    Ok(())
}
