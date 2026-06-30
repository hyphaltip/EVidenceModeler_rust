//! Input partitioning — split genome + GFF3 files into per-contig and per-segment chunks.

use crate::io::fasta::{read_fasta_file, FastaRecord};
use crate::io::partitions::{write_partitions, PartitionEntry};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, Write};

/// Describe one input file to be partitioned.
#[derive(Debug, Clone)]
pub struct InputFile {
    pub file_type: String, // e.g. "gene_predictions"
    pub path: String,
    pub basename: String,
}

impl InputFile {
    pub fn new(file_type: &str, path: &str) -> Self {
        let basename = std::path::Path::new(path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(path)
            .to_string();
        InputFile {
            file_type: file_type.to_string(),
            path: path.to_string(),
            basename,
        }
    }
}

/// Generate overlapping range segments for a sequence of length `seq_len`.
///
/// Mirrors the Perl `get_range_list()` function:
/// windows of `segment_size` with `overlap_size` between them.
pub fn get_range_list(seq_len: u32, segment_size: u32, overlap_size: u32) -> Vec<(u32, u32)> {
    let mut ranges = Vec::new();
    let mut range_lend: u32 = 1;

    loop {
        if range_lend > seq_len.saturating_sub(overlap_size) {
            break;
        }
        let mut range_rend = range_lend + segment_size - 1;
        if range_rend > seq_len {
            range_rend = seq_len;
        }
        ranges.push((range_lend, range_rend));
        range_lend += segment_size - overlap_size;
    }

    if ranges.is_empty() {
        // Short sequence — single segment
        ranges.push((1, seq_len));
    }

    ranges
}

/// Partition all input GFF3 files by contig id, writing per-contig sub-files.
pub fn partition_files_based_on_contig(partition_dir: &str, files: &[InputFile]) -> Result<()> {
    for input in files {
        let mut contig_handles: HashMap<String, Box<dyn Write>> = HashMap::new();
        let f =
            fs::File::open(&input.path).with_context(|| format!("Cannot open {}", input.path))?;
        let reader = std::io::BufReader::new(f);

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let contig_id = line.split('\t').next().unwrap_or("").to_string();
            let contig_adj: String = contig_id
                .chars()
                .map(|c| {
                    if c.is_alphanumeric() || c == '_' {
                        c
                    } else {
                        '_'
                    }
                })
                .collect();
            let dir = format!("{}/{}", partition_dir, contig_adj);
            if !std::path::Path::new(&dir).exists() {
                fs::create_dir_all(&dir).with_context(|| format!("Cannot create dir {}", dir))?;
            }
            let fpath = format!("{}/{}", dir, input.basename);
            let handle = contig_handles.entry(contig_adj).or_insert_with(|| {
                let f = fs::OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&fpath)
                    .expect("Cannot open contig output file");
                Box::new(f) as Box<dyn Write>
            });
            writeln!(handle, "{}", line)?;
        }
    }
    Ok(())
}

/// Extract GFF3 lines for a specific accession and coordinate range from a file.
///
/// If `adjust_to_one` is true, coordinates are shifted so that `range_lend` → 1.
pub fn partition_gff3_range(
    input_path: &str,
    accession: &str,
    range_lend: u32,
    range_rend: u32,
    adjust_to_one: bool,
    output_path: &str,
) -> Result<()> {
    let f = fs::File::open(input_path).with_context(|| format!("Cannot open {}", input_path))?;
    let mut out =
        fs::File::create(output_path).with_context(|| format!("Cannot create {}", output_path))?;

    let reader = std::io::BufReader::new(f);
    let offset: i64 = if adjust_to_one {
        -(range_lend as i64) + 1
    } else {
        0
    };

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.splitn(9, '\t').collect();
        if cols.len() < 5 {
            continue;
        }
        if cols[0] != accession {
            continue;
        }
        let start: u32 = cols[3].parse().unwrap_or(0);
        let end: u32 = cols[4].parse().unwrap_or(0);
        // Keep only features within [range_lend, range_rend]
        if start < range_lend || end > range_rend {
            continue;
        }

        if adjust_to_one {
            let new_start = (start as i64 + offset) as u32;
            let new_end = (end as i64 + offset) as u32;
            let mut new_cols: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
            new_cols[3] = new_start.to_string();
            new_cols[4] = new_end.to_string();
            writeln!(out, "{}", new_cols.join("\t"))?;
        } else {
            writeln!(out, "{}", line)?;
        }
    }
    Ok(())
}

/// Write a genomic sub-sequence as a FASTA file.
pub fn write_genome_partition(
    record: &FastaRecord,
    range_lend: u32,
    range_rend: u32,
    output_path: &str,
) -> Result<()> {
    let seq = &record.sequence;
    let start = (range_lend - 1) as usize;
    let end = range_rend as usize;
    let subseq = if end <= seq.len() {
        &seq[start..end]
    } else {
        &seq[start..]
    };

    let mut out =
        fs::File::create(output_path).with_context(|| format!("Cannot create {}", output_path))?;
    writeln!(out, ">{}", record.accession)?;
    for chunk in subseq.as_bytes().chunks(60) {
        writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
    }
    Ok(())
}

/// Run the full partitioning pipeline.
pub fn run_partition(
    partition_dir: &str,
    genome_path: &str,
    input_files: &[InputFile],
    genome_basename: &str,
    segment_size: u32,
    overlap_size: u32,
    partition_listing_path: &str,
) -> Result<Vec<PartitionEntry>> {
    if !std::path::Path::new(partition_dir).exists() {
        fs::create_dir_all(partition_dir)
            .with_context(|| format!("Cannot create partition dir {}", partition_dir))?;
    }

    // Step 1: partition GFF3 files by contig
    let features_ckpt = format!("{}/features_partitioned.ok", partition_dir);
    if !std::path::Path::new(&features_ckpt).exists() {
        partition_files_based_on_contig(partition_dir, input_files)?;
        fs::write(&features_ckpt, "")?;
    }

    // Step 2: partition sequences
    let seq_ckpt = format!("{}/seqs_partitioned.ok", partition_dir);
    if std::path::Path::new(&seq_ckpt).exists() {
        // Already done — just read back the listing
        return read_existing_listing(partition_listing_path);
    }

    let fasta_records = read_fasta_file(genome_path)?;
    let mut entries: Vec<PartitionEntry> = Vec::new();

    for record in &fasta_records {
        let seq_len = record.sequence.len() as u32;
        let acc_adj: String = record
            .accession
            .chars()
            .map(|c| {
                if c.is_alphanumeric() || c == '_' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let acc_dir = format!("{}/{}", partition_dir, acc_adj);
        if !std::path::Path::new(&acc_dir).exists() {
            log::warn!("Skipping {} — no features partitioned.", record.accession);
            continue;
        }

        // Write full contig genome sequence
        let genome_out = format!("{}/{}", acc_dir, genome_basename);
        let mut gout = fs::File::create(&genome_out)?;
        writeln!(gout, ">{}", record.accession)?;
        for chunk in record.sequence.as_bytes().chunks(60) {
            writeln!(gout, "{}", std::str::from_utf8(chunk).unwrap())?;
        }

        // Ensure input files exist in the acc_dir
        for input in input_files {
            let fp = format!("{}/{}", acc_dir, input.basename);
            if !std::path::Path::new(&fp).exists() {
                fs::File::create(&fp)?; // touch
            }
        }

        let ranges = get_range_list(seq_len, segment_size, overlap_size);

        if ranges.len() == 1 {
            entries.push(PartitionEntry {
                accession: record.accession.clone(),
                base_dir: acc_dir,
                is_partitioned: false,
                partition_dir: None,
            });
        } else {
            for (lend, rend) in ranges {
                // Perl naming: "${accession}_${range_lend}-${range_rend}" — note
                // the hyphen between lend and rend (recombine extracts the range
                // via /(\d+)-(\d+)$/).
                let part_dir = format!("{}/{}_{}-{}", acc_dir, acc_adj, lend, rend);
                if !std::path::Path::new(&part_dir).exists() {
                    fs::create_dir_all(&part_dir)?;
                }
                let ckpt = format!("{}/chunk.ok", part_dir);
                if !std::path::Path::new(&ckpt).exists() {
                    // Write genome partition
                    write_genome_partition(
                        record,
                        lend,
                        rend,
                        &format!("{}/{}", part_dir, genome_basename),
                    )?;
                    // Partition each GFF3 input file
                    for input in input_files {
                        let src = format!("{}/{}", acc_dir, input.basename);
                        let dst = format!("{}/{}", part_dir, input.basename);
                        partition_gff3_range(&src, &record.accession, lend, rend, true, &dst)?;
                    }
                    fs::write(&ckpt, "")?;
                }
                entries.push(PartitionEntry {
                    accession: record.accession.clone(),
                    base_dir: acc_dir.clone(),
                    is_partitioned: true,
                    partition_dir: Some(part_dir),
                });
            }
        }
    }

    // Write listing file
    let mut listing_file = fs::File::create(partition_listing_path)?;
    write_partitions(&mut listing_file, &entries)?;
    fs::write(&seq_ckpt, "")?;

    Ok(entries)
}

fn read_existing_listing(path: &str) -> Result<Vec<PartitionEntry>> {
    crate::io::partitions::read_partitions_file(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_list_no_partition() {
        let ranges = get_range_list(1000, 5000, 500);
        assert_eq!(ranges, vec![(1, 1000)]);
    }

    #[test]
    fn range_list_multiple() {
        // seq_len=10000, seg=5000, overlap=500
        let ranges = get_range_list(10000, 5000, 500);
        // First: 1-5000, second: 4501-9500, third: 9001-10000 (if needed)
        assert!(ranges.len() >= 2);
        assert_eq!(ranges[0], (1, 5000));
        assert_eq!(ranges[1], (4501, 9500));
    }
}
