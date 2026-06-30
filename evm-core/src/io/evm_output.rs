//! Parser and writer for the intermediate EVM text output format.
//!
//! The format produced by evidence_modeler.pl for each partition:
//!
//! ```text
//! !! Predictions spanning range 100-500 [R1]
//! # EVM prediction
//! 100     200     initial  1  ev_type  accession
//! 300     450     terminal 3  ev_type  accession
//!
//! # EVM prediction
//! ...
//! ```
//!
//! Each prediction block is separated by blank lines.
//! Lines starting with `!!` are comment banners.
//! Lines starting with `#` are model headers.
//! Data lines have 6 tab-separated columns: end5, end3, exon_type, phase, ev_type, accession(s).

use anyhow::Result;
use std::io::{BufRead, Write};

/// A single exon row within an EVM prediction block.
#[derive(Debug, Clone)]
pub struct EvmExonRow {
    pub end5: u32,
    pub end3: u32,
    pub exon_type: String,
    pub phase: u8,
    pub ev_info: String, // remaining columns as-is
}

/// A complete EVM prediction block as read from the output file.
#[derive(Debug, Clone)]
pub struct EvmBlock {
    /// Header line text (without the leading `#`).
    pub header: String,
    pub exon_rows: Vec<EvmExonRow>,
    /// True if this block is an ELIMINATED prediction.
    pub is_eliminated: bool,
}

/// Parse an EVM output file produced by evidence_modeler.
pub fn parse_evm_output<R: BufRead>(reader: R) -> Result<Vec<EvmBlock>> {
    let mut blocks: Vec<EvmBlock> = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_rows: Vec<EvmExonRow> = Vec::new();
    let mut current_eliminated = false;

    let flush = |blocks: &mut Vec<EvmBlock>,
                 header: &mut Option<String>,
                 rows: &mut Vec<EvmExonRow>,
                 elim: &mut bool| {
        if let Some(h) = header.take() {
            if !rows.is_empty() {
                blocks.push(EvmBlock {
                    header: h,
                    exon_rows: std::mem::take(rows),
                    is_eliminated: *elim,
                });
            }
            *elim = false;
        }
    };

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("!!") {
            continue;
        } // banner comment

        if line.starts_with('#') {
            // New model header — flush previous block first
            flush(
                &mut blocks,
                &mut current_header,
                &mut current_rows,
                &mut current_eliminated,
            );
            let header_text = line.trim_start_matches('#').trim().to_string();
            current_eliminated = header_text.contains("ELIMINATED");
            current_header = Some(header_text);
            continue;
        }

        if line.trim().is_empty() {
            // Blank line terminates a block
            flush(
                &mut blocks,
                &mut current_header,
                &mut current_rows,
                &mut current_eliminated,
            );
            continue;
        }

        // Data line
        let cols: Vec<&str> = line.splitn(6, '\t').collect();
        if cols.len() >= 4 {
            if let (Ok(end5), Ok(end3)) = (cols[0].parse::<u32>(), cols[1].parse::<u32>()) {
                let exon_type = cols[2].to_string();
                if exon_type == "INTRON" {
                    continue;
                } // skip intron rows
                let phase: u8 = cols[3].parse().unwrap_or(0);
                let ev_info = if cols.len() >= 5 {
                    cols[4].to_string()
                } else {
                    String::new()
                };
                current_rows.push(EvmExonRow {
                    end5,
                    end3,
                    exon_type,
                    phase,
                    ev_info,
                });
            }
        }
    }
    flush(
        &mut blocks,
        &mut current_header,
        &mut current_rows,
        &mut current_eliminated,
    );
    Ok(blocks)
}

/// Write a single EVM prediction block.
pub fn write_evm_block<W: Write>(writer: &mut W, block: &EvmBlock) -> Result<()> {
    let prefix = if block.is_eliminated {
        "#ELIMINATED"
    } else {
        "#"
    };
    writeln!(writer, "{} {}", prefix, block.header)?;
    for row in &block.exon_rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            row.end5, row.end3, row.exon_type, row.phase, row.ev_info
        )?;
    }
    writeln!(writer)?;
    Ok(())
}
