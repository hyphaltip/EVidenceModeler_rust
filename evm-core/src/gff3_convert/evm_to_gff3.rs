//! Convert EVM output format to GFF3.
//!
//! Faithful port of `EVM_to_GFF3.pl` + `Gene_obj::to_GFF3_format`
//! (`convert_EVM_outputs_to_GFF3.pl` just drives `EVM_to_GFF3.pl` per
//! partition; its `join_intronic_preds`/`combine_predictions` subs are dead
//! code — never invoked — so they are intentionally not ported).

use crate::io::partitions::PartitionEntry;
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, Write};

/// EVM exon start-frame [1..=6] → GFF3 CDS phase [0,1,2].
///
/// Composes the two Perl conversions: `EVM_to_GFF3.pl`'s `%phase_conversion`
/// (1→0,2→1,3→2,4→0,5→1,6→2) followed by `Gene_obj::to_GFF3_format`'s GFF3
/// phase swap (1↔2, 0 unchanged). Net: {1,4}→0, {2,5}→2, {3,6}→1.
fn evm_startframe_to_gff3_phase(frame: u8) -> u8 {
    match frame {
        1 | 4 => 0,
        2 | 5 => 2,
        3 | 6 => 1,
        _ => 0,
    }
}

/// URI-escape matching Perl `URI::Escape::uri_escape` defaults: everything
/// except the RFC 3986 unreserved set `A-Za-z0-9-_.~` is percent-encoded.
fn uri_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        let unreserved = b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~');
        if unreserved {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{:02X}", b));
        }
    }
    out
}

/// One parsed EVM model: its exon rows `(end5, end3, start_frame)` in file
/// order, and the evidence-type tag for the GFF3 `source` column.
struct Model {
    coords: Vec<(u32, u32, u8)>,
    ev_type: &'static str,
}

/// Convert a single EVM output file to GFF3 format, writing to `output_path`.
pub fn evm_output_to_gff3(evm_file: &str, contig_id: &str, output_path: &str) -> Result<()> {
    let f =
        fs::File::open(evm_file).with_context(|| format!("Cannot open EVM output {}", evm_file))?;
    let reader = std::io::BufReader::new(f);
    let mut out = fs::File::create(output_path)
        .with_context(|| format!("Cannot create GFF3 output {}", output_path))?;

    // Mirror EVM_to_GFF3.pl's scan: model_id starts at 1 and increments on each
    // blank line; '#' header lines set the current model's ev_type; numeric
    // exon rows (6 tab cols, INTRON excluded) accumulate onto the current model.
    let mut model_id = 1u32;
    let mut models: BTreeMap<u32, Model> = BTreeMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('!') {
            continue;
        } // '^!' comment line

        if line.starts_with('#') {
            // Perl: a '#' line not matching /EVM/ is skipped; otherwise sets ev_type.
            if line.contains("EVM") {
                let ev_type = if line.contains("ELIMINATED") {
                    "EVM_elm"
                } else {
                    "EVM"
                };
                models
                    .entry(model_id)
                    .or_insert_with(|| Model {
                        coords: Vec::new(),
                        ev_type,
                    })
                    .ev_type = ev_type;
            }
            continue;
        }

        // Non-word line (blank) → end of current model.
        if !line.chars().any(|c| c.is_alphanumeric()) {
            model_id += 1;
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        // Perl requires exactly 6 columns with numeric coords.
        if cols.len() == 6 {
            if let (Ok(e5), Ok(e3), Ok(frame)) = (
                cols[0].parse::<u32>(),
                cols[1].parse::<u32>(),
                cols[3].parse::<u8>(),
            ) {
                if cols[2] != "INTRON" {
                    models
                        .entry(model_id)
                        .or_insert_with(|| Model {
                            coords: Vec::new(),
                            ev_type: "EVM",
                        })
                        .coords
                        .push((e5, e3, frame));
                }
            }
        }
    }

    for (mid, model) in &models {
        if model.coords.is_empty() {
            continue;
        }
        emit_model(&mut out, contig_id, *mid, model)?;
    }

    Ok(())
}

/// Emit one gene model (gene + mRNA + interleaved exon/CDS rows) in GFF3.
fn emit_model(out: &mut dyn Write, contig: &str, model_id: u32, model: &Model) -> Result<()> {
    let ev_type = model.ev_type;
    let coords = &model.coords;

    let gene_lend = coords.iter().map(|&(e5, e3, _)| e5.min(e3)).min().unwrap();
    let gene_rend = coords.iter().map(|&(e5, e3, _)| e5.max(e3)).max().unwrap();
    // Strand from the first exon row (Perl Gene_obj orientation).
    let strand = if coords[0].0 <= coords[0].1 { '+' } else { '-' };

    let tu_id = format!("evm.TU.{}.{}", contig, model_id);
    let model_feat = format!("evm.model.{}.{}", contig, model_id);
    let com_name = uri_escape(&format!("EVM prediction {}.{}", contig, model_id));

    // gene + mRNA
    writeln!(
        out,
        "{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={}",
        contig, ev_type, gene_lend, gene_rend, strand, tu_id, com_name
    )?;
    writeln!(
        out,
        "{}\t{}\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};Name={}",
        contig, ev_type, gene_lend, gene_rend, strand, model_feat, tu_id, com_name
    )?;

    // Exons in transcription order: ascending genomic position for '+',
    // descending for '-' (matches Gene_obj::get_exons 5'→3' ordering).
    let mut ordered: Vec<(u32, u32, u8)> = coords.clone();
    ordered.sort_by_key(|&(e5, e3, _)| e5.min(e3));
    if strand == '-' {
        ordered.reverse();
    }

    for (i, &(e5, e3, frame)) in ordered.iter().enumerate() {
        let (cs, ce) = (e5.min(e3), e5.max(e3));
        let exon_num = i + 1;
        writeln!(
            out,
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tID={}.exon{};Parent={}",
            contig, ev_type, cs, ce, strand, model_feat, exon_num, model_feat
        )?;
        let phase = evm_startframe_to_gff3_phase(frame);
        writeln!(
            out,
            "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tID=cds.{};Parent={}",
            contig, ev_type, cs, ce, strand, phase, model_feat, model_feat
        )?;
    }
    writeln!(out)?; // blank spacer between genes
    Ok(())
}

/// Convert EVM outputs for all entries to GFF3.
pub fn convert_all_to_gff3(entries: &[PartitionEntry], output_file_name: &str) -> Result<()> {
    use std::collections::HashMap;
    let mut base_dirs: HashMap<String, String> = HashMap::new();
    for entry in entries {
        base_dirs.insert(entry.accession.clone(), entry.base_dir.clone());
    }
    for (accession, base_dir) in &base_dirs {
        let evm_file = format!("{}/{}", base_dir, output_file_name);
        let gff3_file = format!("{}/{}.gff3", base_dir, output_file_name);
        evm_output_to_gff3(&evm_file, accession, &gff3_file)?;
    }
    Ok(())
}
