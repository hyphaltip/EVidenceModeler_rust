//! Convert GFF3 gene models to BED format.
//!
//! Faithful port of `gene_gff3_to_bed.pl` + `Gene_obj::to_BED_format`. The
//! driver pipes the output through `sort -k1,1 -k2,2g -k3,3g`, which this
//! function reproduces before returning.

use crate::io::gff3::read_gff3_file;
use anyhow::Result;
use std::collections::HashMap;

struct Model {
    gene_id: String,  // TU_feat_name (mRNA Parent)
    model_id: String, // Model_feat_name (mRNA ID)
    com_name: String, // un-escaped Name
    contig: String,
    strand: char,
    exons: Vec<(u32, u32)>, // exon features (block coords)
    cds: Vec<(u32, u32)>,   // CDS features (thick/coding span)
}

/// Convert a GFF3 file containing gene models to BED12, matching
/// `Gene_obj::to_BED_format` and the driver's coordinate sort.
pub fn gff3_to_bed(gff3_path: &str) -> Result<Vec<String>> {
    let records = read_gff3_file(gff3_path)?;

    let mut models: HashMap<String, Model> = HashMap::new();

    for rec in &records {
        match rec.feature.as_str() {
            "mRNA" => {
                let id = rec.attr("ID").unwrap_or("").to_string();
                let gene_id = rec.attr("Parent").unwrap_or("").to_string();
                let com_name = super::uri_unescape(rec.attr("Name").unwrap_or(""));
                let m = models.entry(id.clone()).or_insert_with(|| Model {
                    gene_id: String::new(),
                    model_id: id.clone(),
                    com_name: String::new(),
                    contig: rec.seqid.clone(),
                    strand: rec.strand,
                    exons: Vec::new(),
                    cds: Vec::new(),
                });
                m.gene_id = gene_id;
                m.com_name = com_name;
                m.contig = rec.seqid.clone();
                m.strand = rec.strand;
            }
            "exon" => {
                let parent = rec.attr("Parent").unwrap_or("").to_string();
                models
                    .entry(parent)
                    .or_insert_with(|| Model {
                        gene_id: String::new(),
                        model_id: String::new(),
                        com_name: String::new(),
                        contig: rec.seqid.clone(),
                        strand: rec.strand,
                        exons: Vec::new(),
                        cds: Vec::new(),
                    })
                    .exons
                    .push((rec.start, rec.end));
            }
            "CDS" => {
                let parent = rec.attr("Parent").unwrap_or("").to_string();
                models
                    .entry(parent)
                    .or_insert_with(|| Model {
                        gene_id: String::new(),
                        model_id: String::new(),
                        com_name: String::new(),
                        contig: rec.seqid.clone(),
                        strand: rec.strand,
                        exons: Vec::new(),
                        cds: Vec::new(),
                    })
                    .cds
                    .push((rec.start, rec.end));
            }
            _ => {}
        }
    }

    let mut bed_lines = Vec::new();

    for model in models.values() {
        if model.exons.is_empty() {
            continue;
        }

        // Exons sorted ascending (Perl sorts by end5; for non-overlapping exons
        // this is ascending genomic order on both strands). Blocks are emitted
        // in genomic-ascending order regardless of strand.
        let mut exons = model.exons.clone();
        exons.sort();
        let gene_lend = exons[0].0;
        let gene_rend = exons[exons.len() - 1].1;

        let block_sizes: Vec<String> = exons
            .iter()
            .map(|&(l, r)| (r - l + 1).to_string())
            .collect();
        let block_starts: Vec<String> = exons
            .iter()
            .map(|&(l, _)| (l - gene_lend).to_string())
            .collect();

        // Coding (thick) span from CDS features; falls back to exon span.
        let coding_lend = model.cds.iter().map(|&(l, _)| l).min().unwrap_or(gene_lend);
        let coding_rend = model.cds.iter().map(|&(_, r)| r).max().unwrap_or(gene_rend);

        // name = "ID=<model>;<gene>;<com_name>" with spaces → underscores.
        let mut name = format!("ID={};{};{}", model.model_id, model.gene_id, model.com_name);
        name = name.replace(' ', "_");

        bed_lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            model.contig,
            gene_lend - 1,
            gene_rend,
            name,
            0, // score
            model.strand,
            coding_lend - 1,
            coding_rend,
            "0", // itemRgb
            exons.len(),
            block_sizes.join(","),
            block_starts.join(","),
        ));
    }

    // Reproduce the driver's `sort -k1,1 -k2,2g -k3,3g`: chrom lexical, then
    // start (numeric), then end (numeric).
    bed_lines.sort_by(|a, b| {
        let a: Vec<&str> = a.splitn(4, '\t').collect();
        let b: Vec<&str> = b.splitn(4, '\t').collect();
        a[0].cmp(b[0])
            .then(
                a[1].parse::<i64>()
                    .unwrap_or(0)
                    .cmp(&b[1].parse::<i64>().unwrap_or(0)),
            )
            .then(
                a[2].parse::<i64>()
                    .unwrap_or(0)
                    .cmp(&b[2].parse::<i64>().unwrap_or(0)),
            )
    });

    Ok(bed_lines)
}
