use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process;

#[derive(Clone, Debug)]
struct GeneModel {
    contig: String,
    gene_id: String,
    trans_id: String,
    cds_segments: Vec<(u32, u32)>, // (5' end, 3' end)
    strand: char,
    has_stop_codon: bool,
}

struct AugustusConverter {
    format: InputFormat,
    models: HashMap<String, GeneModel>, // key: contig-gene_id-trans_id
}

#[derive(Clone, Copy, Debug)]
enum InputFormat {
    GFF3,
    GTF,
}

impl AugustusConverter {
    fn new(format: InputFormat) -> Self {
        AugustusConverter {
            format,
            models: HashMap::new(),
        }
    }

    fn parse_input(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        match self.format {
            InputFormat::GFF3 => self.parse_gff3(path),
            InputFormat::GTF => self.parse_gtf(path),
        }
    }

    fn parse_gff3(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;

            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }

            let feat_type = cols[2];
            if feat_type != "CDS" {
                continue;
            }

            let contig = cols[0].to_string();
            let strand = cols[6].chars().next().unwrap_or('.');
            let lend: u32 = cols[3].parse()?;
            let rend: u32 = cols[4].parse()?;
            let info = cols[8];

            // Parse Parent attribute (transcript ID)
            let trans_id = extract_gff3_attribute(info, "Parent")
                .ok_or_else(|| format!("Cannot parse Parent from: {}", info))?;

            let (end5, end3) = if strand == '+' {
                (lend, rend)
            } else {
                (rend, lend)
            };

            let key = format!("{}-{}", contig, trans_id);

            self.models
                .entry(key)
                .or_insert_with(|| GeneModel {
                    contig: contig.clone(),
                    gene_id: format!("gene.{}", trans_id),
                    trans_id: trans_id.clone(),
                    cds_segments: Vec::new(),
                    strand,
                    has_stop_codon: false,
                })
                .cds_segments
                .push((end5, end3));
        }

        Ok(())
    }

    fn parse_gtf(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;

            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }

            let feat_type = cols[2];
            if feat_type != "CDS" && feat_type != "stop_codon" {
                continue;
            }

            let contig = cols[0].to_string();
            let strand = cols[6].chars().next().unwrap_or('.');
            let lend: u32 = cols[3].parse()?;
            let rend: u32 = cols[4].parse()?;
            let info = cols[8];

            // Parse GTF attributes: transcript_id "..."; gene_id "..."
            let trans_id = extract_gtf_attribute(info, "transcript_id")
                .ok_or_else(|| format!("Cannot parse transcript_id from: {}", info))?;
            let gene_id = extract_gtf_attribute(info, "gene_id")
                .ok_or_else(|| format!("Cannot parse gene_id from: {}", info))?;

            let trans_id_full = format!("{}-{}", contig, trans_id);
            let gene_id_full = format!("{}-{}", contig, gene_id);

            let (end5, end3) = if strand == '+' {
                (lend, rend)
            } else {
                (rend, lend)
            };

            let key = format!("{}-{}-{}", contig, gene_id_full, trans_id_full);

            let model = self.models.entry(key).or_insert_with(|| GeneModel {
                contig: contig.clone(),
                gene_id: gene_id_full,
                trans_id: trans_id_full,
                cds_segments: Vec::new(),
                strand,
                has_stop_codon: false,
            });

            if feat_type == "stop_codon" {
                model.has_stop_codon = true;
            } else {
                model.cds_segments.push((end5, end3));
            }
        }

        Ok(())
    }

    fn generate_gff3_output(&self) {
        // Sort by contig for consistent output
        let mut contigs: Vec<String> = self
            .models
            .values()
            .map(|m| m.contig.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        contigs.sort();

        for contig in contigs {
            let models_for_contig: Vec<_> = self
                .models
                .values()
                .filter(|m| m.contig == contig)
                .collect();

            for model in models_for_contig {
                self.output_gene_model(model);
            }
        }
    }

    fn output_gene_model(&self, model: &GeneModel) {
        if model.cds_segments.is_empty() {
            return;
        }

        // Sort CDS segments by 5' position
        let mut sorted_cds = model.cds_segments.clone();
        sorted_cds.sort_by_key(|seg| seg.0);

        // Build exons that encompass CDS segments
        let exons = build_exons_from_cds(&sorted_cds, model.strand);

        let gene_start = exons.iter().map(|e| e.0).min().unwrap_or(0);
        let gene_end = exons.iter().map(|e| e.1).max().unwrap_or(0);

        let strand_str = format!("{}", model.strand);

        // Gene feature
        println!(
            "{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tID={}",
            model.contig, "Augustus", gene_start, gene_end, strand_str, model.gene_id
        );

        // mRNA feature
        println!(
            "{}\t{}\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={}",
            model.contig,
            "Augustus",
            gene_start,
            gene_end,
            strand_str,
            model.trans_id,
            model.gene_id
        );

        // Exon features and CDS features
        for (idx, (exon_start, exon_end)) in exons.iter().enumerate() {
            let exon_num = idx + 1;

            // Exon
            println!(
                "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tID={}.exon{};Parent={}",
                model.contig,
                "Augustus",
                exon_start,
                exon_end,
                strand_str,
                model.trans_id,
                exon_num,
                model.trans_id
            );

            // CDS features that overlap with this exon
            let phase = if model.strand == '+' { 0 } else { 0 };

            for (cds_idx, (cds_end5, cds_end3)) in sorted_cds.iter().enumerate() {
                // Convert end5/end3 to genomic min/max
                let cds_min = (*cds_end5).min(*cds_end3);
                let cds_max = (*cds_end5).max(*cds_end3);

                // Check if CDS overlaps with exon
                if cds_min < *exon_end && cds_max > *exon_start {
                    let cds_num = cds_idx + 1;

                    println!(
                        "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tID={}.cds{};Parent={}",
                        model.contig,
                        "Augustus",
                        cds_min,
                        cds_max,
                        strand_str,
                        phase,
                        model.trans_id,
                        cds_num,
                        model.trans_id
                    );
                }
            }
        }
    }
}

fn build_exons_from_cds(cds_segments: &[(u32, u32)], _strand: char) -> Vec<(u32, u32)> {
    // Exons encompass CDSs
    // Convert end5/end3 format to min/max for exon bounds
    if cds_segments.is_empty() {
        return Vec::new();
    }

    // Convert to actual genomic coordinates (min, max)
    let mut genomic_coords: Vec<(u32, u32)> = cds_segments
        .iter()
        .map(|(a, b)| ((*a).min(*b), (*a).max(*b)))
        .collect();

    // Sort by start position
    genomic_coords.sort_by_key(|seg| seg.0);

    let mut exons = Vec::new();
    let mut current_start = genomic_coords[0].0;
    let mut current_end = genomic_coords[0].1;

    for (seg_start, seg_end) in &genomic_coords[1..] {
        // Check if contiguous (allowing small gaps)
        if *seg_start > current_end + 10000 {
            // Gap > 10kb, new exon
            exons.push((current_start, current_end));
            current_start = *seg_start;
            current_end = *seg_end;
        } else {
            // Extend current exon
            current_end = (*seg_end).max(current_end);
        }
    }

    // Add final exon
    exons.push((current_start, current_end));

    exons
}

fn extract_gff3_attribute(info: &str, key: &str) -> Option<String> {
    for part in info.split(';') {
        let trimmed = part.trim();
        if let Some(eq_pos) = trimmed.find('=') {
            let attr_key = trimmed[..eq_pos].trim();
            if attr_key == key {
                let value = trimmed[eq_pos + 1..].trim();
                return Some(value.to_string());
            }
        }
    }
    None
}

fn extract_gtf_attribute(info: &str, key: &str) -> Option<String> {
    // GTF format: key "value"; ...
    let pattern = format!("{} \"", key);
    if let Some(start_pos) = info.find(&pattern) {
        let start = start_pos + pattern.len();
        if let Some(end_pos) = info[start..].find('"') {
            return Some(info[start..start + end_pos].to_string());
        }
    }
    None
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("usage: {} [--format gff3|gtf] <augustus_output>", args[0]);
        eprintln!();
        eprintln!("Convert Augustus GFF3 or GTF output to EVM-compatible GFF3");
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --format gff3   Input is Augustus GFF3 (default)");
        eprintln!("  --format gtf    Input is Augustus GTF");
        process::exit(1);
    }

    let mut format = InputFormat::GFF3;
    let mut input_file = String::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--format" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: --format requires an argument (gff3 or gtf)");
                    process::exit(1);
                }
                format = match args[i].as_str() {
                    "gff3" => InputFormat::GFF3,
                    "gtf" => InputFormat::GTF,
                    _ => {
                        eprintln!("Error: unknown format '{}' (use 'gff3' or 'gtf')", args[i]);
                        process::exit(1);
                    }
                };
            }
            arg => {
                if input_file.is_empty() {
                    input_file = arg.to_string();
                } else {
                    eprintln!("Error: unexpected argument '{}'", arg);
                    process::exit(1);
                }
            }
        }
        i += 1;
    }

    if input_file.is_empty() {
        eprintln!("Error: no input file specified");
        process::exit(1);
    }

    if !Path::new(&input_file).exists() {
        eprintln!("Error: cannot open file '{}'", input_file);
        process::exit(1);
    }

    let mut converter = AugustusConverter::new(format);

    if let Err(e) = converter.parse_input(&input_file) {
        eprintln!("Error parsing input: {}", e);
        process::exit(1);
    }

    converter.generate_gff3_output();
}
