use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process;

#[derive(Clone, Debug)]
struct Feature {
    feat_type: String,
    contig: String,
    lend: u32,
    rend: u32,
    orient: char,
    parent_id: Option<String>,
    feature_id: String,
}

struct Validator {
    feature_id_to_data: HashMap<String, Feature>,
    parent_to_children: HashMap<String, Vec<String>>,
    child_to_parent: HashMap<String, String>,
}

impl Validator {
    fn new() -> Self {
        Validator {
            feature_id_to_data: HashMap::new(),
            parent_to_children: HashMap::new(),
            child_to_parent: HashMap::new(),
        }
    }

    fn parse_gff3_file(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }

            let contig = cols[0].to_string();
            let feat_type = cols[2].to_string();
            let lend: u32 = cols[3].parse().map_err(|_| {
                format!(
                    "Invalid start coordinate at line {}: {}",
                    line_num + 1,
                    cols[3]
                )
            })?;
            let rend: u32 = cols[4].parse().map_err(|_| {
                format!(
                    "Invalid end coordinate at line {}: {}",
                    line_num + 1,
                    cols[4]
                )
            })?;
            let orient = cols[6].chars().next().unwrap_or('.');

            // Only care about gene, mRNA, exon, and CDS
            if !matches!(feat_type.as_str(), "gene" | "mRNA" | "exon" | "CDS") {
                continue;
            }

            let feature_info = cols[8];

            // Parse ID
            let feature_id = extract_attribute(feature_info, "ID").ok_or_else(|| {
                format!(
                    "Cannot parse ID from entry at line {}\n{}",
                    line_num + 1,
                    line
                )
            })?;

            // Make CDS IDs unique by appending contig and coordinates
            let mut unique_id = feature_id.clone();
            if feat_type == "CDS" {
                unique_id = format!("{}.{}:{}-{}", feature_id, contig, lend, rend);
            }

            // Parse Parent if present
            let parent_id = extract_attribute(feature_info, "Parent");

            let feature = Feature {
                feat_type,
                contig,
                lend,
                rend,
                orient,
                parent_id: parent_id.clone(),
                feature_id: unique_id.clone(),
            };

            // Check for duplicate features
            if let Some(existing) = self.feature_id_to_data.get(&unique_id) {
                if !features_identical(existing, &feature) {
                    eprintln!(
                        "Error, feature: {} is described multiple times with different data values",
                        unique_id
                    );
                    eprintln!("  Existing: {:?}", existing);
                    eprintln!("  New: {:?}", feature);
                }
                continue;
            }

            // Store feature
            self.feature_id_to_data.insert(unique_id.clone(), feature);

            // Build parent-child relationships
            if let Some(pid) = parent_id {
                self.parent_to_children
                    .entry(pid.clone())
                    .or_insert_with(Vec::new)
                    .push(unique_id.clone());
                self.child_to_parent.insert(unique_id, pid);
            }
        }

        Ok(())
    }

    fn check_child_parent_consistency(&self) {
        // Check parent-child relationships according to hierarchy:
        // - exon -> mRNA
        // - CDS -> mRNA
        // - mRNA -> gene

        for (parent_id, children_ids) in &self.parent_to_children {
            let parent = match self.feature_id_to_data.get(parent_id) {
                Some(p) => p,
                None => {
                    eprintln!(
                        "Fatal Error, cannot locate data entry for ID: [{}]",
                        parent_id
                    );
                    continue;
                }
            };

            for child_id in children_ids {
                let child = match self.feature_id_to_data.get(child_id) {
                    Some(c) => c,
                    None => {
                        eprintln!("Fatal error, cannot locate data entry for ID: {}", child_id);
                        continue;
                    }
                };

                // Check contig
                if parent.contig != child.contig {
                    println!(
                        "Error, parent {} ({}, {}) and child {} ({}, {}) are located on different contigs",
                        parent_id, parent.feat_type, parent.contig, child_id, child.feat_type, child.contig
                    );
                }

                // Check orientation
                if parent.orient != child.orient {
                    println!(
                        "Error, parent {} ({}, {}) and child {} ({}, {}) have conflicting orientations.",
                        parent_id, parent.feat_type, parent.orient, child_id, child.feat_type, child.orient
                    );
                }

                // Check coordinate encapsulation
                if !(child.lend >= parent.lend && child.rend <= parent.rend) {
                    println!(
                        "Error, parent {} ({}, {}-{}) does not encapsulate coords of child {} ({}, {}-{})",
                        parent_id, parent.feat_type, parent.lend, parent.rend,
                        child_id, child.feat_type, child.lend, child.rend
                    );
                }

                // Check hierarchy compatibility
                let valid_hierarchy = match child.feat_type.as_str() {
                    "exon" => parent.feat_type == "mRNA",
                    "CDS" => parent.feat_type == "mRNA",
                    "mRNA" => parent.feat_type == "gene",
                    _ => false,
                };

                if !valid_hierarchy {
                    println!(
                        "Error, parent {} ({}) cannot have a child {} of type {}",
                        parent_id, parent.feat_type, child_id, child.feat_type
                    );
                }
            }
        }
    }

    fn check_features_have_parents(&self) {
        for (feature_id, feature) in &self.feature_id_to_data {
            if feature.feat_type != "gene" {
                match self.child_to_parent.get(feature_id) {
                    Some(parent_id) => {
                        // Check if parent has valid type
                        if let Some(parent) = self.feature_id_to_data.get(parent_id) {
                            let valid_parent = match feature.feat_type.as_str() {
                                "exon" => parent.feat_type == "mRNA",
                                "CDS" => parent.feat_type == "mRNA",
                                "mRNA" => parent.feat_type == "gene",
                                _ => false,
                            };

                            if !valid_parent {
                                println!(
                                    "ERROR, feature {} ({}) has a parent {} ({}) and is not allowed!",
                                    feature_id, feature.feat_type, parent_id, parent.feat_type
                                );
                            }
                        }
                    }
                    None => {
                        println!(
                            "ERROR, feature {} ({}) lacks a parent feature",
                            feature_id, feature.feat_type
                        );
                    }
                }
            }
        }
    }

    fn ensure_cds_and_exon_encapsulation(&self) {
        // Every mRNA must have at least one CDS
        // Every CDS must be encapsulated by an exon

        for feature in self.feature_id_to_data.values() {
            if feature.feat_type != "mRNA" {
                continue;
            }

            let mrna_id = &feature.feature_id;
            let children_ids = match self.parent_to_children.get(mrna_id) {
                Some(ids) => ids,
                None => {
                    println!("ERROR, mRNA {} has no child features", mrna_id);
                    continue;
                }
            };

            let mut exons = Vec::new();
            let mut cdss = Vec::new();

            for child_id in children_ids {
                if let Some(child) = self.feature_id_to_data.get(child_id) {
                    match child.feat_type.as_str() {
                        "exon" => exons.push(child.clone()),
                        "CDS" => cdss.push(child.clone()),
                        _ => {
                            eprintln!(
                                "Error, found unexpected child type {} for mRNA {}",
                                child.feat_type, mrna_id
                            );
                        }
                    }
                }
            }

            // Check that mRNA has at least one CDS
            if cdss.is_empty() {
                println!("ERROR, mRNA {} lacks a CDS record", mrna_id);
            }

            // Correlate CDSs with exons
            self.correlate_cdss_with_exons(&exons, &cdss);
        }
    }

    fn correlate_cdss_with_exons(&self, exons: &[Feature], cdss: &[Feature]) {
        let mut exon_used: HashMap<String, String> = HashMap::new();

        for cds in cdss {
            let cds_id = &cds.feature_id;
            let (cds_lend, cds_rend) = (cds.lend, cds.rend);

            let mut found_exon = false;

            for exon in exons {
                let exon_id = &exon.feature_id;
                let (exon_lend, exon_rend) = (exon.lend, exon.rend);

                if cds_lend >= exon_lend && cds_rend <= exon_rend {
                    // Encapsulated by this exon
                    if let Some(other_id) = exon_used.get(exon_id) {
                        println!(
                            "ERROR, CDS {} ({}-{}) maps to exon {} ({}-{}), but this exon already encodes a different CDS record {}",
                            cds_id, cds_lend, cds_rend, exon_id, exon_lend, exon_rend, other_id
                        );
                    }

                    found_exon = true;
                    exon_used.insert(exon_id.clone(), cds_id.clone());
                    break;
                }
            }

            if !found_exon {
                println!(
                    "ERROR, CDS {} does not fully map within an exon record.",
                    cds_id
                );
            }
        }
    }

    fn validate(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        self.parse_gff3_file(path)?;
        self.check_child_parent_consistency();
        self.check_features_have_parents();
        self.ensure_cds_and_exon_encapsulation();
        Ok(())
    }
}

fn extract_attribute(info: &str, key: &str) -> Option<String> {
    for part in info.split(';') {
        let trimmed = part.trim();
        if let Some(value_start) = trimmed.find('=') {
            let attr_key = &trimmed[..value_start].trim();
            if attr_key == &key {
                let value = trimmed[value_start + 1..].trim();
                return Some(value.to_string());
            }
        }
    }
    None
}

fn features_identical(a: &Feature, b: &Feature) -> bool {
    a.feat_type == b.feat_type
        && a.contig == b.contig
        && a.lend == b.lend
        && a.rend == b.rend
        && a.orient == b.orient
        && a.parent_id == b.parent_id
        && a.feature_id == b.feature_id
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("usage: {} gene_annotations.gff3", args[0]);
        process::exit(1);
    }

    let gff3_file = &args[1];

    if !Path::new(gff3_file).exists() {
        eprintln!("Error, cannot open file {}", gff3_file);
        process::exit(1);
    }

    let mut validator = Validator::new();

    if let Err(e) = validator.validate(gff3_file) {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}
