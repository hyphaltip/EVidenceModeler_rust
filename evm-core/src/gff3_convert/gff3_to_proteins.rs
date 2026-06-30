//! Extract protein / CDS / cDNA sequences from GFF3 + genome FASTA.
//!
//! Replaces gff3_file_to_proteins.pl + Gene_obj.pm + Nuc_translator.pm.

use crate::io::fasta::read_fasta_hash;
use crate::io::gff3::read_gff3_file;
use crate::translate::codon_table::translate;
use crate::types::genome::reverse_complement_bytes;
use anyhow::Result;
use std::collections::HashMap;

/// Sequence output type.
pub enum SeqType {
    Prot,
    Cds,
    Cdna,
    Gene,
}

impl SeqType {
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<SeqType> {
        match s {
            "prot" => Some(SeqType::Prot),
            "CDS" => Some(SeqType::Cds),
            "cDNA" => Some(SeqType::Cdna),
            "gene" => Some(SeqType::Gene),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
struct GeneModel {
    gene_id: String,
    model_id: String,
    com_name: String,
    contig: String,
    strand: char,
    cds_exons: Vec<(u32, u32, u8)>, // (start, end, gff3_phase)
}

/// Extract sequences from a GFF3 file and genome FASTA.
pub fn extract_sequences(
    gff3_path: &str,
    genome_path: &str,
    seq_type: &SeqType,
    stop_codons: &[[u8; 3]],
) -> Result<Vec<(String, String)>> {
    let genome = read_fasta_hash(genome_path)?;
    let records = read_gff3_file(gff3_path)?;

    let mut models: HashMap<String, GeneModel> = HashMap::new();
    let mut gene_to_models: HashMap<String, Vec<String>> = HashMap::new();

    for rec in &records {
        match rec.feature.as_str() {
            "mRNA" => {
                let model_id = rec.attr("ID").unwrap_or("").to_string();
                let gene_id = rec.attr("Parent").unwrap_or("").to_string();
                let com_name = super::uri_unescape(rec.attr("Name").unwrap_or(""));
                models.insert(
                    model_id.clone(),
                    GeneModel {
                        gene_id: gene_id.clone(),
                        model_id: model_id.clone(),
                        com_name,
                        contig: rec.seqid.clone(),
                        strand: rec.strand,
                        cds_exons: Vec::new(),
                    },
                );
                gene_to_models.entry(gene_id).or_default().push(model_id);
            }
            "CDS" => {
                let parent = rec.attr("Parent").unwrap_or("").to_string();
                if let Some(m) = models.get_mut(&parent) {
                    let phase = rec.phase.unwrap_or(0);
                    m.cds_exons.push((rec.start, rec.end, phase));
                }
            }
            _ => {}
        }
    }

    let mut results = Vec::new();

    for model in models.values() {
        let genome_seq = match genome.get(&model.contig) {
            Some(s) => s,
            None => {
                log::warn!("No sequence for contig {}", model.contig);
                continue;
            }
        };

        let seq = match seq_type {
            SeqType::Prot | SeqType::Cds | SeqType::Cdna => {
                build_cds_sequence(model, genome_seq, model.strand, stop_codons)?
            }
            SeqType::Gene => {
                let (lend, rend) = model_span(model);
                let sub = &genome_seq.as_bytes()[(lend - 1) as usize..rend as usize];
                if model.strand == '-' {
                    String::from_utf8(reverse_complement_bytes(sub)).unwrap_or_default()
                } else {
                    String::from_utf8(sub.to_vec()).unwrap_or_default()
                }
            }
        };

        // Mirror gff3_file_to_proteins.pl line 128:
        //   ">$isoform_id $gene_id $locus_string $com_name $asmbl:lend-rend(orient)"
        // with an empty $locus_string (no pub_locus), which leaves a double
        // space between $gene_id and $com_name. Perl blanks com_name when it
        // equals the model id.
        let com_name = if model.com_name == model.model_id {
            ""
        } else {
            model.com_name.as_str()
        };
        let header = format!(
            "{} {}  {} {}:{}-{}({})",
            model.model_id,
            model.gene_id,
            com_name,
            model.contig,
            model
                .cds_exons
                .iter()
                .map(|&(s, _, _)| s)
                .min()
                .unwrap_or(0),
            model
                .cds_exons
                .iter()
                .map(|&(_, e, _)| e)
                .max()
                .unwrap_or(0),
            model.strand,
        );

        let final_seq = match seq_type {
            SeqType::Prot => {
                // 5'-partial models start mid-codon: trim the first CDS exon's
                // GFF3 phase (bases to drop to reach the first complete codon)
                // before translating. The CDS/cDNA outputs keep the full
                // sequence. The transcription-first exon is the lowest genomic
                // start on '+' and the highest on '-'.
                let phase = transcription_first_phase(model) as usize;
                let bytes = seq.as_bytes();
                let start = phase.min(bytes.len());
                translate(&bytes[start..], stop_codons)
            }
            _ => seq,
        };

        results.push((header, final_seq));
    }

    Ok(results)
}

fn build_cds_sequence(
    model: &GeneModel,
    genome_seq: &str,
    strand: char,
    _stop_codons: &[[u8; 3]],
) -> Result<String> {
    let mut exons = model.cds_exons.clone();
    // Sort by genomic position
    exons.sort_by_key(|&(s, _, _)| s);

    let mut cds = String::new();
    for (start, end, _) in &exons {
        let s = (start - 1) as usize;
        let e = *end as usize;
        if e <= genome_seq.len() {
            cds.push_str(&genome_seq[s..e].to_ascii_uppercase());
        }
    }

    if strand == '-' {
        let rc = reverse_complement_bytes(cds.as_bytes());
        cds = String::from_utf8(rc).unwrap_or_default();
    }

    Ok(cds)
}

/// GFF3 phase of the transcription-first CDS exon (lowest genomic start on
/// '+', highest on '-'); used to trim a partial leading codon before protein
/// translation.
fn transcription_first_phase(model: &GeneModel) -> u8 {
    if model.strand == '-' {
        model
            .cds_exons
            .iter()
            .max_by_key(|&&(s, _, _)| s)
            .map(|&(_, _, p)| p)
            .unwrap_or(0)
    } else {
        model
            .cds_exons
            .iter()
            .min_by_key(|&&(s, _, _)| s)
            .map(|&(_, _, p)| p)
            .unwrap_or(0)
    }
}

fn model_span(model: &GeneModel) -> (u32, u32) {
    let lend = model
        .cds_exons
        .iter()
        .map(|&(s, _, _)| s)
        .min()
        .unwrap_or(0);
    let rend = model
        .cds_exons
        .iter()
        .map(|&(_, e, _)| e)
        .max()
        .unwrap_or(0);
    (lend, rend)
}
