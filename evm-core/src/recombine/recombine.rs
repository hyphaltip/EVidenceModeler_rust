//! Recombine partial EVM outputs from overlapping partitions.
//!
//! Mirrors the Perl `recombine_EVM_partial_outputs.pl` script.

use crate::io::partitions::PartitionEntry;
use crate::types::prediction::{PartitionPred, PredClass};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, Write};

/// Parse an EVM output file from a single partition, offsetting coordinates
/// by `partition_lend - 1` to map them back to the full-contig space.
pub fn parse_and_add_predictions(
    output_file: &str,
    partition_lend: u32,
    predictions: &mut Vec<PartitionPred>,
) -> Result<()> {
    log::debug!("Parsing {}", output_file);
    let f = match fs::File::open(output_file) {
        Ok(f) => f,
        Err(e) => {
            log::warn!("Cannot open {}: {}", output_file, e);
            return Ok(());
        }
    };
    let reader = std::io::BufReader::new(f);

    let mut current_text = String::new();
    let mut preds: Vec<PartitionPred> = Vec::new();

    let process = |text: &str, offset: u32, out: &mut Vec<PartitionPred>| {
        if text.is_empty() {
            return;
        }
        if let Some(pred) = process_prediction_text(text, offset) {
            out.push(pred);
        }
    };

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("!!") {
            continue;
        }
        if line.starts_with('#') && !line.contains("EVM") {
            continue;
        }

        if line.starts_with(|c: char| c.is_ascii_digit() || c == '#') {
            current_text.push_str(&line);
            current_text.push('\n');
        } else {
            if !current_text.is_empty() {
                process(&current_text, partition_lend, &mut preds);
                current_text.clear();
            }
        }
    }
    if !current_text.is_empty() {
        process(&current_text, partition_lend, &mut preds);
    }

    // Join nested (intronic) predictions
    let final_preds = join_intronic_preds(preds);
    predictions.extend(final_preds);
    Ok(())
}

/// Parse a single prediction text block and adjust coordinates. Faithful port
/// of Perl `process_prediction`: the header's coordspan (whitespace-token index
/// 6, e.g. `842-3150`) and every data row's first two tab columns are offset by
/// `partition_lend - 1`. INTRON rows are retained (Perl keeps them).
fn process_prediction_text(text: &str, partition_lend: u32) -> Option<PartitionPred> {
    let offset = partition_lend.saturating_sub(1);
    let mut lines_iter = text.lines();

    // First line is the header. Offset the coordspan at positional index 6.
    let header_line = lines_iter.next()?;
    let mut header_parts: Vec<String> = header_line
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();
    if header_parts.len() > 6 {
        if let Some((l, r)) = header_parts[6].split_once('-') {
            if let (Ok(l), Ok(r)) = (l.parse::<u32>(), r.parse::<u32>()) {
                header_parts[6] = format!("{}-{}", l + offset, r + offset);
            }
        }
    }
    let mut new_text = format!("{}\n", header_parts.join(" "));

    let mut exon_types: Vec<String> = Vec::new();
    let mut all_coords: Vec<u32> = Vec::new();

    for data_line in lines_iter {
        let cols: Vec<&str> = data_line.split('\t').collect();
        if cols.len() >= 3 {
            if let (Ok(e5), Ok(e3)) = (cols[0].parse::<u32>(), cols[1].parse::<u32>()) {
                let etype = cols[2].to_string();
                let new_e5 = e5 + offset;
                let new_e3 = e3 + offset;
                all_coords.push(new_e5);
                all_coords.push(new_e3);
                exon_types.push(etype);
                // Perl: split on '\t', overwrite x[0]/x[1], rejoin ALL columns —
                // x[2] is the exon type (e.g. `initial+`), which must be kept.
                let rest: String = cols[2..].join("\t");
                new_text.push_str(&format!("{}\t{}\t{}\n", new_e5, new_e3, rest));
            }
        }
    }

    if all_coords.is_empty() {
        return None;
    }

    let gene_lend = *all_coords.iter().min()?;
    let gene_rend = *all_coords.iter().max()?;
    let length = gene_rend - gene_lend + 1;

    let type_str = exon_types.join(",");
    let class = if (type_str.contains("initial") && type_str.contains("terminal"))
        || type_str.contains("single")
    {
        PredClass::Complete
    } else {
        PredClass::Partial
    };

    Some(PartitionPred {
        lend: gene_lend,
        rend: gene_rend,
        class,
        text: new_text,
        length,
        path_score: length,
        prev_link: None,
        intronic_preds: Vec::new(),
        encaps: false,
    })
}

/// Identify predictions encapsulated within the introns of others and nest them.
pub fn join_intronic_preds(mut preds: Vec<PartitionPred>) -> Vec<PartitionPred> {
    preds.sort_by_key(|p| p.lend);

    let n = preds.len();
    let mut encaps = vec![false; n];

    for i in 0..n {
        for j in (i + 1)..n {
            if preds[j].lend > preds[i].lend && preds[j].rend < preds[i].rend {
                encaps[j] = true;
                let nested = preds[j].clone();
                let len_j = nested.length;
                let score_j = nested.path_score;
                preds[i].intronic_preds.push(nested);
                preds[i].length += len_j;
                preds[i].path_score += score_j;
            }
        }
    }

    preds
        .into_iter()
        .zip(encaps)
        .filter_map(|(p, enc)| if enc { None } else { Some(p) })
        .collect()
}

/// Dynamic-programming combination of predictions from multiple overlapping partitions.
/// Selects the maximal set of non-overlapping complete genes.
pub fn combine_predictions(mut preds: Vec<PartitionPred>) -> Vec<PartitionPred> {
    if preds.is_empty() {
        return preds;
    }
    preds.sort_by_key(|p| p.lend);
    let n = preds.len();

    for i in 1..n {
        let lend_i = preds[i].lend;
        let length_i = preds[i].length;
        let mut best_score = preds[i].path_score;
        let mut best_link: Option<usize> = None;

        for j in (0..i).rev() {
            let rend_j = preds[j].rend;
            if rend_j < lend_i {
                let candidate = preds[j].path_score + length_i;
                if candidate > best_score {
                    best_score = candidate;
                    best_link = Some(j);
                }
            }
        }
        preds[i].path_score = best_score;
        preds[i].prev_link = best_link;
    }

    // Find highest-scoring end
    let best_end = preds
        .iter()
        .enumerate()
        .max_by_key(|(_, p)| p.path_score)
        .map(|(i, _)| i);

    let mut result: Vec<PartitionPred> = Vec::new();
    let mut cur = best_end;
    while let Some(idx) = cur {
        result.push(preds[idx].clone());
        cur = preds[idx].prev_link;
    }
    result.reverse();
    result
}

/// Run recombination for all contigs listed in `entries`.
pub fn recombine_outputs(entries: &[PartitionEntry], output_file_name: &str) -> Result<()> {
    // Group entries by base_dir
    let mut base_to_partitions: HashMap<String, Vec<(String, u32)>> = HashMap::new();
    for entry in entries {
        if entry.is_partitioned {
            if let Some(pdir) = &entry.partition_dir {
                // Extract lend from partition dir name "..._LEND-REND". Perl
                // `recombine_EVM_partial_outputs.pl` does
                // `$partition_dir =~ /(\d+)-(\d+)$/ or die` — fail loudly rather
                // than silently defaulting to a wrong (lend=1, offset=0) mapping.
                let lend = extract_partition_lend(pdir).ok_or_else(|| {
                    anyhow::anyhow!("Error, cannot extract coords from partition dir {}", pdir)
                })?;
                base_to_partitions
                    .entry(entry.base_dir.clone())
                    .or_default()
                    .push((pdir.clone(), lend));
            }
        }
    }

    for (base_dir, partition_dirs) in &base_to_partitions {
        let mut all_preds: Vec<PartitionPred> = Vec::new();
        for (pdir, lend) in partition_dirs {
            let output_path = format!("{}/{}", pdir, output_file_name);
            parse_and_add_predictions(&output_path, *lend, &mut all_preds)?;
        }

        let final_preds = combine_predictions(all_preds);
        let out_path = format!("{}/{}", base_dir, output_file_name);
        let mut out =
            fs::File::create(&out_path).with_context(|| format!("Cannot create {}", out_path))?;

        log::debug!("Writing combined output to {}", out_path);
        // Perl prints "$pred_text\n" — pred.text already ends in '\n', so the
        // extra newline yields a blank line separating predictions (which
        // EVM_to_GFF3 relies on to delimit gene models).
        for pred in &final_preds {
            writeln!(out, "{}", pred.text)?;
            for nested in &pred.intronic_preds {
                writeln!(out, "!! Intron-containing prediction")?;
                writeln!(out, "{}", nested.text)?;
            }
        }
    }
    Ok(())
}

fn extract_partition_lend(pdir: &str) -> Option<u32> {
    // Perl: `$partition_dir =~ /(\d+)-(\d+)$/` — the lend is the run of digits
    // immediately before the final `<digits>` group (robust to accessions that
    // themselves contain '-' or '_').
    let name = std::path::Path::new(pdir).file_name()?.to_str()?;
    let (head, rend) = name.rsplit_once('-')?;
    if rend.is_empty() || !rend.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let lend: String = head
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    if lend.is_empty() {
        return None;
    }
    lend.chars().rev().collect::<String>().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_partition_lend_matches_perl_regex() {
        // /(\d+)-(\d+)$/ → lend is the first capture
        assert_eq!(
            extract_partition_lend("partitions/Contig1/Contig1_20001-50000"),
            Some(20001)
        );
        assert_eq!(extract_partition_lend("Contig1_1-30000"), Some(1));
        // accession containing '-' / '_' must not confuse the extractor
        assert_eq!(
            extract_partition_lend("a/scaf-2_b/scaf-2_b_40001-63304"),
            Some(40001)
        );
    }

    #[test]
    fn join_intronic_preds_nests_encapsulated_predictions() {
        // Outer gene spans 100-1000; inner gene sits at 400-500 fully inside.
        let outer = PartitionPred {
            lend: 100,
            rend: 1000,
            class: PredClass::Complete,
            text: "# outer\n100\t200\tinitial+\t1\t1\t\n201\t300\tINTRON\t\t\t\n301\t1000\tterminal+\t1\t3\t\n".to_string(),
            length: 901,
            path_score: 901,
            prev_link: None,
            intronic_preds: Vec::new(),
            encaps: false,
        };
        let inner = PartitionPred {
            lend: 400,
            rend: 500,
            class: PredClass::Complete,
            text: "# inner\n400\t500\tsingle+\t1\t3\t\n".to_string(),
            length: 101,
            path_score: 101,
            prev_link: None,
            intronic_preds: Vec::new(),
            encaps: false,
        };
        let non_overlap = PartitionPred {
            lend: 1100,
            rend: 1200,
            class: PredClass::Complete,
            text: "# non-overlap\n1100\t1200\tsingle+\t1\t3\t\n".to_string(),
            length: 101,
            path_score: 101,
            prev_link: None,
            intronic_preds: Vec::new(),
            encaps: false,
        };

        let joined = join_intronic_preds(vec![outer, inner, non_overlap]);
        // Inner should be removed from top level and nested under outer.
        assert_eq!(
            joined.len(),
            2,
            "outer + non_overlap should remain at top level"
        );
        assert_eq!(joined[0].lend, 100);
        assert_eq!(joined[0].rend, 1000);
        assert_eq!(
            joined[0].intronic_preds.len(),
            1,
            "inner should be nested in outer"
        );
        assert_eq!(joined[0].intronic_preds[0].lend, 400);
        assert_eq!(joined[0].intronic_preds[0].rend, 500);
        // Outer length/path_score should absorb inner's contribution.
        assert_eq!(joined[0].length, 901 + 101);
        assert_eq!(joined[0].path_score, 901 + 101);
        assert_eq!(joined[1].lend, 1100);
        assert_eq!(joined[1].rend, 1200);
    }

    #[test]
    fn process_prediction_offsets_coords_and_keeps_type_column() {
        // partition_lend = 20001 → offset 20000 added to header coordspan (token 6)
        // and to the first two tab columns of each data row. The exon-type column
        // (`initial+` / `INTRON`) must be preserved (Perl rejoins ALL columns).
        let text = "\
# EVM prediction: Mode:STANDARD S-ratio:1.00 g1 100-310 + score(1.0) noncoding(0.0) raw(0.0) offset(0.0)
100\t200\tinitial+\t1\t1\t{src_a;src}
201\t250\tINTRON\t\t\t{src_b;src}
251\t310\tterminal+\t1\t3\t{src_a;src}
";
        let pred = process_prediction_text(text, 20001).expect("parsed");
        // header token 6 (0-based) was the coordspan `100-310`
        let header = pred.text.lines().next().unwrap();
        let tok6 = header.split_whitespace().nth(6).unwrap();
        assert_eq!(tok6, "20100-20310");
        // data rows: coords offset by 20000, type column intact
        let rows: Vec<&str> = pred.text.lines().skip(1).collect();
        assert_eq!(rows[0], "20100\t20200\tinitial+\t1\t1\t{src_a;src}");
        assert_eq!(rows[1], "20201\t20250\tINTRON\t\t\t{src_b;src}");
        assert_eq!(rows[2], "20251\t20310\tterminal+\t1\t3\t{src_a;src}");
        // gene span + completeness (initial + terminal → complete)
        assert_eq!(pred.lend, 20100);
        assert_eq!(pred.rend, 20310);
        assert!(matches!(pred.class, PredClass::Complete));
    }
}
