//! Intergenic region scoring.

use crate::io::gff3::Gff3Record;
use crate::types::evidence::EvWeightMap;
use crate::types::exon::{Exon, ExonType};
use crate::types::genome::MaskVec;
use std::collections::HashMap;

/// Intergenic per-base scores with a cached prefix-sum array for O(1) range
/// queries. `calc_intergenic_score` is called millions of times during the
/// trellis, so summing ranges in a loop became the dominant cost.
#[derive(Clone, Debug)]
pub struct IntergenicScores {
    per_base: Vec<f64>,
    prefix: Vec<f64>,
}

impl IntergenicScores {
    fn from_per_base(per_base: Vec<f64>) -> Self {
        let mut prefix = Vec::with_capacity(per_base.len() + 1);
        prefix.push(0.0);
        let mut sum = 0.0;
        for &v in &per_base {
            sum += v;
            prefix.push(sum);
        }
        Self { per_base, prefix }
    }

    #[inline]
    pub fn calc(&self, lend: u32, rend: u32) -> f64 {
        if lend > rend {
            return 0.0;
        }
        let max_1based = (self.per_base.len().saturating_sub(1)) as u32;
        let l = lend.max(1) as usize;
        let r = rend.min(max_1based) as usize;
        // prefix[k] = sum of per_base[0..k]; sum over [l, r] = prefix[r+1] - prefix[l].
        self.prefix[r + 1] - self.prefix[l]
    }

    pub fn per_base(&self) -> &[f64] {
        &self.per_base
    }

    pub fn per_base_mut(&mut self) -> &mut [f64] {
        &mut self.per_base
    }

    pub fn rebuild_prefix(&mut self) {
        self.prefix.clear();
        self.prefix.push(0.0);
        let mut sum = 0.0;
        for &v in &self.per_base {
            sum += v;
            self.prefix.push(sum);
        }
    }

    pub fn len(&self) -> usize {
        self.per_base.len()
    }

    pub fn is_empty(&self) -> bool {
        self.per_base.is_empty()
    }
}

/// Compute per-base intergenic scores from gene-prediction spans.
///
/// Mirrors the Perl `populate_intergenic_regions`: for each ab-initio
/// prediction type, the per-base score in the gap between two neighbouring
/// genes of that type is incremented by that type's weight. Only
/// `ABINITIO_PREDICTION` types contribute; all genes on both strands are
/// considered. Masked positions are skipped. The genome ends [0,0] and
/// [seq_len,seq_len] are added as sentinels so the flanking regions count.
///
/// The `INTERGENIC_SCORE_ADJUST_FACTOR` is folded into the per-base value
/// here; with the default factor of 1.0 this is identical to applying it in
/// `calc_intergenic_score` as the Perl does.
pub fn populate_intergenic_scores(
    seq_len: usize,
    gene_pred_records: &[Gff3Record],
    ev_weights: &EvWeightMap,
    mask: &MaskVec,
    adjust_factor: f64,
) -> IntergenicScores {
    let mut ig = vec![0.0f64; seq_len + 2];

    // Group CDS spans per (ev_type, model id), forward-coordinate span.
    // model_spans[ev_type][model_id] = (min_coord, max_coord)
    let mut per_type: HashMap<String, HashMap<String, (u32, u32)>> = HashMap::new();

    for rec in gene_pred_records {
        if rec.feature != "CDS" {
            continue;
        }
        let entry = match ev_weights.get(&rec.source) {
            Some(e) => e,
            None => continue,
        };
        if !entry.ev_class.is_abinitio() {
            continue;
        }

        // Perl get_gene_predictions groups CDS by the FULL attribute column
        // ($x[8]), not by Parent. This matters when otherwise-identical CDS of
        // one model carry a differing attribute suffix (e.g. `5_prime_partial=true`):
        // Perl then treats them as separate "genes", creating an intergenic gap
        // in the intron between them. Group by the raw attribute string to match.
        let group_key = rec.raw_attributes.clone();
        if group_key.is_empty() {
            continue;
        }

        let lend = rec.start.min(rec.end);
        let rend = rec.start.max(rec.end);
        let models = per_type.entry(rec.source.clone()).or_default();
        let span = models.entry(group_key).or_insert((lend, rend));
        span.0 = span.0.min(lend);
        span.1 = span.1.max(rend);
    }

    // Iterate over ALL ab-initio types in the weights file, not just those
    // with gene predictions in this partition. Perl's populate_intergenic_regions
    // loops over %PREDICTION_PROGS_CONTRIBUTE_INTERGENIC (populated from the
    // weights file), so a predictor with no predictions in the current partition
    // still contributes its weight to every intergenic position.
    for (ev_type, entry) in ev_weights {
        if !entry.ev_class.is_abinitio() {
            continue;
        }
        let weight = entry.weight * adjust_factor;
        let models = per_type.get(ev_type);

        // Collect gene spans plus genome-boundary sentinels, sorted by lend.
        let mut spans: Vec<(u32, u32)> = match models {
            Some(m) => m.values().copied().collect(),
            None => Vec::new(),
        };
        spans.push((0, 0));
        spans.push((seq_len as u32, seq_len as u32));
        spans.sort_by_key(|&(l, _)| l);

        for w in spans.windows(2) {
            let lend_intergenic = w[0].1;
            let rend_intergenic = w[1].0;
            if lend_intergenic > rend_intergenic {
                continue;
            }
            let mut j = lend_intergenic + 1;
            while j < rend_intergenic {
                if !mask.get(j as usize) {
                    ig[j as usize] += weight;
                }
                j += 1;
            }
        }
    }

    IntergenicScores::from_per_base(ig)
}

/// Compute the sum of intergenic scores over [lend, rend] (inclusive, 1-based).
pub fn calc_intergenic_score(scores: &IntergenicScores, lend: u32, rend: u32) -> f64 {
    scores.calc(lend, rend)
}

/// Faithful port of Perl `augment_intergenic_from_start/stop_peaks`.
///
/// For each start/stop peak, find the closest matching initial/single (start) or
/// terminal/single (stop) exon on the peak's strand within `start_stop_range`,
/// then walk outward to the nearest "relevant" exon of the opposite role and
/// **set** every non-masked intergenic position in that flanking span to
/// `sum_genepred_weights`. Exons are processed in end5-ascending order (the Perl
/// wrapper sorts `@EXONS` before augmenting).
#[allow(clippy::too_many_arguments)]
pub fn augment_intergenic_from_start_stop_peaks(
    ig: &mut IntergenicScores,
    start_peaks: &[(u32, f64, char)],
    end_peaks: &[(u32, f64, char)],
    exons: &[Exon],
    mask: &MaskVec,
    seq_len: u32,
    sum_genepred_weights: f64,
    start_stop_range: u32,
) {
    // Perl sorts @EXONS by end5 ascending before augmenting.
    let mut sorted: Vec<&Exon> = exons.iter().collect();
    sorted.sort_by_key(|e| e.end5);

    let set_range = |ig: &mut IntergenicScores, lo: u32, hi: u32| {
        if lo > hi {
            return;
        }
        let base = ig.per_base_mut();
        for i in lo..=hi {
            let iu = i as usize;
            if iu < base.len() && !mask.get(iu) {
                base[iu] = sum_genepred_weights;
            }
        }
    };

    // ── START peaks ────────────────────────────────────────────────────────
    for &(pos, _score, strand) in start_peaks {
        if strand == '+' {
            // closest initial|single + exon within range
            let closest = match find_closest_exon(
                &sorted,
                &[ExonType::Initial, ExonType::Single],
                '+',
                pos,
                start_stop_range,
            ) {
                Some(e) => e,
                None => continue,
            };
            let position = closest.coords_sorted().0; // lend (5' of initial+)
                                                      // walk right→left: nearest exon with rend < position that ends a gene
            let mut exon_rend = 1u32;
            for &e in sorted.iter().rev() {
                let (_, e_rend) = e.coords_sorted();
                if e_rend >= position {
                    continue;
                }
                let o = e.orientation.as_char();
                if (o == '+' && is_terminal_or_single(e)) || (o == '-' && is_initial_or_single(e)) {
                    exon_rend = e.coords_sorted().1;
                    break;
                }
            }
            set_range(ig, exon_rend, position);
        } else {
            let closest = match find_closest_exon(
                &sorted,
                &[ExonType::Initial, ExonType::Single],
                '-',
                pos,
                start_stop_range,
            ) {
                Some(e) => e,
                None => continue,
            };
            let position = closest.coords_sorted().1; // rend (5' of initial-)
            let mut exon_lend = seq_len;
            for &e in sorted.iter() {
                let (e_lend, _) = e.coords_sorted();
                if e_lend <= position {
                    continue;
                }
                let o = e.orientation.as_char();
                if (o == '+' && is_initial_or_single(e)) || (o == '-' && is_terminal_or_single(e)) {
                    exon_lend = e.coords_sorted().0;
                    break;
                }
            }
            set_range(ig, position, exon_lend);
        }
    }

    // ── STOP peaks ─────────────────────────────────────────────────────────
    for &(pos, _score, strand) in end_peaks {
        if strand == '+' {
            let closest = match find_closest_exon(
                &sorted,
                &[ExonType::Terminal, ExonType::Single],
                '+',
                pos,
                start_stop_range,
            ) {
                Some(e) => e,
                None => continue,
            };
            let position = closest.coords_sorted().1; // rend (3' of terminal+)
            let mut exon_lend = seq_len;
            for &e in sorted.iter() {
                let (e_lend, _) = e.coords_sorted();
                if e_lend <= position {
                    continue;
                }
                let o = e.orientation.as_char();
                if (o == '+' && is_initial_or_single(e)) || (o == '-' && is_terminal_or_single(e)) {
                    exon_lend = e.coords_sorted().0;
                    break;
                }
            }
            set_range(ig, position, exon_lend);
        } else {
            let closest = match find_closest_exon(
                &sorted,
                &[ExonType::Terminal, ExonType::Single],
                '-',
                pos,
                start_stop_range,
            ) {
                Some(e) => e,
                None => continue,
            };
            let position = closest.coords_sorted().0; // lend (3' of terminal-)
            let mut exon_rend = 1u32;
            for &e in sorted.iter().rev() {
                let (_, e_rend) = e.coords_sorted();
                if e_rend >= position {
                    continue;
                }
                let o = e.orientation.as_char();
                if (o == '+' && is_terminal_or_single(e)) || (o == '-' && is_initial_or_single(e)) {
                    exon_rend = e.coords_sorted().1;
                    break;
                }
            }
            set_range(ig, exon_rend, position);
        }
    }

    ig.rebuild_prefix();
}

fn is_initial_or_single(e: &Exon) -> bool {
    matches!(e.exon_type, ExonType::Initial | ExonType::Single)
}
fn is_terminal_or_single(e: &Exon) -> bool {
    matches!(e.exon_type, ExonType::Terminal | ExonType::Single)
}

/// Perl `find_closest_exon_within_range`: among exons matching any `types` and
/// the given `strand`, return the one whose nearest endpoint is closest to
/// `position`, within `range`. Ties keep the first in iteration order.
fn find_closest_exon<'a>(
    exons: &[&'a Exon],
    types: &[ExonType],
    strand: char,
    position: u32,
    range: u32,
) -> Option<&'a Exon> {
    let mut closest: Option<&Exon> = None;
    let mut closest_dist: Option<u32> = None;
    for &e in exons {
        if !types.contains(&e.exon_type) {
            continue;
        }
        if e.orientation.as_char() != strand {
            continue;
        }
        let d5 = e.end5.abs_diff(position);
        let d3 = e.end3.abs_diff(position);
        let delta = d5.min(d3);
        if delta > range {
            continue;
        }
        match closest_dist {
            Some(cd) if cd <= delta => {}
            _ => {
                closest_dist = Some(delta);
                closest = Some(e);
            }
        }
    }
    closest
}

/// Get contiguous intergenic regions between adjacent predictions.
/// Returns Vec<(lend, rend)> pairs.
pub fn get_intergenic_regions(
    predictions_sorted_by_lend: &[(u32, u32)], // (lend, rend)
) -> Vec<(u32, u32)> {
    let mut regions = Vec::new();
    for pair in predictions_sorted_by_lend.windows(2) {
        let rend_prev = pair[0].1;
        let lend_next = pair[1].0;
        if lend_next > rend_prev + 1 {
            regions.push((rend_prev + 1, lend_next - 1));
        }
    }
    regions
}
