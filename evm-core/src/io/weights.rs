//! Weights file parser.
//!
//! File format (whitespace-delimited, 3 columns per line):
//!   PROTEIN         ev_type_name   1
//!   ABINITIO_PREDICTION  genscan  5
//!   # comment lines are ignored

use crate::types::evidence::{EvClass, EvEntry, EvWeightMap};
use anyhow::Result;
use std::io::BufRead;

/// Parse a weights file from any `BufRead` and return an `EvWeightMap`.
pub fn parse_weights<R: BufRead>(reader: R) -> Result<EvWeightMap> {
    let mut map = EvWeightMap::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 3 {
            anyhow::bail!("Weights line has fewer than 3 columns: {}", trimmed);
        }
        let ev_class = EvClass::from_str(parts[0])
            .ok_or_else(|| anyhow::anyhow!("Unknown evidence class: {}", parts[0]))?;
        let ev_type = parts[1].to_string();
        let weight: f64 = parts[2]
            .parse()
            .map_err(|_| anyhow::anyhow!("Non-numeric weight for {}: {}", ev_type, parts[2]))?;
        map.insert(ev_type, EvEntry { ev_class, weight });
    }
    Ok(map)
}

/// Convenience: read a weights file from disk.
pub fn read_weights_file(path: &str) -> Result<EvWeightMap> {
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("Cannot open weights file {}: {}", path, e))?;
    parse_weights(std::io::BufReader::new(f))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::evidence::EvClass;

    #[test]
    fn parse_weights_basic() {
        let data = b"# comment\nPROTEIN\tnap\t1\nABINITIO_PREDICTION\tgenscan\t5\n";
        let w = parse_weights(data.as_ref()).unwrap();
        assert_eq!(w["nap"].weight, 1.0);
        assert_eq!(w["nap"].ev_class, EvClass::Protein);
        assert_eq!(w["genscan"].weight, 5.0);
        assert_eq!(w["genscan"].ev_class, EvClass::AbinitioPrediction);
    }
}
