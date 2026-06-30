//! Streaming GFF3 parser.

use anyhow::Result;
use std::collections::HashMap;
use std::io::BufRead;

/// A single parsed GFF3 feature line.
#[derive(Debug, Clone)]
pub struct Gff3Record {
    pub seqid: String,
    pub source: String,
    pub feature: String,
    pub start: u32, // 1-based, inclusive
    pub end: u32,   // 1-based, inclusive
    pub score: Option<f64>,
    pub strand: char,      // '+' or '-' or '.'
    pub phase: Option<u8>, // 0, 1, or 2
    pub attributes: HashMap<String, String>,
    pub raw_attributes: String,
}

impl Gff3Record {
    /// Parse a GFF3 data line (already confirmed non-comment, non-empty).
    pub fn parse(line: &str) -> Result<Gff3Record> {
        let cols: Vec<&str> = line.splitn(9, '\t').collect();
        if cols.len() < 8 {
            anyhow::bail!("GFF3 line has fewer than 8 columns: {}", line);
        }
        let start: u32 = cols[3]
            .parse()
            .map_err(|_| anyhow::anyhow!("Bad GFF3 start: {}", cols[3]))?;
        let end: u32 = cols[4]
            .parse()
            .map_err(|_| anyhow::anyhow!("Bad GFF3 end: {}", cols[4]))?;
        let score = match cols[5] {
            "." => None,
            s => Some(s.parse::<f64>().unwrap_or(0.0)),
        };
        let strand = cols[6].chars().next().unwrap_or('.');
        let phase = match cols[7] {
            "." => None,
            s => s.parse::<u8>().ok(),
        };
        let raw_attributes = if cols.len() >= 9 {
            cols[8].to_string()
        } else {
            String::new()
        };
        let attributes = parse_attributes(&raw_attributes);

        Ok(Gff3Record {
            seqid: cols[0].to_string(),
            source: cols[1].to_string(),
            feature: cols[2].to_string(),
            start,
            end,
            score,
            strand,
            phase,
            attributes,
            raw_attributes,
        })
    }

    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attributes.get(key).map(|s| s.as_str())
    }
}

/// Parse GFF3 attribute column into a map.
pub fn parse_attributes(attrs: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for kv in attrs.split(';') {
        let kv = kv.trim();
        if let Some(eq) = kv.find('=') {
            let key = kv[..eq].trim().to_string();
            // Take only the first space-delimited token of the value
            // (mimics how the Perl code parses Target/Query).
            let val = kv[eq + 1..]
                .trim()
                .split(' ')
                .next()
                .unwrap_or("")
                .to_string();
            map.insert(key, val);
        }
    }
    map
}

/// Read all GFF3 records from a `BufRead`, skipping comments and blank lines.
pub fn read_gff3<R: BufRead>(reader: R) -> Result<Vec<Gff3Record>> {
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        records.push(Gff3Record::parse(&line)?);
    }
    Ok(records)
}

/// Read all GFF3 records from a file path.
pub fn read_gff3_file(path: &str) -> Result<Vec<Gff3Record>> {
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("Cannot open GFF3 file {}: {}", path, e))?;
    read_gff3(std::io::BufReader::new(f))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cds_line() {
        let line = "chr1\tgenemark\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=mrna1";
        let rec = Gff3Record::parse(line).unwrap();
        assert_eq!(rec.seqid, "chr1");
        assert_eq!(rec.feature, "CDS");
        assert_eq!(rec.start, 100);
        assert_eq!(rec.end, 200);
        assert_eq!(rec.strand, '+');
        assert_eq!(rec.phase, Some(0));
        assert_eq!(rec.attr("Parent"), Some("mrna1"));
    }
}
