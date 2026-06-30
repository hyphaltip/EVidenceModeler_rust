//! Streaming FASTA reader.

use anyhow::Result;
use std::io::{self, BufRead};

/// A single FASTA record.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub header: String,
    /// Accession: first whitespace-delimited token of the header.
    pub accession: String,
    /// Uppercase DNA/protein sequence.
    pub sequence: String,
}

impl FastaRecord {
    /// Return FASTA-formatted string (60-char line wrap).
    pub fn to_fasta(&self) -> String {
        let mut out = format!(">{}\n", self.header);
        for chunk in self.sequence.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push('\n');
        }
        out
    }
}

/// Streaming FASTA reader backed by any `BufRead`.
pub struct FastaReader<R: BufRead> {
    reader: R,
    next_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(mut reader: R) -> Result<Self> {
        let mut line = String::new();
        // skip blank lines / find first header
        loop {
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(FastaReader {
                    reader,
                    next_header: None,
                });
            }
            let trimmed = line.trim_end();
            if let Some(stripped) = trimmed.strip_prefix('>') {
                return Ok(FastaReader {
                    reader,
                    next_header: Some(stripped.to_string()),
                });
            }
        }
    }

    /// Read all records into a Vec.
    pub fn collect_all(reader: R) -> Result<Vec<FastaRecord>> {
        let mut fr = Self::new(reader)?;
        let mut out = Vec::new();
        while let Some(rec) = fr.next()? {
            out.push(rec);
        }
        Ok(out)
    }

    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<FastaRecord>> {
        let header = match self.next_header.take() {
            Some(h) => h,
            None => return Ok(None),
        };
        let accession = header.split_whitespace().next().unwrap_or("").to_string();
        let mut seq = String::new();
        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end();
            if let Some(stripped) = trimmed.strip_prefix('>') {
                self.next_header = Some(stripped.to_string());
                break;
            }
            seq.push_str(trimmed.trim());
        }
        let sequence: String = seq
            .bytes()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        Ok(Some(FastaRecord {
            header,
            accession,
            sequence,
        }))
    }
}

/// Convenience: read all records from a file path.
pub fn read_fasta_file(path: &str) -> Result<Vec<FastaRecord>> {
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("Cannot open FASTA file {}: {}", path, e))?;
    FastaReader::collect_all(io::BufReader::new(f))
}

/// Convenience: read all records into a HashMap keyed by accession.
pub fn read_fasta_hash(path: &str) -> Result<std::collections::HashMap<String, String>> {
    let records = read_fasta_file(path)?;
    let mut map = std::collections::HashMap::new();
    for rec in records {
        map.insert(rec.accession, rec.sequence);
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_fasta() {
        let data = b">seq1 description\nATGCGT\nAAA\n>seq2\nCCCC\n";
        let records = FastaReader::collect_all(data.as_ref()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].accession, "seq1");
        assert_eq!(records[0].sequence, "ATGCGTAAA");
        assert_eq!(records[1].accession, "seq2");
        assert_eq!(records[1].sequence, "CCCC");
    }

    #[test]
    fn to_fasta_wraps_at_60() {
        let rec = FastaRecord {
            header: "seq1".into(),
            accession: "seq1".into(),
            sequence: "A".repeat(120),
        };
        let out = rec.to_fasta();
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], ">seq1");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 60);
    }
}
