pub mod evm_to_gff3;
pub mod gff3_to_bed;
pub mod gff3_to_proteins;

/// Decode percent-encoded octets in a GFF3 attribute value (e.g. `%20` → space),
/// matching the URI un-escaping that Perl `GFF3_utils::index_GFF3_gene_objs`
/// applies when it reads `Name=` into a gene object's `com_name`.
pub(crate) fn uri_unescape(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            let hi = (bytes[i + 1] as char).to_digit(16);
            let lo = (bytes[i + 2] as char).to_digit(16);
            if let (Some(h), Some(l)) = (hi, lo) {
                out.push((h * 16 + l) as u8);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}
