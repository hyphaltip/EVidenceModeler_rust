#!/bin/bash
# End-to-end Rust EVM pipeline test.
# Runs the Rust EVidenceModeler on the testing/ dataset and validates output
# against the Perl golden fixtures in evm/tests/fixtures/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
RUST_BIN="$REPO/target/debug/EVidenceModeler"
FIXTURES="$REPO/tests/fixtures"

if [[ ! -x "$RUST_BIN" ]]; then
    echo "Rust binary not found: $RUST_BIN"
    echo "Run: cargo build"
    exit 1
fi

cd "$SCRIPT_DIR"

# Clean previous Rust outputs (NOT the Perl golden directories)
rm -rf __rusttest-EVM_chckpts rusttest.partitions rusttest.partitions.listing
rm -f rusttest.EVM.gff3 rusttest.EVM.pep rusttest.EVM.cds rusttest.EVM.bed

echo "=== Running Rust EVidenceModeler on testing/ ==="
"$RUST_BIN" \
    --sample_id rusttest \
    --genome genome.fasta \
    --weights ./weights.txt \
    --gene_predictions gene_predictions.gff3 \
    --protein_alignments protein_alignments.gff3 \
    --transcript_alignments transcript_alignments.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000

echo ""
echo "=== Comparing to Perl golden fixtures ==="

# GFF3: byte-identical
diff rusttest.EVM.gff3 "$FIXTURES/Contig1.perl.EVM.gff3" \
    && echo "PASS: GFF3 byte-identical" \
    || { echo "FAIL: GFF3 differs from Perl golden"; exit 1; }

# BED: byte-identical
diff rusttest.EVM.bed "$FIXTURES/Contig1.perl.EVM.bed" \
    && echo "PASS: BED byte-identical" \
    || { echo "FAIL: BED differs from Perl golden"; exit 1; }

# PEP: record-identical (sorted by header; Perl FASTA order is non-deterministic)
sort_fasta() {
    python3 -c "
import sys
recs = {}; header = None; seq = []
with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if header: recs[header] = seq
            header = line; seq = []
        else: seq.append(line)
    if header: recs[header] = seq
for h in sorted(recs):
    print(h)
    for s in recs[h]: print(s)
" "$1"
}

diff <(sort_fasta rusttest.EVM.pep) <(sort_fasta "$FIXTURES/Contig1.perl.EVM.pep") \
    && echo "PASS: PEP records identical (sorted)" \
    || { echo "FAIL: PEP differs from Perl golden"; exit 1; }

diff <(sort_fasta rusttest.EVM.cds) <(sort_fasta "$FIXTURES/Contig1.perl.EVM.cds") \
    && echo "PASS: CDS records identical (sorted)" \
    || { echo "FAIL: CDS differs from Perl golden"; exit 1; }

echo ""
echo "All tests PASSED."
