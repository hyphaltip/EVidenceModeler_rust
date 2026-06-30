#!/bin/sh
set -ex

## run Rust EVM pipeline test
cd "$(dirname "$0")"

RUST_BIN="../target/debug/EVidenceModeler"
if [ ! -x "$RUST_BIN" ]; then
    RUST_BIN="../target/release/EVidenceModeler"
fi
if [ ! -x "$RUST_BIN" ]; then
    echo "Rust binary not found. Run: cargo build"
    exit 1
fi

"$RUST_BIN" \
    --sample_id rusttest \
    --genome genome.fasta \
    --weights ./weights.txt \
    --gene_predictions gene_predictions.gff3 \
    --protein_alignments protein_alignments.gff3 \
    --transcript_alignments transcript_alignments.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000

echo "Done. See rusttest.EVM.* outputs"
