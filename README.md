# EVidenceModeler (Rust)

A Rust port of [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler),
a framework for combining gene predictions and evidence alignments into
consensus gene models. The Rust rewrite is a drop-in replacement for the
Perl `evidence_modeler.pl` and the `EVidenceModeler` wrapper script.

## Building from source

Requirements:
- Rust toolchain (stable, 1.70+)

```sh
git clone https://github.com/hyphaltip/EVidenceModeler_rust.git
cd EVidenceModeler_rust
make            # builds release binaries into target/release/
```

This produces six binaries:

| Binary | Purpose |
|--------|---------|
| `EVidenceModeler` | Main orchestrator (partitions, runs EVM, recombines, converts) |
| `evidence_modeler` | Core EVM engine for a single partition |
| `partition_evm_inputs` | Splits genome + GFF3 into overlapping partitions |
| `recombine_evm_outputs` | Merges per-partition EVM outputs |
| `convert_EVM_outputs_to_GFF3` | Converts raw EVM output to GFF3 |
| `gff3_file_to_proteins` | Translates GFF3 CDS features to protein sequences |

## Installation

### Option 1: Install to a user-local bin directory

```sh
make
cp target/release/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins} ~/bin/
```

Ensure `~/bin` is on your `PATH`:

```sh
export PATH="$HOME/bin:$PATH"
```

### Option 2: System-wide installation

```sh
make
sudo cp target/release/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins} /usr/local/bin/
```

### Option 3: cargo install (from a local checkout)

```sh
cargo install --path evm-cli
cargo install --path evm-utils
```

This installs the binaries into `~/.cargo/bin/`.

## Running the test suite

```sh
make test       # runs cargo test + end-to-end pipeline test
```

Or manually:

```sh
cargo test --all
cd testing && bash runMe.rust.sh
```

## Usage

### As a standalone tool

```sh
EVidenceModeler \
    --sample_id mygenome \
    --genome genome.fa \
    --weights weights.txt \
    --gene_predictions genes.gff3 \
    --protein_alignments proteins.gff3 \
    --transcript_alignments transcripts.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000
```

This produces `mygenome.EVM.{gff3,bed,pep,cds}` in the current directory.

### Integrated with funannotate

The Rust EVM engine integrates with [funannotate](https://github.com/nextgenusfs/funannotate)
as a drop-in replacement for the Perl `evidence_modeler.pl`. See
[INSTALL.funannotate.md](INSTALL.funannotate.md) for detailed setup
instructions.

## Parity with the Perl implementation

The Rust port achieves **99.4% gene-level parity** and **97.0% protein
sequence parity** with the Perl implementation on a real fungal genome
(*Rhodotorula sphaerocarpa*, 24 scaffolds, ~20 Mb). The remaining
differences are due to partition boundary effects and same-score
tie-breaking in the DP trellis. See [PARITY_REPORT.md](PARITY_REPORT.md)
for full details.

## License

MIT — see [LICENSE.txt](LICENSE.txt).
