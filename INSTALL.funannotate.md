# Installing the Rust EVidenceModeler and Integrating with funannotate

This guide covers building the Rust EVM binaries, installing them on
your system, and configuring [funannotate](https://github.com/nextgenusfs/funannotate)
to use them as a drop-in replacement for the Perl `evidence_modeler.pl`.

---

## 1. Prerequisites

- **Rust toolchain** (stable, 1.70+). Install via [rustup](https://rustup.rs/):
  ```sh
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  source "$HOME/.cargo/env"
  ```
- **funannotate** (1.9.0+ with Rust EVM support). The `funannotate-runEVM.py`
  script must contain the `get_evm_engine()` / `get_evm_binpaths()` functions
  that detect and route to the Rust engine.

---

## 2. Build the Rust EVM binaries

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

---

## 3. Install the binaries

### Option A: Copy to a system-wide location

```sh
sudo cp evm/target/release/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins} /usr/local/bin/
```

### Option B: Copy to a user-local bin directory

```sh
mkdir -p ~/bin
cp evm/target/release/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins} ~/bin/
export PATH="$HOME/bin:$PATH"     # add to ~/.bashrc for persistence
```

### Option C: Symlink from the build directory

```sh
mkdir -p ~/bin
for bin in EVidenceModeler evidence_modeler partition_evm_inputs recombine_evm_outputs convert_EVM_outputs_to_GFF3 gff3_file_to_proteins; do
    ln -sf "$(pwd)/evm/target/release/$bin" ~/bin/
done
export PATH="$HOME/bin:$PATH"
```

### Verify installation

```sh
which evidence_modeler
# should print something like: /home/user/bin/evidence_modeler

evidence_modeler --help
# should print the Rust EVM CLI help
```

---

## 4. Configure funannotate to use the Rust EVM engine

funannotate's `runEVM.py` script auto-detects the EVM engine:

1. If the environment variable `FUNANNOTATE_EVM_ENGINE` is set to `rust`
   or `perl`, that engine is used.
2. Otherwise, if `evidence_modeler` is found in `PATH`, the **Rust**
   engine is used.
3. If neither condition is met, funannotate falls back to the **Perl**
   engine (requires `EVM_HOME` to be set).

### Automatic detection (recommended)

Simply ensure the Rust binaries are on your `PATH` before running
funannotate. No environment variables are needed:

```sh
export PATH="$HOME/bin:$PATH"
funannotate predict ...    # will auto-detect and use Rust EVM
```

You should see this in the funannotate log:
```
Using Rust EVM engine (evidence_modeler found in PATH)
```

### Force the Rust engine

If both Perl and Rust EVM are installed, set the environment variable
to force the Rust engine:

```sh
export FUNANNOTATE_EVM_ENGINE=rust
```

### Force the Perl engine (fallback)

If you need to revert to the Perl engine temporarily:

```sh
export FUNANNOTATE_EVM_ENGINE=perl
export EVM_HOME=/path/to/Perl/EVidenceModeler
```

---

## 5. Running funannotate predict with the Rust EVM engine

Once the Rust binaries are on `PATH`, run funannotate predict as usual:

```sh
funannotate predict \
    -i genome.softmasked.fa \
    -o annotation_results \
    -s "Rhodotorula sphaerocarpa" \
    --augustus_species rhodotorula \
    --cpus 16
```

funannotate will internally invoke `funannotate-runEVM.py`, which will
auto-detect the Rust engine and use the Rust `evidence_modeler` binary
for each partition.

### Expected log output

```
[INFO]: Using Rust EVM engine (evidence_modeler found in PATH)
[INFO]: EVM: partitioning input to ~ 35 genes per partition using min 1500 bp interval
[INFO]: Converting to GFF3 and collecting all EVM results
```

---

## 6. Running funannotate update with the Rust EVM engine

The `funannotate update` command also runs EVM (round 2) when
re-annotating with new evidence. The same auto-detection applies:

```sh
funannotate update \
    -i annotation_results \
    -p proteins.gff3 \
    --cpus 16
```

---

## 7. Running the Rust EVM standalone (without funannotate)

The Rust `EVidenceModeler` orchestrator binary can be used directly,
without the funannotate wrapper. This is useful for testing or for
pipelines that don't use funannotate:

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

---

## 8. Troubleshooting

### "Could not find $EVM_HOME and Rust evidence_modeler is not in PATH"

The Rust `evidence_modeler` binary is not on your `PATH`. Verify:

```sh
which evidence_modeler
```

If this prints nothing, install the Rust binaries (see step 3) and
ensure `~/bin` (or wherever you installed them) is on your `PATH`.

### "Rust convert_evm_outputs_to_gff3 binary not found in PATH"

The `convert_EVM_outputs_to_GFF3` binary is missing from your `PATH`.
Rebuild and reinstall the Rust EVM binaries:

```sh
cd EVidenceModeler
make rust
cp evm/target/release/* ~/bin/
```

### funannotate is still using the Perl engine

Check the funannotate log for the engine detection message. If it says
"Using Perl EVM engine", the Rust `evidence_modeler` binary is not on
your `PATH`. Alternatively, set `FUNANNOTATE_EVM_ENGINE=rust` to force
the Rust engine.

### Different gene models between Perl and Rust EVM

The Rust EVM engine achieves 99.4% gene-level parity with the Perl
implementation. The remaining differences are due to:

1. **Partition boundary effects** — genes near partition boundaries may
   be predicted differently depending on the partitioning strategy.
2. **Same-score tie-breaking** — when two paths through the DP trellis
   have the same score, the Perl and Rust implementations may choose
   different paths (182 cases on a real fungal genome).
3. **Small gene model boundary differences** — some genes have slightly
   different 5' or 3' boundaries (different start codon or terminal
   exon selection with equivalent scores).

These differences are algorithmically equivalent and do not represent
bugs. See `PARITY_REPORT.md` for full details.

---

## 9. HPC / SLURM deployment

On UCR HPCC or similar SLURM clusters, the Rust EVM binaries can be
installed as a module or in a shared location:

### Install to a shared project directory

```sh
cd /bigdata/stajichlab/jstajich/projects/EVidenceModeler_rust
make rust
cp evm/target/release/* /bigdata/stajichlab/$USER/bin/
```

### Add to PATH in SLURM scripts

```sh
export PATH="/bigdata/stajichlab/$USER/bin:$PATH"
```

### Or use a module file

Create a module file at `~/.modulefiles/rust-evm/2.1.0`:

```tcl
#%Module1.0
set prefix /bigdata/stajichlab/jstajich/bin
prepend-path PATH $prefix
```

Then in your SLURM script:

```sh
module load rust-evm/2.1.0
```

---

## 10. Uninstalling

To remove the Rust EVM binaries:

```sh
rm ~/bin/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins}
```

Or if installed system-wide:

```sh
sudo rm /usr/local/bin/{EVidenceModeler,evidence_modeler,partition_evm_inputs,recombine_evm_outputs,convert_EVM_outputs_to_GFF3,gff3_file_to_proteins}
```

funannotate will automatically fall back to the Perl engine (if
`EVM_HOME` is set) when the Rust binaries are no longer on `PATH`.
