# FunAnnotate EvmUtils Integration Analysis

**Date:** 2026-06-30  
**Analysis:** All places where EvmUtils Perl scripts are used in funannotate-live, and mapping to Rust toolset

---

## Executive Summary

**Good news:** The three core EVM pipeline scripts have already been ported to Rust and are auto-detected by funannotate.

**To-do:** Two critical script categories remain Perl-only and block full Rust-based deployment:
1. **GFF3 validation** — `gff3_gene_prediction_file_validator.pl`
2. **Format converters** — Augustus/SNAP/etc. converters in `EvmUtils/misc/`

---

## Part 1: Core EVM Pipeline (✅ COMPLETE)

These are the three scripts that drive EVM execution. **All have been ported to Rust** and are fully integrated into funannotate via auto-detection in `funannotate-runEVM.py`.

### 1.1 `evidence_modeler.pl` → `evidence_modeler` (Rust binary)

| Aspect | Details |
|--------|---------|
| **Purpose** | Core EVM scoring engine for a single genome partition |
| **Perl location** | `EvmUtils/evidence_modeler.pl` (main) or `evidence_modeler.pl` (root) |
| **Rust location** | `evm-utils/src/bin/evidence_modeler.rs` |
| **Rust binary** | `evidence_modeler` |
| **Used in** | `funannotate/aux_scripts/funannotate-runEVM.py::build_partition_evm_cmd()` |
| **Integration** | Auto-detected by `get_evm_engine()` and `get_evm_binpaths()` |
| **Status** | ✅ **FULLY WORKING** |

**Command mapping:**

```
Perl:  perl evidence_modeler.pl -G genome.fa -g genes.gff3 -w weights.txt --min_intron_length 10 --exec_dir ./
Rust:  evidence_modeler -G genome.fa -g genes.gff3 -w weights.txt --min_intron_length 10
```

---

### 1.2 `recombine_EVM_partial_outputs.pl` → `recombine_evm_outputs` (Rust binary)

| Aspect | Details |
|--------|---------|
| **Purpose** | Merge EVM outputs from overlapping partitions into a single genome-wide result |
| **Perl location** | `EvmUtils/recombine_EVM_partial_outputs.pl` |
| **Rust location** | `evm-utils/src/bin/recombine_evm_outputs.rs` |
| **Rust binary** | `recombine_evm_outputs` |
| **Used in** | `funannotate/aux_scripts/funannotate-runEVM.py::build_recombine_cmd()` |
| **Integration** | Auto-detected by `get_evm_binpaths()` |
| **Status** | ✅ **FULLY WORKING** |

**Command mapping:**

```
Perl:  perl recombine_EVM_partial_outputs.pl --partitions partitions.txt --output_file_name evm.out
Rust:  recombine_evm_outputs --partitions partitions.txt -O evm.out
```

---

### 1.3 `convert_EVM_outputs_to_GFF3.pl` → `convert_EVM_outputs_to_GFF3` (Rust binary)

| Aspect | Details |
|--------|---------|
| **Purpose** | Convert raw EVM output format to GFF3 gene annotations |
| **Perl location** | `EvmUtils/convert_EVM_outputs_to_GFF3.pl` |
| **Rust location** | `evm-utils/src/bin/convert_EVM_outputs_to_GFF3.rs` |
| **Rust binary** | `convert_EVM_outputs_to_GFF3` |
| **Used in** | `funannotate/aux_scripts/funannotate-runEVM.py::build_convert_cmd()` |
| **Integration** | Auto-detected by `get_evm_binpaths()` |
| **Status** | ✅ **FULLY WORKING** (with check for binary existence at line 87) |

**Command mapping:**

```
Perl:  perl convert_EVM_outputs_to_GFF3.pl --partitions partitions.txt --output evm.out --genome genome.fa
Rust:  convert_EVM_outputs_to_GFF3 --partitions partitions.txt -O evm.out
```

---

## Part 2: Utility Scripts (✅ Ported to Rust, ⚠️ Partially Integrated)

These scripts preprocess or transform data. Rust versions exist but integration varies.

### 2.1 `gff3_file_to_proteins.pl` → `gff3_file_to_proteins` (Rust binary)

| Aspect | Details |
|--------|---------|
| **Purpose** | Translate CDS features in GFF3 to protein FASTA sequences |
| **Perl location** | `EvmUtils/gff3_file_to_proteins.pl` (EvmUtils version) |
| **Actual usage** | funannotate uses **PASA version** (`PASA/misc_utilities/gff3_file_to_proteins.pl`) |
| **Rust location** | `evm-utils/src/bin/gff3_file_to_proteins.rs` |
| **Rust binary** | `gff3_file_to_proteins` |
| **Used in** | `funannotate/train.py` and `funannotate/update.py` |
| **Integration** | ⚠️ **Not integrated** — still uses PASA version as fallback |
| **Status** | ✅ **Rust version exists but not wired into funannotate** |

**Note:** funannotate preferentially uses the PASA version of this script. The EvmUtils version exists but isn't called.

**Recommendation:** Rust `gff3_file_to_proteins` could replace PASA version for better performance, but this is lower priority since it's not a blocker.

---

### 2.2 `gff_range_retriever.pl` (DEPRECATED - Replaced with Python)

| Aspect | Details |
|--------|---------|
| **Purpose** | Extract GFF3 features within a genomic coordinate range |
| **Perl location** | `EvmUtils/gff_range_retriever.pl` |
| **Status** | ✅ **Replaced with Python implementation** |
| **Replacement** | `RangeFinder()` function in `funannotate-runEVM.py` (lines 466–492) |
| **Why replaced** | Python implementation is simpler, faster (no subprocess), and avoids Perl dependency |

The Perl version is no longer called; it's commented out in the source. This is an example of good migration strategy.

---

### 2.3 `partition_EVM_inputs.pl` (DEPRECATED - Replaced with Python)

| Aspect | Details |
|--------|---------|
| **Purpose** | Split genome FASTA + GFF3 into overlapping partitions |
| **Perl location** | `EvmUtils/partition_EVM_inputs.pl` |
| **Status** | ✅ **Replaced with Python implementation** |
| **Replacement** | `create_partitions()` function in `funannotate-runEVM.py` (lines 227–463) |
| **Why replaced** | Python implementation provides better control over gene-aware partitioning |
| **Rust version exists** | Yes, `partition_evm_inputs` binary in `evm-utils/` |
| **Note** | funannotate check.py still references `EvmUtils/partition_EVM_inputs.pl` for existence check, but the script is never actually executed |

---

## Part 3: Scripts NOT Yet Ported (❌ BLOCKING RUST-ONLY DEPLOYMENT)

### 3.1 `gff3_gene_prediction_file_validator.pl` (❌ MISSING)

| Aspect | Details |
|--------|---------|
| **Purpose** | Validate GFF3 gene prediction file format (parent-child consistency, feature hierarchy) |
| **Perl location** | `EvmUtils/gff3_gene_prediction_file_validator.pl` (207 lines) |
| **Used in** | `funannotate/library.py::evmGFFvalidate()` — called before EVM execution |
| **Rust version** | ❌ **Does not exist** |
| **Impact** | **BLOCKS** Rust-only deployment (still requires Perl) |
| **Complexity** | **Low** — Pure validation logic, no I/O or formats to convert |

**Usage in funannotate:**
```python
def evmGFFvalidate(input, evmpath, logfile):
    Validator = os.path.join(evmpath, "EvmUtils", "gff3_gene_prediction_file_validator.pl")
    cmd = ["perl", Validator, os.path.realpath(input)]
    # ... executes Perl script
```

**Priority:** 🔴 **HIGH** — Needs Rust port before full Rust-only deployment

---

### 3.2 Format Converters in `EvmUtils/misc/` (❌ MISSING)

#### Used by funannotate/predict.py:

| Script | Purpose | Used for | Status |
|--------|---------|----------|--------|
| `augustus_GFF3_to_EVM_GFF3.pl` | Convert Augustus GFF3 output to EVM-compatible format | Gene predictions from Augustus GFF3 mode | ❌ **Not ported** |
| `augustus_GTF_to_EVM_GFF3.pl` | Convert Augustus GTF output to EVM-compatible format | Gene predictions from Augustus GTF mode | ❌ **Not ported** |

**Usage locations:**
- Line 1297–1298: Converter variables defined
- Line 1704, 1781, 1824, 1899, 2052, 2289: Called with `["perl", Converter, ...]`

**Complexity:** **Medium** — Format understanding required, but mostly string parsing/reformatting

**Alternatives to porting:**
1. **Port to Rust:** Full Rust binary for each converter
2. **Implement in Python:** Add Python module to funannotate itself (simpler, no separate binary needed)
3. **Keep Perl as fallback:** Accept that this requires Perl unless user opts into pure Rust

**Priority:** 🟠 **MEDIUM-HIGH** — Blocks Rust-only deployment for Augustus gene predictions

---

#### Other `misc/` scripts (not currently used by funannotate):

Many other converters exist for SNAP, MAKER, Exonerate, etc., but **are not called by funannotate**:
- `SNAP_CDS_to_gff3.pl` — Not used (SNAP output format conversion)
- `braker_GTF_to_EVM_GFF3.pl` — Not used
- 25+ other format converters — Not used by current funannotate

**Recommendation:** Leave these in EvmUtils as reference but don't port unless funannotate needs them.

---

## Part 4: Legacy/Unused Scripts (ℹ️ Not called by funannotate)

The following scripts in `EvmUtils/` are **not referenced anywhere in funannotate**:

| Script | Probable purpose | Status |
|--------|------------------|--------|
| `create_weights_file.pl` | Generate EVM weights file template | Not called |
| `write_EVM_commands.pl` | Generate EVM command list (pre-funannotate) | Not called |
| `execute_EVM_commands.pl` | Execute EVM command list (pre-funannotate) | Not called |
| `summarize_btab_tophits.pl` | BLAST tab summary | Not called |
| `BPbtab.pl` | Protein alignment format conversion | Not called |
| `gff3_file_fix_CDS_phases.pl` | CDS phase correction | Not called |
| `gff3_file_show_phase_based_translations.pl` | Phase-based translation display | Not called |
| `transcript_gff3_to_bed.pl` | GFF3→BED conversion for transcripts | Not called |
| `gene_gff3_to_bed.pl` | GFF3→BED conversion for genes | Not called |
| `extract_complete_proteins.pl` | Extract complete ORF proteins | Not called |
| `EVM_to_GFF3.pl` | EVM output to GFF3 (older version?) | Not called |

**Recommendation:** These are legacy scripts from the original EVidenceModeler pipeline (pre-funannotate). Keep them in `EvmUtils/` as reference, but mark as deprecated. Do not port.

---

## Part 5: Integration Status Summary

| Component | Perl script | Rust binary | Status | funannotate integration |
|-----------|------------|-------------|--------|------------------------|
| **Core EVM** | evidence_modeler.pl | ✅ evidence_modeler | ✅ Ported | ✅ Auto-detected |
| **Recombine** | recombine_EVM_partial_outputs.pl | ✅ recombine_evm_outputs | ✅ Ported | ✅ Auto-detected |
| **Convert** | convert_EVM_outputs_to_GFF3.pl | ✅ convert_EVM_outputs_to_GFF3 | ✅ Ported | ✅ Auto-detected |
| **Proteins** | gff3_file_to_proteins.pl | ✅ gff3_file_to_proteins | ✅ Ported | ⚠️ Uses PASA version |
| **Range finder** | gff_range_retriever.pl | — | ✅ Replaced (Python) | ✅ No Perl call |
| **Partitioner** | partition_EVM_inputs.pl | ✅ partition_evm_inputs | ✅ Replaced (Python) | ✅ No execution |
| **Validator** | gff3_gene_prediction_file_validator.pl | ❌ Missing | ❌ Not ported | ❌ **BLOCKS** |
| **Augustus→EVM** | augustus_GFF3_to_EVM_GFF3.pl | ❌ Missing | ❌ Not ported | ❌ **BLOCKS** |
| **Augustus→EVM** | augustus_GTF_to_EVM_GFF3.pl | ❌ Missing | ❌ Not ported | ❌ **BLOCKS** |

---

## Part 6: Action Items

### PRIORITY 1: Unblock Rust-only deployment 🔴

These must be addressed before claiming "Rust version ready":

#### Task 1.1: Port `gff3_gene_prediction_file_validator.pl` to Rust
- **File:** `evm-utils/src/bin/gff3_gene_prediction_file_validator.rs`
- **Input:** GFF3 file
- **Output:** Stdout/stderr with validation errors, exit code
- **Effort:** ~200 lines Rust code
- **Logic to port:**
  - Parse GFF3 file
  - Build parent-child map from `Parent=` attributes
  - Validate hierarchy (mRNA→gene, exon/CDS→mRNA)
  - Check feature completeness (CDS/exons within mRNA bounds)
  - Report errors with line numbers

#### Task 1.2: Port/replace Augustus format converters
- **Option A:** Port `augustus_GFF3_to_EVM_GFF3.pl` and `augustus_GTF_to_EVM_GFF3.pl` to Rust
  - Create: `evm-utils/src/bin/augustus_gff3_to_evm_gff3.rs`
  - Create: `evm-utils/src/bin/augustus_gtf_to_evm_gff3.rs`
  - Effort: ~100 lines each (format string manipulation)

- **Option B:** Implement in Python in funannotate
  - Create: `funannotate/aux_scripts/augustus_converters.py`
  - More maintainable if funannotate-only feature
  - Effort: ~150 lines Python

**Recommendation:** Option A (Rust) keeps all tools cohesive in evm-utils; Option B (Python) is faster if you want to unblock now.

---

### PRIORITY 2: Optimization (nice-to-have) 🟡

#### Task 2.1: Integrate `gff3_file_to_proteins` Rust version
- Currently: `train.py` and `update.py` use PASA version
- Goal: Use Rust version for better performance
- Work: Update two function calls in funannotate to prefer Rust binary if available
- Effort: ~20 lines Python code (fallback logic)

---

### PRIORITY 3: Cleanup (optional) 🟢

#### Task 3.1: Mark legacy Perl scripts as deprecated
- Add README in `EvmUtils/` documenting which scripts are kept for backward compatibility only
- Move rarely-used converters to `EvmUtils/deprecated/` or mark with `.deprecated` suffix
- No code changes needed, just documentation

#### Task 3.2: Update `funannotate/check.py`
- Remove the check for `EvmUtils/partition_EVM_inputs.pl` (line ~195)
- It's no longer executed; the check is misleading

---

## Part 7: Testing Recommendations

Once Task 1.1 and 1.2 are complete, test with:

```bash
# Test case 1: Run funannotate predict with Rust tools only
export FUNANNOTATE_EVM_ENGINE=rust
funannotate predict -i genome.fa -o output --cpus 8

# Test case 2: Verify all Rust binaries are in PATH
which evidence_modeler recombine_evm_outputs convert_EVM_outputs_to_GFF3 \
      gff3_gene_prediction_file_validator

# Test case 3: Augustus gene prediction (uses converters)
funannotate predict -i genome.fa -o output --augustus_species myspecies --cpus 8
```

---

## Summary

**Current status:** The core EVM pipeline is fully ported to Rust and working. Two script categories remain as Perl-only blockers:

1. **Validator** (high impact, easy port)
2. **Format converters** (medium impact, medium effort)

**Effort to full Rust readiness:** ~8-10 hours of development + testing

**Recommended path:** Port the validator first (quick win), then decide on converters (Rust vs Python).
