# EVidenceModeler Rust - Conda Recipe

This directory contains the conda recipe for building and distributing the Rust-optimized EVidenceModeler package.

## Building the Package

### Prerequisites

- conda-build (or mamba build)
- Git
- Rust toolchain (will be installed by conda)

### Build Instructions

```bash
# Option 1: Build locally from this directory
cd EVidenceModeler_rust/conda-recipe
conda build . -c conda-forge

# Option 2: Build using mamba (faster)
mamba build . -c conda-forge

# Option 3: Build with specific output directory
conda build . -c conda-forge --croot /path/to/build/output
```

### Testing the Build

```bash
# Test the built package
conda install --use-local evidencemodeler-rust

# Verify installation
evidence_modeler --version
```

## Publishing to Conda-Forge

To publish this package to [bioconda](https://bioconda.github.io/), follow the bioconda contribution guide:

1. Fork the bioconda-recipes repository
2. Create a new recipe directory: `recipes/evidencemodeler-rust/`
3. Copy this recipe's contents
4. Submit a pull request

### Example bioconda recipe structure:
```
recipes/evidencemodeler-rust/
├── meta.yaml
├── build.sh
└── LICENSE.txt
```

## Environment Variables

When installed, this package sets no special environment variables. Binaries are available directly from `$PATH`:

- `evidence_modeler` - Main EVM binary
- `EVidenceModeler` - Alias for evidence_modeler
- `partition_evm_inputs` - Utility for partitioning inputs
- `recombine_evm_outputs` - Utility for combining outputs
- `convert_EVM_outputs_to_GFF3` - Output format converter
- `gff3_file_to_proteins` - Protein extraction tool
- `gff3_gene_prediction_file_validator` - Validation utility
- `augustus_to_evm_gff3` - Augustus format converter

## Notes

- This recipe clones from GitHub during build, so internet access is required
- Build time: ~10-20 minutes depending on system
- The package uses the `install.sh` script from the EVidenceModeler_rust repository
- Architecture: Linux 64-bit only (Windows and macOS support would require additional work)
