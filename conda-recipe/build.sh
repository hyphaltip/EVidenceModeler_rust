#!/usr/bin/env bash
set -euo pipefail

# Use the install.sh script to build and install
mkdir -p "$PREFIX/opt/evm-rust"
CONDA_PREFIX="$PREFIX" "$SRC_DIR/scripts/install.sh" --install-prefix "$PREFIX/opt/evm-rust"

# Move binaries to $PREFIX/bin for easy access
mkdir -p "$PREFIX/bin"
cp "$PREFIX/opt/evm-rust/bin"/* "$PREFIX/bin/"

echo "EVidenceModeler installed to $PREFIX"
ls -la "$PREFIX/bin" | grep evidence_modeler
