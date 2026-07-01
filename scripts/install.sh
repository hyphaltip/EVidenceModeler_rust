#!/usr/bin/env bash
# Install EVidenceModeler (Rust) into a conda/pixi environment.
#
# Usage:
#   scripts/install.sh [--install-prefix PATH]
#
# If --install-prefix is not provided, uses $CONDA_PREFIX if set.
# Otherwise defaults to /opt/evm-rust
#
# This script:
# - Builds all EVidenceModeler binaries using cargo
# - Installs binaries to $INSTALL_PREFIX/bin
# - Is idempotent: safe to run multiple times

set -euo pipefail

# Determine install prefix
INSTALL_PREFIX="${CONDA_PREFIX:-/opt/evm-rust}"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        *)
            echo "Error: unknown option $1" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVM_ROOT="$(dirname "${SCRIPT_DIR}")"
BIN_DIR="${INSTALL_PREFIX}/bin"

# List of EVidenceModeler binaries produced by cargo build
EVM_BINS=(
    EVidenceModeler
    evidence_modeler
    partition_evm_inputs
    recombine_evm_outputs
    convert_EVM_outputs_to_GFF3
    gff3_file_to_proteins
    gff3_gene_prediction_file_validator
    augustus_to_evm_gff3
)

# Check if binaries are already installed
if [ -x "${BIN_DIR}/evidence_modeler" ]; then
    echo "[EVidenceModeler install] Already installed at ${INSTALL_PREFIX}"
    exit 0
fi

echo "[EVidenceModeler install] Installing to ${INSTALL_PREFIX}"
mkdir -p "${BIN_DIR}"

# Build Rust binaries
echo "[EVidenceModeler install] Building EVidenceModeler (Rust)..."
if ! (cd "${EVM_ROOT}" && cargo build --release); then
    echo "[EVidenceModeler install] ERROR: Build failed" >&2
    exit 1
fi

# Copy binaries
echo "[EVidenceModeler install] Installing binaries..."
for bin in "${EVM_BINS[@]}"; do
    src_bin="${EVM_ROOT}/target/release/${bin}"
    if [ -x "${src_bin}" ]; then
        cp "${src_bin}" "${BIN_DIR}/${bin}"
        chmod +x "${BIN_DIR}/${bin}"
    else
        echo "[EVidenceModeler install] WARNING: Binary not found: ${bin}" >&2
    fi
done

echo "[EVidenceModeler install] Installation complete at ${INSTALL_PREFIX}"
exit 0
