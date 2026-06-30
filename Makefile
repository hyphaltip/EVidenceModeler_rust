all: rust

# --- Rust build (default) ---
rust:
	cargo build --release

rust-debug:
	cargo build

test:
	cargo test --all
	cd testing && bash runMe.rust.sh

large_sample_data:
	git clone https://github.com/EVidenceModeler/EVM_sample_data.git
