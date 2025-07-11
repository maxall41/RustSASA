# RustSASA
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/maxall41/RustSASA/rust.yml)
![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/rust-sasa)
![Crates.io License](https://img.shields.io/crates/l/rust-sasa)

RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm[1]. It can be used in Rust and Python!

## Features:
- 🦀 Written in Pure Rust
- ⚡️ Ludicrously fast. *<>X* Faster than Biopython and *<>%* faster than Freesasa
- 🧪 Full test coverage
- 🐍 Python support
- 🤖 Command line interface

## Using in Rust 🦀

```rust
use pdbtbx::StrictnessLevel;
use rust_sasa::{Atom, calculate_sasa, calculate_sasa_internal, SASALevel};
let (mut pdb, _errors) = pdbtbx::open("./example.cif").unwrap();
let result = calculate_sasa(&pdb,None,None,SASALevel::Residue);
```
Full documentation can be found [here](https://docs.rs/rust-sasa/latest/rust_sasa/)

## Using in Python 🐍

You can now utilize RustSasa within Python to speed up your scripts! Take a look at [rust-sasa-python](https://github.com/maxall41/rust-sasa-python)!

Installation:
```
pip install rust-sasa-python
```
Example:
```python
from rust_sasa_python import calculate_sasa_at_residue_level
# Also supports mmCIF files!
residue_sasa_values = calculate_sasa_at_residue_level("path_to_pdb_file.pdb")
```
See full docs [here](https://github.com/maxall41/rust-sasa-python/blob/main/DOCS.md)

## Using CLI

```shell
# Single file
# Also supports .xml, .pdb, and .cif!
rust-sasa path_to_pdb_file.pdb output.json
# Entire directory
rust-sasa input_directory/ output_directory/ --format json
```

## Benchmarking


Methodology:

We computed residue level SASA values for the entire AlphaFold E. coli proteome structure database using RustSASA, Freesasa, and Biopython. Benchmarks were run with Hyperfine with options: --warmup 3 --runs 3. All three methods ran across 8 cores on an Apple M3 Macbook.

- RustSasa: 11.498 s ±  0.340 s

- Freesasa: 54.914 s ±  0.455 s

- Biopython:

## Citations:
1: Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol. 1973 Sep 15;79(2):351-71. doi: 10.1016/0022-2836(73)90011-9. PMID: 4760134.
