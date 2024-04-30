# RustSASA
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/maxall41/RustSASA/rust.yml)
![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/rust-sasa)
![Crates.io License](https://img.shields.io/crates/l/rust-sasa)

RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm[1]. It can be used in Rust and Python!

## Features:
- 🦀 Written in Pure Rust
- ⚡️ 3X Faster than Biopython and ~120% faster than Freesasa
- 🧪 Full test coverage

## Using in Python 🐍

You can now utilize RustSasa within Python to speed up your scripts! Take a look at [rust-sasa-python](https://github.com/maxall41/rust-sasa-python)!

Installation:
```
pip install rust-sasa-python
```
Example:
```python
from rust_sasa_python import calculate_sasa_at_residue_level
residue_sasa_values = calculate_sasa_at_residue_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```
See full docs [here](https://github.com/maxall41/rust-sasa-python/blob/main/DOCS.md)

## Using in Rust 🦀

```rust
use pdbtbx::StrictnessLevel;
use rust_sasa::{Atom, calculate_sasa, calculate_sasa_internal, SASALevel};
let (mut pdb, _errors) = pdbtbx::open(
             "./example.cif",
             StrictnessLevel::Medium
).unwrap();
let result = calculate_sasa(&pdb,None,None,SASALevel::Residue);
```
Full documentation can be found [here](https://docs.rs/rust-sasa/latest/rust_sasa/)

## Benchmarking
Benchmarks were performed on an M2 Apple Mac with 8GB of RAM and 8 Cores with the protein AF-A0A2K5XT84-F1 (AlphaFold).

- Biopython: ~150ms

- Freesasa: ~90ms

- RustSASA: ~40ms

## Citations:
1: Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol. 1973 Sep 15;79(2):351-71. doi: 10.1016/0022-2836(73)90011-9. PMID: 4760134.
