# RustSASA
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/maxall41/RustSASA/rust.yml)
![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/rust-sasa)
![Crates.io License](https://img.shields.io/crates/l/rust-sasa)

RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm[1].
## Features:
- 🦀 Written in Pure Rust
- ⚡️ 3X Faster than Biopython
- 🧪 Full test coverage
  
## RustSASA Implementation vs Biopython Implementation
Benchmarks were performed on an M2 Apple Mac with 8GB of RAM and 8 Cores with the protein AF-A0A2K5XT84-F1 (AlphaFold).

Biopython: ~150ms

RustSasa: ~50ms

## Example Usage with `pdbtbx`:
```rust
use nalgebra::{Point3, Vector3};
use pdbtbx::StrictnessLevel;
use rust_sasa::{Atom, calculate_sasa};
let (mut pdb, _errors) = pdbtbx::open(
          "./example.cif",
          StrictnessLevel::Medium
).unwrap();
let mut atoms = vec![];
for atom in pdb.atoms() {
  atoms.push(Atom {
              position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
              radius: atom.element().unwrap().atomic_radius().van_der_waals.unwrap() as f32,
              id: atom.serial_number()
  })
}
let sasa = calculate_sasa(&atoms, None, None);
```

## Documentation:
See https://docs.rs/rust-sasa/latest/rust_sasa/

## Citations:
1: Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol. 1973 Sep 15;79(2):351-71. doi: 10.1016/0022-2836(73)90011-9. PMID: 4760134.
