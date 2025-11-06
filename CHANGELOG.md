# CHANGELOG

## Version 0.4.0 (Latest Update)

- Replace custom SIMD implementation with `pulp` theoretically improving performance on x86 CPUs with AVX SIMD extensions.
- Better support for PDB files with bad records (i.e: old space group like `H 3` or bad `SEQADV` record).
- Update package dependencies.

## Version 0.3.1

- ⚡️ Slightly faster due to memory allocation optimization
- PGO Builds

## Version 0.3.0

- 🔥 ~12X Faster
- 🤖 Implemented Command-line interface (CLI)!
- 🧪 Better tests
- 🆙 Upgraded packages
- 🦀 Upgraded to Rust Edition 2024
- ⛓️‍💥 BREAKING: New SASAOptions interface replaces calculate_sasa!

## Version 0.2.0 & 0.2.1
- RustSASA now includes convenience methods for working directly with pdbtbx
- RustSASA now supports specifying return resolution (e.g: Protein, Chain, Residue, Atom)
- Minor restructuring

## VERSION 0.2.2

- RustSASA now computes the total polar & non-polar values when outputing a protein result
- `is_polar` boolean now included in RustSASA Residue struct.

## VERSION 0.2.3

???

## VERSION 0.2.4

- Updated to pdbtbx 0.12.0
