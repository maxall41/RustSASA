# CHANGELOG

## Version 0.8.0 (Latest update)

* Removed `nalgebra` dependancy.
* Updated `criterion` to version 0.7.0.
* Paper updated
* CLI dependencies gated behind `cli` feature to make library more lightweight. Library dependencies removed: clap, serde, serde_json, quick-xml, indicatif, and mimalloc.
* Improved SASA calculation accuracy by expanding neighbor search radius.
* Improved performance of SpatialGrid data structure (using Verlet lists) to compensate for performance penalty of increased neighbor search radius.

## Version 0.7.0 

* RustSASA now excludes HETATM records (non-standard amino acids and ligands) by default. This can be overridden by passing `--include-hetatoms` (or API equivalent).
* Improved CLI error reporting and handling.
* Added support for customizing number of threads used for parallelization. This can be done via `--threads` (or API equivalent).
* RustSASA now throws an error if an atom is not found in the provided/default radii config. This error can be ignored via `--allow-vdw-fallback`  (or API equivalent).
* Paper updated to include single-threaded RustSASA benchmark among other changes.

## Version 0.6.0

* RustSASA now excludes hydrogens by default and uses ProtOr radii for improved accuracy. Hydrogens can be included by passing `--include-hydrogens` (or API equivalent). If you include hydrogen atoms, you should also provide a custom atomic radii file designed to work with the hydrogens you are including. See README and documentation for more information.
* Improved error handling and reduced code duplication.
* Expose probe radius in cli via `--probe-radius`.

## Version 0.5.0

* Implemented new SIMD SASA calculation kernel and improved spatial grid data structure, thus improving performance by ~35% compared to v0.4.0.
* Fixed benchmarking issue that exaggerated RustSASA's performance (see issue #40); paper and readme have been updated accordingly.
* Added `--n-points` option to CLI to enable customization of the Shrake-Rupley algorithm.  

## Version 0.4.0

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
