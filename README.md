# RustSASA
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/maxall41/RustSASA/rust.yml)
![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/rust-sasa)
![Crates.io License](https://img.shields.io/crates/l/rust-sasa)
![rustc 1.85+](https://img.shields.io/badge/msrv-rustc_1.85+-red.svg)
[![codecov](https://codecov.io/github/maxall41/rustsasa/graph/badge.svg?token=SHM6RRMKSL)](https://codecov.io/github/maxall41/rustsasa)

⚡ Ludicrously fast **Rust crate** for protein solvent accessible surface area (SASA) calculations - **46x faster** than Biopython, **7x faster** than FreeSASA. Pure Rust with Python bindings & CLI. Implements Shrake-Rupley algorithm [1].

# Features:
- 🦀 Written in Pure Rust.
- ⚡️ Ludicrously fast. **46X faster** than Biopython, **14X faster** than mdakit_sasa, and **7X faster** than Freesasa.
- 🧪 Full test coverage.
- 🐍 Python support.
- 🤖 Command line interface.

# Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
  - [Rust](#using-in-rust-)
  - [Python](#using-in-python-)
  - [CLI](#using-cli-)
  - [MDAnalysis](#using-with-mdanalysis)
- [Benchmarks](#benchmarking)
- [Validation](#validation-against-freesasa)
- [Contributing](#contributing)
- [License](#license)

# Installation

## Rust 🦀

```
cargo add rust-sasa
```

## Python 🐍
```
pip install rust-sasa-python
```

## MDAnalysis package <img src="https://github.com/maxall41/RustSASA/blob/radical/imgs/mdanalysis-logo.png" width="25" height="25">

```
pip install mdsasa-bolt
```

## Command-line interface 🤖

**1. Install Cargo Bin Install**

```
curl -L --proto '=https' --tlsv1.2 -sSf https://raw.githubusercontent.com/cargo-bins/cargo-binstall/main/install-from-binstall-release.sh | bash
```

**2. Install rust-sasa**

```
cargo binstall rust-sasa
```

# Quick start

## Using in Rust 🦀

```rust
use pdbtbx::StrictnessLevel;
use rust_sasa::options::{SASAOptions, ResidueLevel};

let (mut pdb, _errors) = pdbtbx::open("./example.cif").unwrap();
let result = SASAOptions::<ResidueLevel>::new().process(&pdb);

```
Full documentation can be found [here](https://docs.rs/rust-sasa/latest/rust_sasa/).

## Using in Python 🐍

You can now utilize RustSasa within Python to speed up your scripts! Take a look at [rust-sasa-python](https://github.com/maxall41/rust-sasa-python)!

```python
import rust_sasa_python as sasa

# Simple calculation - use convenience function
result = sasa.calculate_protein_sasa("protein.pdb")
print(f"Total SASA: {result.total:.2f}")
```

See full docs [here](https://github.com/maxall41/rust-sasa-python/blob/main/DOCS.md).

## Using CLI 🤖

**Processing single file**

```
rust-sasa path_to_pdb_file.pdb output.json # Also supports .xml, .pdb, and .cif!
```

**Processing an entire directory**

```
rust-sasa input_directory/ output_directory/ --format json # Also supports .xml, .pdb, and .cif!
```

## Using with MDAnalysis <img src="https://github.com/maxall41/RustSASA/blob/radical/imgs/mdanalysis-logo.png" width="25" height="25">

RustSASA can be used with MDAnalysis to calculate SASA for a protein in a trajectory. RustSASA is **17x faster** than mdakit_sasa.

```python
import MDAnalysis as mda
from mdsasa_bolt import SASAAnalysis

# Load your trajectory
u = mda.Universe("topology.pdb", "trajectory.dcd")

# Create SASA analysis
sasa_analysis = SASAAnalysis(u, select="protein")

# Run the analysis
sasa_analysis.run()

# Access results
print(f"Mean total SASA: {sasa_analysis.results.mean_total_area:.2f} Ų")
print(f"SASA per frame: {sasa_analysis.results.total_area}")
print(f"SASA per residue: {sasa_analysis.results.residue_area}")
```

See the [mdsasa-bolt](https://github.com/maxall41/mdsasa-bolt) package for more information.

# Benchmarking

## Results:

- RustSasa: *8.071 s ±  0.361 s*

- Freesasa: *54.914 s ±  0.455 s*

- Biopython: *368.025 s ± 51.156 s*

## Methodology:

We computed residue level SASA values for the entire AlphaFold E. coli proteome structure database using RustSASA, Freesasa, and Biopython. Benchmarks were run with Hyperfine with options: --warmup 3 --runs 3. All three methods ran across 8 cores on an Apple M3 Macbook with 24GB of unified memory. The RustSASA CLI was used to take advantage of profile guided optimization. GNU Parallel was used to run Freesasa and Biopython in parallel.


# Validation against Freesasa

![Comparing Freesasa and RustSasa on E. coli proteome](https://github.com/maxall41/RustSASA/blob/main/imgs/sasa_chain_comparison_E_coli.svg)


![Comparing Freesasa and RustSasa on Freesasa comparison dataset](https://github.com/maxall41/RustSASA/blob/main/imgs/sasa_chain_comparison_freesasa_ds.svg)

# Other

## License
MIT

## Latest update (0.3.1)

- ⚡️ Slightly faster due to memory allocation optimization
- PGO Builds

Also see [changelog](https://github.com/maxall41/rustsasa/blob/master/CHANGELOG.md).

## Contributing

Contributions are welcome! Please feel free to submit pull requests and open issues. As this is an actively developed library, I encourage sharing your thoughts, ideas, suggestions, and feedback.

## 🗺️ Roadmap

- [ ] Use mimalloc V3: Once out of beta, if it provides better performance.
- [ ] Automated MacOS PGO builds: See https://github.com/Kobzol/cargo-pgo/issues/68.
- [ ] Automated BOLT for Linux builds.
- [ ] Use minimal parser.
- [ ] Implement PowerSasa algorithm. Faster?

# Citations:

1: Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol. 1973 Sep 15;79(2):351-71. doi: 10.1016/0022-2836(73)90011-9. PMID: 4760134.
