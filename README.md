# RustSASA
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/maxall41/RustSASA/rust.yml)
![Crates.io Downloads (recent)](https://img.shields.io/crates/dr/rust-sasa)
![Crates.io License](https://img.shields.io/crates/l/rust-sasa)
![rustc 1.85+](https://img.shields.io/badge/msrv-rustc_1.85+-red.svg)
[![codecov](https://codecov.io/github/maxall41/rustsasa/graph/badge.svg?token=SHM6RRMKSL)](https://codecov.io/github/maxall41/rustsasa)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.09537/status.svg)](https://doi.org/10.21105/joss.09537)

⚡ Ludicrously fast **Rust crate** for protein solvent accessible surface area (SASA) calculations - **63x faster** than Biopython, **5x faster** than FreeSASA. Pure Rust with Python bindings & CLI. Implements Shrake-Rupley algorithm [1].

# Features:
- 🦀 Written in Pure Rust.
- ⚡️ Ludicrously fast. **63X faster** than Biopython, **14X faster** than mdakit_sasa, and **5X faster** than FreeSASA.
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
- [How to cite](#how-to-cite)

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

- RustSasa: *5.237 s ± 0.049 s*

- FreeSASA: *28.042 s ± 2.269 s*

- Biopython: *368.025 s ± 51.156 s*

## Methodology:

We computed residue level SASA values for the entire AlphaFold E. coli proteome structure database using RustSASA, FreeSASA, and Biopython. Benchmarks were run with Hyperfine with options: --warmup 3 --runs 3. All three methods ran across 8 cores on an Apple M3 Macbook with 24GB of unified memory. The RustSASA CLI was used to take advantage of profile guided optimization. GNU Parallel was used to run FreeSASA and Biopython in parallel.


# Validation against FreeSASA

![Comparing FreeSASA and RustSasa](https://github.com/maxall41/RustSASA/blob/main/imgs/sasa_chain_comparison_combined.svg)

# Other

## How to cite

If you use the RustSASA library in your publication please cite our paper in The Journal of Open Source Software. To cite RustSASA scroll up to the top of this page, and then click on the "Cite this repository" button in the right hand GitHub side bar. This will give you a citation in your desired format (e.g., BiBTeX, APA).

## License
MIT

## Building from source

First, make sure you have the Rust compiler installed. See https://rust-lang.org/tools/install/ for installation instructions.

To build RustSASA from source start by initializing git submodules with the following command:
```
git submodule update --init
```
Then build the binary with:
```
cargo build --release
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests and open issues. As this is an actively developed library, I encourage sharing your thoughts, ideas, suggestions, and feedback.

## Hydrogen & Non-standard amino acid handling 

By default, RustSASA strips hydrogen atoms and uses the ProtOr radii config. If you want to include hydrogens, you can use the CLI argument `--include-hydrogens`. If you do so, you should provide your own atom radii config designed to work with hydrogens. Custom radii configs can be provided with `--radii-file`. `--radii-file` accepts a FreeSASA style `.config` file see configs [here](https://github.com/mittinatten/freesasa/tree/master/share).

Additionally, RustSASA filters HETATM records by default. If you want to include HETATM records, you can use the CLI argument `--include-hetatoms` or the API equivalent (see docs). If you include HETATMs you will need to provide a custom radii file that specifies radii for the non-standard amino acids/ligands in your input. 

# Citations:

1: Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol. 1973 Sep 15;79(2):351-71. doi: 10.1016/0022-2836(73)90011-9. PMID: 4760134.
