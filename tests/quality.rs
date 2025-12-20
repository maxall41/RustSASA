// Copyright (c) 2024 Maxwell Campbell. Licensed under the MIT License.

// Quality regression tests that runs RustSASA over full FreeSASA set and ensures mean RMSE is below cutoff.
mod common;

#[cfg(test)]
mod tests {
    use assert_cmd::{Command, cargo::*};
    use std::collections::HashMap;
    use std::env;
    use std::fs;
    use std::path::{Path, PathBuf};

    use crate::common::io::FreeSASAOutput;
    use rust_sasa::SASAResult;

    const RMSE_BASELINE: f64 = 614.26; //RustSASA RMSE as of v0.8.0
    const TOLERANCE: f64 = 0.01; // 1%

    /// Load FreeSASA reference data and extract chain totals or file-level total
    /// For use_file_total=true, sums all chains and returns a single value keyed by filename
    /// For use_file_total=false, returns individual chain totals keyed by chain label
    fn load_freesasa_chains(
        file_path: &Path,
        use_file_total: bool,
    ) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(file_path)?;
        let data: FreeSASAOutput = serde_json::from_str(&content)?;

        let mut chain_totals = HashMap::new();

        if use_file_total {
            // Sum all chains for file-level comparison (atom/protein depths)
            let mut total = 0.0;
            for result in data.results {
                for structure in result.structure {
                    for chain in structure.chains {
                        total += chain.area.total;
                    }
                }
            }
            let filename = file_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown");
            chain_totals.insert(filename.to_string(), total);
        } else {
            // Extract individual chain totals (residue/chain depths)
            for result in data.results {
                for structure in result.structure {
                    for chain in structure.chains {
                        chain_totals.insert(chain.label, chain.area.total);
                    }
                }
            }
        }

        Ok(chain_totals)
    }

    /// Load RustSASA output and aggregate for comparison (handles all output depths)
    /// Returns either chain-level totals or file-level total depending on output depth
    fn load_rustsasa_chains(
        file_path: &Path,
    ) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(file_path)?;
        let data: SASAResult = serde_json::from_str(&content)?;

        let mut chain_totals: HashMap<String, f64> = HashMap::new();

        match data {
            SASAResult::Atom(atoms) => {
                // Atom level has no chain information, so we sum all atoms
                // and use filename as key for comparison
                let total: f64 = atoms.iter().map(|&v| v as f64).sum();
                let filename = file_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                chain_totals.insert(filename.to_string(), total);
            }
            SASAResult::Residue(residues) => {
                for residue in residues {
                    *chain_totals.entry(residue.chain_id.clone()).or_insert(0.0) +=
                        residue.value as f64;
                }
            }
            SASAResult::Chain(chains) => {
                for chain in chains {
                    chain_totals.insert(chain.name.clone(), chain.value as f64);
                }
            }
            SASAResult::Protein(protein) => {
                // Protein level returns a single total, use filename as key
                let filename = file_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                chain_totals.insert(filename.to_string(), protein.global_total as f64);
            }
        }

        Ok(chain_totals)
    }

    /// Compare two sets of chain data and return matching values
    fn extract_matching_values(
        freesasa: &HashMap<String, f64>,
        rustsasa: &HashMap<String, f64>,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut freesasa_values = Vec::new();
        let mut rustsasa_values = Vec::new();

        for (chain_id, &freesasa_val) in freesasa {
            if let Some(&rustsasa_val) = rustsasa.get(chain_id) {
                freesasa_values.push(freesasa_val);
                rustsasa_values.push(rustsasa_val);
            }
        }

        (freesasa_values, rustsasa_values)
    }

    /// Calculate Root Mean Square Error
    fn calculate_rmse(values1: &[f64], values2: &[f64]) -> f64 {
        assert_eq!(values1.len(), values2.len());

        let sum_squared_diff: f64 = values1
            .iter()
            .zip(values2.iter())
            .map(|(v1, v2)| (v1 - v2).powi(2))
            .sum();

        (sum_squared_diff / values1.len() as f64).sqrt()
    }

    /// Generic quality test that runs RustSASA with specified output depth and validates RMSE
    fn run_quality_test(output_depth: &str) {
        // Determine if we should use file-level totals (for atom/protein) or chain-level (for residue/chain)
        let use_file_total = matches!(output_depth, "atom" | "protein");

        // Setup directories
        let reference_dir = PathBuf::from("./tests/data/freesasa_reference/");
        let pdb_dir = PathBuf::from("./tests/data/freesasa_pdbs/");
        let output_dir = env::temp_dir().join(format!("rustsasa_quality_test_{}", output_depth));

        // Clean and create output directory
        if output_dir.exists() {
            fs::remove_dir_all(&output_dir).expect("Failed to clean output directory");
        }
        fs::create_dir_all(&output_dir).expect("Failed to create output directory");

        // Run RustSASA on the PDB directory with specified output depth
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        cmd.arg(pdb_dir.to_str().unwrap())
            .arg(output_dir.to_str().unwrap())
            .arg("--format")
            .arg("json")
            .arg("--output-depth")
            .arg(output_depth);

        cmd.assert().success();

        // Find matching JSON files
        let reference_files: Vec<_> = fs::read_dir(&reference_dir)
            .expect("Failed to read reference directory")
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("json"))
            .map(|e| e.file_name().to_string_lossy().to_string())
            .collect();

        let output_files: Vec<_> = fs::read_dir(&output_dir)
            .expect("Failed to read output directory")
            .filter_map(|e| e.ok())
            .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("json"))
            .map(|e| e.file_name().to_string_lossy().to_string())
            .collect();

        // Get intersection of files
        let matching_files: Vec<String> = reference_files
            .into_iter()
            .filter(|f| output_files.contains(f))
            .collect();

        assert!(
            !matching_files.is_empty(),
            "No matching JSON files found between reference and output directories"
        );

        // Collect all matching chain values
        let mut all_freesasa_values = Vec::new();
        let mut all_rustsasa_values = Vec::new();

        for filename in &matching_files {
            let reference_path = reference_dir.join(filename);
            let output_path = output_dir.join(filename);

            match (
                load_freesasa_chains(&reference_path, use_file_total),
                load_rustsasa_chains(&output_path),
            ) {
                (Ok(freesasa_chains), Ok(rustsasa_chains)) => {
                    let (freesasa_vals, rustsasa_vals) =
                        extract_matching_values(&freesasa_chains, &rustsasa_chains);

                    all_freesasa_values.extend(freesasa_vals);
                    all_rustsasa_values.extend(rustsasa_vals);
                }
                (Err(e), _) => {
                    eprintln!("Error loading FreeSASA data from {}: {}", filename, e);
                }
                (_, Err(e)) => {
                    eprintln!("Error loading RustSASA data from {}: {}", filename, e);
                }
            }
        }

        assert!(
            !all_freesasa_values.is_empty(),
            "No matching chains found across all files"
        );

        // Calculate RMSE
        let rmse = calculate_rmse(&all_freesasa_values, &all_rustsasa_values);
        let max_rmse = RMSE_BASELINE * (1.0 + TOLERANCE);

        assert!(
            rmse <= max_rmse,
            "[{}] RMSE ({:.2}) exceeds maximum threshold ({:.2})",
            output_depth,
            rmse,
            max_rmse
        );
        println!(
            "[{}] Quality test passed! RMSE: {:.2} (threshold: {:.2})",
            output_depth, rmse, max_rmse
        );
    }

    #[test]
    fn test_quality_atom_depth() {
        run_quality_test("atom");
    }

    #[test]
    fn test_quality_residue_depth() {
        run_quality_test("residue");
    }

    #[test]
    fn test_quality_chain_depth() {
        run_quality_test("chain");
    }

    #[test]
    fn test_quality_protein_depth() {
        run_quality_test("protein");
    }
}
