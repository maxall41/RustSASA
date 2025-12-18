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
    const MAX_RMSE: f64 = RMSE_BASELINE * (1.0 + TOLERANCE);

    /// Load FreeSASA reference data and extract chain totals
    fn load_freesasa_chains(
        file_path: &Path,
    ) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(file_path)?;
        let data: FreeSASAOutput = serde_json::from_str(&content)?;

        let mut chain_totals = HashMap::new();
        for result in data.results {
            for structure in result.structure {
                for chain in structure.chains {
                    chain_totals.insert(chain.label, chain.area.total);
                }
            }
        }

        Ok(chain_totals)
    }

    /// Load RustSASA output and aggregate by chain
    fn load_rustsasa_chains(
        file_path: &Path,
    ) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(file_path)?;
        let data: SASAResult = serde_json::from_str(&content)?;

        let mut chain_totals: HashMap<String, f64> = HashMap::new();
        if let SASAResult::Residue(residues) = data {
            for residue in residues {
                *chain_totals.entry(residue.chain_id).or_insert(0.0) += residue.value as f64;
            }
        } else {
            return Err("Expected Residue level SASAResult".into());
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

    #[test]
    fn test_quality_against_freesasa() {
        // Setup directories
        let reference_dir = PathBuf::from("./tests/data/freesasa_reference/");
        let pdb_dir = PathBuf::from("./tests/data/freesasa_pdbs/");
        let output_dir = env::temp_dir().join("rustsasa_quality_test");

        // Clean and create output directory
        if output_dir.exists() {
            fs::remove_dir_all(&output_dir).expect("Failed to clean output directory");
        }
        fs::create_dir_all(&output_dir).expect("Failed to create output directory");

        // Run RustSASA on the PDB directory
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        cmd.arg(pdb_dir.to_str().unwrap())
            .arg(output_dir.to_str().unwrap())
            .arg("--format")
            .arg("json")
            .arg("--output-depth")
            .arg("residue");

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
                load_freesasa_chains(&reference_path),
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

        assert!(
            rmse <= MAX_RMSE,
            "RMSE ({:.2}) exceeds maximum threshold ({:.2})",
            rmse,
            MAX_RMSE
        );
        println!("Quality test passed! RMSE within acceptable threshold.");
    }
}
