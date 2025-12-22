// Copyright (c) 2024 Maxwell Campbell. Licensed under the MIT License.
#[cfg(test)]
mod tests {
    use pdbtbx::{Format, ReadOptions};
    use rust_sasa::structures::atomic::{ChainResult, ProteinResult, ResidueResult};
    use rust_sasa::{SASAResult, sasa_result_to_protein_object};
    use std::io::{BufReader, Cursor};

    #[test]
    fn test_sasa_result_to_protein_object_atom() {
        // Create a simple PDB with 3 atoms
        let pdb_content = r#"ATOM      1  N   ALA A   1      20.154  16.967  25.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030  16.155  25.000  1.00 15.00           C
ATOM      3  C   ALA A   1      17.948  16.712  25.000  1.00 20.00           C
END
"#;
        let mut pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .read_raw(BufReader::new(Cursor::new(pdb_content.as_bytes())))
            .unwrap()
            .0;

        // Create atom-level SASA reResult::Atom(vec![5.0, 10.0, 15.0]);
        let sasa_result = SASAResult::Atom(vec![5.0, 10.0, 15.0]);
        // Apply the result to the PDB
        sasa_result_to_protein_object(&mut pdb, &sasa_result).unwrap();

        // Check that b-factors were set correctly
        let atoms: Vec<_> = pdb.atoms().collect();
        assert_eq!(atoms.len(), 3);
        assert!((atoms[0].b_factor() - 5.0).abs() < 0.001);
        assert!((atoms[1].b_factor() - 10.0).abs() < 0.001);
        assert!((atoms[2].b_factor() - 15.0).abs() < 0.001);
    }

    #[test]
    fn test_sasa_result_to_protein_object_residue() {
        // Create a simple PDB with 2 residues
        let pdb_content = r#"ATOM      1  N   ALA A   1      20.154  16.967  25.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030  16.155  25.000  1.00 15.00           C
ATOM      3  N   GLY A   2      17.948  16.712  25.000  1.00 20.00           N
ATOM      4  CA  GLY A   2      16.500  17.000  25.000  1.00 25.00           C
END
"#;
        let mut pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .read_raw(BufReader::new(Cursor::new(pdb_content.as_bytes())))
            .unwrap()
            .0;

        // Create residue-level SASA result
        let sasa_result = SASAResult::Residue(vec![
            ResidueResult {
                serial_number: 1,
                value: 100.0,
                name: "ALA".to_string(),
                is_polar: false,
                chain_id: "A".to_string(),
            },
            ResidueResult {
                serial_number: 2,
                value: 200.0,
                name: "GLY".to_string(),
                is_polar: false,
                chain_id: "A".to_string(),
            },
        ]);

        // Apply the result to the PDB
        sasa_result_to_protein_object(&mut pdb, &sasa_result).unwrap();

        // Check that b-factors were set correctly for each residue
        let residues: Vec<_> = pdb.residues().collect();
        assert_eq!(residues.len(), 2);

        // First residue atoms should have b-factor 100.0
        for atom in residues[0].atoms() {
            assert!((atom.b_factor() - 100.0).abs() < 0.001);
        }

        // Second residue atoms should have b-factor 200.0
        for atom in residues[1].atoms() {
            assert!((atom.b_factor() - 200.0).abs() < 0.001);
        }
    }

    #[test]
    fn test_sasa_result_to_protein_object_chain() {
        // Create a simple PDB with 2 chains
        let pdb_content = r#"ATOM      1  N   ALA A   1      20.154  16.967  25.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030  16.155  25.000  1.00 15.00           C
ATOM      3  N   GLY B   1      17.948  16.712  25.000  1.00 20.00           N
ATOM      4  CA  GLY B   1      16.500  17.000  25.000  1.00 25.00           C
END
"#;
        let mut pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .read_raw(BufReader::new(Cursor::new(pdb_content.as_bytes())))
            .unwrap()
            .0;

        // Create chain-level SASA result
        let sasa_result = SASAResult::Chain(vec![
            ChainResult {
                name: "A".to_string(),
                value: 300.0,
            },
            ChainResult {
                name: "B".to_string(),
                value: 400.0,
            },
        ]);

        // Apply the result to the PDB
        sasa_result_to_protein_object(&mut pdb, &sasa_result).unwrap();

        // Check that b-factors were set correctly for each chain
        let chains: Vec<_> = pdb.chains().collect();
        assert_eq!(chains.len(), 2);

        // Chain A atoms should have b-factor 300.0
        for atom in chains[0].atoms() {
            assert!((atom.b_factor() - 300.0).abs() < 0.001);
        }

        // Chain B atoms should have b-factor 400.0
        for atom in chains[1].atoms() {
            assert!((atom.b_factor() - 400.0).abs() < 0.001);
        }
    }

    #[test]
    fn test_sasa_result_to_protein_object_protein() {
        // Create a simple PDB with multiple atoms
        let pdb_content = r#"ATOM      1  N   ALA A   1      20.154  16.967  25.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030  16.155  25.000  1.00 15.00           C
ATOM      3  N   GLY A   2      17.948  16.712  25.000  1.00 20.00           N
END
"#;
        let mut pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .read_raw(BufReader::new(Cursor::new(pdb_content.as_bytes())))
            .unwrap()
            .0;

        // Create protein-level SASA result
        let sasa_result = SASAResult::Protein(ProteinResult {
            global_total: 500.0,
            polar_total: 200.0,
            non_polar_total: 300.0,
        });

        // Apply the result to the PDB
        sasa_result_to_protein_object(&mut pdb, &sasa_result).unwrap();

        // Check that all atoms have the same b-factor (global_total)
        for atom in pdb.atoms() {
            assert!((atom.b_factor() - 500.0).abs() < 0.001);
        }
    }

    #[test]
    fn test_hetatm_excluded_gets_zero_sasa_with_multi_chain() {
        use rust_sasa::{ResidueLevel, SASAOptions};

        // PDB with Chain A (residues 1,2) and Chain B (residue 1 = same serial number!)
        // Plus HETATM water molecules in both chains.
        // This tests that:
        // 1. Residue serial numbers don't collide across chains
        // 2. HETATM residues get SASA=0.0 when include_hetatms=false
        let pdb_content = r#"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ATOM      3  N   GLY A   2       3.000   0.000   0.000  1.00 10.00           N
ATOM      4  CA  GLY A   2       4.458   0.000   0.000  1.00 10.00           C
ATOM      5  N   ALA B   1      10.000   0.000   0.000  1.00 10.00           N
ATOM      6  CA  ALA B   1      11.458   0.000   0.000  1.00 10.00           C
HETATM    7  O   HOH A   3       6.000   0.000   0.000  1.00 10.00           O
HETATM    8  O   HOH B   2      13.000   0.000   0.000  1.00 10.00           O
END
"#;
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .read_raw(BufReader::new(Cursor::new(pdb_content.as_bytes())))
            .unwrap();

        // Run with include_hetatms=false (default)
        let result = SASAOptions::<ResidueLevel>::new()
            .with_allow_vdw_fallback(true)
            .process(&pdb)
            .unwrap();

        // Verify: ATOM residues have non-zero SASA, HETATM residues have 0.0
        for res in &result {
            if res.name == "HOH" {
                assert_eq!(
                    res.value, 0.0,
                    "HETATM residue {} in chain {} should have SASA 0.0, got {}",
                    res.serial_number, res.chain_id, res.value
                );
            } else {
                assert!(
                    res.value > 0.0,
                    "ATOM residue {} in chain {} should have non-zero SASA",
                    res.serial_number, res.chain_id
                );
            }
        }

        // Verify no chain A/B value leakage (both chain A res 1 and chain B res 1 exist)
        let chain_a_res1: Vec<_> = result
            .iter()
            .filter(|r| r.chain_id == "A" && r.serial_number == 1)
            .collect();
        let chain_b_res1: Vec<_> = result
            .iter()
            .filter(|r| r.chain_id == "B" && r.serial_number == 1)
            .collect();

        assert_eq!(chain_a_res1.len(), 1, "Should have exactly one A:1 residue");
        assert_eq!(chain_b_res1.len(), 1, "Should have exactly one B:1 residue");
    }
}
