use pdbtbx::PDB;
use quick_xml::SeError as XmlError;
use serde_json::Error as JsonError;

use crate::structures::atomic::SASAResult;

pub fn sasa_result_to_json(result: &SASAResult) -> Result<String, JsonError> {
    serde_json::to_string(result)
}

pub fn sasa_result_to_xml(result: &SASAResult) -> Result<String, XmlError> {
    quick_xml::se::to_string(result)
}

pub fn sasa_result_to_protein_object(
    original_pdb: &mut PDB,
    result: &SASAResult,
) -> Result<(), String> {
    match result {
        SASAResult::Atom(v) => {
            for (i, atom) in original_pdb.atoms_mut().enumerate() {
                let item = v[i];
                atom.set_b_factor(item as f64)?;
            }
        }
        SASAResult::Residue(v) => {
            for (i, residue) in original_pdb.residues_mut().enumerate() {
                let item = &v[i];
                assert!(residue.serial_number() == item.serial_number);
                for atom in residue.atoms_mut() {
                    atom.set_b_factor(item.value as f64)?;
                }
            }
        }
        SASAResult::Chain(v) => {
            for (i, chain) in original_pdb.chains_mut().enumerate() {
                let id = chain.id().to_string();
                for residue in chain.residues_mut() {
                    for atom in residue.atoms_mut() {
                        assert!(v[i].name == id);
                        atom.set_b_factor(v[i].value as f64)?;
                    }
                }
            }
        }
        SASAResult::Protein(v) => {
            for residue in original_pdb.residues_mut() {
                for atom in residue.atoms_mut() {
                    atom.set_b_factor(v.global_total as f64)?;
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structures::atomic::{ChainResult, ProteinResult, ResidueResult};
    use pdbtbx::{Format, ReadOptions};
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

        // Create atom-level SASA result
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
}
