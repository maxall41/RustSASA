use crate::SASAResult;
use pdbtbx::PDB;
use quick_xml::SeError as XmlError;
use serde_json::Error as JsonError;

pub fn sasa_result_to_json(result: &SASAResult) -> Result<String, JsonError> {
    Ok(serde_json::to_string(result)?)
}

pub fn sasa_result_to_xml(result: &SASAResult) -> Result<String, XmlError> {
    Ok(quick_xml::se::to_string(result)?)
}

pub fn sasa_result_to_protein_object(
    original_pdb: &mut PDB,
    result: &SASAResult,
) -> Result<(), String> {
    match result {
        SASAResult::Atom(v) => {
            let mut i = 0;
            for atom in original_pdb.atoms_mut() {
                let item = v[i];
                atom.set_b_factor(item as f64)?;
                i += 1;
            }
        }
        SASAResult::Residue(v) => {
            let mut i = 0;
            for residue in original_pdb.residues_mut() {
                let item = &v[i];
                assert!(residue.serial_number() == item.serial_number);
                for atom in residue.atoms_mut() {
                    atom.set_b_factor(item.value as f64)?;
                }
                i += 1;
            }
        }
        SASAResult::Chain(v) => {
            let mut i = 0;
            for chain in original_pdb.chains_mut() {
                let id = chain.id().to_string();
                for residue in chain.residues_mut() {
                    for atom in residue.atoms_mut() {
                        assert!(v[i].name == id);
                        atom.set_b_factor(v[i].value as f64)?;
                    }
                }
                i += 1;
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
