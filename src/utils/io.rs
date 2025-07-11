use crate::SASAResult;
use pdbtbx::PDB;
use quick_xml::SeError as XmlError;
use serde_json::Error as JsonError;

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
