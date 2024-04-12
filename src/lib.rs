//! RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm[1].
mod utils;
mod test;

use std::collections::HashMap;
use nalgebra::{Point3, Vector3};
use pdbtbx::{PDB};
use rayon::prelude::*;
use rstar::{PointDistance, RTree, RTreeObject, AABB};
use std::sync::Arc;
use snafu::{OptionExt};
use snafu::prelude::*;
use crate::utils::{simd_sum};

/// This struct represents an individual Atom
#[derive(Clone)]
pub struct Atom {
    /// The 3D position of the atom
    pub position: Point3<f32>,
    /// The Van Der Walls radius of the atom
    pub radius: f32,
    /// A unique Id for the atom
    pub id: usize,
    /// Parent Id
    pub parent_id: Option<isize>
}

/// Can be used to specify output resolution of SASA computation for convenience.
pub enum SASALevel {
    Atom,
    Residue,
    Chain,
    Protein,
}

#[derive(Debug,PartialEq)]
pub struct ChainResult {
    pub name: String,
    pub value: f32
}

#[derive(Debug,PartialEq)]
pub struct ResidueResult {
    pub serial_number: isize,
    pub value: f32
}

#[derive(Debug,PartialEq)]
pub enum SASAResult {
    Atom(Vec<f32>),
    Residue(Vec<ResidueResult>),
    Chain(Vec<ChainResult>),
    Protein(f32)
}

#[derive(Debug, Snafu)]
pub enum SASACalcError {
    #[snafu(display("Element missing for atom"))]
    ElementMissing,

    #[snafu(display("Van der Waals radius missing for element"))]
    VanDerWaalsMissing,

    #[snafu(display("Failed to map atoms back to level element"))]
    AtomMapToLevelElementFailed,
}

fn serialize_chain_id(s: &str) -> isize {
    let mut result = 0;
    for c in s.chars() {
        if c.is_ascii_alphabetic() {
            let position = c.to_ascii_uppercase() as isize - 64;
            result = result * 10 + position;
        }
    }
    result
}


impl RTreeObject for Atom {
    type Envelope = AABB<[f32; 3]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(<[f32; 3]>::from(self.position))
    }
}

impl PointDistance for Atom {
    fn distance_2(&self, other: &[f32; 3]) -> f32 {
        let xyz = self.position.coords.xyz();
        let z = xyz[2];
        let y = xyz[1];
        let x = xyz[0];
        // No square root as that is required by the package
        (other[2] - z).mul_add(
            other[2] - z,
            (other[1] - y).mul_add(other[1] - y, (other[0] - x).powi(2)),
        )
    }
}

/// Generates points on a sphere using the Golden Section Spiral algorithm
fn generate_sphere_points(n_points: usize) -> Vec<Vector3<f32>> {
    let mut points = Vec::with_capacity(n_points);
    let golden_ratio = (1.0 + 5f32.sqrt()) / 2.0;
    let angle_increment = 2.0 * std::f32::consts::PI * golden_ratio;

    for i in 0..n_points {
        let t = i as f32 / n_points as f32;
        let inclination = (1.0 - 2.0 * t).acos();
        let azimuth = angle_increment * i as f32;

        let x = inclination.sin() * azimuth.cos();
        let y = inclination.sin() * azimuth.sin();
        let z = inclination.cos();

        points.push(Vector3::new(x, y, z));
    }

    points
}

fn is_accessible_rstar(
    test_point: &Point3<f32>,
    atom: &Atom,
    atoms: &RTree<Atom>,
    probe_radius: f32,
    max_radii: f32,
) -> bool {
    let xyz = test_point.coords.xyz();
    let sr = probe_radius + (max_radii * 2.0);
    let candidates = atoms.locate_within_distance([xyz[0], xyz[1], xyz[2]], sr * sr);
    for candidate in candidates {
        if atom.id != candidate.id
            && (test_point - candidate.position).norm() < (candidate.radius + probe_radius)
        {
            return false;
        }
    }
    true
}

/// Takes the probe radius and number of points to use along with a list of Atoms as inputs and returns a Vec with SASA values for each atom.
/// For most users it is recommend that you use `calculate_sasa` instead. This method can be used directly if you do not want to use pdbtbx to load PDB/mmCIF files or want to load them from a different source.
/// Probe Radius Default: 1.4
/// Point Count Default: 100
/// ## Example using pdbtbx:
/// ```
/// use nalgebra::{Point3, Vector3};
/// use pdbtbx::StrictnessLevel;
/// use rust_sasa::{Atom, calculate_sasa_internal};
/// let (mut pdb, _errors) = pdbtbx::open(
///             "./example.cif",
///             StrictnessLevel::Medium
///         ).unwrap();
/// let mut atoms = vec![];
/// for atom in pdb.atoms() {
///     atoms.push(Atom {
///                 position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
///                 radius: atom.element().unwrap().atomic_radius().van_der_waals.unwrap() as f32,
///                 id: atom.serial_number(),
///                 parent_id: None
///     })
///  }
///  let sasa = calculate_sasa_internal(&atoms, None, None);
/// ```
pub fn calculate_sasa_internal(
    atoms: &[Atom],
    in_probe_radius: Option<f32>,
    in_n_points: Option<usize>,
) -> Vec<f32> {
    // Load defaults if not specified
    let mut probe_radius = 1.4;
    let mut n_points = 100;
    if let Some(in_probe_radius) = in_probe_radius {
        probe_radius = in_probe_radius;
    }
    if let Some(in_n_points) = in_n_points {
        n_points = in_n_points;
    }
    //
    let sphere_points = generate_sphere_points(n_points);

    // Create R*-tree from atoms for spatial lookup
    let tree = RTree::bulk_load(atoms.to_vec());
    let tree_arc = Arc::new(tree); // Use Arc for safe sharing among threads
    let mut max_radii = 0.0;
    for atom in atoms {
        if atom.radius > max_radii {
            max_radii = atom.radius;
        }
    }
    atoms
        .par_iter()
        .map(|atom| {
            let mut accessible_points = 0;

            for sphere_point in &sphere_points {
                let test_point = atom.position + sphere_point * (atom.radius + probe_radius);
                if is_accessible_rstar(&test_point, atom, &tree_arc, probe_radius, max_radii) {
                    accessible_points += 1;
                }
            }
            4.0 * std::f32::consts::PI
                * (atom.radius + probe_radius).powi(2)
                * (accessible_points as f32)
                / (n_points as f32)
        })
        .collect()
}

/// This function calculates the SASA for a given protein. The output level can be specified with the level attribute e.g: (SASALevel::Atom,SASALevel::Residue,etc...).
/// Probe radius and n_points can be customized if not customized will default to 1.4, and 100 respectively.
/// If you want more fine-grained control you may want to use [calculate_sasa_internal] instead.
/// ## Example
/// ```
/// use pdbtbx::StrictnessLevel;
/// use rust_sasa::{Atom, calculate_sasa, calculate_sasa_internal, SASALevel};
/// let (mut pdb, _errors) = pdbtbx::open(
///             "./example.cif",
///             StrictnessLevel::Medium
/// ).unwrap();
/// let result = calculate_sasa(&pdb,None,None,SASALevel::Residue);
/// ```
pub fn calculate_sasa(pdb: &PDB,probe_radius: Option<f32>, n_points: Option<usize>,level: SASALevel) -> Result<SASAResult, SASACalcError> {
    let mut atoms = vec![];
    let mut parent_to_atoms = HashMap::new();
    match level {
        SASALevel::Atom | SASALevel::Protein => {
            for atom in pdb.atoms() {
                atoms.push(Atom {
                    position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
                    radius: atom.element().context(ElementMissingSnafu)?.atomic_radius().van_der_waals.context(VanDerWaalsMissingSnafu)? as f32,
                    id: atom.serial_number(),
                    parent_id: None
                })
            }
        }
        SASALevel::Residue => {
            let mut i = 0;
            for residue in pdb.residues() {
                let mut temp = vec![];
                for atom in residue.atoms() {
                    atoms.push(Atom {
                        position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
                        radius: atom.element().context(ElementMissingSnafu)?.atomic_radius().van_der_waals.context(VanDerWaalsMissingSnafu)? as f32,
                        id: atom.serial_number(),
                        parent_id: Some(residue.serial_number())
                    });
                    temp.push(i);
                    i += 1;
                }
                parent_to_atoms.insert(residue.serial_number(),temp);
            }
        }
        SASALevel::Chain => {
            let mut i = 0;
            for chain in pdb.chains() {
                let mut temp = vec![];
                let chain_id = serialize_chain_id(chain.id());
                for atom in chain.atoms() {
                    atoms.push(Atom {
                        position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
                        radius: atom.element().context(ElementMissingSnafu)?.atomic_radius().van_der_waals.context(VanDerWaalsMissingSnafu)? as f32,
                        id: atom.serial_number(),
                        parent_id: Some(chain_id)
                    });
                    temp.push(i);
                    i += 1
                }
                parent_to_atoms.insert(chain_id,temp);
            }
        }
    }
    let atom_sasa = calculate_sasa_internal(&atoms, probe_radius, n_points);
    return match level {
        SASALevel::Atom => {
            Ok(SASAResult::Atom(atom_sasa))
        },
        SASALevel::Chain => {
            let mut chain_sasa = vec![];
            for chain in pdb.chains() {
                let chain_id = serialize_chain_id(chain.id());
                let chain_atom_index = parent_to_atoms.get(&chain_id).context(AtomMapToLevelElementFailedSnafu)?;
                let chain_atoms: Vec<_> = chain_atom_index.iter().map(|&index| atom_sasa[index]).collect();
                let sum = simd_sum(chain_atoms.as_slice());
                chain_sasa.push(ChainResult {
                    name: chain.id().to_string(),
                    value: sum,
                })
            }
            Ok(SASAResult::Chain(chain_sasa))
        },
        SASALevel::Residue => {
            let mut residue_sasa = vec![];
            for residue in pdb.residues() {
                let residue_atom_index = parent_to_atoms.get(&residue.serial_number()).context(AtomMapToLevelElementFailedSnafu)?;
                let residue_atoms: Vec<_> = residue_atom_index.iter().map(|&index| atom_sasa[index]).collect();
                let sum = simd_sum(residue_atoms.as_slice());
                residue_sasa.push(ResidueResult {
                    serial_number: residue.serial_number(),
                    value: sum,
                })
            }
            Ok(SASAResult::Residue(residue_sasa))
        }
        SASALevel::Protein => {
            let sum = simd_sum(atom_sasa.as_slice());
           Ok(SASAResult::Protein(sum))
        }
    }
}