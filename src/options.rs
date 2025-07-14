use crate::structures::atomic::{ChainResult, ProteinResult, ResidueResult};
use crate::utils::consts::POLAR_AMINO_ACIDS;
use crate::utils::{serialize_chain_id, simd_sum};
use crate::{Atom, calculate_sasa_internal};
use nalgebra::Point3;
use pdbtbx::PDB;
use snafu::OptionExt;
use snafu::prelude::*;
use std::collections::HashMap;
use std::marker::PhantomData;

/// Options for configuring SASA (Solvent Accessible Surface Area) calculations.
///
/// This struct provides configuration options for SASA calculations at different levels
/// of granularity (atom, residue, chain, or protein level). The type parameter `T`
/// determines the output type and processing behavior.
///
/// # Type Parameters
///
/// * `T` - The processing level, which must implement [`SASAProcessor`]. Available levels:
///   - [`AtomLevel`] - Returns SASA values for individual atoms
///   - [`ResidueLevel`] - Returns SASA values aggregated by residue
///   - [`ChainLevel`] - Returns SASA values aggregated by chain
///   - [`ProteinLevel`] - Returns SASA values aggregated for the entire protein
///
/// # Fields
///
/// * `probe_radius` - Radius of the solvent probe sphere in Angstroms (default: 1.4)
/// * `n_points` - Number of points on the sphere surface for sampling (default: 100)
/// * `parallel` - Whether to use parallel processing (default: true)
///
/// # Examples
///
/// ```rust
/// use rust_sasa::options::{SASAOptions, ResidueLevel};
/// use pdbtbx::PDB;
///
/// // Create options with default settings
/// let options = SASAOptions::<ResidueLevel>::new();
///
/// // Customize the configuration
/// let custom_options = SASAOptions::<ResidueLevel>::new()
///     .with_probe_radius(1.2)
///     .with_n_points(200)
///     .with_parallel(true);
///
/// // Process a PDB structure
/// # let pdb = PDB::new();
/// let result = custom_options.process(&pdb)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct SASAOptions<T> {
    probe_radius: f32,
    n_points: usize,
    parallel: bool,
    _marker: PhantomData<T>,
}

// Zero-sized marker types for each level
pub struct AtomLevel;
pub struct ResidueLevel;
pub struct ChainLevel;
pub struct ProteinLevel;

pub type AtomsMappingResult = Result<(Vec<Atom>, HashMap<isize, Vec<usize>>), SASACalcError>;

// Trait that defines the processing behavior for each level
pub trait SASAProcessor {
    type Output;

    fn process_atoms(
        atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &HashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError>;

    fn build_atoms_and_mapping(pdb: &PDB) -> AtomsMappingResult;
}

impl SASAProcessor for AtomLevel {
    type Output = Vec<f32>;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        _pdb: &PDB,
        _parent_to_atoms: &HashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError> {
        Ok(atom_sasa.to_vec())
    }

    fn build_atoms_and_mapping(
        pdb: &PDB,
    ) -> Result<(Vec<Atom>, HashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        for atom in pdb.atoms() {
            atoms.push(Atom {
                position: Point3::new(
                    atom.pos().0 as f32,
                    atom.pos().1 as f32,
                    atom.pos().2 as f32,
                ),
                radius: atom
                    .element()
                    .context(ElementMissingSnafu)?
                    .atomic_radius()
                    .van_der_waals
                    .context(VanDerWaalsMissingSnafu)? as f32,
                id: atom.serial_number(),
                parent_id: None,
            })
        }
        Ok((atoms, HashMap::new()))
    }
}

impl SASAProcessor for ResidueLevel {
    type Output = Vec<ResidueResult>;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &HashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError> {
        let mut residue_sasa = vec![];
        for chain in pdb.chains() {
            for residue in chain.residues() {
                let residue_atom_index = parent_to_atoms
                    .get(&residue.serial_number())
                    .context(AtomMapToLevelElementFailedSnafu)?;
                let residue_atoms: Vec<_> = residue_atom_index
                    .iter()
                    .map(|&index| atom_sasa[index])
                    .collect();
                let sum = simd_sum(residue_atoms.as_slice());
                let name = residue
                    .name()
                    .context(FailedToGetResidueNameSnafu)?
                    .to_string();
                residue_sasa.push(ResidueResult {
                    serial_number: residue.serial_number(),
                    value: sum,
                    is_polar: POLAR_AMINO_ACIDS.contains(&name),
                    chain_id: chain.id().to_string(),
                    name,
                })
            }
        }
        Ok(residue_sasa)
    }

    fn build_atoms_and_mapping(
        pdb: &PDB,
    ) -> Result<(Vec<Atom>, HashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = HashMap::new();
        let mut i = 0;
        for residue in pdb.residues() {
            let mut temp = vec![];
            for atom in residue.atoms() {
                atoms.push(Atom {
                    position: Point3::new(
                        atom.pos().0 as f32,
                        atom.pos().1 as f32,
                        atom.pos().2 as f32,
                    ),
                    radius: atom
                        .element()
                        .context(ElementMissingSnafu)?
                        .atomic_radius()
                        .van_der_waals
                        .context(VanDerWaalsMissingSnafu)? as f32,
                    id: atom.serial_number(),
                    parent_id: Some(residue.serial_number()),
                });
                temp.push(i);
                i += 1;
            }
            parent_to_atoms.insert(residue.serial_number(), temp);
        }
        Ok((atoms, parent_to_atoms))
    }
}

impl SASAProcessor for ChainLevel {
    type Output = Vec<ChainResult>;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &HashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError> {
        let mut chain_sasa = vec![];
        for chain in pdb.chains() {
            let chain_id = serialize_chain_id(chain.id());
            let chain_atom_index = parent_to_atoms
                .get(&chain_id)
                .context(AtomMapToLevelElementFailedSnafu)?;
            let chain_atoms: Vec<_> = chain_atom_index
                .iter()
                .map(|&index| atom_sasa[index])
                .collect();
            let sum = simd_sum(chain_atoms.as_slice());
            chain_sasa.push(ChainResult {
                name: chain.id().to_string(),
                value: sum,
            })
        }
        Ok(chain_sasa)
    }

    fn build_atoms_and_mapping(
        pdb: &PDB,
    ) -> Result<(Vec<Atom>, HashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = HashMap::new();
        let mut i = 0;
        for chain in pdb.chains() {
            let mut temp = vec![];
            let chain_id = serialize_chain_id(chain.id());
            for atom in chain.atoms() {
                atoms.push(Atom {
                    position: Point3::new(
                        atom.pos().0 as f32,
                        atom.pos().1 as f32,
                        atom.pos().2 as f32,
                    ),
                    radius: atom
                        .element()
                        .context(ElementMissingSnafu)?
                        .atomic_radius()
                        .van_der_waals
                        .context(VanDerWaalsMissingSnafu)? as f32,
                    id: atom.serial_number(),
                    parent_id: Some(chain_id),
                });
                temp.push(i);
                i += 1
            }
            parent_to_atoms.insert(chain_id, temp);
        }
        Ok((atoms, parent_to_atoms))
    }
}

impl SASAProcessor for ProteinLevel {
    type Output = ProteinResult;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &HashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError> {
        let mut polar_total: f32 = 0.0;
        let mut non_polar_total: f32 = 0.0;
        for residue in pdb.residues() {
            let residue_atom_index = parent_to_atoms
                .get(&residue.serial_number())
                .context(AtomMapToLevelElementFailedSnafu)?;
            let residue_atoms: Vec<_> = residue_atom_index
                .iter()
                .map(|&index| atom_sasa[index])
                .collect();
            let sum = simd_sum(residue_atoms.as_slice());
            let name = residue
                .name()
                .context(FailedToGetResidueNameSnafu)?
                .to_string();
            if POLAR_AMINO_ACIDS.contains(&name) {
                polar_total += sum
            } else {
                non_polar_total += sum
            }
        }
        let global_sum = simd_sum(atom_sasa);
        Ok(ProteinResult {
            global_total: global_sum,
            polar_total,
            non_polar_total,
        })
    }

    fn build_atoms_and_mapping(
        pdb: &PDB,
    ) -> Result<(Vec<Atom>, HashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = HashMap::new();
        let mut i = 0;
        for residue in pdb.residues() {
            let mut temp = vec![];
            for atom in residue.atoms() {
                atoms.push(Atom {
                    position: Point3::new(
                        atom.pos().0 as f32,
                        atom.pos().1 as f32,
                        atom.pos().2 as f32,
                    ),
                    radius: atom
                        .element()
                        .context(ElementMissingSnafu)?
                        .atomic_radius()
                        .van_der_waals
                        .context(VanDerWaalsMissingSnafu)? as f32,
                    id: atom.serial_number(),
                    parent_id: Some(residue.serial_number()),
                });
                temp.push(i);
                i += 1;
            }
            parent_to_atoms.insert(residue.serial_number(), temp);
        }
        Ok((atoms, parent_to_atoms))
    }
}

#[derive(Debug, Snafu)]
pub enum SASACalcError {
    #[snafu(display("Element missing for atom"))]
    ElementMissing,

    #[snafu(display("Van der Waals radius missing for element"))]
    VanDerWaalsMissing,

    #[snafu(display("Failed to map atoms back to level element"))]
    AtomMapToLevelElementFailed,

    #[snafu(display("Failed to get residue name"))]
    FailedToGetResidueName,
}

impl Default for SASAOptions<ResidueLevel> {
    fn default() -> Self {
        Self {
            probe_radius: 1.4, // Standard water probe radius in Angstroms
            n_points: 100,     // Number of points on sphere for sampling
            parallel: true,    // Parallel processing by default
            _marker: PhantomData,
        }
    }
}

impl<T> SASAOptions<T> {
    /// Create a new SASAOptions with the specified level type
    pub fn new() -> SASAOptions<T> {
        SASAOptions {
            probe_radius: 1.4,
            n_points: 100,
            parallel: false,
            _marker: PhantomData,
        }
    }

    /// Set the probe radius (default: 1.4 Angstroms)
    pub fn with_probe_radius(mut self, radius: f32) -> Self {
        self.probe_radius = radius;
        self
    }

    /// Set the number of points on the sphere for sampling (default: 100)
    pub fn with_n_points(mut self, points: usize) -> Self {
        self.n_points = points;
        self
    }

    /// Enable or disable parallel processing (default: true)
    pub fn with_parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }
}

// Convenience constructors for each level
impl SASAOptions<AtomLevel> {
    pub fn atom_level() -> Self {
        Self::new()
    }
}

impl SASAOptions<ResidueLevel> {
    pub fn residue_level() -> Self {
        Self::new()
    }
}

impl SASAOptions<ChainLevel> {
    pub fn chain_level() -> Self {
        Self::new()
    }
}

impl SASAOptions<ProteinLevel> {
    pub fn protein_level() -> Self {
        Self::new()
    }
}

impl<T: SASAProcessor> SASAOptions<T> {
    /// This function calculates the SASA for a given protein. The output type is determined by the level type parameter.
    /// Probe radius and n_points can be customized, defaulting to 1.4 and 100 respectively.
    /// If you want more fine-grained control you may want to use [calculate_sasa_internal] instead.
    /// ## Example
    /// ```
    /// use pdbtbx::StrictnessLevel;
    /// use rust_sasa::options::{SASAOptions, ResidueLevel};
    /// let (mut pdb, _errors) = pdbtbx::open("./example.cif").unwrap();
    /// let result = SASAOptions::<ResidueLevel>::new().process(&pdb);
    /// ```
    pub fn process(&self, pdb: &PDB) -> Result<T::Output, SASACalcError> {
        let (atoms, parent_to_atoms) = T::build_atoms_and_mapping(pdb)?;
        let atom_sasa =
            calculate_sasa_internal(&atoms, self.probe_radius, self.n_points, self.parallel);
        T::process_atoms(&atoms, &atom_sasa, pdb, &parent_to_atoms)
    }
}
