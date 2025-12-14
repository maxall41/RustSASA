use crate::structures::atomic::{ChainResult, ProteinResult, ResidueResult};
use crate::utils::consts::{POLAR_AMINO_ACIDS, load_radii_from_file};
use crate::utils::{get_radius, serialize_chain_id, simd_sum};
use crate::{Atom, calculate_sasa_internal};
use fnv::FnvHashMap;
use nalgebra::Point3;
use pdbtbx::PDB;
use snafu::OptionExt;
use snafu::prelude::*;
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
/// * `threads` - Number of threads to use for parallel processing (default: -1 for all cores)
/// * `include_hydrogens` - Whether to include hydrogen atoms in calculations (default: false)
/// * `radii_config` - Optional custom radii configuration (default: uses embedded protor.config)
/// * `allow_vdw_fallback` - Allow fallback to PDBTBX van der Waals radii when radius is not found in radii file (default: false)
/// * `include_hetatms` - Whether to include HETATM records (e.g. non-standard amino acids) in calculations (default: false)
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
///     .with_threads(-1)
///     .with_include_hydrogens(false)
///     .with_allow_vdw_fallback(true)
///     .with_include_hetatms(false);
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
    threads: isize,
    include_hydrogens: bool,
    radii_config: Option<FnvHashMap<String, FnvHashMap<String, f32>>>,
    allow_vdw_fallback: bool,
    include_hetatms: bool,
    _marker: PhantomData<T>,
}

// Zero-sized marker types for each level
pub struct AtomLevel;
pub struct ResidueLevel;
pub struct ChainLevel;
pub struct ProteinLevel;

pub type AtomsMappingResult = Result<(Vec<Atom>, FnvHashMap<isize, Vec<usize>>), SASACalcError>;

/// Macro to reduce duplication in atom building logic
macro_rules! build_atom {
    ($atoms:expr, $atom:expr, $element:expr, $residue_name:expr, $atom_name:expr, $parent_id:expr, $radii_config:expr, $allow_vdw_fallback:expr) => {{
        let radius = match get_radius($residue_name, $atom_name, $radii_config) {
            Some(r) => r,
            None => {
                if $allow_vdw_fallback {
                    $element
                        .atomic_radius()
                        .van_der_waals
                        .context(VanDerWaalsMissingSnafu)? as f32
                } else {
                    return Err(SASACalcError::RadiusMissing {
                        residue_name: $residue_name.to_string(),
                        atom_name: $atom_name.to_string(),
                        element: $element.to_string(),
                    });
                }
            }
        };

        $atoms.push(Atom {
            position: Point3::new(
                $atom.pos().0 as f32,
                $atom.pos().1 as f32,
                $atom.pos().2 as f32,
            ),
            radius,
            id: $atom.serial_number(),
            parent_id: $parent_id,
        });
    }};
}

// Trait that defines the processing behavior for each level
pub trait SASAProcessor {
    type Output;

    fn process_atoms(
        atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &FnvHashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError>;

    fn build_atoms_and_mapping(
        pdb: &PDB,
        radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
        allow_vdw_fallback: bool,
        include_hydrogens: bool,
        include_hetatms: bool,
    ) -> AtomsMappingResult;
}

impl SASAProcessor for AtomLevel {
    type Output = Vec<f32>;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        _pdb: &PDB,
        _parent_to_atoms: &FnvHashMap<isize, Vec<usize>>,
    ) -> Result<Self::Output, SASACalcError> {
        Ok(atom_sasa.to_vec())
    }

    fn build_atoms_and_mapping(
        pdb: &PDB,
        radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
        allow_vdw_fallback: bool,
        include_hydrogens: bool,
        include_hetatms: bool,
    ) -> Result<(Vec<Atom>, FnvHashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        for residue in pdb.residues() {
            let residue_name = residue.name().context(FailedToGetResidueNameSnafu)?;
            for atom in residue.atoms() {
                let element = atom.element().context(ElementMissingSnafu)?;
                let atom_name = atom.name();
                if element == &pdbtbx::Element::H && !include_hydrogens {
                    continue;
                };
                if atom.hetero() && !include_hetatms {
                    continue;
                }
                build_atom!(
                    atoms,
                    atom,
                    element,
                    residue_name,
                    atom_name,
                    None,
                    radii_config,
                    allow_vdw_fallback
                );
            }
        }
        Ok((atoms, FnvHashMap::default()))
    }
}

impl SASAProcessor for ResidueLevel {
    type Output = Vec<ResidueResult>;

    fn process_atoms(
        _atoms: &[Atom],
        atom_sasa: &[f32],
        pdb: &PDB,
        parent_to_atoms: &FnvHashMap<isize, Vec<usize>>,
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
        radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
        allow_vdw_fallback: bool,
        include_hydrogens: bool,
        include_hetatms: bool,
    ) -> Result<(Vec<Atom>, FnvHashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = FnvHashMap::default();
        let mut i = 0;
        for residue in pdb.residues() {
            let residue_name = residue.name().context(FailedToGetResidueNameSnafu)?;
            let mut temp = vec![];
            for atom in residue.atoms() {
                let element = atom.element().context(ElementMissingSnafu)?;
                let atom_name = atom.name();
                if element == &pdbtbx::Element::H && !include_hydrogens {
                    continue;
                };
                if atom.hetero() && !include_hetatms {
                    continue;
                }
                build_atom!(
                    atoms,
                    atom,
                    element,
                    residue_name,
                    atom_name,
                    Some(residue.serial_number()),
                    radii_config,
                    allow_vdw_fallback
                );
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
        parent_to_atoms: &FnvHashMap<isize, Vec<usize>>,
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
        radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
        allow_vdw_fallback: bool,
        include_hydrogens: bool,
        include_hetatms: bool,
    ) -> Result<(Vec<Atom>, FnvHashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = FnvHashMap::default();
        let mut i = 0;
        for chain in pdb.chains() {
            let mut temp = vec![];
            let chain_id = serialize_chain_id(chain.id());
            for residue in chain.residues() {
                let residue_name = residue.name().context(FailedToGetResidueNameSnafu)?;
                for atom in residue.atoms() {
                    let element = atom.element().context(ElementMissingSnafu)?;
                    let atom_name = atom.name();
                    if element == &pdbtbx::Element::H && !include_hydrogens {
                        continue;
                    };
                    if atom.hetero() && !include_hetatms {
                        continue;
                    }
                    build_atom!(
                        atoms,
                        atom,
                        element,
                        residue_name,
                        atom_name,
                        Some(chain_id),
                        radii_config,
                        allow_vdw_fallback
                    );
                    temp.push(i);
                    i += 1
                }
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
        parent_to_atoms: &FnvHashMap<isize, Vec<usize>>,
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
        radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
        allow_vdw_fallback: bool,
        include_hydrogens: bool,
        include_hetatms: bool,
    ) -> Result<(Vec<Atom>, FnvHashMap<isize, Vec<usize>>), SASACalcError> {
        let mut atoms = vec![];
        let mut parent_to_atoms = FnvHashMap::default();
        let mut i = 0;
        for residue in pdb.residues() {
            let residue_name = residue.name().context(FailedToGetResidueNameSnafu)?;
            let mut temp = vec![];
            for atom in residue.atoms() {
                let element = atom.element().context(ElementMissingSnafu)?;
                let atom_name = atom.name();
                if element == &pdbtbx::Element::H && !include_hydrogens {
                    continue;
                };
                if atom.hetero() && !include_hetatms {
                    continue;
                }
                build_atom!(
                    atoms,
                    atom,
                    element,
                    residue_name,
                    atom_name,
                    Some(residue.serial_number()),
                    radii_config,
                    allow_vdw_fallback
                );
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

    #[snafu(display(
        "Radius not found for residue '{}' atom '{}' of type '{}'. This error can can be ignored, if you are using the CLI pass --allow-vdw-fallback or use with_allow_vdw_fallback if you are using the API.",
        residue_name,
        atom_name,
        element
    ))]
    RadiusMissing {
        residue_name: String,
        atom_name: String,
        element: String,
    },

    #[snafu(display("Failed to map atoms back to level element"))]
    AtomMapToLevelElementFailed,

    #[snafu(display("Failed to get residue name"))]
    FailedToGetResidueName,

    #[snafu(display("Failed to load radii file: {source}"))]
    RadiiFileLoad { source: std::io::Error },
}

impl<T> SASAOptions<T> {
    /// Create a new SASAOptions with the specified level type
    pub fn new() -> SASAOptions<T> {
        SASAOptions {
            probe_radius: 1.4,
            n_points: 100,
            threads: -1,
            include_hydrogens: false,
            radii_config: None,
            allow_vdw_fallback: false,
            include_hetatms: false,
            _marker: PhantomData,
        }
    }

    /// Set the probe radius (default: 1.4 Angstroms)
    pub fn with_probe_radius(mut self, radius: f32) -> Self {
        self.probe_radius = radius;
        self
    }

    /// Include or exclude HETATM records in protein.
    pub fn with_include_hetatms(mut self, include_hetatms: bool) -> Self {
        self.include_hetatms = include_hetatms;
        self
    }

    /// Set the number of points on the sphere for sampling (default: 100)
    pub fn with_n_points(mut self, points: usize) -> Self {
        self.n_points = points;
        self
    }

    /// Configure the number of threads to use for parallel processing
    ///   - `-1`: Use all available CPU cores (default)
    ///   - `1`: Single-threaded execution (disables parallelism)
    ///   - `> 1`: Use specified number of threads
    pub fn with_threads(mut self, threads: isize) -> Self {
        self.threads = threads;
        self
    }

    /// Include or exclude hydrogen atoms in calculations (default: false)
    pub fn with_include_hydrogens(mut self, include_hydrogens: bool) -> Self {
        self.include_hydrogens = include_hydrogens;
        self
    }

    /// Load custom radii configuration from a file (default: uses embedded protor.config)
    pub fn with_radii_file(mut self, path: &str) -> Result<Self, std::io::Error> {
        self.radii_config = Some(load_radii_from_file(path)?);
        Ok(self)
    }

    /// Allow fallback to PDBTBX van der Waals radii when radius is not found in radii config file (default: false)
    pub fn with_allow_vdw_fallback(mut self, allow: bool) -> Self {
        self.allow_vdw_fallback = allow;
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

impl<T> Default for SASAOptions<T> {
    fn default() -> Self {
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
    /// let (mut pdb, _errors) = pdbtbx::open("./pdbs/example.cif").unwrap();
    /// let result = SASAOptions::<ResidueLevel>::new().process(&pdb);
    /// ```
    pub fn process(&self, pdb: &PDB) -> Result<T::Output, SASACalcError> {
        let (atoms, parent_to_atoms) = T::build_atoms_and_mapping(
            pdb,
            self.radii_config.as_ref(),
            self.allow_vdw_fallback,
            self.include_hydrogens,
            self.include_hetatms,
        )?;
        let atom_sasa =
            calculate_sasa_internal(&atoms, self.probe_radius, self.n_points, self.threads);
        T::process_atoms(&atoms, &atom_sasa, pdb, &parent_to_atoms)
    }
}
