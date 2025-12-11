use nalgebra::Point3;
use serde::Serialize;

#[repr(C)]
pub(crate) struct NeighborData {
    pub(crate) threshold_squared: f32,
    pub(crate) idx: u32,
}

/// This struct represents an individual Atom
#[derive(Clone)]
#[repr(C)]
pub struct Atom {
    /// The 3D position of the atom (12 bytes)
    pub position: Point3<f32>,
    /// The Van Der Walls radius of the atom (4 bytes)
    pub radius: f32,
    /// A unique Id for the atom (8 bytes)
    pub id: usize,
    /// Parent Id (8 bytes)
    pub parent_id: Option<isize>,
    /// Whether this atom is a hydrogen atom (1 byte + padding)
    pub is_hydrogen: bool,
}

/// Can be used to specify output resolution of SASA computation for convenience.
#[derive(clap::ValueEnum, Clone, Default, Debug)]
pub enum SASALevel {
    Atom,
    #[default]
    Residue,
    Chain,
    Protein,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct ChainResult {
    /// Chain name
    pub name: String,
    /// Chain SASA value
    pub value: f32,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct ResidueResult {
    /// Residue serial number
    pub serial_number: isize,
    /// SASA value for residue
    pub value: f32,
    //// The name of the residue
    pub name: String,
    /// Whether the residue is polar
    pub is_polar: bool,
    /// Chain ID
    pub chain_id: String,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct ProteinResult {
    /// The total SASA value for the entire protein
    pub global_total: f32,
    /// The total polar SASA value for the entire protein
    pub polar_total: f32,
    /// The total *non*-polar SASA value for the entire protein
    pub non_polar_total: f32,
}

#[derive(Debug, PartialEq, Serialize)]
pub enum SASAResult {
    Atom(Vec<f32>),
    Residue(Vec<ResidueResult>),
    Chain(Vec<ChainResult>),
    Protein(ProteinResult),
}
