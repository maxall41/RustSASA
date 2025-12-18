#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy)]
#[repr(C)]
pub struct NeighborData {
    pub threshold_squared: f32,
    pub idx: u32,
}

/// This struct represents an individual Atom
#[derive(Clone)]
#[repr(C)]
pub struct Atom {
    /// The 3D position of the atom (12 bytes)
    pub position: [f32; 3],
    /// The Van Der Walls radius of the atom (4 bytes)
    pub radius: f32,
    /// A unique Id for the atom (8 bytes)
    pub id: usize,
    /// Parent Id (8 bytes)
    pub parent_id: Option<isize>,
}

#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChainResult {
    /// Chain name
    pub name: String,
    /// Chain SASA value
    pub value: f32,
}

#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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

#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProteinResult {
    /// The total SASA value for the entire protein
    pub global_total: f32,
    /// The total polar SASA value for the entire protein
    pub polar_total: f32,
    /// The total *non*-polar SASA value for the entire protein
    pub non_polar_total: f32,
}

#[derive(Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SASAResult {
    Atom(Vec<f32>),
    Residue(Vec<ResidueResult>),
    Chain(Vec<ChainResult>),
    Protein(ProteinResult),
}
