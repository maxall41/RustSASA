//! RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm[1].
pub mod options;

// Re-export the new level types and processor trait
pub use options::{AtomLevel, ChainLevel, ProteinLevel, ResidueLevel, SASAProcessor};
mod test;
mod utils;

// Re-export io functions for use in the binary crate
pub use utils::io::{sasa_result_to_json, sasa_result_to_protein_object, sasa_result_to_xml};

use nalgebra::{Point3, Vector3};
use rayon::prelude::*;
use serde::Serialize;

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
    /// Wether the residue is polar
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

/// Fast 3D grid for spatial queries
#[repr(C)]
struct SpatialGrid {
    grid: Vec<Vec<Vec<Vec<usize>>>>, // 3D grid of atom indices
    min_bounds: [f32; 3],            // 12 bytes
    cell_size: f32,                  // 4 bytes
    grid_size: [usize; 3],           // 24 bytes
}

impl SpatialGrid {
    fn new(atoms: &[Atom], cell_size: f32) -> Self {
        // Calculate bounds
        let mut min_bounds = [f32::INFINITY; 3];
        let mut max_bounds = [f32::NEG_INFINITY; 3];

        for atom in atoms {
            let pos = atom.position.coords.xyz();
            for i in 0..3 {
                min_bounds[i] = min_bounds[i].min(pos[i]);
                max_bounds[i] = max_bounds[i].max(pos[i]);
            }
        }

        // Add padding
        for i in 0..3 {
            min_bounds[i] -= cell_size;
            max_bounds[i] += cell_size;
        }

        // Calculate grid size
        let grid_size = [
            ((max_bounds[0] - min_bounds[0]) / cell_size).ceil() as usize + 1,
            ((max_bounds[1] - min_bounds[1]) / cell_size).ceil() as usize + 1,
            ((max_bounds[2] - min_bounds[2]) / cell_size).ceil() as usize + 1,
        ];

        // Initialize grid
        let mut grid = vec![vec![vec![Vec::new(); grid_size[2]]; grid_size[1]]; grid_size[0]];

        // Add atoms to grid
        for (atom_idx, atom) in atoms.iter().enumerate() {
            let pos = atom.position.coords.xyz();
            let cell_x = ((pos[0] - min_bounds[0]) / cell_size) as usize;
            let cell_y = ((pos[1] - min_bounds[1]) / cell_size) as usize;
            let cell_z = ((pos[2] - min_bounds[2]) / cell_size) as usize;

            if cell_x < grid_size[0] && cell_y < grid_size[1] && cell_z < grid_size[2] {
                grid[cell_x][cell_y][cell_z].push(atom_idx);
            }
        }

        SpatialGrid {
            grid,
            cell_size,
            min_bounds,
            grid_size,
        }
    }

    #[inline(never)]
    fn locate_within_distance(
        &self,
        point: [f32; 3],
        radius_squared: f32,
        atoms: &[Atom],
        result: &mut Vec<usize>,
    ) {
        result.clear();
        let radius = radius_squared.sqrt();

        // Calculate cell range to search - use integer math where possible
        let inv_cell_size = 1.0 / self.cell_size;
        let search_radius = radius + self.cell_size;

        let min_cell = [
            (((point[0] - search_radius) - self.min_bounds[0]) * inv_cell_size).max(0.0) as usize,
            (((point[1] - search_radius) - self.min_bounds[1]) * inv_cell_size).max(0.0) as usize,
            (((point[2] - search_radius) - self.min_bounds[2]) * inv_cell_size).max(0.0) as usize,
        ];
        let max_cell = [
            (((point[0] + search_radius) - self.min_bounds[0]) * inv_cell_size)
                .min(self.grid_size[0] as f32 - 1.0) as usize,
            (((point[1] + search_radius) - self.min_bounds[1]) * inv_cell_size)
                .min(self.grid_size[1] as f32 - 1.0) as usize,
            (((point[2] + search_radius) - self.min_bounds[2]) * inv_cell_size)
                .min(self.grid_size[2] as f32 - 1.0) as usize,
        ];

        // Search cells with manual loop unrolling for better performance
        for x in min_cell[0]..=max_cell[0] {
            for y in min_cell[1]..=max_cell[1] {
                for z in min_cell[2]..=max_cell[2] {
                    let cell = &self.grid[x][y][z];
                    for &atom_idx in cell {
                        let atom = &atoms[atom_idx];
                        let pos = atom.position.coords.xyz();
                        // Use fused multiply-add for better performance
                        let dx = point[0] - pos[0];
                        let dy = point[1] - pos[1];
                        let dz = point[2] - pos[2];
                        let dist_sq = dx.mul_add(dx, dy.mul_add(dy, dz * dz));
                        if dist_sq <= radius_squared {
                            result.push(atom_idx);
                        }
                    }
                }
            }
        }
    }
}

/// Generates points on a sphere using the Golden Section Spiral algorithm
#[inline(never)]
fn generate_sphere_points(n_points: usize) -> Vec<Vector3<f32>> {
    let mut points = Vec::with_capacity(n_points);

    // Precompute constants
    const GOLDEN_RATIO: f32 = 1.618033988749894;
    const ANGLE_INCREMENT: f32 = 2.0 * std::f32::consts::PI * GOLDEN_RATIO;
    let inv_n_points = 1.0 / n_points as f32;

    for i in 0..n_points {
        let i_f32 = i as f32;
        let t = i_f32 * inv_n_points;
        let inclination = (1.0 - 2.0 * t).acos();
        let azimuth = ANGLE_INCREMENT * i_f32;

        // Use sin_cos for better performance
        let (sin_azimuth, cos_azimuth) = azimuth.sin_cos();
        let sin_inclination = inclination.sin();

        let x = sin_inclination * cos_azimuth;
        let y = sin_inclination * sin_azimuth;
        let z = inclination.cos();

        points.push(Vector3::new(x, y, z));
    }

    points
}

#[repr(C)]
struct NeighborData {
    threshold_squared: f32, // 4 bytes - put smaller field first
    idx: usize,             // 8 bytes
}

#[inline(never)]
fn is_accessible_precomputed(
    test_point: &Point3<f32>,
    atom: &Atom,
    neighbors: &[NeighborData],
    atoms: &[Atom],
) -> bool {
    // Optimize neighbor checks by processing in pairs for potential SIMD benefits
    // and better cache utilization.
    // let mut neighbors_iter = neighbors.chunks_exact(2);
    for neighbor_data in neighbors {
        let neighbor = &atoms[neighbor_data.idx];
        if atom.id != neighbor.id {
            // Manual distance calculation is faster than nalgebra's norm_squared
            let pos = neighbor.position.coords.xyz();
            let test_pos = test_point.coords.xyz();
            let dx = test_pos[0] - pos[0];
            let dy = test_pos[1] - pos[1];
            let dz = test_pos[2] - pos[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < neighbor_data.threshold_squared {
                return false;
            }
        }
    }
    true
}

#[inline(never)]
fn precompute_neighbors(
    atoms: &[Atom],
    probe_radius: f32,
    max_radii: f32,
) -> Vec<Vec<NeighborData>> {
    let cell_size = probe_radius + max_radii;
    let grid = SpatialGrid::new(atoms, cell_size);

    let mut neighbors = Vec::with_capacity(atoms.len());
    let mut temp_candidates = Vec::with_capacity(64); // Reuse buffer to avoid allocations
    let sr = probe_radius + (max_radii * 2.0);
    let sr_squared = sr * sr;

    for atom in atoms.iter() {
        let xyz = atom.position.coords.xyz();
        grid.locate_within_distance(
            [xyz[0], xyz[1], xyz[2]],
            sr_squared,
            atoms,
            &mut temp_candidates,
        );

        // Precompute squared thresholds for each neighbor
        let mut neighbor_data = Vec::with_capacity(temp_candidates.len());
        for &idx in &temp_candidates {
            let neighbor = &atoms[idx];
            let threshold = neighbor.radius + probe_radius;
            neighbor_data.push(NeighborData {
                idx,
                threshold_squared: threshold * threshold,
            });
        }

        neighbors.push(neighbor_data);
    }

    neighbors
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
///  let sasa = calculate_sasa_internal(&atoms, 1.4, 100,true);
/// ```
#[inline(never)]
pub fn calculate_sasa_internal(
    atoms: &[Atom],
    probe_radius: f32,
    n_points: usize,
    parallel: bool,
) -> Vec<f32> {
    let sphere_points = generate_sphere_points(n_points);

    let mut max_radii = 0.0;
    for atom in atoms {
        if atom.radius > max_radii {
            max_radii = atom.radius;
        }
    }

    // Precompute constants for better performance
    let inv_n_points = 1.0 / n_points as f32;

    // Use precomputed neighbors for better performance
    let neighbor_indices = precompute_neighbors(atoms, probe_radius, max_radii);

    let calculate_atom_sasa_optimized = |atom_idx: usize| {
        let atom = &atoms[atom_idx];
        let mut accessible_points = 0u32;
        let expanded_radius = atom.radius + probe_radius;

        // Get precomputed neighbors for this atom
        let neighbors = &neighbor_indices[atom_idx];

        // Cache atom position for better memory access
        let atom_pos = atom.position.coords.xyz();

        // Use iterator for better performance
        for sphere_point in sphere_points.iter() {
            // Manual vector math is faster than nalgebra operations
            let test_point = Point3::new(
                atom_pos[0] + sphere_point.x * expanded_radius,
                atom_pos[1] + sphere_point.y * expanded_radius,
                atom_pos[2] + sphere_point.z * expanded_radius,
            );
            if is_accessible_precomputed(&test_point, atom, neighbors, atoms) {
                accessible_points += 1;
            }
        }

        // Precompute constants outside the loop
        let surface_area = 4.0 * std::f32::consts::PI * expanded_radius * expanded_radius;
        surface_area * (accessible_points as f32 * inv_n_points)
    };

    if parallel {
        return (0..atoms.len())
            .into_par_iter()
            .map(calculate_atom_sasa_optimized)
            .collect();
    } else {
        return (0..atoms.len())
            .into_iter()
            .map(calculate_atom_sasa_optimized)
            .collect();
    }
}
