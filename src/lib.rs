//! RustSASA is a Rust library for computing the absolute solvent accessible surface area (ASA/SASA) of each atom in a given protein structure using the Shrake-Rupley algorithm.
//! Example:
//! ```rust
//! use pdbtbx::StrictnessLevel;
//! use rust_sasa::{SASAOptions, ResidueLevel};
//!
//! let (mut pdb, _errors) = pdbtbx::open("./pdbs/example.cif").unwrap();
//! let result = SASAOptions::<ResidueLevel>::new().process(&pdb);
//! ```
pub mod options;
// Re-export the new level types and processor trait
pub use options::{AtomLevel, ChainLevel, ProteinLevel, ResidueLevel, SASAProcessor};
use utils::consts::ANGLE_INCREMENT;
pub mod structures;
mod test;
mod utils;

pub use crate::options::*;
pub use crate::structures::atomic::*;

use structures::spatial_grid::SpatialGrid;
// Re-export io functions for use in the binary crate
pub use utils::io::{sasa_result_to_json, sasa_result_to_protein_object, sasa_result_to_xml};

use nalgebra::{Point3, Vector3};
use rayon::prelude::*;

/// Generates points on a sphere using the Golden Section Spiral algorithm
fn generate_sphere_points(n_points: usize) -> Vec<Vector3<f32>> {
    let mut points = Vec::with_capacity(n_points);
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

fn is_accessible_precomputed(
    test_point: &Point3<f32>,
    atom: &Atom,
    neighbors: &[NeighborData],
    atoms: &[Atom],
) -> bool {
    for neighbor_data in neighbors {
        let neighbor = &atoms[neighbor_data.idx as usize];
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
                idx: idx as u32,
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
///             "./pdbs/example.cif",
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
        (0..atoms.len())
            .into_par_iter()
            .map(calculate_atom_sasa_optimized)
            .collect()
    } else {
        (0..atoms.len())
            .map(calculate_atom_sasa_optimized)
            .collect()
    }
}
