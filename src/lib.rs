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
use crate::utils::ARCH;
use structures::spatial_grid::SpatialGrid;
// Re-export io functions for use in the binary crate
use rayon::prelude::*;
pub use utils::io::{sasa_result_to_json, sasa_result_to_protein_object, sasa_result_to_xml};

struct SpherePointsSoA {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
}

impl SpherePointsSoA {
    fn len(&self) -> usize {
        self.x.len()
    }
}

/// Generates points on a sphere using the Golden Section Spiral algorithm
fn generate_sphere_points(n_points: usize) -> SpherePointsSoA {
    let mut x = Vec::with_capacity(n_points);
    let mut y = Vec::with_capacity(n_points);
    let mut z = Vec::with_capacity(n_points);

    let inv_n_points = 1.0 / n_points as f32;

    for i in 0..n_points {
        let i_f32 = i as f32;
        let t = i_f32 * inv_n_points;
        let inclination = (1.0 - 2.0 * t).acos();
        let azimuth = ANGLE_INCREMENT * i_f32;

        // Use sin_cos for better performance
        let (sin_azimuth, cos_azimuth) = azimuth.sin_cos();
        let sin_inclination = inclination.sin();

        x.push(sin_inclination * cos_azimuth);
        y.push(sin_inclination * sin_azimuth);
        z.push(inclination.cos());
    }

    SpherePointsSoA { x, y, z }
}

#[inline(never)]
fn precompute_neighbors(
    atoms: &[Atom],
    active_indices: &[usize],
    probe_radius: f32,
    max_radii: f32,
) -> Vec<Vec<NeighborData>> {
    let cell_size = probe_radius + max_radii;
    let grid = SpatialGrid::new(atoms, active_indices, cell_size);

    let mut neighbors = Vec::with_capacity(active_indices.len());
    let mut temp_candidates = Vec::with_capacity(64);
    let sr = probe_radius + (max_radii * 2.0);
    let sr_squared = sr * sr;

    for &orig_idx in active_indices {
        let atom = &atoms[orig_idx];
        let xyz = atom.position.coords.xyz();
        grid.locate_within_distance([xyz[0], xyz[1], xyz[2]], sr_squared, &mut temp_candidates);

        // Sort candidates by distance (closest first for early exit optimization)
        let center_pos = atom.position;
        temp_candidates.sort_unstable_by(|&a_idx, &b_idx| {
            let pos_a = atoms[a_idx].position;
            let pos_b = atoms[b_idx].position;

            let dx_a = center_pos.x - pos_a.x;
            let dy_a = center_pos.y - pos_a.y;
            let dz_a = center_pos.z - pos_a.z;
            let dist_sq_a = dx_a * dx_a + dy_a * dy_a + dz_a * dz_a;

            let dx_b = center_pos.x - pos_b.x;
            let dy_b = center_pos.y - pos_b.y;
            let dz_b = center_pos.z - pos_b.z;
            let dist_sq_b = dx_b * dx_b + dy_b * dy_b + dz_b * dz_b;

            dist_sq_a
                .partial_cmp(&dist_sq_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

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

struct AtomSasaKernel<'a> {
    atom_index: usize,
    atoms: &'a [Atom],
    neighbors: &'a [NeighborData],
    sphere_points: &'a SpherePointsSoA,
    probe_radius: f32,
}

impl<'a> pulp::WithSimd for AtomSasaKernel<'a> {
    type Output = f32;

    #[inline(always)]
    fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
        let atom = &self.atoms[self.atom_index];
        let center_pos = atom.position;
        let r = atom.radius + self.probe_radius;
        let r2 = r * r;

        let (sx_chunks, sx_rem) = S::as_simd_f32s(&self.sphere_points.x);
        let (sy_chunks, sy_rem) = S::as_simd_f32s(&self.sphere_points.y);
        let (sz_chunks, sz_rem) = S::as_simd_f32s(&self.sphere_points.z);

        let mut accessible_points = 0.0;

        // Precompute true mask for NOT operation
        let zero = simd.splat_f32s(0.0);
        let true_mask = simd.equal_f32s(zero, zero);

        // Process chunks
        for i in 0..sx_chunks.len() {
            let sx = sx_chunks[i];
            let sy = sy_chunks[i];
            let sz = sz_chunks[i];

            // Initialize with all false (0.0).
            let mut chunk_mask = simd.less_than_f32s(simd.splat_f32s(1.0), simd.splat_f32s(0.0));

            for neighbor in self.neighbors {
                if self.atoms[neighbor.idx as usize].id == atom.id {
                    continue;
                }

                let neighbor_pos = self.atoms[neighbor.idx as usize].position;
                let vx_scalar = center_pos.x - neighbor_pos.x;
                let vy_scalar = center_pos.y - neighbor_pos.y;
                let vz_scalar = center_pos.z - neighbor_pos.z;
                let v_mag_sq =
                    vx_scalar * vx_scalar + vy_scalar * vy_scalar + vz_scalar * vz_scalar;

                let t = neighbor.threshold_squared;
                let limit_scalar = (t - v_mag_sq - r2) / (2.0 * r);

                let vx = simd.splat_f32s(vx_scalar);
                let vy = simd.splat_f32s(vy_scalar);
                let vz = simd.splat_f32s(vz_scalar);
                let limit = simd.splat_f32s(limit_scalar);

                let dot =
                    simd.mul_add_f32s(sx, vx, simd.mul_add_f32s(sy, vy, simd.mul_f32s(sz, vz)));

                let occ = simd.less_than_f32s(dot, limit);
                chunk_mask = simd.or_m32s(chunk_mask, occ);

                let not_occ = simd.xor_m32s(chunk_mask, true_mask);
                if simd.first_true_m32s(not_occ) == S::F32_LANES {
                    break;
                }
            }

            // accessible points are !chunk_mask
            let not_occ = simd.xor_m32s(chunk_mask, true_mask);
            let contribution =
                simd.select_f32s(not_occ, simd.splat_f32s(1.0), simd.splat_f32s(0.0));
            accessible_points += simd.reduce_sum_f32s(contribution);
        }

        // Process remainder
        for i in 0..sx_rem.len() {
            let sx = sx_rem[i];
            let sy = sy_rem[i];
            let sz = sz_rem[i];
            let mut occluded = false;

            for neighbor in self.neighbors {
                if self.atoms[neighbor.idx as usize].id == atom.id {
                    continue;
                }
                let n_pos = self.atoms[neighbor.idx as usize].position;
                let vx = center_pos.x - n_pos.x;
                let vy = center_pos.y - n_pos.y;
                let vz = center_pos.z - n_pos.z;
                let v_mag_sq = vx * vx + vy * vy + vz * vz;

                let t = neighbor.threshold_squared;
                let limit = (t - v_mag_sq - r2) / (2.0 * r);

                let dot = sx * vx + sy * vy + sz * vz;
                if dot < limit {
                    occluded = true;
                    break;
                }
            }

            if !occluded {
                accessible_points += 1.0;
            }
        }

        let surface_area = 4.0 * std::f32::consts::PI * r2;
        let inv_n_points = 1.0 / (self.sphere_points.len() as f32);
        surface_area * accessible_points * inv_n_points
    }
}

/// Takes the probe radius and number of points to use along with a list of Atoms as inputs and returns a Vec with SASA values for each atom.
/// For most users it is recommend that you use `calculate_sasa` instead. This method can be used directly if you do not want to use pdbtbx to load PDB/mmCIF files or want to load them from a different source.
/// Probe Radius Default: 1.4
/// Point Count Default: 100
/// Include Hydrogens Default: false
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
///                 parent_id: None,
///                 is_hydrogen: atom.element() == Some(&pdbtbx::Element::H)
///     })
///  }
///  let sasa = calculate_sasa_internal(&atoms, 1.4, 100, true, false);
/// ```
pub fn calculate_sasa_internal(
    atoms: &[Atom],
    probe_radius: f32,
    n_points: usize,
    parallel: bool,
    include_hydrogens: bool,
) -> Vec<f32> {
    let active_indices: Vec<usize> = if include_hydrogens {
        (0..atoms.len()).collect()
    } else {
        atoms
            .iter()
            .enumerate()
            .filter_map(|(i, atom)| (!atom.is_hydrogen).then_some(i))
            .collect()
    };

    let sphere_points = generate_sphere_points(n_points);

    let max_radii = active_indices
        .iter()
        .map(|&i| atoms[i].radius)
        .fold(0.0f32, f32::max);

    let neighbor_lists = precompute_neighbors(atoms, &active_indices, probe_radius, max_radii);

    let process_atom = |(list_idx, neighbors): (usize, &Vec<NeighborData>)| {
        let orig_idx = active_indices[list_idx];
        ARCH.dispatch(AtomSasaKernel {
            atom_index: orig_idx,
            atoms,
            neighbors,
            sphere_points: &sphere_points,
            probe_radius,
        })
    };

    let active_results: Vec<f32> = if parallel {
        neighbor_lists
            .par_iter()
            .enumerate()
            .map(process_atom)
            .collect()
    } else {
        neighbor_lists
            .iter()
            .enumerate()
            .map(process_atom)
            .collect()
    };

    // Map results back to original indices (hydrogens get 0.0)
    let mut results = vec![0.0; atoms.len()];
    for (list_idx, &orig_idx) in active_indices.iter().enumerate() {
        results[orig_idx] = active_results[list_idx];
    }
    results
}
