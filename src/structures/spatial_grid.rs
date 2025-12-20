// Copyright (c) 2024 Maxwell Campbell. Licensed under the MIT License.
use crate::{Atom, NeighborData};

pub struct SpatialGrid {
    /// Atom indices sorted by cell (contiguous per cell)
    atom_indices: Vec<u32>,

    /// Positions in SoA layout
    positions_x: Vec<f32>,
    positions_y: Vec<f32>,
    positions_z: Vec<f32>,

    /// Radii for each sorted atom (parallel to positions)
    radii: Vec<f32>,

    /// Start index in atom_indices for each cell
    cell_starts: Vec<u32>,

    /// Grid parameters
    grid_dims: [u32; 3],
    num_cells: usize,

    /// Precomputed half-shell offsets for the search extent
    half_shell_offsets: Vec<(i32, i32, i32)>,
}

impl SpatialGrid {
    pub fn new(
        atoms: &[Atom],
        active_indices: &[usize],
        cell_size: f32,
        max_search_radius: f32,
    ) -> Self {
        // Calculate bounds
        let (min_bounds, max_bounds) = Self::calculate_bounds(atoms, active_indices, cell_size);
        let inv_cell_size = 1.0 / cell_size;

        // Grid dimensions
        let grid_dims = [
            ((max_bounds[0] - min_bounds[0]) * inv_cell_size).ceil() as u32 + 1,
            ((max_bounds[1] - min_bounds[1]) * inv_cell_size).ceil() as u32 + 1,
            ((max_bounds[2] - min_bounds[2]) * inv_cell_size).ceil() as u32 + 1,
        ];
        let num_cells = (grid_dims[0] * grid_dims[1] * grid_dims[2]) as usize;

        // Calculate search extent
        let search_extent = (max_search_radius / cell_size).ceil() as i32;

        // Precompute half-shell offsets for this search extent
        let half_shell_offsets = Self::compute_half_shell_offsets(search_extent);

        // Count atoms per cell
        let mut cell_counts = vec![0u32; num_cells];
        for &idx in active_indices {
            let cell = Self::get_cell_index_static(
                &atoms[idx].position,
                &min_bounds,
                inv_cell_size,
                &grid_dims,
            );
            cell_counts[cell] += 1;
        }

        // Build cell_starts (exclusive prefix sum)
        let mut cell_starts = vec![0u32; num_cells + 1];
        for i in 0..num_cells {
            cell_starts[i + 1] = cell_starts[i] + cell_counts[i];
        }

        // Allocate and fill sorted arrays
        let n_active = active_indices.len();
        let mut atom_indices = vec![0u32; n_active];
        let mut positions_x = vec![0.0f32; n_active];
        let mut positions_y = vec![0.0f32; n_active];
        let mut positions_z = vec![0.0f32; n_active];
        let mut radii = vec![0.0f32; n_active];

        let mut write_pos = cell_starts[..num_cells].to_vec();

        for &orig_idx in active_indices {
            let atom = &atoms[orig_idx];
            let pos = &atom.position;
            let cell = Self::get_cell_index_static(pos, &min_bounds, inv_cell_size, &grid_dims);

            let wp = write_pos[cell] as usize;
            atom_indices[wp] = orig_idx as u32;
            positions_x[wp] = pos[0];
            positions_y[wp] = pos[1];
            positions_z[wp] = pos[2];
            radii[wp] = atom.radius;

            write_pos[cell] += 1;
        }

        SpatialGrid {
            atom_indices,
            positions_x,
            positions_y,
            positions_z,
            radii,
            cell_starts,
            grid_dims,
            num_cells,
            half_shell_offsets,
        }
    }

    fn calculate_bounds(
        atoms: &[Atom],
        active_indices: &[usize],
        padding: f32,
    ) -> ([f32; 3], [f32; 3]) {
        let mut min_b = [f32::INFINITY; 3];
        let mut max_b = [f32::NEG_INFINITY; 3];

        for &idx in active_indices {
            let pos = &atoms[idx].position;
            for i in 0..3 {
                min_b[i] = min_b[i].min(pos[i]);
                max_b[i] = max_b[i].max(pos[i]);
            }
        }

        for i in 0..3 {
            min_b[i] -= padding;
            max_b[i] += padding;
        }

        (min_b, max_b)
    }

    #[inline(always)]
    fn get_cell_index_static(
        pos: &[f32; 3],
        min_bounds: &[f32; 3],
        inv_cell_size: f32,
        grid_dims: &[u32; 3],
    ) -> usize {
        let x = ((pos[0] - min_bounds[0]) * inv_cell_size) as u32;
        let y = ((pos[1] - min_bounds[1]) * inv_cell_size) as u32;
        let z = ((pos[2] - min_bounds[2]) * inv_cell_size) as u32;
        (x + y * grid_dims[0] + z * grid_dims[0] * grid_dims[1]) as usize
    }

    #[inline(always)]
    fn cell_coords_to_index(&self, cx: i32, cy: i32, cz: i32) -> Option<usize> {
        if cx < 0 || cy < 0 || cz < 0 {
            return None;
        }
        let cx = cx as u32;
        let cy = cy as u32;
        let cz = cz as u32;
        if cx >= self.grid_dims[0] || cy >= self.grid_dims[1] || cz >= self.grid_dims[2] {
            return None;
        }
        Some((cx + cy * self.grid_dims[0] + cz * self.grid_dims[0] * self.grid_dims[1]) as usize)
    }

    #[inline(always)]
    fn index_to_cell_coords(&self, idx: usize) -> (i32, i32, i32) {
        let idx = idx as u32;
        let cz = idx / (self.grid_dims[0] * self.grid_dims[1]);
        let remainder = idx % (self.grid_dims[0] * self.grid_dims[1]);
        let cy = remainder / self.grid_dims[0];
        let cx = remainder % self.grid_dims[0];
        (cx as i32, cy as i32, cz as i32)
    }

    /// Compute half-shell offsets for a given search extent [See http://doi.acm.org/10.1145/1862648.1862653]
    ///
    /// Half-shell means: for each pair of cells, only one cell "owns" the check.
    /// We select cells where (dz > 0) OR (dz == 0 && dy > 0) OR (dz == 0 && dy == 0 && dx >= 0)
    /// Note: dx >= 0 includes self (0,0,0) which we handle specially
    fn compute_half_shell_offsets(extent: i32) -> Vec<(i32, i32, i32)> {
        let mut offsets = Vec::new();

        for dz in -extent..=extent {
            for dy in -extent..=extent {
                for dx in -extent..=extent {
                    // Half-shell condition
                    let include =
                        (dz > 0) || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx >= 0);

                    if include {
                        offsets.push((dx, dy, dz));
                    }
                }
            }
        }

        offsets
    }

    /// Build neighbor lists for all active atoms
    pub fn build_all_neighbor_lists(
        &self,
        atoms: &[Atom],
        active_indices: &[usize],
        probe_radius: f32,
        max_radius: f32,
    ) -> Vec<Vec<NeighborData>> {
        let n_atoms = atoms.len();
        let n_active = active_indices.len();

        // Map original index -> active index
        let mut orig_to_active = vec![u32::MAX; n_atoms];
        for (active_idx, &orig_idx) in active_indices.iter().enumerate() {
            orig_to_active[orig_idx] = active_idx as u32;
        }

        // Preallocate neighbor lists
        let mut neighbors: Vec<Vec<NeighborData>> =
            (0..n_active).map(|_| Vec::with_capacity(80)).collect();

        // Maximum possible search radius (for quick rejection)
        let max_search_radius = max_radius + max_radius + 2.0 * probe_radius;
        let max_search_radius_sq = max_search_radius * max_search_radius;

        // Iterate through all cells using half-shell pattern
        for cell_a in 0..self.num_cells {
            let start_a = self.cell_starts[cell_a] as usize;
            let end_a = self.cell_starts[cell_a + 1] as usize;

            if start_a == end_a {
                continue;
            }

            let (cx, cy, cz) = self.index_to_cell_coords(cell_a);

            // Process this cell against all cells in half-shell
            for &(dx, dy, dz) in &self.half_shell_offsets {
                let cell_b = match self.cell_coords_to_index(cx + dx, cy + dy, cz + dz) {
                    Some(c) => c,
                    None => continue,
                };

                let start_b = self.cell_starts[cell_b] as usize;
                let end_b = self.cell_starts[cell_b + 1] as usize;

                if start_b == end_b {
                    continue;
                }

                let is_self = dx == 0 && dy == 0 && dz == 0;

                if is_self {
                    self.process_self_cell(
                        atoms,
                        &orig_to_active,
                        start_a,
                        end_a,
                        probe_radius,
                        max_radius,
                        max_search_radius_sq,
                        &mut neighbors,
                    );
                } else {
                    self.process_neighbor_cells(
                        atoms,
                        &orig_to_active,
                        start_a,
                        end_a,
                        start_b,
                        end_b,
                        probe_radius,
                        max_radius,
                        max_search_radius_sq,
                        &mut neighbors,
                    );
                }
            }
        }

        // Sort neighbors by distance for early-exit optimization
        self.sort_neighbors_by_distance(atoms, active_indices, &mut neighbors);

        neighbors
    }

    #[inline(always)]
    fn process_self_cell(
        &self,
        atoms: &[Atom],
        orig_to_active: &[u32],
        start: usize,
        end: usize,
        probe_radius: f32,
        max_radius: f32,
        max_search_radius_sq: f32,
        neighbors: &mut [Vec<NeighborData>],
    ) {
        for i in start..end {
            let orig_i = self.atom_indices[i] as usize;
            let active_i = orig_to_active[orig_i];
            if active_i == u32::MAX {
                continue;
            }

            let xi = self.positions_x[i];
            let yi = self.positions_y[i];
            let zi = self.positions_z[i];
            let ri = self.radii[i];
            let id_i = atoms[orig_i].id;

            // Search radius for atom i
            let sr_i = ri + max_radius + 2.0 * probe_radius;
            let sr_i_sq = sr_i * sr_i;

            for j in (i + 1)..end {
                let orig_j = self.atom_indices[j] as usize;

                // Skip if same atom ID
                if atoms[orig_j].id == id_i {
                    continue;
                }

                let dx = xi - self.positions_x[j];
                let dy = yi - self.positions_y[j];
                let dz = zi - self.positions_z[j];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                // Quick rejection
                if dist_sq > max_search_radius_sq {
                    continue;
                }

                let rj = self.radii[j];

                // Search radius for atom j
                let sr_j = rj + max_radius + 2.0 * probe_radius;
                let sr_j_sq = sr_j * sr_j;

                // Check if i finds j
                if dist_sq <= sr_i_sq {
                    let thresh_j = rj + probe_radius;
                    neighbors[active_i as usize].push(NeighborData {
                        idx: orig_j as u32,
                        threshold_squared: thresh_j * thresh_j,
                    });
                }

                // Check if j finds i
                if dist_sq <= sr_j_sq {
                    let active_j = orig_to_active[orig_j];
                    if active_j != u32::MAX {
                        let thresh_i = ri + probe_radius;
                        neighbors[active_j as usize].push(NeighborData {
                            idx: orig_i as u32,
                            threshold_squared: thresh_i * thresh_i,
                        });
                    }
                }
            }
        }
    }

    #[inline(always)]
    fn process_neighbor_cells(
        &self,
        atoms: &[Atom],
        orig_to_active: &[u32],
        start_a: usize,
        end_a: usize,
        start_b: usize,
        end_b: usize,
        probe_radius: f32,
        max_radius: f32,
        max_search_radius_sq: f32,
        neighbors: &mut [Vec<NeighborData>],
    ) {
        for i in start_a..end_a {
            let orig_i = self.atom_indices[i] as usize;
            let active_i = orig_to_active[orig_i];
            if active_i == u32::MAX {
                continue;
            }

            let xi = self.positions_x[i];
            let yi = self.positions_y[i];
            let zi = self.positions_z[i];
            let ri = self.radii[i];
            let id_i = atoms[orig_i].id;

            // Search radius for atom i
            let sr_i = ri + max_radius + 2.0 * probe_radius;
            let sr_i_sq = sr_i * sr_i;

            for j in start_b..end_b {
                let orig_j = self.atom_indices[j] as usize;

                // Skip if same atom ID
                if atoms[orig_j].id == id_i {
                    continue;
                }

                let dx = xi - self.positions_x[j];
                let dy = yi - self.positions_y[j];
                let dz = zi - self.positions_z[j];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                // Quick rejection
                if dist_sq > max_search_radius_sq {
                    continue;
                }

                let rj = self.radii[j];

                // Search radius for atom j
                let sr_j = rj + max_radius + 2.0 * probe_radius;
                let sr_j_sq = sr_j * sr_j;

                // Check if i finds j
                if dist_sq <= sr_i_sq {
                    let thresh_j = rj + probe_radius;
                    neighbors[active_i as usize].push(NeighborData {
                        idx: orig_j as u32,
                        threshold_squared: thresh_j * thresh_j,
                    });
                }

                // Check if j finds i
                if dist_sq <= sr_j_sq {
                    let active_j = orig_to_active[orig_j];
                    if active_j != u32::MAX {
                        let thresh_i = ri + probe_radius;
                        neighbors[active_j as usize].push(NeighborData {
                            idx: orig_i as u32,
                            threshold_squared: thresh_i * thresh_i,
                        });
                    }
                }
            }
        }
    }

    fn sort_neighbors_by_distance(
        &self,
        atoms: &[Atom],
        active_indices: &[usize],
        neighbors: &mut [Vec<NeighborData>],
    ) {
        for (active_idx, neighbor_list) in neighbors.iter_mut().enumerate() {
            if neighbor_list.len() <= 1 {
                continue;
            }

            let center = atoms[active_indices[active_idx]].position;

            neighbor_list.sort_unstable_by(|a, b| {
                let pa = atoms[a.idx as usize].position;
                let pb = atoms[b.idx as usize].position;

                let da = (center[0] - pa[0]).powi(2)
                    + (center[1] - pa[1]).powi(2)
                    + (center[2] - pa[2]).powi(2);
                let db = (center[0] - pb[0]).powi(2)
                    + (center[1] - pb[1]).powi(2)
                    + (center[2] - pb[2]).powi(2);

                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            });
        }
    }
}
