use crate::Atom;

/// Spatial grid for optimized search
#[repr(C)]
pub(crate) struct SpatialGrid {
    /// A single contiguous buffer of atom indices, sorted by grid cell.
    sorted_atom_indices: Vec<usize>,
    /// A parallel buffer of positions (x, y, z) corresponding to `sorted_atom_indices`.
    sorted_positions: Vec<f32>,
    /// The start index in `sorted_atom_indices` for each flat cell index.
    cell_starts: Vec<u32>,
    /// min bounds for grid
    min_bounds: [f32; 3],
    /// Size of each cell
    cell_size: f32,
    /// Dimensions of grid
    grid_dims: [usize; 3],
    /// Dimension strides for flattened atom vec
    strides: [usize; 3],
}

impl SpatialGrid {
    pub(crate) fn new(atoms: &[Atom], cell_size: f32) -> Self {
        // Calculate Bounds
        let mut min_bounds = [f32::INFINITY; 3];
        let mut max_bounds = [f32::NEG_INFINITY; 3];

        for atom in atoms {
            let pos = atom.position.coords.xyz();
            for i in 0..3 {
                min_bounds[i] = min_bounds[i].min(pos[i]);
                max_bounds[i] = max_bounds[i].max(pos[i]);
            }
        }

        // Add padding to bounds
        for i in 0..3 {
            min_bounds[i] -= cell_size;
            max_bounds[i] += cell_size;
        }

        // Calculate grid size
        let grid_dims = [
            ((max_bounds[0] - min_bounds[0]) / cell_size).ceil() as usize + 1,
            ((max_bounds[1] - min_bounds[1]) / cell_size).ceil() as usize + 1,
            ((max_bounds[2] - min_bounds[2]) / cell_size).ceil() as usize + 1,
        ];

        // Use stride to index into 3D grid as 1D vec
        let strides = [1, grid_dims[0], grid_dims[0] * grid_dims[1]];
        let num_cells = grid_dims[0] * grid_dims[1] * grid_dims[2];

        // Count number of atoms in each cell for future indexing into cells
        let mut cell_counts = vec![0u32; num_cells];
        let inv_cell_size = 1.0 / cell_size;

        // Helper to get flat index using strides
        let get_cell_idx = |pos: &[f32; 3]| -> usize {
            let x = ((pos[0] - min_bounds[0]) * inv_cell_size) as usize;
            let y = ((pos[1] - min_bounds[1]) * inv_cell_size) as usize;
            let z = ((pos[2] - min_bounds[2]) * inv_cell_size) as usize;
            x * strides[0] + y * strides[1] + z * strides[2]
        };

        for atom in atoms {
            let pos = [atom.position.x, atom.position.y, atom.position.z];
            let idx = get_cell_idx(&pos);
            cell_counts[idx] += 1
        }

        // Create cell start index using cell_counts
        let mut cell_starts = vec![0u32; num_cells + 1];
        let mut current_offset = 0;
        for i in 0..num_cells {
            cell_starts[i] = current_offset;
            current_offset += cell_counts[i];
        }
        cell_starts[num_cells] = current_offset;

        // Fill sorted_positions and sorted_atom_indices
        let mut sorted_atom_indices = vec![0usize; atoms.len()];
        let mut sorted_positions = vec![0.0f32; atoms.len() * 3];
        let mut write_heads = cell_starts.clone();

        for (atom_idx, atom) in atoms.iter().enumerate() {
            let pos = [atom.position.x, atom.position.y, atom.position.z];
            let cell_idx = get_cell_idx(&pos);

            let write_pos = write_heads[cell_idx] as usize;

            // Store original index
            sorted_atom_indices[write_pos] = atom_idx;

            // Store position contiguously (Structure of Arrays-ish)
            // Storing as flat floats [x,y,z, x,y,z] improves prefetching
            let pos_offset = write_pos * 3;
            sorted_positions[pos_offset] = pos[0];
            sorted_positions[pos_offset + 1] = pos[1];
            sorted_positions[pos_offset + 2] = pos[2];

            // Increment head
            write_heads[cell_idx] += 1;
        }

        SpatialGrid {
            sorted_atom_indices,
            sorted_positions,
            cell_starts,
            min_bounds,
            cell_size,
            grid_dims,
            strides,
        }
    }

    #[inline(always)]
    pub(crate) fn locate_within_distance(
        &self,
        point: [f32; 3],
        radius_squared: f32,
        result: &mut Vec<usize>,
    ) {
        result.clear();
        let radius = radius_squared.sqrt();
        let inv_cell_size = 1.0 / self.cell_size;

        // Global bounds for clamping
        let max_grid_x = self.grid_dims[0] as isize - 1;
        let max_grid_y = self.grid_dims[1] as isize - 1;
        let max_grid_z = self.grid_dims[2] as isize - 1;

        // Initial Z range
        let min_z = (((point[2] - radius) - self.min_bounds[2]) * inv_cell_size).floor() as isize;
        let max_z = (((point[2] + radius) - self.min_bounds[2]) * inv_cell_size).floor() as isize;
        let min_z = min_z.max(0) as usize;
        let max_z = max_z.min(max_grid_z) as usize;

        let stride_y = self.strides[1];
        let stride_z = self.strides[2];

        let px = point[0];
        let py = point[1];
        let pz = point[2];

        for z in min_z..=max_z {
            let cell_min_z = self.min_bounds[2] + (z as f32) * self.cell_size;
            let cell_max_z = cell_min_z + self.cell_size;
            let dz = if pz < cell_min_z {
                cell_min_z - pz
            } else if pz > cell_max_z {
                pz - cell_max_z
            } else {
                0.0
            };
            let dz_sq = dz * dz;
            if dz_sq > radius_squared {
                continue;
            }

            // Calculate remaining radius for XY plane
            let r_xy_sq = radius_squared - dz_sq;
            let r_xy = r_xy_sq.sqrt();

            // Calculate Y range for this Z-slice
            let min_y = (((py - r_xy) - self.min_bounds[1]) * inv_cell_size).floor() as isize;
            let max_y = (((py + r_xy) - self.min_bounds[1]) * inv_cell_size).floor() as isize;
            let min_y = min_y.max(0) as usize;
            let max_y = max_y.min(max_grid_y) as usize;

            let z_offset = z * stride_z;

            for y in min_y..=max_y {
                let cell_min_y = self.min_bounds[1] + (y as f32) * self.cell_size;
                let cell_max_y = cell_min_y + self.cell_size;
                let dy = if py < cell_min_y {
                    cell_min_y - py
                } else if py > cell_max_y {
                    py - cell_max_y
                } else {
                    0.0
                };
                let dy_sq = dy * dy;

                if dy_sq > r_xy_sq {
                    continue;
                }

                // Calculate remaining radius for X line
                let r_x_sq = r_xy_sq - dy_sq;
                let r_x = r_x_sq.sqrt();

                // Calculate X range for this row
                let min_x = (((px - r_x) - self.min_bounds[0]) * inv_cell_size).floor() as isize;
                let max_x = (((px + r_x) - self.min_bounds[0]) * inv_cell_size).floor() as isize;
                let min_x = min_x.max(0) as usize;
                let max_x = max_x.min(max_grid_x) as usize;

                let y_offset = z_offset + y * stride_y;

                // Iterate row cells
                for x in min_x..=max_x {
                    let cell_idx = y_offset + x;
                    let start = self.cell_starts[cell_idx] as usize;
                    let end = self.cell_starts[cell_idx + 1] as usize;

                    if start == end {
                        continue;
                    }

                    for i in start..end {
                        let pos_idx = i * 3;
                        let ax = self.sorted_positions[pos_idx];
                        let ay = self.sorted_positions[pos_idx + 1];
                        let az = self.sorted_positions[pos_idx + 2];

                        let dx = px - ax;
                        let dy = py - ay;
                        let dz = pz - az;

                        if dx * dx + dy * dy + dz * dz <= radius_squared {
                            result.push(self.sorted_atom_indices[i]);
                        }
                    }
                }
            }
        }
    }
}
