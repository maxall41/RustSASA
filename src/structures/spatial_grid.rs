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
        let search_radius = radius + self.cell_size;
        let inv_cell_size = 1.0 / self.cell_size;

        // Calculate index ranges
        // Using min/max and clamping to grid dimensions
        let min_x =
            (((point[0] - search_radius) - self.min_bounds[0]) * inv_cell_size).max(0.0) as usize;
        let min_y =
            (((point[1] - search_radius) - self.min_bounds[1]) * inv_cell_size).max(0.0) as usize;
        let min_z =
            (((point[2] - search_radius) - self.min_bounds[2]) * inv_cell_size).max(0.0) as usize;

        let max_x = (((point[0] + search_radius) - self.min_bounds[0]) * inv_cell_size)
            .min(self.grid_dims[0] as f32 - 1.0) as usize;
        let max_y = (((point[1] + search_radius) - self.min_bounds[1]) * inv_cell_size)
            .min(self.grid_dims[1] as f32 - 1.0) as usize;
        let max_z = (((point[2] + search_radius) - self.min_bounds[2]) * inv_cell_size)
            .min(self.grid_dims[2] as f32 - 1.0) as usize;

        let stride_y = self.strides[1];
        let stride_z = self.strides[2];

        // Pre-calculate query point components for registers
        let px = point[0];
        let py = point[1];
        let pz = point[2];

        // Iterate flat indices
        for z in min_z..=max_z {
            let z_offset = z * stride_z;
            for y in min_y..=max_y {
                let y_offset = z_offset + y * stride_y;

                let row_start_idx = y_offset + min_x;
                let row_end_idx = y_offset + max_x;

                for cell_idx in row_start_idx..=row_end_idx {
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

                        let dist_sq = dx * dx + dy * dy + dz * dz;

                        if dist_sq <= radius_squared {
                            result.push(self.sorted_atom_indices[i]);
                        }
                    }
                }
            }
        }
    }
}
