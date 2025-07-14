use crate::Atom;

/// Fast 3D grid for spatial queries
#[repr(C)]
pub(crate) struct SpatialGrid {
    grid: Vec<Vec<Vec<Vec<usize>>>>, // 3D grid of atom indices
    min_bounds: [f32; 3],            // 12 bytes
    cell_size: f32,                  // 4 bytes
    grid_size: [usize; 3],           // 24 bytes
}

impl SpatialGrid {
    pub(crate) fn new(atoms: &[Atom], cell_size: f32) -> Self {
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

    pub(crate) fn locate_within_distance(
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
