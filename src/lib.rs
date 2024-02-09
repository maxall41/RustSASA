use nalgebra::{Point3, Vector3};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use rstar::{RTree, RTreeObject, AABB, PointDistance};

// A simple structure to model an Atom
#[derive(Clone)]
struct Atom {
    position: Point3<f32>, // You can adjust types based on precision requirements
    radius: f32,
    id: usize
}

impl RTreeObject for Atom {
    type Envelope = AABB<[f32; 3]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(<[f32; 3]>::from(self.position))
    }
}

impl PointDistance for Atom {
    fn distance_2(&self, other: &[f32; 3]) -> f32 {
        let xyz = self.position.coords.xyz();
        let z = xyz[2];
        let y = xyz[1];
        let x = xyz[0];
        // No square root as that is required by the package
        (other[2] - z).mul_add(
            other[2] - z,
            (other[1] - y).mul_add(other[1] - y, (other[0] - x).powi(2)),
        )
    }
}

// The core ASA calculation function
fn calculate_asa(atoms: &[Atom], probe_radius: f32, n_points: usize) -> Vec<f32> {
    let sphere_points = generate_sphere_points(n_points);

    // Create R*-tree from atoms for efficient spatial lookup
    let tree = RTree::bulk_load(atoms.to_vec());
    let tree_arc = Arc::new(tree); // Use Arc for safe sharing among threads
    let mut max_radii = 0.0;
    for atom in atoms {
        if atom.radius > max_radii {
            max_radii = atom.radius;
        }
    }
    atoms.par_iter().map(|atom| {
        let mut accessible_points = 0;

        for sphere_point in &sphere_points {
            let test_point = atom.position + sphere_point * (atom.radius + probe_radius);
            if is_accessible_rstar(&test_point, atom, &tree_arc, probe_radius, max_radii) {
                accessible_points += 1;
            }
        }
        4.0 * std::f32::consts::PI * (atom.radius + probe_radius).powi(2) * (accessible_points as f32) / (n_points as f32)
    }).collect()
}

// Generates points on a sphere using the Golden Section Spiral algorithm
fn generate_sphere_points(n_points: usize) -> Vec<Vector3<f32>> {
    let mut points = Vec::with_capacity(n_points);
    let golden_ratio = (1.0 + 5f32.sqrt()) / 2.0;
    let angle_increment = 2.0 * std::f32::consts::PI * golden_ratio;

    for i in 0..n_points {
        let t = i as f32 / n_points as f32;
        let inclination = (1.0 - 2.0 * t).acos();
        let azimuth = angle_increment * i as f32;

        let x = inclination.sin() * azimuth.cos();
        let y = inclination.sin() * azimuth.sin();
        let z = inclination.cos();

        points.push(Vector3::new(x, y, z));
    }

    points
}

// Determines if a given point is accessible
fn is_accessible(test_point: &Point3<f32>, atom: &Atom, atoms: &[Atom], probe_radius: f32) -> bool {
    for other_atom in atoms {
        if !std::ptr::eq(atom, other_atom) && (test_point - other_atom.position).norm() < (other_atom.radius + probe_radius) {
            return false;
        }
    }
    true
}

fn is_accessible_rstar(test_point: &Point3<f32>, atom: &Atom, atoms: &RTree<Atom>, probe_radius: f32, max_radii: f32) -> bool {
    let xyz = test_point.coords.xyz();
    let sr = probe_radius + (max_radii * 2.0);
    let candidates = atoms.locate_within_distance([xyz[0],xyz[1],xyz[2]],sr * sr);
    for candidate in candidates {
        if atom.id != candidate.id && (test_point - candidate.position).norm() < (candidate.radius + probe_radius) {
            return false;
        }
    }
    true
}


#[cfg(test)]
mod tests {
    use std::time::Instant;
    use pdbtbx::StrictnessLevel;
    use super::*;

    #[test]
    fn testing() {
        let (mut pdb, _errors) = pdbtbx::open(
            "./AF-A0A2K5XT84-F1-model_v4.cif",
            StrictnessLevel::Medium
        ).unwrap();
        let mut atoms = vec![];
        for atom in pdb.atoms() {
            atoms.push(Atom {
                position: Point3::new(atom.pos().0 as f32, atom.pos().1 as f32, atom.pos().2 as f32),
                radius: atom.element().unwrap().atomic_radius().van_der_waals.unwrap() as f32,
                id: atom.serial_number()
            })
        }
        let probe_radius = 1.4;
        let n_points = 100;
        let start = Instant::now();
        let sasa = calculate_asa(&atoms, probe_radius, n_points);
        let duration = start.elapsed();
        println!("{:?}", sasa);
        println!("Time elapsed: {:?}", duration);
    }
}