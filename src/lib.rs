// use std::collections::{HashMap, HashSet};
// use std::f64::consts::PI;
// use ndarray::{Array2, Ix};
// use ndarray::prelude::*;
// use pdbtbx::{Element, PDB, StrictnessLevel};
// use ndarray_stats::QuantileExt;
// use ordered_float::OrderedFloat;
// use rstar::RTree;
//
// #[derive(Hash, Eq, PartialEq)]
// struct HashIndex {
//     p1: OrderedFloat<f64>,
//     p2: OrderedFloat<f64>,
//     p3: OrderedFloat<f64>
// }
//
// // Calculate euclidean distance between points
// fn distance(p1: &(f64, f64, f64), p2: &(f64, f64, f64)) -> f64 {
//     let (x1, y1, z1) = p1;
//     let (x2, y2, z2) = p2;
//
//     let dx = x1 - x2;
//     let dy = y1 - y2;
//     let dz = z1 - z2;
//
//     (dx*dx + dy*dy + dz*dz).sqrt()
// }
//
// /// Implementation is based off of biopython. Biopython implementation can be found here: https://github.com/biopython/biopython/blob/master/Bio/PDB/SASA.py
// struct ShrakeRupley {
//     probe_radius: f64,
//     n_points: u64,
//     sphere: Array2<f64>
// }
//
// impl ShrakeRupley {
//     fn compute(&mut self,pdb: &PDB) {
//         let num_atoms = pdb.atom_count();
//
//         let mut atom_coords : Vec<Array1<f64>> = vec![];
//         let mut radii = vec![];
//         let mut atom_serial_numbers = vec![];
//         for atom in pdb.atoms() {
//             let pos = atom.pos();
//             atom_coords.push(array![pos.0,pos.1,pos.2]);
//             atom_serial_numbers.push(atom.serial_number());
//             radii.push(atom.element().unwrap().atomic_radius().van_der_waals.unwrap())
//         }
//         let mut radii_arr = Array::from_vec(radii);
//         radii_arr = radii_arr + self.probe_radius;
//         let twice_maxradii = radii_arr.max().unwrap() * 2.0;
//
//
//
//         let tree = pdb.create_atom_rtree();
//
//
//         let mut asa_array = Array::<f64, Ix2>::zeros((num_atoms as Ix, 1).f());
//
//         for i in 0..num_atoms {
//             let r_i = radii_arr[i];
//             // Move sphere to atom
//             let s_on_i : Array::<f64, Ix2> = (self.sphere.clone() * r_i) + &atom_coords[i];
//             let mut avset : HashMap<HashIndex,u64>  = HashMap::new();
//             let mut sub_tree = RTree::new();
//             let mut j = 0;
//             for x in s_on_i.rows() {
//                let mut tp = [0.0, 0.0, 0.0];
//                 let mut i = 0;
//                 for vl in x {
//                     tp[i] = *vl;
//                     i +=  1;
//                 }
//                 sub_tree.insert(tp);
//                 avset.insert(HashIndex {
//                     p1: OrderedFloat(tp[0]),
//                     p2: OrderedFloat(tp[1]),
//                     p3: OrderedFloat(tp[2]),
//                 },j);
//                 j += 1
//             }
//
//             let qp = &atom_coords[i];
//             let query_point = (qp[0],qp[1],qp[2]);
//
//             // locate_within_distance takes squared distance
//             for atom in tree.locate_within_distance(query_point,twice_maxradii * twice_maxradii) {
//                 // Skip current atom
//                 if atom.serial_number() == atom_serial_numbers[i] {
//                     continue
//                 }
//                 let current_radii = atom.element().unwrap().atomic_radius().van_der_waals.unwrap();
//                 if distance(&atom.pos(),&query_point) < (r_i + current_radii) {
//                     // sub_tree.locate_within_distance((), ());
//                     for ix in sub_tree.locate_within_distance([atom.pos().0,atom.pos().1,atom.pos().2],current_radii * current_radii) {
//                         avset.remove(&HashIndex {
//                             p1: OrderedFloat(ix[0]),
//                             p2: OrderedFloat(ix[1]),
//                             p3: OrderedFloat(ix[2]),
//                         });
//                     }
//                 }
//             }
//         }
//
//
//         //pdb.atoms()
//     }
// }
//
// fn setup_sphere(n: f64) -> Array<f64, Ix2> {
//     // Uses the golden spiral algorithm to place points 'evenly' on the sphere
//     // surface. We compute this once and then move the sphere to the centroid
//     // of each atom as we compute the ASAs.
//
//     let dl = PI * (3.0 - 5.0_f64.powf(0.5));
//     let dz = 2.0 / n;
//
//     let mut longitude: f64 = 0.0;
//     let mut z = 1.0 - dz / 2.0;
//
//     let mut coords = Array::<f64, Ix2>::zeros((n as Ix, 3).f());
//
//     for k in 0..n as u64 {
//         let r = (1.0 - z * z).powf(0.5);
//         coords.slice_mut(s![k as i32, 0]).fill(f64::cos(longitude) * r);
//         coords.slice_mut(s![k as i32, 1]).fill(f64::sin(longitude) * r);
//         coords.slice_mut(s![k as i32, 2]).fill(z);
//         z = z - dz;
//         longitude += dl
//     }
//
//     return coords
// }
//
// fn init_shrake_rupley(probe_radius: Option<f64>,n_points: Option<u64>) -> ShrakeRupley {
//     let mut probe_radius_default = 1.4;
//     let mut n_points_default = 100;
//     if probe_radius.is_some() {
//         probe_radius_default = probe_radius.unwrap();
//     }
//     if n_points.is_some() {
//         n_points_default = n_points.unwrap();
//     }
//     let sphere = setup_sphere(n_points_default as f64);
//     let mut new = ShrakeRupley {
//         probe_radius: probe_radius_default,
//         n_points: n_points_default,
//         sphere
//     };
//     return new;
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn testing() {
//         let mut result = init_shrake_rupley(Some(1.4), Some(100));
//         let (mut pdb, _errors) = pdbtbx::open(
//             "./AF-A0A2K5XT84-F1-model_v4.cif",
//             StrictnessLevel::Medium
//         ).unwrap();
//         result.compute(&pdb);
//     }
// }

use nalgebra::{Point3, Vector3};
use rayon::prelude::*;
use std::sync::Mutex;

// A simple structure to model an Atom
struct Atom {
    position: Point3<f32>, // You can adjust types based on precision requirements
    radius: f32,
}

// The core ASA calculation function
fn calculate_asa(atoms: &[Atom], probe_radius: f32, n_points: usize) -> Vec<f32> {
    let sphere_points = generate_sphere_points(n_points);

    atoms.par_iter().map(|atom| {
        let mut accessible_points = 0;

        for sphere_point in &sphere_points {
            let test_point = atom.position + sphere_point * (atom.radius + probe_radius);
            if is_accessible(&test_point, atom, atoms, probe_radius) {
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn testing() {
        // Example usage
        let atoms = vec![
            Atom {
                position: Point3::new(0.0, 0.0, 0.0),
                radius: 1.5,
            },
            Atom {
                position: Point3::new(3.0, 0.0, 0.0),
                radius: 1.5,
            },
        ];
        let probe_radius = 1.4;
        let n_points = 1000;
        let sasa = calculate_asa(&atoms, probe_radius, n_points);
        println!("{:?}", sasa);
    }
}