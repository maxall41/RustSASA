use approx::assert_relative_eq;
use rust_sasa::{Atom, calculate_sasa_internal};
use std::f32::consts::PI;

const PROBE_RADIUS: f32 = 1.4;
const HIGH_PRECISION_N_POINTS: usize = 50000;
const RELATIVE_TOLERANCE: f32 = 0.005; // 0.5% tolerance

fn create_atom(x: f32, y: f32, z: f32, radius: f32, id: usize) -> Atom {
    Atom {
        position: [x, y, z],
        radius,
        id,
        parent_id: None,
    }
}

#[test]
fn test_single_sphere() {
    let radius = 2.0;
    let atoms = vec![create_atom(0.0, 0.0, 0.0, radius, 1)];

    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    let effective_radius = radius + PROBE_RADIUS;
    let expected_area = 4.0 * PI * effective_radius.powi(2);

    assert_relative_eq!(sasa[0], expected_area, max_relative = RELATIVE_TOLERANCE);
}

#[test]
fn test_two_non_overlapping_spheres() {
    let radius = 2.0;
    // Distance 10.0 is enough to be non-overlapping
    let atoms = vec![
        create_atom(0.0, 0.0, 0.0, radius, 1),
        create_atom(10.0, 0.0, 0.0, radius, 2),
    ];

    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    let effective_radius = radius + PROBE_RADIUS;
    let expected_single_area = 4.0 * PI * effective_radius.powi(2);

    assert_relative_eq!(
        sasa[0],
        expected_single_area,
        max_relative = RELATIVE_TOLERANCE
    );
    assert_relative_eq!(
        sasa[1],
        expected_single_area,
        max_relative = RELATIVE_TOLERANCE
    );
    assert_relative_eq!(
        sasa.iter().sum::<f32>(),
        2.0 * expected_single_area,
        max_relative = RELATIVE_TOLERANCE
    );
}

#[test]
fn test_two_overlapping_spheres() {
    let radius = 2.0;
    let dist = 4.0;
    // Effective radius r = 3.4. Dist = 4.0.
    // 3.4 + 3.4 = 6.8 > 4.0 -> Overlap.
    let atoms = vec![
        create_atom(0.0, 0.0, 0.0, radius, 1),
        create_atom(dist, 0.0, 0.0, radius, 2),
    ];

    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    let r = radius + PROBE_RADIUS;
    // Distance from center to radical plane
    let x = (r.powi(2) - r.powi(2) + dist.powi(2)) / (2.0 * dist); // = dist / 2
    let h_buried = r - x;
    let area_buried = 2.0 * PI * r * h_buried;
    let expected_exposed = 4.0 * PI * r.powi(2) - area_buried;

    assert_relative_eq!(sasa[0], expected_exposed, max_relative = RELATIVE_TOLERANCE);
    assert_relative_eq!(sasa[1], expected_exposed, max_relative = RELATIVE_TOLERANCE);
}

#[test]
fn test_contained_sphere() {
    let r_large = 10.0;
    let r_small = 2.0;
    let dist = 2.0;
    // r_eff_large = 11.4
    // r_eff_small = 3.4
    // dist + r_eff_small = 5.4 < 11.4 -> Contained

    let atoms = vec![
        create_atom(0.0, 0.0, 0.0, r_large, 1),
        create_atom(dist, 0.0, 0.0, r_small, 2),
    ];

    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    let r_eff_large = r_large + PROBE_RADIUS;
    let expected_large_area = 4.0 * PI * r_eff_large.powi(2);

    assert_relative_eq!(
        sasa[0],
        expected_large_area,
        max_relative = RELATIVE_TOLERANCE
    );
    assert_relative_eq!(sasa[1], 0.0, epsilon = RELATIVE_TOLERANCE);
}

#[test]
fn test_three_spheres_linear_chain() {
    let radius = 2.0;
    let dist = 5.0;
    // 0 -- 5 -- 10
    // r_eff = 3.4
    // Neighbor overlap at dist 5. (6.8 > 5)
    // Non-neighbor dist 10. (6.8 < 10) -> No overlap between 1 and 3.

    let atoms = vec![
        create_atom(0.0, 0.0, 0.0, radius, 1),
        create_atom(dist, 0.0, 0.0, radius, 2),
        create_atom(dist * 2.0, 0.0, 0.0, radius, 3),
    ];

    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    let r = radius + PROBE_RADIUS;
    // Overlap with neighbor at dist 5
    let x = (r.powi(2) - r.powi(2) + dist.powi(2)) / (2.0 * dist); // 2.5
    let h_buried = r - x; // 3.4 - 2.5 = 0.9
    let area_buried_one_side = 2.0 * PI * r * h_buried;

    // End spheres have 1 neighbor
    let expected_end = 4.0 * PI * r.powi(2) - area_buried_one_side;

    // Middle sphere has 2 neighbors on opposite sides
    let expected_middle = 4.0 * PI * r.powi(2) - 2.0 * area_buried_one_side;

    assert_relative_eq!(sasa[0], expected_end, max_relative = RELATIVE_TOLERANCE);
    assert_relative_eq!(sasa[2], expected_end, max_relative = RELATIVE_TOLERANCE);
    assert_relative_eq!(sasa[1], expected_middle, max_relative = RELATIVE_TOLERANCE);
}

#[test]
fn test_empty_atom_list() {
    let atoms: Vec<Atom> = vec![];
    let sasa = calculate_sasa_internal(&atoms, PROBE_RADIUS, HIGH_PRECISION_N_POINTS, 1);

    assert!(
        sasa.is_empty(),
        "Empty atom list should return empty result"
    );
}
