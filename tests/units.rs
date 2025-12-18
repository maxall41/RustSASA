mod common;

#[cfg(test)]
mod tests {
    use super::common::data::FIXED_LOW_RES_ATOMS;
    use approx::assert_abs_diff_eq;
    use rust_sasa::options::SASAOptions;
    use rust_sasa::structures::spatial_grid::SpatialGrid;
    use rust_sasa::utils::consts::PROTOR_RADII;
    use rust_sasa::utils::get_protor_radius;
    use rust_sasa::{
        Atom, AtomLevel, ChainLevel, ProteinLevel, ResidueLevel, calculate_sasa_internal,
    };
    use std::time::Instant;

    #[test]
    fn internal_test() {
        let (pdb, _errors) = pdbtbx::open("./pdbs/example.cif").unwrap();
        let mut atoms = vec![];
        for atom in pdb.atoms() {
            let element = atom.element().unwrap();
            atoms.push(Atom {
                position: [
                    atom.pos().0 as f32,
                    atom.pos().1 as f32,
                    atom.pos().2 as f32,
                ],
                radius: element.atomic_radius().van_der_waals.unwrap() as f32,
                id: atom.serial_number(),
                parent_id: None,
            })
        }
        let start = Instant::now();
        let sasa = calculate_sasa_internal(&atoms, 1.4, 100, -1);
        let duration = start.elapsed();
        // Compare element by element since Vec<f32> doesn't implement AbsDiffEq
        assert_eq!(sasa.len(), FIXED_LOW_RES_ATOMS.len());
        for (actual, expected) in sasa.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 25.0);
        }
        println!("Time elapsed (INTERNAL): {duration:?}");
    }

    #[test]
    fn external_test() {
        let (pdb, _errors) = pdbtbx::open("./pdbs/example.cif").unwrap();
        let protein_sasa = SASAOptions::<ProteinLevel>::new().process(&pdb).unwrap();
        let chain_sasa = SASAOptions::<ChainLevel>::new().process(&pdb).unwrap();

        let _residue_sasa = SASAOptions::<ResidueLevel>::new().process(&pdb).unwrap();

        let start = Instant::now();
        let atom_sasa = SASAOptions::<AtomLevel>::new().process(&pdb).unwrap();
        let duration = start.elapsed();
        println!("Time elapsed (ATOM): {duration:?}");

        assert_abs_diff_eq!(protein_sasa.global_total, 20268.004, epsilon = 1500.0);
        // Compare element by element since Vec<f32> doesn't implement AbsDiffEq
        for (actual, expected) in atom_sasa.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 25.0);
        }
        println!("PROTEIN SASA {protein_sasa:?}");
        println!("CHAIN SASA {chain_sasa:?}");
    }

    #[test]
    fn check_pdb_w_bad_seqadv_record() {
        let (pdb, _errors) = pdbtbx::open("./pdbs/bad_seqadv_1A06.pdb").unwrap();
        let protein_sasa = SASAOptions::<ProteinLevel>::new().process(&pdb).unwrap();

        let start = Instant::now();
        let duration = start.elapsed();
        println!("Time elapsed (ATOM): {duration:?}");

        assert_abs_diff_eq!(protein_sasa.global_total, 14466.709, epsilon = 1500.0);
        println!("PROTEIN SASA {protein_sasa:?}");
    }

    #[test]
    fn check_pdb_w_atypical_spacegroup() {
        let (pdb, _errors) = pdbtbx::open("./pdbs/151L_H3.pdb").unwrap();
        let protein_sasa = SASAOptions::<ProteinLevel>::new().process(&pdb).unwrap();

        let start = Instant::now();
        let duration = start.elapsed();
        println!("Time elapsed (ATOM): {duration:?}");

        assert_abs_diff_eq!(protein_sasa.global_total, 9558.812, epsilon = 1500.0);
        println!("PROTEIN SASA {protein_sasa:?}");
    }

    #[test]
    fn external_test_high_res() {
        let (pdb, _errors) = pdbtbx::open("./pdbs/example.cif").unwrap();
        let protein_sasa = SASAOptions::<ProteinLevel>::new()
            .with_n_points(960)
            .process(&pdb)
            .unwrap();
        let chain_sasa = SASAOptions::<ChainLevel>::new()
            .with_n_points(960)
            .process(&pdb)
            .unwrap();
        let _residue_sasa = SASAOptions::<ResidueLevel>::new()
            .with_n_points(960)
            .process(&pdb)
            .unwrap();

        let start = Instant::now();
        let atom_sasa = SASAOptions::<AtomLevel>::new()
            .with_n_points(960)
            .process(&pdb)
            .unwrap();
        let duration = start.elapsed();
        println!("Time elapsed (ATOM): {duration:?}");

        assert_abs_diff_eq!(protein_sasa.global_total, 20131.227, epsilon = 1500.0);
        assert_abs_diff_eq!(protein_sasa.polar_total, 4279.8906, epsilon = 1500.0);
        assert_abs_diff_eq!(protein_sasa.non_polar_total, 15999.43, epsilon = 1500.0);

        assert_abs_diff_eq!(chain_sasa[0].value, 20131.227, epsilon = 1500.0);
        assert_eq!(chain_sasa[0].name, "A");

        // Compare element by element since Vec<f32> doesn't implement AbsDiffEq
        assert_eq!(atom_sasa.len(), FIXED_LOW_RES_ATOMS.len());
        for (actual, expected) in atom_sasa.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 25.0);
        }
    }

    #[test]
    fn spatial_grid() {
        // Create test atoms with varying distances
        let atoms = vec![
            Atom {
                position: [0.0, 0.0, 0.0],
                radius: 1.5,
                id: 1,
                parent_id: None,
            },
            Atom {
                position: [3.0, 0.0, 0.0],
                radius: 1.5,
                id: 2,
                parent_id: None,
            },
            Atom {
                position: [0.0, 3.0, 0.0],
                radius: 1.5,
                id: 3,
                parent_id: None,
            },
            Atom {
                position: [20.0, 20.0, 20.0],
                radius: 1.5,
                id: 4,
                parent_id: None,
            },
        ];

        // Create spatial grid with cell size of 5.0 and probe radius
        let active_indices: Vec<usize> = (0..atoms.len()).collect();
        let probe_radius = 1.4;
        let max_radius = 1.5;
        let cell_size = 5.0;
        let max_search_radius = max_radius + max_radius + 2.0 * probe_radius;

        let grid = SpatialGrid::new(&atoms, &active_indices, cell_size, max_search_radius);

        // Build neighbor lists for all atoms
        let neighbors =
            grid.build_all_neighbor_lists(&atoms, &active_indices, probe_radius, max_radius);

        // Test 1: Atom 0 at (0,0,0) should have atoms 1 and 2 as neighbors (distance ~3.0)
        // but not atom 3 at (20,20,20) which is too far away
        assert!(
            neighbors[0].len() >= 2,
            "Atom 0 should have at least 2 neighbors"
        );
        let neighbor_indices: Vec<u32> = neighbors[0].iter().map(|n| n.idx).collect();
        assert!(
            neighbor_indices.contains(&1),
            "Atom 0 should have atom 1 as neighbor"
        );
        assert!(
            neighbor_indices.contains(&2),
            "Atom 0 should have atom 2 as neighbor"
        );
        assert!(
            !neighbor_indices.contains(&3),
            "Atom 0 should not have atom 3 (too far)"
        );

        // Test 2: Atom 3 at (20,20,20) should have no neighbors (it's isolated)
        assert_eq!(neighbors[3].len(), 0, "Atom 3 should have no neighbors");

        // Test 3: Atoms 1 and 2 should each have at least atom 0 as a neighbor
        let neighbor_indices_1: Vec<u32> = neighbors[1].iter().map(|n| n.idx).collect();
        assert!(
            neighbor_indices_1.contains(&0),
            "Atom 1 should have atom 0 as neighbor"
        );

        let neighbor_indices_2: Vec<u32> = neighbors[2].iter().map(|n| n.idx).collect();
        assert!(
            neighbor_indices_2.contains(&0),
            "Atom 2 should have atom 0 as neighbor"
        );
    }

    #[test]
    fn test_protor_radii_loading() {
        assert_eq!(
            get_protor_radius("ASN", "CA"),
            Some(1.88),
            "ASN CA should map to C4H1 with radius 1.88"
        );
        assert_eq!(
            get_protor_radius("ASN", "N"),
            Some(1.64),
            "ASN N should map to N3H2 with radius 1.64"
        );
        assert_eq!(
            get_protor_radius("ASN", "CB"),
            Some(1.88),
            "ASN CB should map to C4H2 with radius 1.88"
        );
        assert_eq!(
            get_protor_radius("CYS", "SG"),
            Some(1.77),
            "CYS SG should map to S2H1 with radius 1.77"
        );
        assert_eq!(
            get_protor_radius("XXX", "YY"),
            None,
            "Unknown combinations should return None"
        );
        assert_eq!(
            get_protor_radius("ALA", "CA"),
            Some(1.88),
            "ALA CA should map to C4H1 with radius 1.88"
        );
        assert_eq!(
            get_protor_radius("GLY", "CA"),
            Some(1.88),
            "GLY CA should map to C4H2 with radius 1.88"
        );
        assert_eq!(
            get_protor_radius("TYR", "OH"),
            Some(1.46),
            "TYR OH should map to O2H1 with radius 1.46"
        );
        assert_eq!(
            get_protor_radius("ASN", "CA"),
            PROTOR_RADII
                .get("ASN")
                .and_then(|inner| inner.get("CA"))
                .copied(),
            "Helper function should match direct HashMap access"
        );
        assert_eq!(
            get_protor_radius("CYS", "SG"),
            PROTOR_RADII
                .get("CYS")
                .and_then(|inner| inner.get("SG"))
                .copied(),
            "Helper function should match direct HashMap access"
        );
        assert_eq!(
            get_protor_radius("XXX", "YY"),
            PROTOR_RADII
                .get("XXX")
                .and_then(|inner| inner.get("YY"))
                .copied(),
            "Helper function should return None for unknown combinations"
        );
    }
}
