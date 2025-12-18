mod common;

#[cfg(test)]
mod tests {
    use super::common::data::FIXED_LOW_RES_ATOMS;
    use super::common::io::{read_json_result, read_xml_result};
    use approx::assert_abs_diff_eq;
    use assert_cmd::{Command, cargo::*};
    use rust_sasa::SASAResult;
    use std::env;

    #[test]
    fn test_cli_sasa_calc_output_cif() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.cif");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("atom");
        cmd.assert().success();
        assert!(dir.exists());

        let (pdb, _) = pdbtbx::open(dir.to_str().unwrap()).unwrap();
        let atoms: Vec<_> = pdb.atoms().collect();
        assert_eq!(atoms.len(), FIXED_LOW_RES_ATOMS.len());
        for (atom, expected) in atoms.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
            assert_abs_diff_eq!(atom.b_factor() as f32, expected, epsilon = 25.0);
        }
        Ok(())
    }
    #[test]
    fn test_cli_sasa_calc_output_pdb() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.pdb");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("atom");
        cmd.assert().success();
        assert!(dir.exists());

        let (pdb, _) = pdbtbx::open(dir.to_str().unwrap()).unwrap();
        let atoms: Vec<_> = pdb.atoms().collect();
        assert_eq!(atoms.len(), FIXED_LOW_RES_ATOMS.len());
        for (atom, expected) in atoms.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
            assert_abs_diff_eq!(atom.b_factor() as f32, expected, epsilon = 25.0);
        }
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_json() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.json");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("atom");
        cmd.assert().success();
        assert!(dir.exists());

        let result = read_json_result(&dir);
        if let SASAResult::Atom(sasa) = result {
            assert_eq!(sasa.len(), FIXED_LOW_RES_ATOMS.len());
            for (actual, expected) in sasa.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
                assert_abs_diff_eq!(actual, expected, epsilon = 25.0);
            }
        } else {
            panic!("Expected SASAResult::Atom");
        }
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("atom");
        cmd.assert().success();
        assert!(dir.exists());

        let result = read_xml_result(&dir);
        if let SASAResult::Atom(sasa) = result {
            assert_eq!(sasa.len(), FIXED_LOW_RES_ATOMS.len());
            for (actual, expected) in sasa.iter().zip(FIXED_LOW_RES_ATOMS.iter()) {
                assert_abs_diff_eq!(actual, expected, epsilon = 25.0);
            }
        } else {
            panic!("Expected SASAResult::Atom");
        }

        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--probe-radius")
            .arg("1.5")
            .arg("--n-points")
            .arg("200")
            .arg("--include-hetatms")
            .arg("--include-hydrogens")
            .arg("--allow-vdw-fallback");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom_radii_file() -> Result<(), Box<dyn std::error::Error>>
    {
        let mut dir = env::temp_dir();
        dir.push("out.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--probe-radius")
            .arg("1.5")
            .arg("--radii-file")
            .arg("./radii/protor.config")
            .arg("--include-hetatms")
            .arg("--include-hydrogens")
            .arg("--allow-vdw-fallback");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom_threads_and_output_depth_residue()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--threads")
            .arg("1")
            .arg("--output-depth")
            .arg("residue");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom_output_depth_atom()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("atom");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_folder() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("xml");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_json_folder() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("json");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_pdb_folder() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("pdb");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_cif_folder() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("cif");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_folder_customized() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("xml")
            .arg("--probe-radius")
            .arg("1.5")
            .arg("--threads")
            .arg("1")
            .arg("--allow-vdw-fallback")
            .arg("--n-points")
            .arg("200");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom_output_depth_chain()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_chain.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("chain");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_sasa_calc_output_xml_w_custom_output_depth_protein()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_protein.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.clone().into_os_string().into_string().unwrap())
            .arg("--output-depth")
            .arg("protein");
        cmd.assert().success();
        assert!(dir.exists());
        Ok(())
    }

    #[test]
    fn test_cli_input_file_not_found() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_fail.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/non_existent.cif")
            .arg(dir.into_os_string().into_string().unwrap());
        cmd.assert().failure();
        Ok(())
    }

    #[test]
    fn test_cli_input_directory_not_found() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_fail_dir.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/non_existent_dir/")
            .arg(dir.into_os_string().into_string().unwrap())
            .arg("--format")
            .arg("xml");
        cmd.assert().failure();
        Ok(())
    }

    #[test]
    fn test_cli_directory_missing_format() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_fail_format.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/")
            .arg(dir.into_os_string().into_string().unwrap());
        cmd.assert().failure();
        Ok(())
    }

    #[test]
    fn test_cli_invalid_radii_file_path() -> Result<(), Box<dyn std::error::Error>> {
        let mut dir = env::temp_dir();
        dir.push("out_fail_radii.xml");
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.into_os_string().into_string().unwrap())
            .arg("--radii-file")
            .arg("./radii/non_existent.config");
        cmd.assert().failure();
        Ok(())
    }

    #[test]
    fn test_cli_output_is_directory_but_input_is_file() -> Result<(), Box<dyn std::error::Error>> {
        let dir = env::temp_dir();
        // Use the temp dir itself as the output file path, which is a directory
        let mut cmd = Command::new(cargo_bin!("rust-sasa"));
        let cmd = cmd
            .arg("./pdbs/example.cif")
            .arg(dir.into_os_string().into_string().unwrap());
        cmd.assert().failure();
        Ok(())
    }
}
