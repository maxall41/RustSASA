use std::collections::HashSet;
use std::sync::LazyLock;

use fnv::FnvHashMap;

pub(crate) static POLAR_AMINO_ACIDS: LazyLock<HashSet<String>> = LazyLock::new(|| {
    let mut m = HashSet::new();
    m.insert("SER".to_string());
    m.insert("THR".to_string());
    m.insert("CYS".to_string());
    m.insert("ASN".to_string());
    m.insert("GLN".to_string());
    m.insert("TYR".to_string());
    m
});

pub(crate) const GOLDEN_RATIO: f32 = 1.618_034;
pub(crate) const ANGLE_INCREMENT: f32 = 2.0 * std::f32::consts::PI * GOLDEN_RATIO;

/// Macro to load ProtOr radii config at compile time and parse lazily
macro_rules! load_protor_radii {
    () => {
        LazyLock::new(|| {
            let config = include_str!("../../radii/protor.config");
            parse_radii_config(config)
        })
    };
}

pub(crate) fn parse_radii_config(content: &str) -> FnvHashMap<String, FnvHashMap<String, f32>> {
    let mut types: FnvHashMap<String, f32> = FnvHashMap::default();
    let mut atoms: FnvHashMap<String, FnvHashMap<String, f32>> = FnvHashMap::default();

    let mut in_types = false;
    let mut in_atoms = false;

    for line in content.lines() {
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') || line.starts_with("name:") {
            continue;
        }

        // Track which section we're in
        if line == "types:" {
            in_types = true;
            in_atoms = false;
            continue;
        }
        if line == "atoms:" {
            in_types = false;
            in_atoms = true;
            continue;
        }

        if in_types {
            // Format: TYPE_NAME RADIUS POLARITY
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let Ok(radius) = parts[1].parse::<f32>() {
                    types.insert(parts[0].to_string(), radius);
                }
            }
        } else if in_atoms {
            // Format: RESIDUE ATOM_NAME TYPE_NAME
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                if let Some(&radius) = types.get(parts[2]) {
                    #[allow(clippy::unwrap_or_default)]
                    atoms
                        .entry(parts[0].to_string())
                        .or_insert_with(FnvHashMap::default)
                        .insert(parts[1].to_string(), radius);
                }
            }
        }
    }

    atoms
}

pub(crate) fn load_radii_from_file(
    path: &str,
) -> Result<FnvHashMap<String, FnvHashMap<String, f32>>, std::io::Error> {
    let content = std::fs::read_to_string(path)?;
    Ok(parse_radii_config(&content))
}

pub(crate) static PROTOR_RADII: LazyLock<FnvHashMap<String, FnvHashMap<String, f32>>> =
    load_protor_radii!();
