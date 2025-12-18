pub mod consts;
pub mod io;
use std::sync::LazyLock;

use crate::utils::consts::PROTOR_RADII;
use fnv::{FnvHashMap, FnvHasher};
use pulp::Arch;
use rayon::ThreadPoolBuilder;
use std::hash::{Hash, Hasher};

pub(crate) static ARCH: LazyLock<Arch> = LazyLock::new(Arch::new);

pub(crate) fn simd_sum(values: &[f32]) -> f32 {
    let mut total = 0f32;
    ARCH.dispatch(|| {
        for x in values {
            total += x;
        }
    });
    total
}

pub fn serialize_chain_id(s: &str) -> isize {
    let mut result = 0;
    for c in s.chars() {
        if c.is_ascii_alphabetic() {
            let position = c.to_ascii_uppercase() as isize - 64;
            result = result * 10 + position;
        }
    }
    result
}

pub fn get_protor_radius(residue: &str, atom: &str) -> Option<f32> {
    PROTOR_RADII.get(residue)?.get(atom).copied()
}

/// Helper function to get atomic radius from custom config or default protor config
pub fn get_radius(
    residue_name: &str,
    atom_name: &str,
    radii_config: Option<&FnvHashMap<String, FnvHashMap<String, f32>>>,
) -> Option<f32> {
    // Check custom config first
    if let Some(config) = radii_config {
        if let Some(radius) = config
            .get(residue_name)
            .and_then(|inner| inner.get(atom_name))
        {
            return Some(*radius);
        }
    }
    // Fall back to default protor config
    get_protor_radius(residue_name, atom_name)
}

/// Configure the global rayon thread pool based on the threads argument
///   - `-1`: Use all available CPU cores (default rayon behavior)
///   - `1`: Single-threaded execution
///   - `> 1`: Use specified number of threads
///   - `0`: Invalid returns error
pub(crate) fn configure_thread_pool(threads: isize) -> Result<(), std::io::Error> {
    if threads == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Thread count must be -1 (all cores) or a positive number",
        ));
    }

    // Configures global thread pool
    if threads > 0 {
        ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build_global()
            .map_err(|e| std::io::Error::other(format!("Failed to configure thread pool: {e}")))?;
    }
    // If threads == -1, use default rayon behavior (all cores)

    Ok(())
}

pub(crate) fn combine_hash(s: &str, n: usize) -> usize {
    let mut hasher = FnvHasher::default();

    // Hash a tuple containing both values
    (s, n).hash(&mut hasher);

    hasher.finish() as usize
}
