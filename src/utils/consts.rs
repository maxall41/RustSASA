use std::collections::HashSet;
use std::sync::LazyLock;

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
