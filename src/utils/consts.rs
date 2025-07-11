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
