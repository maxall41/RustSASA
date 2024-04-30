use std::collections::HashSet;
use lazy_static::lazy_static;

lazy_static! {
    pub(crate) static ref POLAR_AMINO_ACIDS: HashSet<String> = {
        let mut m = HashSet::new();
        m.insert("SER".to_string());
        m.insert("THR".to_string());
        m.insert("CYS".to_string());
        m.insert("ASN".to_string());
        m.insert("GLN".to_string());
        m.insert("TYR".to_string());
        m
    };
}