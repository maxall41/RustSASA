use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use pdbtbx::ReadOptions;
use rust_sasa::options::{ResidueLevel, SASAOptions, SASAProcessor};
use rust_sasa::utils::get_radius;

fn load_pdb() -> pdbtbx::PDB {
    let (pdb, _) = ReadOptions::default()
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read("pdbs/example.cif")
        .expect("Failed to load PDB file");
    pdb
}

fn bench_get_radius(c: &mut Criterion) {
    c.bench_function("get_radius", |b| {
        b.iter(|| {
            // Benchmark common amino acid and atom combinations
            black_box(get_radius("ALA", "CA", None));
            black_box(get_radius("GLY", "N", None));
            black_box(get_radius("TRP", "CD2", None));
            black_box(get_radius("ARG", "CZ", None));
            black_box(get_radius("PHE", "CE1", None));
        })
    });
}

fn bench_build_atoms_and_mapping(c: &mut Criterion) {
    let pdb = load_pdb();

    c.bench_function("build_atoms_and_mapping", |b| {
        b.iter(|| {
            black_box(ResidueLevel::build_atoms_and_mapping(
                &pdb, None, false, false, false,
            ))
        })
    });
}

fn bench_full_sasa_residue_level(c: &mut Criterion) {
    let pdb = load_pdb();

    c.bench_function("full_sasa_residue_level", |b| {
        b.iter(|| {
            let options = SASAOptions::<ResidueLevel>::new();
            black_box(options.process(&pdb))
        })
    });
}

criterion_group!(
    benches,
    bench_get_radius,
    bench_build_atoms_and_mapping,
    bench_full_sasa_residue_level
);
criterion_main!(benches);
