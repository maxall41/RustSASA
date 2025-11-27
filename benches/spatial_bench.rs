use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rust_sasa::{calculate_sasa_internal, Atom};
use nalgebra::Point3;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

fn generate_atoms(n: usize) -> Vec<Atom> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut atoms = Vec::with_capacity(n);
    for i in 0..n {
        atoms.push(Atom {
            position: Point3::new(
                rng.gen_range(0.0..100.0),
                rng.gen_range(0.0..100.0),
                rng.gen_range(0.0..100.0),
            ),
            radius: 1.5,
            id: i,
            parent_id: None,
        });
    }
    atoms
}

fn bench_sasa(c: &mut Criterion) {
    let atoms = generate_atoms(5000);
    // Low n_points to emphasize neighbor search cost
    c.bench_function("sasa_5000_atoms_10pts", |b| {
        b.iter(|| {
            calculate_sasa_internal(black_box(&atoms), 1.4, 10, false)
        })
    });
}

criterion_group!(benches, bench_sasa);
criterion_main!(benches);
