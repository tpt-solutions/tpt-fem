//! Benchmarks for the native 3-D tet-mesh generator.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_fem_mesh_gen::{box_mesh, delaunay_3d, Point3};

/// Random points in the unit cube, deterministically seeded.
fn random_points(n: usize) -> Vec<Point3> {
    // Simple LCG so the benchmark is reproducible without pulling in `rand`.
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut pts = Vec::with_capacity(n);
    for _ in 0..n {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let x = (state >> 33) as f64 / (u64::MAX as f64);
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let y = (state >> 33) as f64 / (u64::MAX as f64);
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let z = (state >> 33) as f64 / (u64::MAX as f64);
        pts.push([x, y, z]);
    }
    pts
}

fn bench_delaunay(c: &mut Criterion) {
    let mut group = c.benchmark_group("delaunay_3d");
    for &n in &[200usize, 500, 1000] {
        let pts = random_points(n);
        group.bench_function(format!("points={n}"), |b| {
            b.iter(|| {
                let mesh = delaunay_3d(black_box(&pts));
                black_box(mesh);
            });
        });
    }
    group.finish();
}

fn bench_box_mesh(c: &mut Criterion) {
    let mut group = c.benchmark_group("box_mesh");
    for &n in &[10usize, 20, 40] {
        group.bench_function(format!("n={n}^3"), |b| {
            b.iter(|| {
                let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [n, n, n]);
                black_box(mesh);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_delaunay, bench_box_mesh);
criterion_main!(benches);
