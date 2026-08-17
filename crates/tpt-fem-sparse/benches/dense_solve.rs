//! Benchmarks for the default (dense LU) solve backend.
//!
//! The `russell` sparse-direct backend is unavailable in CI (it needs the
//! SuiteSparse/MUMPS Fortran/C toolchain), so these benchmarks cover the
//! in-house dense `solve` at a few problem sizes as a baseline. Swap in
//! `solve_russell` once that toolchain is present to compare sparse vs dense.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_fem_sparse::Coo;

/// 2-D 5-point Laplacian on an `n × n` grid, assembled as a [`Coo`].
fn laplacian_2d(n: usize) -> Coo {
    let mut c = Coo::new();
    let idx = |i: usize, j: usize| i * n + j;
    for i in 0..n {
        for j in 0..n {
            let k = idx(i, j);
            c.push(k, k, 4.0);
            if i > 0 {
                c.push(k, idx(i - 1, j), -1.0);
            }
            if i + 1 < n {
                c.push(k, idx(i + 1, j), -1.0);
            }
            if j > 0 {
                c.push(k, idx(i, j - 1), -1.0);
            }
            if j + 1 < n {
                c.push(k, idx(i, j + 1), -1.0);
            }
        }
    }
    c
}

fn bench_dense_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("dense_solve");
    for &n in &[50usize, 100, 200] {
        let coo = laplacian_2d(n);
        let rhs = vec![1.0; n * n];
        group.bench_function(format!("n={n}"), |b| {
            b.iter(|| {
                let sol = tpt_fem_sparse::solve(black_box(&coo), black_box(&rhs)).unwrap();
                black_box(sol);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_dense_solve);
criterion_main!(benches);
