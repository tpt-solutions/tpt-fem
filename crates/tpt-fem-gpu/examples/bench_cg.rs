//! CPU vs GPU conjugate-gradient timing on a 2-D Laplacian.
//!
//! `cargo run --release --example bench_cg -- 300` (grid size m, n = m*m).

use std::time::Instant;
use tpt_fem_gpu::{solve_cg_gpu, GpuCgOptions, GpuContext};
use tpt_fem_sparse::{solve_cg, solve_skyline, CgOptions, Coo};

fn laplacian_2d(m: usize) -> Coo {
    let mut c = Coo::new();
    let id = |i: usize, j: usize| i * m + j;
    for i in 0..m {
        for j in 0..m {
            c.push(id(i, j), id(i, j), 4.0);
            if i + 1 < m {
                c.push(id(i, j), id(i + 1, j), -1.0);
                c.push(id(i + 1, j), id(i, j), -1.0);
            }
            if j + 1 < m {
                c.push(id(i, j), id(i, j + 1), -1.0);
                c.push(id(i, j + 1), id(i, j), -1.0);
            }
        }
    }
    c
}

fn main() {
    let m: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let n = m * m;
    let a = laplacian_2d(m);
    let b = vec![1.0; n];
    println!("2-D Laplacian, n = {n}");

    let gpu = GpuContext::new().expect("GPU adapter");
    println!("GPU: {} ({})", gpu.adapter_name(), gpu.backend());

    let t = Instant::now();
    let s = solve_cg_gpu(&gpu, &a, &b, &GpuCgOptions::default()).expect("gpu cg");
    println!(
        "GPU PCG (f32 + f64 refinement): {:>9.3?}  inner its {:>5}  refinements {}  rel.res {:.2e}",
        t.elapsed(),
        s.gpu_iterations,
        s.refinements,
        s.relative_residual
    );

    let t = Instant::now();
    let c = solve_cg(
        &a,
        &b,
        &CgOptions {
            tol: 1e-10,
            max_iter: 0,
        },
    )
    .expect("cpu cg");
    println!(
        "CPU PCG (f64):                  {:>9.3?}  its {:>5}  rel.res {:.2e}",
        t.elapsed(),
        c.iterations,
        c.relative_residual
    );

    let t = Instant::now();
    let d = solve_skyline(&a, &b).expect("skyline");
    println!("CPU skyline Cholesky (direct):  {:>9.3?}", t.elapsed());

    let err =
        s.x.iter()
            .zip(&d)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
    println!("max |x_gpu - x_direct| = {err:.2e}");
}
