//! GPU tests. They skip (pass with a notice) when no adapter is available, so
//! they stay green on GPU-less CI runners.

use tpt_fem_gpu::{solve_cg_gpu, GpuCgOptions, GpuContext, GpuError};
use tpt_fem_sparse::{solve_skyline, Coo};

fn context() -> Option<GpuContext> {
    match GpuContext::new() {
        Ok(c) => {
            eprintln!("GPU: {} ({})", c.adapter_name(), c.backend());
            Some(c)
        }
        Err(GpuError::NoAdapter) | Err(GpuError::Device(_)) => {
            eprintln!("no GPU adapter: skipping");
            None
        }
        Err(e) => panic!("{e}"),
    }
}

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

#[test]
fn spmv_matches_cpu() {
    let Some(gpu) = context() else { return };
    let a = laplacian_2d(37).to_csr();
    let x: Vec<f64> = (0..a.nrows).map(|i| (i as f64 * 0.11).sin()).collect();
    let cpu = a.matvec(&x);
    let got = gpu.spmv(&a, &x).unwrap();
    for (c, g) in cpu.iter().zip(&got) {
        assert!((c - g).abs() < 1e-4, "{c} vs {g}");
    }
}

#[test]
fn cg_reaches_double_precision_and_matches_direct_solve() {
    let Some(gpu) = context() else { return };
    let a = laplacian_2d(40);
    let b: Vec<f64> = (0..1600).map(|i| (i as f64 * 0.07).cos() + 2.0).collect();
    let sol = solve_cg_gpu(&gpu, &a, &b, &GpuCgOptions::default()).unwrap();
    assert!(sol.relative_residual <= 1e-10, "{}", sol.relative_residual);
    let direct = solve_skyline(&a, &b).unwrap();
    for (x, y) in sol.x.iter().zip(&direct) {
        assert!((x - y).abs() < 1e-8, "{x} vs {y}");
    }
}

#[test]
fn cg_rejects_bad_input() {
    let Some(gpu) = context() else { return };
    let a = laplacian_2d(4);
    assert!(matches!(
        solve_cg_gpu(&gpu, &a, &[1.0], &GpuCgOptions::default()),
        Err(GpuError::InvalidInput(_))
    ));
    let mut neg = Coo::new();
    neg.push(0, 0, -1.0);
    assert!(matches!(
        solve_cg_gpu(&gpu, &neg, &[1.0], &GpuCgOptions::default()),
        Err(GpuError::InvalidInput(_))
    ));
    let zero = solve_cg_gpu(&gpu, &a, &[0.0; 16], &GpuCgOptions::default()).unwrap();
    assert_eq!(zero.gpu_iterations, 0);
}
