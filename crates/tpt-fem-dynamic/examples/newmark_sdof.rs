//! SDOF oscillator integrated with `newmark`, vs the closed form `cos(ωt)`.
//!
//! Run with: `cargo run -p tpt-fem-dynamic --example newmark_sdof`
//!
//! Solves M·ü + K·u = 0 with M=1, K=4 (so ω=2) from u(0)=1, v(0)=0. The exact
//! free response is u(t) = cos(2t); the numerical result is asserted to match
//! it within a stated tolerance, so this example doubles as a verification.

use tpt_fem_dynamic::{newmark, NewmarkOptions};
use tpt_fem_sparse::Coo;

fn main() {
    let m = Coo {
        rows: vec![0],
        cols: vec![0],
        vals: vec![1.0],
    };
    let k = Coo {
        rows: vec![0],
        cols: vec![0],
        vals: vec![4.0],
    };
    let c = Coo::new(); // no damping

    let dt = 0.01;
    let nsteps = 200; // t = 2.0
    let opts = NewmarkOptions {
        dt,
        beta: 0.25,
        gamma: 0.5,
    };
    let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);

    let (t, u) = hist[nsteps].clone();
    let want = (2.0 * t).cos(); // ω = 2
    let err = (u[0] - want).abs();
    println!("newmark (implicit, average-acceleration, β=0.25, γ=0.5)");
    println!("  t            = {t}");
    println!("  u(t)         = {:.6}", u[0]);
    println!("  cos(2t)      = {want:.6}");
    println!("  |error|      = {err:.2e}");
    assert!(err < 1e-3, "got {} want {}", u[0], want);
    println!("\nverified: |u - cos(2t)| < 1e-3");
}
