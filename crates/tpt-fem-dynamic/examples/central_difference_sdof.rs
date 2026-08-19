//! Explicit `central_difference` on the same SDOF oscillator, with the
//! `Δt < 2/ω` stability limit noted.
//!
//! Run with: `cargo run -p tpt-fem-dynamic --example central_difference_sdof`
//!
//! Solves M·ü + K·u = 0 with M=1, K=4 (ω=2) from u(0)=1, v(0)=0. The explicit
//! central-difference scheme is only stable when the time step respects
//!
//! ```text
//! Δt < 2/ω = 2/2 = 1.0
//! ```
//!
//! so we use `dt = 0.002`. The result is checked against `cos(2t)`.

use tpt_fem_dynamic::{central_difference, CentralOptions};
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

    let omega = 2.0;
    let dt = 0.002; // must satisfy dt < 2/omega = 1.0
    let nsteps = 1000; // t = 2.0
    println!("central-difference (explicit), mass lumped internally");
    println!("  ω            = {omega}");
    println!("  stability Δt < {:.3}", 2.0 / omega);
    println!("  chosen Δt    = {dt}");

    let opts = CentralOptions { dt };
    let hist = central_difference(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);

    let (t, u) = hist[nsteps].clone();
    let want = (omega * t).cos();
    let err = (u[0] - want).abs();
    println!("  t            = {t}");
    println!("  u(t)         = {:.6}", u[0]);
    println!("  cos(2t)      = {want:.6}");
    println!("  |error|      = {err:.2e}");
    assert!(err < 5e-3, "got {} want {}", u[0], want);
    println!(
        "\nverified: |u - cos(2t)| < 5e-3 (explicit scheme stayed within the stability limit)"
    );
}
