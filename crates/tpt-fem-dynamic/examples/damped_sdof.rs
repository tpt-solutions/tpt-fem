//! Damped SDOF response compared to the analytical decay.
//!
//! Run with: `cargo run -p tpt-fem-dynamic --example damped_sdof`
//!
//! Solves M·ü + C·v + K·u = 0 with M=1, K=4 (ω_n = 2), C = 0.4. The damping
//! ratio is ζ = C / (2√(KM)) = 0.4 / 4 = 0.1, giving the closed-form underdamped
//! response
//!
//! ```text
//! u(t) = e^{-ζ·ωₙ·t} · ( cos(ω_d·t) + (ζ·ωₙ/ω_d)·sin(ω_d·t) ),
//!   ω_d = ωₙ·√(1 - ζ²).
//! ```
//!
//! The `newmark` result is asserted against this expression within a tolerance.

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
    let c = Coo {
        rows: vec![0],
        cols: vec![0],
        vals: vec![0.4],
    }; // damping C

    let wn: f64 = 2.0;
    let zeta = 0.4 / (2.0 * wn); // 0.1
    let wd = wn * (1.0 - zeta * zeta).sqrt();

    let dt = 0.005;
    let nsteps = 400; // t = 2.0
    let opts = NewmarkOptions {
        dt,
        beta: 0.25,
        gamma: 0.5,
    };
    let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);

    let (t, u) = hist[nsteps].clone();
    let want = (-zeta * wn * t).exp() * ((wd * t).cos() + zeta * wn / wd * (wd * t).sin());
    let err = (u[0] - want).abs();
    println!("damped SDOF: newmark with C = 0.4 (ζ = {zeta})");
    println!("  t            = {t}");
    println!("  u(t)         = {:.6}", u[0]);
    println!("  closed form  = {want:.6}");
    println!("  |error|      = {err:.2e}");
    assert!(err < 1e-2, "got {} want {}", u[0], want);
    println!("\nverified: |u - analytical damped response| < 1e-2");
}
