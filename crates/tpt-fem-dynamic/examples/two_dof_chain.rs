//! A 2-DOF spring chain, integrated and checked against the two-mode closed form.
//!
//! Run with: `cargo run -p tpt-fem-dynamic --example two_dof_chain`
//!
//! Masses and springs give
//!
//! ```text
//! M = [[1, 0], [0, 1]],   K = [[3, -1], [-1, 3]].
//! ```
//!
//! The two natural frequencies are ω = √2 and ω = 2, with mode shapes [1,1]
//! and [1,-1]. From u(0)=[1,0], v(0)=[0,0] the exact response is the modal sum
//!
//! ```text
//! u₁(t) = ½·(cos(√2·t) + cos(2·t)),   u₂(t) = ½·(cos(√2·t) - cos(2·t)).
//! ```
//!
//! `newmark` is asserted against this at the final time.

use tpt_fem_dynamic::{newmark, NewmarkOptions};
use tpt_fem_sparse::Coo;

fn main() {
    let mut m = Coo::new();
    let mut k = Coo::new();
    for i in 0..2 {
        m.push(i, i, 1.0);
    }
    k.push(0, 0, 3.0);
    k.push(0, 1, -1.0);
    k.push(1, 0, -1.0);
    k.push(1, 1, 3.0);
    let c = Coo::new();

    let dt = 0.01;
    let nsteps = 200; // t = 2.0
    let opts = NewmarkOptions {
        dt,
        beta: 0.25,
        gamma: 0.5,
    };
    let hist = newmark(
        &m,
        &c,
        &k,
        &[1.0, 0.0],
        &[0.0, 0.0],
        |_| vec![0.0, 0.0],
        &opts,
        nsteps,
    );

    let (t, u) = hist[nsteps].clone();
    let sq2 = 2.0_f64.sqrt();
    let want1 = 0.5 * ((sq2 * t).cos() + (2.0 * t).cos());
    let want2 = 0.5 * ((sq2 * t).cos() - (2.0 * t).cos());
    let err1 = (u[0] - want1).abs();
    let err2 = (u[1] - want2).abs();
    println!("2-DOF chain: M = I, K = [[3,-1],[-1,3]], newmark");
    println!("  t            = {t}");
    println!("  u₁(t), u₂(t) = [{:.6}, {:.6}]", u[0], u[1]);
    println!("  closed form  = [{want1:.6}, {want2:.6}]");
    println!("  |errors|     = [{err1:.2e}, {err2:.2e}]");
    assert!(err1 < 1e-2, "u1 got {} want {}", u[0], want1);
    assert!(err2 < 1e-2, "u2 got {} want {}", u[1], want2);
    println!("\nverified: both DOFs match the two-mode closed form within 1e-2");
}
