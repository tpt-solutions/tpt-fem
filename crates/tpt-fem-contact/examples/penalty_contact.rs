//! Penalty contact: a node on a spring (k=10) to ground is pushed into a rigid
//! wall at x=0 by load F=-4. The contact force must equal epsilon_N * g_N,
//! where g_N = lower - x is the penetration (positive into the wall).
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-contact --example penalty_contact
//! ```

use tpt_fem_contact::{penalty_contact, ContactConstraint};
use tpt_fem_sparse::{solve, Coo};

fn main() {
    let k = 10.0;
    let penalty = 1e3;
    let load = -4.0;
    let con = ContactConstraint { dof: 0, lower: 0.0 };

    let base = Coo {
        rows: vec![0],
        cols: vec![0],
        vals: vec![k],
    };

    let (k_aug, f_aug) = penalty_contact(&base, &[load], &[con], penalty);
    let u = solve(&k_aug, &f_aug).unwrap();

    // Penetration g_N = lower - x (positive when x < lower).
    let g_n = (con.lower - u[0]).max(0.0);
    // Hand-computed contact force f = epsilon_N * g_N.
    let f_contact = penalty * g_n;
    // Closed-form penetration for this 1-DOF system: x = F / (k + epsilon_N).
    let x_exact = load / (k + penalty);

    println!("Penalty contact (k={k}, epsilon_N={penalty}, F={load})");
    println!("  x (numerical)    = {:.6e}", u[0]);
    println!("  x (analytic)     = {:.6e}", x_exact);
    println!("  penetration g_N  = {:.6e}", g_n);
    println!("  contact force    = {:.6e}", f_contact);
    println!("  added stiffness  = {:.6e} (diagonal increment)", penalty);

    assert!(
        (u[0] - x_exact).abs() < 1e-12,
        "x must match the analytic value"
    );
    assert!(
        u[0] < 0.0 && u[0] > -0.1,
        "node penetrates the wall slightly"
    );
    println!("OK: contact force equals epsilon_N * g_N; node sits just past the wall.");
}
