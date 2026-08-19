//! Augmented-Lagrangian contact vs pure penalty: the same 1-DOF spring problem
//! is solved with penalty (soft, leaves a penetration) and with the augmented
//! Lagrangian iteration (hard, drives penetration to ~0). The multiplier
//! converges to the wall reaction.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-contact --example augmented_lagrangian
//! ```

use tpt_fem_contact::{penalty_contact, ContactConstraint};
use tpt_fem_sparse::{solve, Coo};

fn main() {
    let k = 10.0;
    let penalty = 1e4;
    let load = -4.0;
    let con = ContactConstraint { dof: 0, lower: 0.0 };

    let base = Coo {
        rows: vec![0],
        cols: vec![0],
        vals: vec![k],
    };

    // --- Pure penalty reference ---
    let (k_p, f_p) = penalty_contact(&base, &[load], &[con], penalty);
    let u_pen = solve(&k_p, &f_p).unwrap();
    let pen_penetration = (con.lower - u_pen[0]).max(0.0);

    // --- Augmented Lagrangian iteration (reimplemented to expose history) ---
    let max_iter = 100;
    let tol = 1e-9;
    let mut lambda = 0.0_f64;
    let mut u = 0.0_f64;
    println!("Augmented-Lagrangian iteration (penalty = {penalty})");
    println!("  {:>4} {:>14} {:>14}", "iter", "lambda", "x");
    for it in 0..max_iter {
        // Single DOF: fold the multiplier into the load.
        let f_eff = load + lambda;
        let (k_aug, f_aug) = penalty_contact(&base, &[f_eff], &[con], penalty);
        let u_new = solve(&k_aug, &f_aug).unwrap()[0];
        let viol = con.lower - u_new;
        lambda = (lambda + penalty * viol).max(0.0);
        u = u_new;
        if it < 3 || it % 10 == 0 {
            println!("  {:4} {:14.6e} {:14.6e}", it, lambda, u);
        }
        if viol.abs() < tol {
            break;
        }
    }
    let al_penetration = (con.lower - u).max(0.0);

    println!();
    println!("penalty penetration             = {:.6e}", pen_penetration);
    println!("augmented-Lagrange penetration  = {:.6e}", al_penetration);
    println!(
        "final multiplier lambda         = {:.6e}  (wall reaction ~ 4)",
        lambda
    );

    assert!(
        al_penetration < pen_penetration / 10.0,
        "AL must reduce penetration far below the penalty result"
    );
    assert!(
        (lambda - 4.0).abs() < 1e-2,
        "multiplier must balance the applied load"
    );
    println!("OK: AL enforces hard contact; penalty leaves a small penetration.");
}
