//! FE-scale one-sided contact: a bar compressed against a rigid stop.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example contact_column_showcase --features contact
//! ```
//!
//! A bar is discretized with `Line` elements (`BarAxial`, Young interpreted as
//! `EA`), clamped at `x = 0`, assembled with `tpt-fem-assembly`, and pressed
//! against a rigid stop modelled with `tpt-fem-contact` inequality constraints
//! `u_i >= 0` at every free node (the Signorini condition, here in 1-D).
//!
//! Two load cases show the one-sidedness of contact:
//!
//! * **Compression** (axial load pushing toward the clamp): every node *wants*
//!   to move through the stop, so contact engages along the whole bar. The
//!   augmented-Lagrangian iteration drives penetration to ~0 while the
//!   multipliers converge to the total applied load (each becomes the local
//!   stop reaction).
//! * **Tension**: the bar lifts off the stop, no constraint activates, and the
//!   solution must equal the plain elastic solve.

use tpt_fem_assembly::assemble;
use tpt_fem_contact::{augmented_lagrangian, penalty_contact, ContactConstraint};
use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel};
use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_sparse::{solve, Coo};

/// Drop all triplets touching DOF 0 (the clamp) and put 1 on its diagonal.
fn clamp_first_dof(coo: &Coo) -> Coo {
    let mut out = Coo {
        rows: Vec::new(),
        cols: Vec::new(),
        vals: Vec::new(),
    };
    for i in 0..coo.rows.len() {
        if coo.rows[i] != 0 && coo.cols[i] != 0 {
            out.rows.push(coo.rows[i]);
            out.cols.push(coo.cols[i]);
            out.vals.push(coo.vals[i]);
        }
    }
    out.rows.push(0);
    out.cols.push(0);
    out.vals.push(1.0);
    out
}


fn main() {
    let l = 1.0; // bar length, m
    let n_el = 20; // Line elements
    let ea = 1.0e7; // axial stiffness EA, N
    let weight_per_length = 2.0e4; // N/m pushing the bar INTO the floor
    let penalty = 1.0e9;

    // 1-D mesh along x in [0, l].
    let mut b = MeshBuilder::new();
    let nodes: Vec<usize> = (0..=n_el)
        .map(|i| b.add_node(vec![l * i as f64 / n_el as f64]))
        .collect();
    for w in nodes.windows(2) {
        b.add_element(CellType::Line, vec![w[0], w[1]]);
    }
    let mesh = b.build();

    // Global stiffness via the standard assembly operator.
    let k: Coo = assemble(&mesh, 1, |eid, _m| {
        elasticity_element_matrix(&mesh, eid, ElasticModel::BarAxial, ea, 0.0, 2)
            .expect("bar element stiffness")
    });
    let ndof = mesh.node_count();

    // Consistent nodal loads from the distributed load (trapezoid rule):
    // interior nodes carry w*l/n_el, ends carry half. The clamp at dof 0 takes
    // its share into the support, so the stop only feels the free-node load.
    let dl = weight_per_length * l / n_el as f64;
    let compressive: Vec<f64> = (0..ndof)
        .map(|i| -dl * if i == 0 || i == ndof - 1 { 0.5 } else { 1.0 })
        .collect();
    let tensile: Vec<f64> = compressive.iter().map(|f| -f).collect();
    let total: f64 = compressive[1..].iter().sum();

    let cons: Vec<ContactConstraint> = (1..ndof)
        .map(|dof| ContactConstraint { dof, lower: 0.0 })
        .collect();

    // Eliminate the clamp DOF so every solve below is on a nonsingular system.
    let k = clamp_first_dof(&k);

    // ---- Case 1: compression -> contact engages ------------------------
    println!("bar on rigid stop: EA = {ea:.1e} N, w = {weight_per_length:.1e} N/m");
    println!("total load pressing down = {:.4e} N\n", -total);

    let (k_p, f_p) = penalty_contact(&k, &compressive, &cons, penalty);
    let u_pen = solve(&k_p, &f_p).expect("penalty solve");

    // Hard-contact solve via the library's augmented-Lagrangian iteration.
    let tol = 1e-10;
    let (u_al, lambda) =
        augmented_lagrangian(&k, &compressive, &cons, penalty, 200, tol);
    let al_max = u_al[1..].iter().cloned().fold(0.0_f64, f64::min);
    let pen_max = u_pen[1..].iter().cloned().fold(0.0_f64, f64::min);
    let sum_lambda: f64 = lambda.iter().sum();
    println!("  AL converged: max violation below {tol:.0e}\n");
    println!("\n  penalty max penetration  = {pen_max:.6e} m");
    println!("  AL       max penetration = {al_max:.6e} m");
    println!(
        "  sum of multipliers       = {sum_lambda:.6e} N  (total load {:.6e} N)",
        -total
    );
    assert!(pen_max < 0.0, "compression must penetrate under pure penalty");
    assert!(
        al_max.abs() < pen_max.abs() / 10.0,
        "AL penetration must be much smaller than penalty penetration"
    );
    assert!(
        (sum_lambda + total).abs() / (-total) < 1e-3,
        "multipliers must carry the full applied load"
    );

    // ---- Case 2: tension -> lift-off, contact stays silent -------------
    // `penalty_contact` adds its springs unconditionally on constrained DOFs,
    // so the clean way to show one-sidedness is through the multipliers: with
    // the bar pulling AWAY from the stop every violation is negative, the
    // update max(0, ...) keeps every multiplier at zero (no adhesion), and
    // non-penetration holds automatically.
    let (u_lift, lambda_t) = augmented_lagrangian(&k, &tensile, &cons, penalty, 200, tol);
    let sum_lambda_t: f64 = lambda_t.iter().sum();
    let min_gap = u_lift[1..].iter().cloned().fold(f64::INFINITY, f64::min);
    println!("\ntension case (lift-off):");
    println!("  sum of multipliers = {sum_lambda_t:.6e} N (must be ~0: stop carries nothing)");
    println!("  min displacement   = {min_gap:.6e} m (non-penetration u >= 0)");
    assert!(
        sum_lambda_t.abs() < 1e-6 * (-total),
        "inactive contact must carry no load"
    );
    assert!(min_gap > -1e-9, "no penetration may occur under tension");

    // The lifted bar behaves as a plain elastic bar; verify against the closed
    // form u(l) = w l^2 / (2 EA) for internal force N(x) = w (l - x).
    let u_free = solve(&k, &tensile).expect("unconstrained solve");
    let tip_exact = weight_per_length * l * l / (2.0 * ea);
    println!(
        "  tip u = {:.6e} vs closed form wl^2/(2EA) = {tip_exact:.6e}",
        u_free[ndof - 1]
    );
    assert!(
        (u_free[ndof - 1] - tip_exact).abs() / tip_exact < 1e-6,
        "tensile tip deflection must match the closed form"
    );
    println!("\nOK: one-sided contact engages under compression, stays silent in tension.");
}

