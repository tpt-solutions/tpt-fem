//! Low-Reynolds-number Navier-Stokes: Poiseuille consistency and a cavity sweep.
//!
//! Run with: `cargo run -p tpt-fem-fluid --example navier_stokes_low_re`
//!
//! `transient_navier_stokes` adds a Picard-linearised convective term to the
//! penalty-Stokes operator and marches the backward-Euler system
//!
//! ```text
//! (M/dt + K + C(u*)) u^{n+1} = f(t^{n+1}) + (M/dt) u^n
//! ```
//!
//! to steady state. Two checks are printed:
//!
//! 1. *Poiseuille consistency* — unidirectional channel flow has `u . grad u = 0`,
//!    so the parabola `u_x(y) = (G/2mu) y (1 - y)` is an exact solution of the
//!    full Navier-Stokes equations at every Reynolds number. The convective term
//!    must therefore leave it untouched; the example asserts agreement with the
//!    closed form to round-off.
//! 2. *Lid-driven cavity sweep* — the same solver on a recirculating flow, where
//!    `u . grad u != 0`, at `Re = U L rho / mu` of 1, 10 and 100. The primary
//!    vortex migrates and the return flow weakens as `Re` grows, which the table
//!    shows against the `steady_stokes` (`Re -> 0`) reference.
//!
//! There is no upwinding or SUPG stabilisation, so the sweep is limited to
//! moderate `Re` on a coarse mesh.

use tpt_fem_fluid::{steady_stokes, transient_navier_stokes};
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Structured `n` x `n` Quad4 grid on the unit square.
fn unit_square(n: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=n {
        let y = j as f64 / n as f64;
        let mut r = Vec::new();
        for i in 0..=n {
            r.push(b.add_node(vec![i as f64 / n as f64, y]));
        }
        rows.push(r);
    }
    for j in 0..n {
        for i in 0..n {
            b.add_element(
                CellType::Quad,
                vec![
                    rows[j][i],
                    rows[j][i + 1],
                    rows[j + 1][i + 1],
                    rows[j + 1][i],
                ],
            );
        }
    }
    b.build()
}

fn is_boundary(mesh: &Mesh, nd: usize) -> bool {
    let c = mesh.node_coords(nd);
    const TOL: f64 = 1e-9;
    c[0] < TOL || c[0] > 1.0 - TOL || c[1] < TOL || c[1] > 1.0 - TOL
}

/// Index of the node at `(x, y)`.
fn node_at(mesh: &Mesh, x: f64, y: f64) -> usize {
    (0..mesh.node_count())
        .find(|&nd| {
            (mesh.node_coords(nd)[0] - x).abs() < 1e-9 && (mesh.node_coords(nd)[1] - y).abs() < 1e-9
        })
        .expect("probe node must exist")
}

fn poiseuille_consistency() -> f64 {
    let n = 6;
    let mesh = unit_square(n);
    let mu = 1.0;
    let penalty = 1e3;
    let probe = node_at(&mesh, 0.5, 0.5);

    println!("1) Poiseuille consistency (u . grad u = 0, so the parabola is exact)");
    println!("      G       Re      u_x N-S      u_x Stokes    G/(8 mu)     |N-S - exact|");
    println!("   ------  -------  -----------   -----------   ----------   -------------");
    let mut worst = 0.0_f64;
    for g in [0.01, 0.1, 1.0] {
        let mut bc = Vec::new();
        for nd in 0..mesh.node_count() {
            if is_boundary(&mesh, nd) {
                let c = mesh.node_coords(nd);
                bc.push((nd * 2, (g / (2.0 * mu)) * c[1] * (1.0 - c[1])));
                bc.push((nd * 2 + 1, 0.0));
            }
        }
        let u_ns =
            transient_navier_stokes(&mesh, mu, |_, _| vec![g, 0.0], &bc, penalty, 0.01, 60, 3);
        let (u_st, _p) = steady_stokes(&mesh, mu, |_| vec![g, 0.0], &bc, penalty);
        let exact = g / (8.0 * mu);
        let re = exact / mu; // rho = L = 1, U = peak velocity
        let err = (u_ns[probe * 2] - exact).abs();
        println!(
            "   {g:6.2}  {re:7.4}  {:11.6}   {:11.6}   {exact:10.6}   {err:13.3e}",
            u_ns[probe * 2],
            u_st[probe * 2]
        );
        worst = worst.max(err / exact);
    }
    worst
}

fn cavity_sweep() {
    let n = 8;
    let mesh = unit_square(n);
    let penalty = 1e3;
    let u_lid = 1.0;

    // Lid at y = 1 sliding at u_x = u_lid, no-slip elsewhere.
    let mut bc = Vec::new();
    for nd in 0..mesh.node_count() {
        if is_boundary(&mesh, nd) {
            let lid = mesh.node_coords(nd)[1] > 1.0 - 1e-9;
            bc.push((nd * 2, if lid { u_lid } else { 0.0 }));
            bc.push((nd * 2 + 1, 0.0));
        }
    }

    // Centreline column x = 0.5, sorted by height.
    let mut column: Vec<usize> = (0..mesh.node_count())
        .filter(|&nd| (mesh.node_coords(nd)[0] - 0.5).abs() < 1e-9)
        .collect();
    column.sort_by(|&a, &b| {
        mesh.node_coords(a)[1]
            .partial_cmp(&mesh.node_coords(b)[1])
            .unwrap()
    });
    let extremum = |u: &[f64]| -> (f64, f64) {
        let mut best = (0.0, 0.0);
        for &nd in &column {
            if u[nd * 2] < best.1 {
                best = (mesh.node_coords(nd)[1], u[nd * 2]);
            }
        }
        best
    };

    println!("\n2) Lid-driven cavity sweep, {n} x {n} Quad4, u_lid = {u_lid}");
    println!("   (return-flow extremum on the centreline x = 0.5)");
    println!("      Re      mu       y_min    u_x,min     vs Stokes");
    println!("   ------  -------   -------  ---------   -----------");
    let (u_stokes, _p) = steady_stokes(&mesh, 1.0, |_| vec![0.0, 0.0], &bc, penalty);
    let (y_st, ux_st) = extremum(&u_stokes);
    println!("      ->0  1.0e0    {y_st:7.3}  {ux_st:9.6}   (reference)");
    let mut last_delta = 0.0;
    for mu in [1.0, 0.1, 0.01] {
        let re = u_lid / mu; // rho = L = 1
        let u_ns =
            transient_navier_stokes(&mesh, mu, |_, _| vec![0.0, 0.0], &bc, penalty, 0.05, 40, 3);
        let (y_ns, ux_ns) = extremum(&u_ns);
        last_delta = (ux_ns - ux_st).abs();
        println!(
            "   {re:6.0}  {mu:7.1e}   {y_ns:7.3}  {ux_ns:9.6}   {:+11.6}",
            ux_ns - ux_st
        );
    }
    // At Re = 100 convection must have measurably distorted the Stokes field.
    assert!(
        last_delta > 1e-3,
        "convection should alter the cavity flow at Re = 100, delta = {last_delta:.3e}"
    );
    println!("\n   convective distortion at Re = 100: {last_delta:.3e} (non-zero)");
}

fn main() {
    let worst = poiseuille_consistency();
    println!("\n   worst relative deviation from the closed form: {worst:.3e}");
    assert!(
        worst < 1e-8,
        "Poiseuille must be reproduced exactly by the N-S solver, got {worst:.3e}"
    );
    cavity_sweep();
    println!("\nverified: Poiseuille reproduced to better than 1e-8 relative, and the");
    println!("cavity flow responds to convection as Re grows");
}
