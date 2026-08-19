//! Steady Stokes channel flow versus the analytic Poiseuille parabola.
//!
//! Run with: `cargo run -p tpt-fem-fluid --example steady_stokes_poiseuille`
//!
//! A unit square channel is driven by a constant body force `G` in `x` with
//! no-slip walls at `y = 0` and `y = 1`. The exact solution of
//! `−μ∇²u + ∇p = f`, `∇·u = 0` is the parabola
//!
//! ```text
//! u_x(y) = (G / 2μ)·y·(1 − y),   u_y = 0,   ∇p = 0
//! ```
//!
//! The parabola (plus `u_y = 0`) is imposed on the boundary only; every interior
//! velocity DOF is free, so recovering the parabola in the interior is a genuine
//! patch test of the divergence-penalty operator. The example asserts the
//! interior profile, the no-slip walls, and that the recovered pressure stays at
//! its (zero) gauge value, then repeats the solve on a refined mesh.

use tpt_fem_fluid::steady_stokes;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Structured `nx` x `ny` Quad4 grid on the unit square.
fn channel(nx: usize, ny: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=ny {
        let y = j as f64 / ny as f64;
        let mut r = Vec::new();
        for i in 0..=nx {
            r.push(b.add_node(vec![i as f64 / nx as f64, y]));
        }
        rows.push(r);
    }
    for j in 0..ny {
        for i in 0..nx {
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

fn on_boundary(mesh: &Mesh, n: usize) -> bool {
    let c = mesh.node_coords(n);
    const TOL: f64 = 1e-9;
    c[0] < TOL || c[0] > 1.0 - TOL || c[1] < TOL || c[1] > 1.0 - TOL
}

/// Analytic profile `u_x(y) = (G/2μ)·y·(1 − y)`.
fn parabola(y: f64, g: f64, mu: f64) -> f64 {
    (g / (2.0 * mu)) * y * (1.0 - y)
}

/// Solve one mesh and return `(max profile error, max |u_y|, max |p|)`.
fn solve(nx: usize, ny: usize, g: f64, mu: f64, penalty: f64, verbose: bool) -> (f64, f64, f64) {
    let mesh = channel(nx, ny);
    let mut bc = Vec::new();
    for n in 0..mesh.node_count() {
        if on_boundary(&mesh, n) {
            let c = mesh.node_coords(n);
            bc.push((n * 2, parabola(c[1], g, mu))); // u_x
            bc.push((n * 2 + 1, 0.0)); // u_y
        }
    }
    let (u, p) = steady_stokes(&mesh, mu, |_| vec![g, 0.0], &bc, penalty);

    // Interior column at x = 0.5 (present whenever nx is even).
    let x_probe = 0.5;
    let mut err = 0.0_f64;
    let mut max_uy = 0.0_f64;
    let mut max_p = 0.0_f64;
    if verbose {
        println!("     y      u_x (FEM)   u_x (exact)     error      u_y        p");
        println!("  ------   ----------   -----------   ---------  ---------  ---------");
    }
    let mut column: Vec<usize> = (0..mesh.node_count())
        .filter(|&n| (mesh.node_coords(n)[0] - x_probe).abs() < 1e-9)
        .collect();
    column.sort_by(|&a, &b| {
        mesh.node_coords(a)[1]
            .partial_cmp(&mesh.node_coords(b)[1])
            .unwrap()
    });
    for n in column {
        let c = mesh.node_coords(n);
        let want = parabola(c[1], g, mu);
        let got = u[n * 2];
        if verbose {
            println!(
                "  {:6.3}   {got:10.6}   {want:11.6}   {:9.2e}  {:9.2e}  {:9.2e}",
                c[1],
                (got - want).abs(),
                u[n * 2 + 1],
                p[n]
            );
        }
        err = err.max((got - want).abs());
    }
    for n in 0..mesh.node_count() {
        max_uy = max_uy.max(u[n * 2 + 1].abs());
        max_p = max_p.max(p[n].abs());
        // No-slip walls must be reproduced exactly (they are prescribed).
        let y = mesh.node_coords(n)[1];
        if !(1e-9..=1.0 - 1e-9).contains(&y) {
            assert!(u[n * 2].abs() < 1e-12, "wall slip u_x = {}", u[n * 2]);
        }
    }
    (err, max_uy, max_p)
}

fn main() {
    let g = 1.0; // body force per unit volume, +x
    let mu = 1.0; // dynamic viscosity
    let penalty = 1e6; // divergence-penalty weight 1/eps

    println!("Poiseuille channel: G = {g}, mu = {mu}, penalty = {penalty:.0e}");
    println!(
        "centreline value u_x(0.5) = G/(8 mu) = {:.6}\n",
        g / (8.0 * mu)
    );

    println!("8 x 8 Quad4 mesh, vertical profile at x = 0.5");
    let (e8, uy8, p8) = solve(8, 8, g, mu, penalty, true);
    println!("\n  max |u_x - exact| = {e8:.3e}   max |u_y| = {uy8:.3e}   max |p| = {p8:.3e}");

    println!("\nmesh refinement (the patch test is exact on every mesh, so the");
    println!("error stays at round-off instead of decaying)");
    println!("   mesh     max error");
    println!("  ------   -----------");
    for n in [4usize, 8, 16] {
        let (e, _, _) = solve(n, n, g, mu, penalty, false);
        println!("  {n:>2} x {n:<2}   {e:11.3e}");
    }

    // The interpolant of the exact parabola is divergence free, so the penalty
    // solve reproduces it to solver accuracy and the pressure keeps its zero
    // gauge value (the flow is driven by the body force, not a pressure drop).
    assert!(e8 < 1e-6, "profile error too large: {e8:.3e}");
    assert!(uy8 < 1e-9, "transverse velocity too large: {uy8:.3e}");
    assert!(p8 < 1e-6, "pressure gauge drifted: {p8:.3e}");
    println!("\nverified: profile error < 1e-6, |u_y| < 1e-9, |p| < 1e-6");
}
