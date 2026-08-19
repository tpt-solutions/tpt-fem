//! Transient creeping flow relaxing onto the steady Stokes solution.
//!
//! Run with: `cargo run -p tpt-fem-fluid --example transient_stokes_startup`
//!
//! Unsteady Stokes start-up flow: a channel at rest is loaded at `t = 0` by a
//! constant body force `G` in `x`. The momentum balance
//!
//! ```text
//! rho du/dt = mu grad^2 u - grad p + G,   div u = 0
//! ```
//!
//! is first order in time, so the velocity relaxes monotonically onto the steady
//! Poiseuille parabola `u_x(y) = (G/2mu)*y*(1 - y)` with a viscous time scale of
//! order `rho L^2 / mu`. `transient_stokes` integrates it with
//! `tpt-fem-dynamic`'s implicit Newmark scheme (the viscous + penalty operator
//! acts as the damping matrix) and returns the velocity history.
//!
//! The example prints the centre-node time history against the steady value and
//! asserts that (a) the history is monotone while it rises and (b) the final
//! velocity matches the analytic parabola to better than 1%.

use tpt_fem_dynamic::NewmarkOptions;
use tpt_fem_fluid::{steady_stokes, transient_stokes};
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

fn main() {
    let n = 6;
    let mesh = unit_square(n);
    let g = 1.0; // body force per unit volume, +x
    let mu = 0.25; // viscosity (lower mu -> longer viscous time scale)
    let penalty = 1e4; // divergence-penalty weight 1/eps
    let exact_centre = g / (8.0 * mu); // (G/2mu) * 0.5 * 0.5

    // Steady Poiseuille parabola on the boundary, u_y = 0 everywhere on it.
    let mut bc = Vec::new();
    for nd in 0..mesh.node_count() {
        let c = mesh.node_coords(nd);
        const TOL: f64 = 1e-9;
        if c[0] < TOL || c[0] > 1.0 - TOL || c[1] < TOL || c[1] > 1.0 - TOL {
            bc.push((nd * 2, (g / (2.0 * mu)) * c[1] * (1.0 - c[1])));
            bc.push((nd * 2 + 1, 0.0));
        }
    }

    let probe = (0..mesh.node_count())
        .find(|&nd| {
            let c = mesh.node_coords(nd);
            (c[0] - 0.5).abs() < 1e-9 && (c[1] - 0.5).abs() < 1e-9
        })
        .expect("centre node");

    let (u_steady, _p) = steady_stokes(&mesh, mu, |_| vec![g, 0.0], &bc, penalty);
    let steady_centre = u_steady[probe * 2];

    let dt = 0.01;
    let nsteps = 40;
    let opts = NewmarkOptions {
        dt,
        beta: 0.25,
        gamma: 0.5,
    };
    let hist = transient_stokes(&mesh, mu, |_, _| vec![g, 0.0], &bc, penalty, &opts, nsteps);

    println!("Stokes start-up flow: {n} x {n} Quad4, G = {g}, mu = {mu}, penalty = {penalty:.0e}");
    println!("steady_stokes centre u_x = {steady_centre:.6}");
    println!("analytic  G/(8 mu)       = {exact_centre:.6}\n");
    println!("   step      t        u_x(0.5, 0.5)   u_x / steady   |u_x - steady|");
    println!("  ------  -------   --------------   ------------   --------------");
    for (i, (t, u)) in hist.iter().enumerate() {
        if i % 4 == 0 || i == hist.len() - 1 {
            let ux = u[probe * 2];
            println!(
                "  {i:>6}  {t:7.3}   {ux:14.6}   {:12.4}   {:14.3e}",
                ux / steady_centre,
                (ux - steady_centre).abs()
            );
        }
    }

    // The response rises monotonically while it is still climbing towards the
    // steady value (the tail settles onto the discrete steady state from above,
    // ~0.2% higher than the continuum value, so only the rise is monotone).
    for i in 1..=nsteps {
        let prev = hist[i - 1].1[probe * 2];
        let cur = hist[i].1[probe * 2];
        if prev < 0.98 * steady_centre {
            assert!(
                cur >= prev - 1e-12,
                "history must rise monotonically at step {i}: {prev} -> {cur}"
            );
        }
    }

    let final_ux = hist[nsteps].1[probe * 2];
    let rel = (final_ux - exact_centre).abs() / exact_centre;
    println!("\nfinal t = {:.3}", hist[nsteps].0);
    println!("  u_x(0.5, 0.5)        = {final_ux:.6}");
    println!("  analytic steady      = {exact_centre:.6}");
    println!("  relative difference  = {:.3e}", rel);
    assert!(
        rel < 1e-2,
        "transient did not reach the steady value: {rel:.3e}"
    );
    println!("\nverified: monotone rise and final velocity within 1% of G/(8 mu)");
}
