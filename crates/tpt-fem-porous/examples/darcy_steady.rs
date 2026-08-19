//! Steady Darcy seepage through a 1-D column with fixed heads at both ends.
//! The analytical solution for constant permeability k and zero source is a
//! linear pressure profile p(x) = p0 + (p1 - p0) x / L.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-porous --example darcy_steady
//! ```

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_fem_porous::solve_darcy;

fn line_mesh(n: usize, length: f64) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut prev = b.add_node(vec![0.0]);
    for i in 1..=n {
        let node = b.add_node(vec![length * i as f64 / n as f64]);
        b.add_element(CellType::Line, vec![prev, node]);
        prev = node;
    }
    b.build()
}

fn main() {
    let n = 8;
    let l = 1.0;
    let mesh = line_mesh(n, l);
    let k = 2.0;
    // Heads: p = 1.0 at x=0, p = 0.0 at x=L; zero source.
    let p = solve_darcy(&mesh, k, &[], &[(0, 1.0), (n, 0.0)]).unwrap();

    println!("Steady 1-D Darcy flow (k = {k}, L = {l})");
    println!(
        "  {:>6} {:>12} {:>12} {:>8}",
        "x", "p(num)", "p(analyt)", "err"
    );
    let mut max_err = 0.0_f64;
    for node in 0..=n {
        let x = l * node as f64 / n as f64;
        let pa = 1.0 - x; // analytic linear profile
        let err = (p[node] - pa).abs();
        max_err = max_err.max(err);
        println!("  {:6.3} {:12.6} {:12.6} {:8.2e}", x, p[node], pa, err);
    }
    println!();
    println!("max error = {:.3e}", max_err);
    assert!(
        max_err < 1e-6,
        "numerical profile must match the analytic one"
    );
    println!("OK: numerical pressure matches the linear analytic profile.");
}
