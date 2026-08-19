//! Steady 1-D Poisson solve with a known analytic solution.
//!
//! `-u'' = 1` on `[0, 1]` with `u(0) = u(1) = 0` has exact solution
//! `u(x) = ½(x - x²)`, whose midpoint value is `u(½) = 0.125`. The Galerkin
//! solution on a fine line mesh converges to it.

use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_thermal::solve_poisson;

fn main() {
    let nx = 20;
    let mut b = MeshBuilder::new();
    let mut nodes = Vec::new();
    for i in 0..=nx {
        nodes.push(b.add_node(vec![i as f64 / nx as f64]));
    }
    for w in nodes.windows(2) {
        b.add_element(CellType::Line, vec![w[0], w[1]]);
    }
    let mesh = b.build();

    let u = solve_poisson(
        &mesh,
        1.0,
        2,
        |_| 1.0,
        &[(nodes[0], 0.0), (*nodes.last().unwrap(), 0.0)],
        None,
        None,
    )
    .expect("poisson solve");

    let mid = nx / 2;
    let exact = 0.5 * (0.5 - 0.25);
    println!("u(0.5) = {} (exact = {})", u[mid], exact);
    assert!(
        (u[mid] - exact).abs() < 1e-2,
        "got {} expected {}",
        u[mid],
        exact
    );
    println!("OK: converges to the analytic midpoint value");
}
