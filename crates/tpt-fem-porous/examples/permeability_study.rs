//! Parameter study for steady Darcy flow: the pressure field for a fixed head
//! difference is independent of permeability k, while the Darcy flux
//! q = -k ∇p scales linearly with k. This example confirms that behaviour.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-porous --example permeability_study
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
    let ks = [0.5_f64, 1.0, 2.0, 4.0];

    println!("Darcy parameter study (head p0=1 at x=0, p1=0 at x=L)");
    println!("  {:>8} {:>12} {:>10}", "k", "flux q", "q/k");
    for &k in &ks {
        let p = solve_darcy(&mesh, k, &[], &[(0, 1.0), (n, 0.0)]).unwrap();
        // Gradient via the first element (linear solution => constant).
        let dx = l / n as f64;
        let grad = (p[1] - p[0]) / dx;
        let q = -k * grad; // Darcy flux per unit area
        println!("  {:8.3} {:12.6} {:10.4}", k, q, q / k);
        // Flux must scale linearly with k: q/k == 1 (since |grad| = 1).
        assert!((q / k - 1.0).abs() < 1e-6, "flux must equal k*|grad|");
    }
    println!();
    println!("OK: pressure profile is k-independent; flux q = k|∇p| ∝ k.");
}
