//! Verifies and prints the defining properties of Lagrange shape functions:
//! partition of unity everywhere, Kronecker delta at the reference nodes, and
//! quadratic (P2) enrichment via `Tri6`.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-element --example shape_functions
//! ```

use tpt_fem_element::{ReferenceElement, Tri3, Tri6};

fn main() {
    // --- Linear triangle (Tri3) -------------------------------------------
    println!("Tri3 shape functions at the centroid ξ = (1/3, 1/3):");
    let s3 = Tri3::shape(&[1.0 / 3.0, 1.0 / 3.0]);
    for (i, v) in s3.iter().enumerate() {
        println!("  N_{i} = {v:.6}");
    }
    assert!((s3.iter().sum::<f64>() - 1.0).abs() < 1e-12);

    // Partition of unity at a set of interior points.
    for p in [[0.2_f64, 0.3], [0.5, 0.1], [0.25, 0.25]] {
        let sum: f64 = Tri3::shape(&p).iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "partition of unity failed at {p:?}");
    }

    // Kronecker delta: N_j(node i) = δ_ij.
    println!("\nKronecker-delta property (N_j at node i):");
    let nodes3 = Tri3::nodes();
    for (i, n) in nodes3.iter().enumerate() {
        let s = Tri3::shape(n);
        print!("  node {i}: ");
        for (j, v) in s.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((v - expected).abs() < 1e-12);
            print!("{v:.1} ");
        }
        println!();
    }

    // --- Quadratic triangle (Tri6) -----------------------------------------
    println!("\nTri6 (P2) properties:");
    let centre = [1.0_f64 / 3.0, 1.0 / 3.0];
    let s6 = Tri6::shape(&centre);
    println!("  Σ N_i(centroid) = {:.12}", s6.iter().sum::<f64>());
    assert!((s6.iter().sum::<f64>() - 1.0).abs() < 1e-12);

    // A P2 triangle reproduces linear fields exactly: interpolating x with the
    // six nodal values must return x anywhere inside the element.
    let nodes6 = Tri6::nodes();
    let x_nodal: Vec<f64> = nodes6.iter().map(|n| n[0]).collect();
    for p in [[0.15_f64, 0.4], [0.4, 0.2], [1.0 / 3.0, 1.0 / 3.0]] {
        let interp: f64 = Tri6::shape(&p).iter().zip(&x_nodal).map(|(n, x)| n * x).sum();
        assert!(
            (interp - p[0]).abs() < 1e-12,
            "P2 must reproduce linear fields: got {interp}, want {}",
            p[0]
        );
    }
    println!("  P2 interpolation reproduces a linear field exactly.");

    println!("\nAll shape-function checks passed.");
}
