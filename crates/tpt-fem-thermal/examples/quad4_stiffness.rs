//! Single Quad4 conductivity (element stiffness) matrix: symmetry and rigid-body
//! check.
//!
//! For `-∇·(k∇u) = 0` with constant `k`, the element matrix
//! `Kᵢⱼ = ∫ k ∇Nᵢ·∇Nⱼ dΩ`. Because `Σⱼ ∇Nⱼ = ∇(Σⱼ Nⱼ) = ∇(1) = 0`, every row
//! of `K` sums to zero (rigid-body translation). The matrix is also symmetric —
//! both are hand-derived, exact self-checks.

use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_thermal::poisson_element_matrix;

fn main() {
    let mut b = MeshBuilder::new();
    let n00 = b.add_node(vec![0.0, 0.0]);
    let n10 = b.add_node(vec![1.0, 0.0]);
    let n11 = b.add_node(vec![1.0, 1.0]);
    let n01 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Quad, vec![n00, n10, n11, n01]);
    let mesh = b.build();

    let k = poisson_element_matrix(&mesh, 0, 1.0, 2);

    // Symmetry.
    for i in 0..4 {
        for j in 0..4 {
            assert!((k[i][j] - k[j][i]).abs() < 1e-12, "k[{i}][{j}] asymmetric");
        }
    }
    // Zero row sums (rigid-body translation).
    for i in 0..4 {
        let row: f64 = k[i].iter().sum();
        assert!(row.abs() < 1e-12, "row {i} sum = {row}");
    }

    println!("4x4 Quad4 conductivity matrix (k=1):");
    for row in &k {
        println!(
            "  {:8.4} {:8.4} {:8.4} {:8.4}",
            row[0], row[1], row[2], row[3]
        );
    }
    println!("OK: symmetric and rigid-body (zero row-sum) consistent");
}
