//! Tri3 plane-stress element matrix: symmetry and rigid-body check.
//!
//! The element stiffness `K_e = ∫ Bᵀ D B dΩ` must be symmetric. It also admits
//! the two rigid-body translation modes `u = (1,0)` and `u = (0,1)` (zero
//! strain), so `K_e · [1,0,1,0,1,0] = 0` and `K_e · [0,1,0,1,0,1] = 0`. Hence
//! every row of `K_e` sums to zero — a cheap, hand-derived self-check.

use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    let k = elasticity_element_matrix(&mesh, 0, ElasticModel::PlaneStress, 2.0, 0.25, 2)
        .expect("tri3 element matrix");

    // Symmetry.
    for i in 0..6 {
        for j in 0..6 {
            assert!((k[i][j] - k[j][i]).abs() < 1e-12, "k[{i}][{j}] asymmetric");
        }
    }

    // Zero row sums (rigid-body translation).
    for i in 0..6 {
        let row: f64 = k[i].iter().sum();
        assert!(row.abs() < 1e-12, "row {i} sum = {row}");
    }

    println!("6x6 plane-stress Tri3 stiffness (E=2, nu=0.25):");
    for row in &k {
        println!(
            "  {:8.4} {:8.4} {:8.4} {:8.4} {:8.4} {:8.4}",
            row[0], row[1], row[2], row[3], row[4], row[5]
        );
    }
    println!("OK: symmetric and rigid-body (zero row-sum) consistent");
}
