//! Single Tri3 conductivity (element stiffness) matrix — hand-computed check.
//!
//! For `-∇·(k∇u) = 0` on the unit right triangle with vertices `(0,0)`,
//! `(1,0)`, `(0,1)` and `k = 1`, the constant shape-function gradients are
//! `∇N₁=(-1,-1)`, `∇N₂=(1,0)`, `∇N₃=(0,1)` and the area is `½`. The element
//! matrix is `Kᵢⱼ = k·A·(∇Nᵢ·∇Nⱼ)`, giving the closed-form matrix asserted here.

use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_thermal::poisson_element_matrix;

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    let k = poisson_element_matrix(&mesh, 0, 1.0, 2);

    println!("Tri3 conductivity matrix (k=1):");
    for row in &k {
        println!("  {:8.4} {:8.4} {:8.4}", row[0], row[1], row[2]);
    }

    let expected = [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]];
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (k[i][j] - expected[i][j]).abs() < 1e-12,
                "k[{i}][{j}] = {}",
                k[i][j]
            );
        }
    }
    println!("OK: matches hand-computed k·A·(∇Nᵢ·∇Nⱼ)");
}
