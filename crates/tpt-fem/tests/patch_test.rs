//! End-to-end patch test for the `tpt-fem` core pipeline.
//!
//! Builds a 2-element, 3-node 1-D mesh, computes the element stiffness matrices
//! for the steady 1-D diffusion operator `-u'' = 0` using the reference-element
//! shape functions, isoparametric Jacobian, and Gauss quadrature, assembles the
//! global matrix as triplets, applies Dirichlet boundary conditions, and solves
//! with the sparse LU backend.
//!
//! The exact solution `u(x) = x` is linear, so the P1 finite-element
//! discretisation reproduces it to machine precision — the classic patch test.
#![allow(clippy::needless_range_loop)]

use tpt_fem::{
    gauss_legendre, solve_with_dirichlet, CellType, Coo, Line2, Map, MeshBuilder, ReferenceElement,
};

/// Assemble the global stiffness matrix for `-u'' = 0` on the mesh and solve
/// with the given Dirichlet conditions `(dof, value)`.
fn patch_solve(mesh: &tpt_fem::Mesh, bcs: &[(usize, f64)]) -> Vec<f64> {
    let ndof = mesh.node_count();
    let quad = gauss_legendre(1); // single point at xi = 0, weight 2
    let local_grad = Line2::grad(&[0.0]); // constant: [[-0.5], [0.5]]

    let mut coo = Coo::new();
    for elem in &mesh.elements {
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        // Element stiffness: Ke_ij = sum_q w_q * (dNi/dx)(dNj/dx) * |J|.
        let map = Map::from_nodes_and_grad(&phys, &local_grad);
        let detj = map.determinant;
        let mut ke = [[0.0_f64; 2]; 2];
        for (qi, xi) in quad.points.iter().enumerate() {
            let _ = xi;
            let w = quad.weights[qi];
            let b: Vec<f64> = (0..2)
                .map(|n| map.physical_grad(&local_grad[n])[0])
                .collect();
            for i in 0..2 {
                for j in 0..2 {
                    ke[i][j] += w * b[i] * b[j] * detj;
                }
            }
        }
        for i in 0..2 {
            for j in 0..2 {
                coo.push(elem.nodes[i], elem.nodes[j], ke[i][j]);
            }
        }
    }

    solve_with_dirichlet(&coo, &vec![0.0; ndof], bcs).expect("solve")
}

#[test]
fn patch_test_linear_solution_is_exact() {
    // Two linear elements on [0, 1]: nodes at 0.0, 0.5, 1.0.
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0]);
    let n1 = b.add_node(vec![0.5]);
    let n2 = b.add_node(vec![1.0]);
    b.add_element(CellType::Line, vec![n0, n1]);
    b.add_element(CellType::Line, vec![n1, n2]);
    let mesh = b.build();

    let u = patch_solve(&mesh, &[(n0, 0.0), (n2, 1.0)]);

    // Exact solution u(x) = x.
    assert!((u[n0] - 0.0).abs() < 1e-12);
    assert!((u[n1] - 0.5).abs() < 1e-12);
    assert!((u[n2] - 1.0).abs() < 1e-12);
}

#[test]
fn patch_test_three_element_gradient() {
    // Three linear elements on [0, 1]: nodes at 0, 1/3, 2/3, 1.
    let mut b = MeshBuilder::new();
    let nodes: Vec<usize> = (0..4).map(|i| b.add_node(vec![i as f64 / 3.0])).collect();
    for w in nodes.windows(2) {
        b.add_element(CellType::Line, vec![w[0], w[1]]);
    }
    let mesh = b.build();

    let u = patch_solve(&mesh, &[(nodes[0], 0.0), (nodes[3], 1.0)]);

    for (i, &nid) in nodes.iter().enumerate() {
        assert!(
            (u[nid] - i as f64 / 3.0).abs() < 1e-12,
            "node {i}: got {}",
            u[nid]
        );
    }
}
