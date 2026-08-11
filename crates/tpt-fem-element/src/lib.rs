//! Reference elements, Lagrange shape functions, and the isoparametric
//! Jacobian mapping from reference to physical coordinates.
//!
//! Five linear (`P1`) Lagrange elements are provided:
//!
//! | Element | Spatial dim | Reference domain            |
//! |---------|-------------|----------------------------|
//! | `Line2` | 1           | `[-1, 1]`                  |
//! | `Tri3`  | 2           | `(0,0),(1,0),(0,1)`        |
//! | `Quad4` | 2           | `[-1, 1]²`                 |
//! | `Tet4`  | 3           | `(0,0,0),(1,0,0),(0,1,0),(0,0,1)` |
//! | `Hex8`  | 3           | `[-1, 1]³`                 |
//!
//! Each element exposes its reference-node coordinates, its shape-function
//! values, and the shape-function gradients with respect to the reference
//! coordinates. The [`Map`] type assembles the isoparametric Jacobian from
//! physical node coordinates and the local gradients, and maps local gradients
//! to their physical-space counterparts.
//!
//! Quadratic (`P2`) elements are a tracked follow-up and are intentionally not
//! implemented in this pass.
#![allow(clippy::needless_range_loop)]

use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};

/// A reference finite element.
///
/// Node positions, shape-function values, and gradients are stored in
/// flat `Vec`s (a node coordinate and a gradient are each a length-`DIM`
/// vector) so the trait works uniformly across spatial dimensions without
/// depending on unstable const-generic array features.
pub trait ReferenceElement {
    /// Spatial dimension of the element.
    const DIM: usize;
    /// Number of nodes.
    const NUM_NODES: usize;
    /// Reference-coordinate node positions, one `Vec` per node.
    fn nodes() -> Vec<Vec<f64>>;
    /// Shape-function values `N_0..N_{n-1}` at local coordinates `xi`.
    fn shape(xi: &[f64]) -> Vec<f64>;
    /// Gradient of each shape function w.r.t. the reference coordinates.
    ///
    /// Entry `n` is `[dN_n/dxi_0, ..., dN_n/dxi_{DIM-1}]`.
    fn grad(xi: &[f64]) -> Vec<Vec<f64>>;
}

/// Linear 2-node line element on `[-1, 1]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Line2;
/// Linear 3-node triangle on `(0,0),(1,0),(0,1)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tri3;
/// Linear 4-node quadrilateral on `[-1, 1]²`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Quad4;
/// Linear 4-node tetrahedron on `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tet4;
/// Linear 8-node hexahedron on `[-1, 1]³`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hex8;

impl ReferenceElement for Line2 {
    const DIM: usize = 1;
    const NUM_NODES: usize = 2;
    fn nodes() -> Vec<Vec<f64>> {
        vec![vec![-1.0], vec![1.0]]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let x = xi[0];
        vec![(1.0 - x) / 2.0, (1.0 + x) / 2.0]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let _ = xi;
        vec![vec![-0.5], vec![0.5]]
    }
}

impl ReferenceElement for Tri3 {
    const DIM: usize = 2;
    const NUM_NODES: usize = 3;
    fn nodes() -> Vec<Vec<f64>> {
        vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y) = (xi[0], xi[1]);
        vec![1.0 - x - y, x, y]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let _ = xi;
        vec![vec![-1.0, -1.0], vec![1.0, 0.0], vec![0.0, 1.0]]
    }
}

impl ReferenceElement for Quad4 {
    const DIM: usize = 2;
    const NUM_NODES: usize = 4;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![-1.0, -1.0],
            vec![1.0, -1.0],
            vec![1.0, 1.0],
            vec![-1.0, 1.0],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y) = (xi[0], xi[1]);
        let s = [-1.0, 1.0, 1.0, -1.0];
        let t = [-1.0, -1.0, 1.0, 1.0];
        (0..4)
            .map(|i| 0.25 * (1.0 + s[i] * x) * (1.0 + t[i] * y))
            .collect()
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let (x, y) = (xi[0], xi[1]);
        let s = [-1.0, 1.0, 1.0, -1.0];
        let t = [-1.0, -1.0, 1.0, 1.0];
        (0..4)
            .map(|i| {
                vec![
                    0.25 * s[i] * (1.0 + t[i] * y),
                    0.25 * (1.0 + s[i] * x) * t[i],
                ]
            })
            .collect()
    }
}

impl ReferenceElement for Tet4 {
    const DIM: usize = 3;
    const NUM_NODES: usize = 4;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y, z) = (xi[0], xi[1], xi[2]);
        vec![1.0 - x - y - z, x, y, z]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let _ = xi;
        vec![
            vec![-1.0, -1.0, -1.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ]
    }
}

impl ReferenceElement for Hex8 {
    const DIM: usize = 3;
    const NUM_NODES: usize = 8;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![-1.0, -1.0, -1.0],
            vec![1.0, -1.0, -1.0],
            vec![1.0, 1.0, -1.0],
            vec![-1.0, 1.0, -1.0],
            vec![-1.0, -1.0, 1.0],
            vec![1.0, -1.0, 1.0],
            vec![1.0, 1.0, 1.0],
            vec![-1.0, 1.0, 1.0],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y, z) = (xi[0], xi[1], xi[2]);
        let s = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
        let t = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
        let u = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
        (0..8)
            .map(|i| 0.125 * (1.0 + s[i] * x) * (1.0 + t[i] * y) * (1.0 + u[i] * z))
            .collect()
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let (x, y, z) = (xi[0], xi[1], xi[2]);
        let s = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
        let t = [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0];
        let u = [-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0];
        (0..8)
            .map(|i| {
                vec![
                    0.125 * s[i] * (1.0 + t[i] * y) * (1.0 + u[i] * z),
                    0.125 * (1.0 + s[i] * x) * t[i] * (1.0 + u[i] * z),
                    0.125 * (1.0 + s[i] * x) * (1.0 + t[i] * y) * u[i],
                ]
            })
            .collect()
    }
}

/// Determinant of an `n×n` matrix stored row-major in a flat slice.
fn mat_det(n: usize, m: &[f64]) -> f64 {
    match n {
        1 => m[0],
        2 => m[0] * m[3] - m[1] * m[2],
        3 => {
            m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6])
        }
        _ => panic!("mat_det: unsupported dimension {n}"),
    }
}

/// Inverse of an `n×n` matrix stored row-major, returned row-major.
fn mat_inv(n: usize, m: &[f64]) -> Vec<f64> {
    match n {
        1 => vec![1.0 / m[0]],
        2 => {
            let det = m[0] * m[3] - m[1] * m[2];
            let inv = 1.0 / det;
            vec![m[3] * inv, -m[1] * inv, -m[2] * inv, m[0] * inv]
        }
        3 => {
            let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6]);
            let inv = 1.0 / det;
            vec![
                (m[4] * m[8] - m[5] * m[7]) * inv,
                (m[5] * m[6] - m[3] * m[8]) * inv,
                (m[3] * m[7] - m[4] * m[6]) * inv,
                (m[2] * m[7] - m[1] * m[8]) * inv,
                (m[0] * m[8] - m[2] * m[6]) * inv,
                (m[1] * m[6] - m[0] * m[7]) * inv,
                (m[1] * m[5] - m[2] * m[4]) * inv,
                (m[2] * m[3] - m[0] * m[5]) * inv,
                (m[0] * m[4] - m[1] * m[3]) * inv,
            ]
        }
        _ => panic!("mat_inv: unsupported dimension {n}"),
    }
}

/// The isoparametric mapping from reference to physical coordinates.
///
/// The Jacobian `J` satisfies `J[p][r] = dx_p / dxi_r = Σ_n X_n[p] · (dN_n/dxi_r)`,
/// where `X_n` are the physical node coordinates. [`Map::physical_grad`] maps a
/// shape-function gradient expressed in reference coordinates to physical
/// coordinates via `dN/dx = (dN/dxi) · J⁻¹`.
pub struct Map {
    /// Dimension of the element.
    pub dim: usize,
    /// Jacobian `J` stored row-major (row = physical axis, column = reference axis).
    pub jacobian: Vec<f64>,
    /// `|det J|`, the integration measure.
    pub determinant: f64,
    /// Inverse of `J`, stored row-major.
    pub inverse: Vec<f64>,
}

impl Map {
    /// Assemble the isoparametric Jacobian from physical node coordinates and
    /// the corresponding local shape-function gradients.
    pub fn from_nodes_and_grad(physical: &[Vec<f64>], local_grad: &[Vec<f64>]) -> Self {
        let n = physical[0].len();
        let mut jac = vec![0.0; n * n];
        for node in physical.iter().zip(local_grad) {
            let (x, g) = node;
            for p in 0..n {
                for r in 0..n {
                    jac[p * n + r] += x[p] * g[r];
                }
            }
        }
        let det = mat_det(n, &jac);
        let inverse = mat_inv(n, &jac);
        Map {
            dim: n,
            jacobian: jac,
            determinant: det,
            inverse,
        }
    }

    /// Map a single reference-coordinate gradient to physical coordinates.
    ///
    /// `physical[i] = Σ_r local[r] · inverse[r][i]`.
    pub fn physical_grad(&self, local: &[f64]) -> Vec<f64> {
        let n = self.dim;
        let mut out = vec![0.0; n];
        for i in 0..n {
            let mut s = 0.0;
            for r in 0..n {
                s += local[r] * self.inverse[r * n + i];
            }
            out[i] = s;
        }
        out
    }

    /// Map every local gradient in `local_grads` to physical coordinates.
    pub fn map_gradients(&self, local_grads: &[Vec<f64>]) -> Vec<Vec<f64>> {
        local_grads.iter().map(|g| self.physical_grad(g)).collect()
    }
}

/// 1-D Gauss–Legendre rule on `[-1, 1]` for a line element.
pub fn line_rule(order: usize) -> tpt_fem_quadrature::Quad1D {
    gauss_legendre(order)
}
/// Tensor-product rule on `[-1, 1]²` for a quadrilateral element.
pub fn quad_rule(order: usize) -> tpt_fem_quadrature::Quad2D {
    tensor_square(&gauss_legendre(order))
}
/// Tensor-product rule on `[-1, 1]³` for a hexahedral element.
pub fn hex_rule(order: usize) -> tpt_fem_quadrature::Quad3D {
    tensor_cube(&gauss_legendre(order))
}
/// Rule on the reference triangle.
pub fn tri_rule(rule: TriangleRule) -> tpt_fem_quadrature::Quad2D {
    triangle(rule)
}
/// Rule on the reference tetrahedron.
pub fn tet_rule(rule: TetrahedronRule) -> tpt_fem_quadrature::Quad3D {
    tetrahedron(rule)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < tol)
    }

    #[test]
    fn partition_of_unity_tri() {
        for p in [[0.2_f64, 0.3], [0.5, 0.1], [0.0, 0.0], [0.7, 0.25]] {
            let s = Tri3::shape(&p);
            assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn partition_of_unity_quad_tet_hex() {
        for p in [[0.1_f64, -0.3], [0.4, 0.6], [-0.2, 0.9]] {
            let s = Quad4::shape(&p);
            assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        for p in [[0.1_f64, -0.3, 0.2], [0.4, 0.1, -0.2]] {
            let s = Hex8::shape(&p);
            assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn kronecker_at_nodes() {
        let nodes = Tri3::nodes();
        for (j, n) in nodes.iter().enumerate() {
            let s = Tri3::shape(n);
            for (i, v) in s.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((v - expected).abs() < 1e-12, "tri node {j} shape {i}");
            }
        }
        let nodes = Hex8::nodes();
        for (j, n) in nodes.iter().enumerate() {
            let s = Hex8::shape(n);
            for (i, v) in s.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((v - expected).abs() < 1e-12, "hex node {j} shape {i}");
            }
        }
    }

    #[test]
    fn jacobian_identity_on_reference_tri() {
        let nodes = Tri3::nodes();
        let local = Tri3::grad(&[0.3, 0.2]);
        let m = Map::from_nodes_and_grad(&nodes, &local);
        assert!((m.determinant - 1.0).abs() < 1e-12);
        let id = [1.0, 0.0, 0.0, 1.0];
        assert!(all_close(&m.jacobian, &id, 1e-12));
    }

    #[test]
    fn jacobian_identity_at_center_quad_hex() {
        let qnodes = Quad4::nodes();
        let ql = Quad4::grad(&[0.0, 0.0]);
        let qm = Map::from_nodes_and_grad(&qnodes, &ql);
        assert!(
            (qm.determinant - 1.0).abs() < 1e-12,
            "det={}",
            qm.determinant
        );
        assert!(all_close(&qm.jacobian, &[1.0, 0.0, 0.0, 1.0], 1e-12));

        let hnodes = Hex8::nodes();
        let hl = Hex8::grad(&[0.0, 0.0, 0.0]);
        let hm = Map::from_nodes_and_grad(&hnodes, &hl);
        assert!(
            (hm.determinant - 1.0).abs() < 1e-12,
            "det={}",
            hm.determinant
        );
    }

    #[test]
    fn jacobian_line_scaling() {
        let nodes = [vec![2.0_f64], vec![5.0]];
        let local = Line2::grad(&[0.0]);
        let m = Map::from_nodes_and_grad(&nodes, &local);
        assert!((m.determinant - 1.5).abs() < 1e-12);
    }

    #[test]
    fn jacobian_sheared_tri() {
        // Nodes (1,1),(3,1),(1,4): J = [[2,0],[0,3]], det = 6.
        let nodes = [vec![1.0_f64, 1.0], vec![3.0, 1.0], vec![1.0, 4.0]];
        let local = Tri3::grad(&[0.2, 0.3]);
        let m = Map::from_nodes_and_grad(&nodes, &local);
        assert!((m.determinant - 6.0).abs() < 1e-12, "det={}", m.determinant);
        assert!(all_close(&m.jacobian, &[2.0, 0.0, 0.0, 3.0], 1e-12));
        let g = m.physical_grad(&[1.0, 0.0]);
        assert!(all_close(&g, &[0.5, 0.0], 1e-12));
    }

    #[test]
    fn physical_grad_inverse_consistency() {
        let nodes = [vec![0.0_f64, 0.0], vec![2.0, 0.5], vec![0.3, 1.7]];
        let local = Tri3::grad(&[0.4, 0.2]);
        let m = Map::from_nodes_and_grad(&nodes, &local);
        for p in 0..2 {
            for r in 0..2 {
                let mut s = 0.0;
                for k in 0..2 {
                    s += m.jacobian[p * 2 + k] * m.inverse[k * 2 + r];
                }
                let exp = if p == r { 1.0 } else { 0.0 };
                assert!((s - exp).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn jacobian_3d_sheared_tet() {
        // Map reference tet by scaling axes: X = (2x, 3y, 4z) plus offset.
        let nodes = [
            vec![1.0_f64, 1.0, 1.0],
            vec![3.0, 1.0, 1.0],
            vec![1.0, 4.0, 1.0],
            vec![1.0, 1.0, 5.0],
        ];
        let local = Tet4::grad(&[0.2, 0.2, 0.2]);
        let m = Map::from_nodes_and_grad(&nodes, &local);
        assert!(
            (m.determinant - 24.0).abs() < 1e-12,
            "det={}",
            m.determinant
        );
        assert!(all_close(
            &m.jacobian,
            &[2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0],
            1e-12
        ));
    }
}
