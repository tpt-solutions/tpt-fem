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
//! Quadratic (`P2`) elements are provided as a follow-up: `Tri6`, `Quad8`
//! (serendipity), `Quad9` (biquadratic Lagrange), `Tet10`, `Hex20`
//! (serendipity) and `Hex27` (triquadratic Lagrange).
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

/// Quadratic 6-node triangle (P2 Lagrange) on the same reference domain as
/// `Tri3`: vertices `(0,0)`, `(1,0)`, `(0,1)` plus edge midpoints
/// `(0.5,0)`, `(0,0.5)`, `(0.5,0.5)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tri6;
/// Quadratic 9-node quadrilateral (biquadratic Lagrange) on `[-1, 1]²`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Quad9;
/// Quadratic 8-node serendipity quadrilateral on `[-1, 1]²` (no centre node).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Quad8;
/// Quadratic 10-node tetrahedron (P2 Lagrange) on the same reference domain as
/// `Tet4`: the four vertices plus the six edge midpoints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tet10;
/// Quadratic 27-node hexahedron (triquadratic Lagrange) on `[-1, 1]³`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hex27;
/// Quadratic 20-node serendipity hexahedron on `[-1, 1]³` (corners + edge
/// midpoints; no face or centre nodes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hex20;

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

/// 1-D quadratic Lagrange basis value at node position `at` ∈ {-1, 0, 1},
/// evaluated at `x`.
fn lag2(x: f64, at: f64) -> f64 {
    if at == -1.0 {
        x * (x - 1.0) / 2.0
    } else if at == 0.0 {
        1.0 - x * x
    } else {
        debug_assert_eq!(at, 1.0);
        x * (x + 1.0) / 2.0
    }
}

/// Derivative of [`lag2`] with respect to `x`.
fn dlag2(x: f64, at: f64) -> f64 {
    if at == -1.0 {
        (2.0 * x - 1.0) / 2.0
    } else if at == 0.0 {
        -2.0 * x
    } else {
        debug_assert_eq!(at, 1.0);
        (2.0 * x + 1.0) / 2.0
    }
}

/// Reference node list for `Quad9`/`Hex27`: the tensor product of
/// `{-1, 0, 1}` over the element's spatial dimension.
///
/// # Panics
///
/// Panics (in debug and release builds) if `dim` is not `2` or `3`. This is a
/// contracted precondition: the only reference elements that call
/// `tensor_nodes` are `Quad9`/`Hex27`, whose spatial dimension is always 2/3.
/// See also the Phase 9b / 10a hardening notes in `todo.md`; converting this
/// to a `Result` was intentionally deferred because the public `Map` /
/// `from_nodes_and_grad` API would otherwise ripple for negligible gain.
fn tensor_nodes(dim: usize) -> Vec<Vec<f64>> {
    debug_assert!(
        dim == 2 || dim == 3,
        "tensor_nodes: unsupported dimension {dim} (only 2 or 3 are valid)"
    );
    let levels = [-1.0_f64, 0.0, 1.0];
    let mut v = Vec::new();
    match dim {
        2 => {
            for &a in &levels {
                for &b in &levels {
                    v.push(vec![a, b]);
                }
            }
        }
        3 => {
            for &a in &levels {
                for &b in &levels {
                    for &c in &levels {
                        v.push(vec![a, b, c]);
                    }
                }
            }
        }
        _ => unreachable!("tensor_nodes: unsupported dimension {dim}"),
    }
    v
}

impl ReferenceElement for Tri6 {
    const DIM: usize = 2;
    const NUM_NODES: usize = 6;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![0.5, 0.0],
            vec![0.5, 0.5],
            vec![0.0, 0.5],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        // Barycentric (area) coordinates on the right-triangle reference domain.
        let l1 = 1.0 - xi[0] - xi[1];
        let l2 = xi[0];
        let l3 = xi[1];
        vec![
            l1 * (2.0 * l1 - 1.0), // vertex 1
            l2 * (2.0 * l2 - 1.0), // vertex 2
            l3 * (2.0 * l3 - 1.0), // vertex 3
            4.0 * l1 * l2,         // mid 1-2
            4.0 * l2 * l3,         // mid 2-3
            4.0 * l1 * l3,         // mid 1-3
        ]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let l1 = 1.0 - xi[0] - xi[1];
        let l2 = xi[0];
        let l3 = xi[1];
        // dN/dL_i for each node: (dL1, dL2, dL3).
        let dl = [
            [4.0 * l1 - 1.0, 0.0, 0.0],
            [0.0, 4.0 * l2 - 1.0, 0.0],
            [0.0, 0.0, 4.0 * l3 - 1.0],
            [4.0 * l2, 4.0 * l1, 0.0],
            [0.0, 4.0 * l3, 4.0 * l2],
            [4.0 * l3, 0.0, 4.0 * l1],
        ];
        dl.iter()
            .map(|d| {
                // dL1/dx = -1, dL2/dx = 1, dL3/dx = 0; dN/dx = -dL1 + dL2.
                // dL1/dy = -1, dL2/dy = 0, dL3/dy = 1; dN/dy = -dL1 + dL3.
                vec![-d[0] + d[1], -d[0] + d[2]]
            })
            .collect()
    }
}

impl ReferenceElement for Quad9 {
    const DIM: usize = 2;
    const NUM_NODES: usize = 9;
    fn nodes() -> Vec<Vec<f64>> {
        tensor_nodes(2)
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let ns = Self::nodes();
        ns.iter()
            .map(|n| lag2(xi[0], n[0]) * lag2(xi[1], n[1]))
            .collect()
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let ns = Self::nodes();
        ns.iter()
            .map(|n| {
                vec![
                    dlag2(xi[0], n[0]) * lag2(xi[1], n[1]),
                    lag2(xi[0], n[0]) * dlag2(xi[1], n[1]),
                ]
            })
            .collect()
    }
}

impl ReferenceElement for Quad8 {
    const DIM: usize = 2;
    const NUM_NODES: usize = 8;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![-1.0, -1.0],
            vec![1.0, -1.0],
            vec![1.0, 1.0],
            vec![-1.0, 1.0],
            vec![0.0, -1.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![-1.0, 0.0],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y) = (xi[0], xi[1]);
        // Corners 0-3: 1/4 (1+a x)(1+b y)(a x + b y - 1).
        // Mids 4-7: bottom/right/top/left serendipity bubbles.
        let corner = |a: f64, b: f64| -> f64 {
            0.25 * (1.0 + a * x) * (1.0 + b * y) * (a * x + b * y - 1.0)
        };
        vec![
            corner(-1.0, -1.0),
            corner(1.0, -1.0),
            corner(1.0, 1.0),
            corner(-1.0, 1.0),
            0.5 * (1.0 - x * x) * (1.0 - y), // bottom mid
            0.5 * (1.0 - y * y) * (1.0 + x), // right mid
            0.5 * (1.0 - x * x) * (1.0 + y), // top mid
            0.5 * (1.0 - y * y) * (1.0 - x), // left mid
        ]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let (x, y) = (xi[0], xi[1]);
        let corner = |a: f64, b: f64| -> [f64; 2] {
            let aa = 1.0 + a * x;
            let bb = 1.0 + b * y;
            let cc = a * x + b * y - 1.0;
            // dN/dx = a/4 * bb * (cc + aa); dN/dy = b/4 * aa * (cc + bb).
            [0.25 * a * bb * (cc + aa), 0.25 * b * aa * (cc + bb)]
        };
        vec![
            corner(-1.0, -1.0).to_vec(),
            corner(1.0, -1.0).to_vec(),
            corner(1.0, 1.0).to_vec(),
            corner(-1.0, 1.0).to_vec(),
            vec![-x * (1.0 - y), -0.5 * (1.0 - x * x)], // bottom mid
            vec![0.5 * (1.0 - y * y), -y * (1.0 + x)],  // right mid
            vec![-x * (1.0 + y), 0.5 * (1.0 - x * x)],  // top mid
            vec![-0.5 * (1.0 - y * y), -y * (1.0 - x)], // left mid
        ]
    }
}

impl ReferenceElement for Tet10 {
    const DIM: usize = 3;
    const NUM_NODES: usize = 10;
    fn nodes() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![0.5, 0.0, 0.0], // mid 1-2
            vec![0.5, 0.5, 0.0], // mid 2-3
            vec![0.0, 0.5, 0.0], // mid 1-3
            vec![0.0, 0.0, 0.5], // mid 1-4
            vec![0.5, 0.0, 0.5], // mid 2-4
            vec![0.0, 0.5, 0.5], // mid 3-4
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let l1 = 1.0 - xi[0] - xi[1] - xi[2];
        let l2 = xi[0];
        let l3 = xi[1];
        let l4 = xi[2];
        vec![
            l1 * (2.0 * l1 - 1.0), // vertex 1
            l2 * (2.0 * l2 - 1.0), // vertex 2
            l3 * (2.0 * l3 - 1.0), // vertex 3
            l4 * (2.0 * l4 - 1.0), // vertex 4
            4.0 * l1 * l2,         // mid 1-2
            4.0 * l2 * l3,         // mid 2-3
            4.0 * l1 * l3,         // mid 1-3
            4.0 * l1 * l4,         // mid 1-4
            4.0 * l2 * l4,         // mid 2-4
            4.0 * l3 * l4,         // mid 3-4
        ]
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let l1 = 1.0 - xi[0] - xi[1] - xi[2];
        let l2 = xi[0];
        let l3 = xi[1];
        let l4 = xi[2];
        let dl = [
            [4.0 * l1 - 1.0, 0.0, 0.0, 0.0],
            [0.0, 4.0 * l2 - 1.0, 0.0, 0.0],
            [0.0, 0.0, 4.0 * l3 - 1.0, 0.0],
            [0.0, 0.0, 0.0, 4.0 * l4 - 1.0],
            [4.0 * l2, 4.0 * l1, 0.0, 0.0],
            [0.0, 4.0 * l3, 4.0 * l2, 0.0],
            [4.0 * l3, 0.0, 4.0 * l1, 0.0],
            [4.0 * l4, 0.0, 0.0, 4.0 * l1],
            [0.0, 4.0 * l4, 0.0, 4.0 * l2],
            [0.0, 0.0, 4.0 * l4, 4.0 * l3],
        ];
        dl.iter()
            .map(|d| {
                // dL1 = (-1,-1,-1), dL2 = (1,0,0), dL3 = (0,1,0), dL4 = (0,0,1).
                // dN/dx = -dL1 + dL2; dN/dy = -dL1 + dL3; dN/dz = -dL1 + dL4.
                vec![-d[0] + d[1], -d[0] + d[2], -d[0] + d[3]]
            })
            .collect()
    }
}

impl ReferenceElement for Hex27 {
    const DIM: usize = 3;
    const NUM_NODES: usize = 27;
    fn nodes() -> Vec<Vec<f64>> {
        tensor_nodes(3)
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let ns = Self::nodes();
        ns.iter()
            .map(|n| lag2(xi[0], n[0]) * lag2(xi[1], n[1]) * lag2(xi[2], n[2]))
            .collect()
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let ns = Self::nodes();
        ns.iter()
            .map(|n| {
                vec![
                    dlag2(xi[0], n[0]) * lag2(xi[1], n[1]) * lag2(xi[2], n[2]),
                    lag2(xi[0], n[0]) * dlag2(xi[1], n[1]) * lag2(xi[2], n[2]),
                    lag2(xi[0], n[0]) * lag2(xi[1], n[1]) * dlag2(xi[2], n[2]),
                ]
            })
            .collect()
    }
}

impl ReferenceElement for Hex20 {
    const DIM: usize = 3;
    const NUM_NODES: usize = 20;
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
            vec![0.0, -1.0, -1.0],
            vec![1.0, 0.0, -1.0],
            vec![0.0, 1.0, -1.0],
            vec![-1.0, 0.0, -1.0],
            vec![0.0, -1.0, 1.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![-1.0, 0.0, 1.0],
            vec![-1.0, -1.0, 0.0],
            vec![1.0, -1.0, 0.0],
            vec![1.0, 1.0, 0.0],
            vec![-1.0, 1.0, 0.0],
        ]
    }
    fn shape(xi: &[f64]) -> Vec<f64> {
        let (x, y, z) = (xi[0], xi[1], xi[2]);
        Self::nodes()
            .iter()
            .map(|n| hex20_shape_at(n, x, y, z))
            .collect()
    }
    fn grad(xi: &[f64]) -> Vec<Vec<f64>> {
        let (x, y, z) = (xi[0], xi[1], xi[2]);
        Self::nodes()
            .iter()
            .map(|n| hex20_grad_at(n, x, y, z))
            .collect()
    }
}

/// `Hex20` shape-function value for a node with reference coords `n` at
/// `(x, y, z)`.
fn hex20_shape_at(n: &[f64], x: f64, y: f64, z: f64) -> f64 {
    let zeros = n.iter().filter(|v| **v == 0.0).count();
    if zeros == 0 {
        // Corner: 1/8 (1+a x)(1+b y)(1+c z)(a x + b y + c z - 2).
        let (a, b, c) = (n[0], n[1], n[2]);
        (1.0 + a * x) * (1.0 + b * y) * (1.0 + c * z) * (a * x + b * y + c * z - 2.0) / 8.0
    } else {
        // Edge midpoint: the single zero coordinate selects the edge axis.
        // `zeros != 0` guarantees exactly one coordinate is `0.0`, so this
        // `position` always finds a match. Contracted precondition (mirrors the
        // `mat_det`/`mat_inv` treatment, Phase 9b/10a; see `todo.md` 11b):
        // conversion to `Result` was deferred as the invariant is provable.
        debug_assert!(n.contains(&0.0));
        match n.iter().position(|v| *v == 0.0).unwrap() {
            0 => 0.25 * (1.0 - x * x) * (1.0 + n[1] * y) * (1.0 + n[2] * z), // ξ = 0
            1 => 0.25 * (1.0 + n[0] * x) * (1.0 - y * y) * (1.0 + n[2] * z), // η = 0
            _ => 0.25 * (1.0 + n[0] * x) * (1.0 + n[1] * y) * (1.0 - z * z), // ζ = 0
        }
    }
}

/// `Hex20` shape-function gradient for a node with reference coords `n` at
/// `(x, y, z)`.
fn hex20_grad_at(n: &[f64], x: f64, y: f64, z: f64) -> Vec<f64> {
    let zeros = n.iter().filter(|v| **v == 0.0).count();
    if zeros == 0 {
        let (a, b, c) = (n[0], n[1], n[2]);
        let aa = 1.0 + a * x;
        let bb = 1.0 + b * y;
        let cc = 1.0 + c * z;
        let dd = a * x + b * y + c * z - 2.0;
        // dN/dx = a/8 * bb * cc * (dd + aa); etc.
        vec![
            0.125 * a * bb * cc * (dd + aa),
            0.125 * b * aa * cc * (dd + bb),
            0.125 * c * aa * bb * (dd + cc),
        ]
    } else {
        // Edge midpoint: the single zero coordinate selects the edge axis.
        // `zeros != 0` guarantees exactly one coordinate is `0.0`; see the
        // `mat_det`/`mat_inv` contracted-precondition note (Phase 9b/10a,
        // `todo.md` 11b).
        debug_assert!(n.contains(&0.0));
        match n.iter().position(|v| *v == 0.0).unwrap() {
            0 => vec![
                -0.5 * x * (1.0 + n[1] * y) * (1.0 + n[2] * z),
                0.25 * (1.0 - x * x) * n[1] * (1.0 + n[2] * z),
                0.25 * (1.0 - x * x) * (1.0 + n[1] * y) * n[2],
            ],
            1 => vec![
                0.25 * n[0] * (1.0 - y * y) * (1.0 + n[2] * z),
                -0.5 * y * (1.0 + n[0] * x) * (1.0 + n[2] * z),
                0.25 * (1.0 + n[0] * x) * (1.0 - y * y) * n[2],
            ],
            _ => vec![
                0.25 * n[0] * (1.0 + n[1] * y) * (1.0 - z * z),
                0.25 * (1.0 + n[0] * x) * n[1] * (1.0 - z * z),
                -0.5 * z * (1.0 + n[0] * x) * (1.0 + n[1] * y),
            ],
        }
    }
}

/// Determinant of an `n×n` matrix stored row-major in a flat slice.
///
/// # Panics
///
/// Panics if `n` is not `1`, `2`, or `3`. This is a contracted precondition:
/// every `Map`/`from_nodes_and_grad` instance in this crate has a known spatial
/// dimension of 1/2/3, so an out-of-range `n` is unreachable in practice. See
/// the Phase 9b / 10a notes in `todo.md` — converting this (and `mat_inv`) to a
/// `Result` was deferred because it would ripple the public Jacobian API for
/// negligible gain.
fn mat_det(n: usize, m: &[f64]) -> f64 {
    debug_assert!(
        (1..=3).contains(&n),
        "mat_det: unsupported dimension {n} (only 1, 2, or 3 are valid)"
    );
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
///
/// # Panics
///
/// Panics if `n` is not `1`, `2`, or `3` (a contracted precondition; see the
/// `mat_det` doc note for the rationale and the deferred Phase 9b / 10a
/// conversion-to-`Result` discussion in `todo.md`).
fn mat_inv(n: usize, m: &[f64]) -> Vec<f64> {
    debug_assert!(
        (1..=3).contains(&n),
        "mat_inv: unsupported dimension {n} (only 1, 2, or 3 are valid)"
    );
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
                // Row 0 of the true inverse (cofactor matrix transposed).
                (m[4] * m[8] - m[5] * m[7]) * inv,
                (m[2] * m[7] - m[1] * m[8]) * inv,
                (m[1] * m[5] - m[2] * m[4]) * inv,
                // Row 1.
                (m[5] * m[6] - m[3] * m[8]) * inv,
                (m[0] * m[8] - m[2] * m[6]) * inv,
                (m[2] * m[3] - m[0] * m[5]) * inv,
                // Row 2.
                (m[3] * m[7] - m[4] * m[6]) * inv,
                (m[1] * m[6] - m[0] * m[7]) * inv,
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
        // The element (reference) dimension is the width of a local gradient, not
        // the physical coordinate count: a planar 2-D element imported from Gmsh
        // carries 3-component coordinates, but its Jacobian is the 2×2 planar one.
        let n = local_grad[0].len();
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

#[cfg(test)]
mod p2_tests {
    use super::*;

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < tol)
    }

    /// A Lagrange or serendipity element reproduces affine maps exactly, so
    /// mapping the reference nodes to themselves yields the identity Jacobian
    /// (determinant 1) at every point in the reference domain.
    fn check_jacobian_identity<T: ReferenceElement>() {
        let nodes = T::nodes();
        let pt = vec![0.2_f64; T::DIM]
            .iter()
            .map(|v| v * 0.5 + 0.1)
            .collect::<Vec<_>>();
        let grad = T::grad(&pt);
        let m = Map::from_nodes_and_grad(&nodes, &grad);
        assert!(
            (m.determinant - 1.0).abs() < 1e-9,
            "det={} ({} nodes)",
            m.determinant,
            T::NUM_NODES
        );
        let id = (0..T::DIM)
            .flat_map(|r| (0..T::DIM).map(move |c| if r == c { 1.0 } else { 0.0 }))
            .collect::<Vec<_>>();
        assert!(
            close(&m.jacobian, &id, 1e-9),
            "jac={:?} ({} nodes)",
            m.jacobian,
            T::NUM_NODES
        );
    }

    fn check_partition_and_kronecker<T: ReferenceElement>() {
        let samples: Vec<Vec<f64>> = match T::DIM {
            2 => vec![
                vec![0.25, 0.15],
                vec![0.5, 0.1],
                vec![-0.2, 0.4],
                vec![0.0, 0.0],
            ],
            3 => vec![
                vec![0.2, 0.1, 0.15],
                vec![0.1, -0.2, 0.3],
                vec![0.0, 0.0, 0.0],
            ],
            _ => vec![vec![0.0; T::DIM]],
        };
        for p in &samples {
            let s = T::shape(p);
            let sum = s.iter().sum::<f64>();
            assert!((sum - 1.0).abs() < 1e-9, "partition sum={sum}");
        }
        let nodes = T::nodes();
        for (j, n) in nodes.iter().enumerate() {
            let s = T::shape(n);
            for (i, v) in s.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (v - expected).abs() < 1e-9,
                    "{} node {j} shape {i} = {v}",
                    T::NUM_NODES
                );
            }
        }
    }

    #[test]
    fn tri6_props() {
        check_partition_and_kronecker::<Tri6>();
        check_jacobian_identity::<Tri6>();
    }

    #[test]
    fn quad9_props() {
        check_partition_and_kronecker::<Quad9>();
        check_jacobian_identity::<Quad9>();
    }

    #[test]
    fn quad8_props() {
        check_partition_and_kronecker::<Quad8>();
        check_jacobian_identity::<Quad8>();
    }

    #[test]
    fn tet10_props() {
        check_partition_and_kronecker::<Tet10>();
        check_jacobian_identity::<Tet10>();
    }

    #[test]
    fn hex27_props() {
        check_partition_and_kronecker::<Hex27>();
        check_jacobian_identity::<Hex27>();
    }

    #[test]
    fn hex20_props() {
        check_partition_and_kronecker::<Hex20>();
        check_jacobian_identity::<Hex20>();
    }
}
