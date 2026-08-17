//! Heat-conduction and Poisson element formulations for `tpt-fem`.
//!
//! The steady-state scalar problem solved here is
//!
//! ```text
//! -∇·(k ∇u) = f   in Ω,
//! ```
//!
//! equipped with Dirichlet (essential), Neumann (outward-flux), and Robin
//! (convective) boundary conditions. The element stiffness and load vectors are
//! integrated with the reference-element quadrature from `tpt-fem-quadrature`
//! and the isoparametric Jacobian from `tpt-fem-element`, then scattered into
//! a global system and solved by `tpt-fem-assembly` + `tpt-fem-sparse`.
//!
//! Scalar fields (one degree of freedom per node) are assumed.

use tpt_fem_assembly::{apply_neumann, apply_robin, assemble, solve_with_dirichlet};
use tpt_fem_element::{
    Hex20, Hex27, Hex8, Line2, Map, Quad4, Quad8, Quad9, ReferenceElement, Tet10, Tet4, Tri3, Tri6,
};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};
use tpt_fem_sparse::SparseError;

/// Quadrature points (reference coordinates) and weights for a cell type.
fn cell_quad(cell: CellType, order: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    match cell {
        CellType::Line => {
            let r = gauss_legendre(order);
            (r.points.iter().map(|x| vec![*x]).collect(), r.weights)
        }
        CellType::Tri => {
            let r = triangle(TriangleRule::Degree2);
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Quad => {
            let r = tensor_square(&gauss_legendre(order));
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Tet => {
            let r = tetrahedron(TetrahedronRule::Degree2);
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
        CellType::Hex => {
            let r = tensor_cube(&gauss_legendre(order));
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
        CellType::Tri6 => {
            let r = triangle(TriangleRule::HammerStroud);
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Quad8 | CellType::Quad9 => {
            let r = tensor_square(&gauss_legendre(order + 1));
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Tet10 => {
            let r = tetrahedron(TetrahedronRule::Keast4);
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
        CellType::Hex20 | CellType::Hex27 => {
            let r = tensor_cube(&gauss_legendre(order + 1));
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
    }
}

fn ref_shape(cell: CellType, xi: &[f64]) -> Vec<f64> {
    match cell {
        CellType::Line => Line2::shape(xi),
        CellType::Tri => Tri3::shape(xi),
        CellType::Quad => Quad4::shape(xi),
        CellType::Tet => Tet4::shape(xi),
        CellType::Hex => Hex8::shape(xi),
        CellType::Tri6 => Tri6::shape(xi),
        CellType::Quad8 => Quad8::shape(xi),
        CellType::Quad9 => Quad9::shape(xi),
        CellType::Tet10 => Tet10::shape(xi),
        CellType::Hex20 => Hex20::shape(xi),
        CellType::Hex27 => Hex27::shape(xi),
    }
}

fn ref_grad(cell: CellType, xi: &[f64]) -> Vec<Vec<f64>> {
    match cell {
        CellType::Line => Line2::grad(xi),
        CellType::Tri => Tri3::grad(xi),
        CellType::Quad => Quad4::grad(xi),
        CellType::Tet => Tet4::grad(xi),
        CellType::Hex => Hex8::grad(xi),
        CellType::Tri6 => Tri6::grad(xi),
        CellType::Quad8 => Quad8::grad(xi),
        CellType::Quad9 => Quad9::grad(xi),
        CellType::Tet10 => Tet10::grad(xi),
        CellType::Hex20 => Hex20::grad(xi),
        CellType::Hex27 => Hex27::grad(xi),
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Element stiffness matrix `K_e` for `-∇·(k ∇u) = f` with constant
/// conductivity `k`, returned in node order (1 DOF/node).
pub fn poisson_element_matrix(
    mesh: &Mesh,
    eid: usize,
    conductivity: f64,
    quad_order: usize,
) -> Vec<Vec<f64>> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let (qpts, qw) = cell_quad(cell, quad_order);
    let n = elem.nodes.len();
    let mut k = vec![vec![0.0; n]; n];
    for (qp, w) in qpts.iter().zip(&qw) {
        let local = ref_grad(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        let gphys: Vec<Vec<f64>> = local.iter().map(|g| map.physical_grad(g)).collect();
        for i in 0..n {
            for j in 0..n {
                k[i][j] += w * conductivity * dot(&gphys[i], &gphys[j]) * det;
            }
        }
    }
    k
}

/// Element load vector from a source term `f(x)`, returned in node order.
pub fn poisson_source_vector(
    mesh: &Mesh,
    eid: usize,
    source: impl Fn(&[f64]) -> f64,
    quad_order: usize,
) -> Vec<f64> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let (qpts, qw) = cell_quad(cell, quad_order);
    let n = elem.nodes.len();
    let mut f = vec![0.0; n];
    for (qp, w) in qpts.iter().zip(&qw) {
        let ns = ref_shape(cell, qp);
        let local = ref_grad(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        let mut x = vec![0.0; phys[0].len()];
        for k in 0..n {
            for d in 0..x.len() {
                x[d] += ns[k] * phys[k][d];
            }
        }
        let s = source(&x);
        for i in 0..n {
            f[i] += w * s * ns[i] * det;
        }
    }
    f
}

/// Solve the steady Poisson/heat-conduction problem.
///
/// * `conductivity` — constant scalar `k`.
/// * `source` — volumetric source `f(x)`.
/// * `dirichlet` — `(dof, value)` essential conditions (one DOF per node).
/// * `neumann` / `robin` — optional natural boundary fluxes `g(x, n)` and
///   Robin coefficients `h(x, n)`; applied via `tpt-fem-assembly`.
#[allow(clippy::type_complexity)]
pub fn solve_poisson<S>(
    mesh: &Mesh,
    conductivity: f64,
    quad_order: usize,
    source: S,
    dirichlet: &[(usize, f64)],
    neumann: Option<&dyn Fn(&[f64], &[f64]) -> f64>,
    robin: Option<&dyn Fn(&[f64], &[f64]) -> f64>,
) -> Result<Vec<f64>, SparseError>
where
    S: Fn(&[f64]) -> f64,
{
    let ndof = mesh.node_count();
    let mut coo = assemble(mesh, 1, |eid, m| {
        poisson_element_matrix(m, eid, conductivity, quad_order)
    });

    let mut rhs = vec![0.0; ndof];
    for eid in 0..mesh.elements.len() {
        let f = poisson_source_vector(mesh, eid, &source, quad_order);
        let elem = &mesh.elements[eid];
        for (i, &node) in elem.nodes.iter().enumerate() {
            rhs[node] += f[i];
        }
    }

    if let Some(nf) = neumann {
        apply_neumann(mesh, 1, nf, &mut rhs);
    }
    if let Some(cf) = robin {
        apply_robin(mesh, 1, cf, &mut coo);
    }

    solve_with_dirichlet(&coo, &rhs, dirichlet)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    #[test]
    fn poisson_1d_quadratic_source() {
        // -u'' = 1 on [0,1], u(0)=u(1)=0. Exact u = 0.5*(x - x^2).
        let nx = 16;
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
        .unwrap();
        let mid = nx / 2;
        let exact = 0.5 * (0.5 - 0.25);
        assert!(
            (u[mid] - exact).abs() < 2e-3,
            "got {} expected {}",
            u[mid],
            exact
        );
    }

    #[test]
    fn poisson_p2_tri6_linear_is_exact() {
        // Single Tri6 element on the reference triangle. `u = x` is linear, so
        // the P2 interpolant is exact; with Dirichlet on the three vertices the
        // free mid-nodes must recover the linear field (verifies the full
        // P2 assembly + solve path end to end).
        use tpt_fem_element::Tri6;
        let phys = Tri6::nodes();
        let mut b = MeshBuilder::new();
        let mut ids = Vec::new();
        for p in &phys {
            ids.push(b.add_node(vec![p[0], p[1]]));
        }
        b.add_element(CellType::Tri6, ids.clone());
        let mesh = b.build();
        let u = solve_poisson(
            &mesh,
            1.0,
            4,
            |x| x[0],
            &[(ids[0], 0.0), (ids[1], 1.0), (ids[2], 0.0)],
            None,
            None,
        )
        .unwrap();
        // Mid nodes: (0.5,0) -> 0.5, (0.5,0.5) -> 0.5, (0,0.5) -> 0.0.
        assert!((u[ids[3]] - 0.5).abs() < 1e-9, "got {}", u[ids[3]]);
        assert!((u[ids[4]] - 0.5).abs() < 1e-9, "got {}", u[ids[4]]);
        assert!((u[ids[5]] - 0.0).abs() < 1e-9, "got {}", u[ids[5]]);
    }

    #[test]
    fn poisson_2d_linear_is_exact() {
        // u = x + y is harmonic; P1 interpolant is exact. Dirichlet on the four
        // corners, free centre node must recover u = 1.
        let mut b = MeshBuilder::new();
        let c00 = b.add_node(vec![0.0, 0.0]);
        let c10 = b.add_node(vec![1.0, 0.0]);
        let c01 = b.add_node(vec![0.0, 1.0]);
        let c11 = b.add_node(vec![1.0, 1.0]);
        let mid = b.add_node(vec![0.5, 0.5]);
        b.add_element(CellType::Tri, vec![c00, c10, mid]);
        b.add_element(CellType::Tri, vec![c10, c11, mid]);
        b.add_element(CellType::Tri, vec![c11, c01, mid]);
        b.add_element(CellType::Tri, vec![c01, c00, mid]);
        let mesh = b.build();
        let u = solve_poisson(
            &mesh,
            1.0,
            2,
            |_| 0.0,
            &[(c00, 0.0), (c10, 1.0), (c01, 1.0), (c11, 2.0)],
            None,
            None,
        )
        .unwrap();
        assert!((u[mid] - 1.0).abs() < 1e-10, "got {}", u[mid]);
        assert!((u[c00] - 0.0).abs() < 1e-12);
        assert!((u[c11] - 2.0).abs() < 1e-12);
    }
}
