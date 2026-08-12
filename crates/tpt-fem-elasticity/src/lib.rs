//! Linear-elasticity element formulations for `tpt-fem`.
//!
//! Supports:
//!
//! * `Line2` — 1-D axial bar (`EA` stiffness),
//! * `Tri3` / `Quad4` — 2-D plane-stress and plane-strain continua,
//! * `Tet4` / `Hex8` — 3-D isotropic continua.
//!
//! Each element stiffness is `K_e = ∫_Ω Bᵀ D B dΩ`, with `B` the
//! strain–displacement matrix and `D` the isotropic constitutive matrix. The
//! per-element matrices are scattered by `tpt-fem-assembly` and solved with
//! `tpt-fem-sparse`.

use tpt_fem_assembly::{assemble, solve_with_dirichlet};
use tpt_fem_element::{Hex8, Line2, Map, Quad4, ReferenceElement, Tet4, Tri3};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};
use tpt_fem_sparse::SparseError;

/// The elasticity model (selects the constitutive matrix `D`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElasticModel {
    /// 1-D axial bar (`young` interpreted as `EA`).
    BarAxial,
    /// 2-D plane stress.
    PlaneStress,
    /// 2-D plane strain.
    PlaneStrain,
    /// 3-D isotropic continuum.
    Continuum3D,
}

fn ref_dim(cell: CellType) -> usize {
    match cell {
        CellType::Line => Line2::DIM,
        CellType::Tri => Tri3::DIM,
        CellType::Quad => Quad4::DIM,
        CellType::Tet => Tet4::DIM,
        CellType::Hex => Hex8::DIM,
    }
}

fn ref_shape(cell: CellType, xi: &[f64]) -> Vec<f64> {
    match cell {
        CellType::Line => Line2::shape(xi),
        CellType::Tri => Tri3::shape(xi),
        CellType::Quad => Quad4::shape(xi),
        CellType::Tet => Tet4::shape(xi),
        CellType::Hex => Hex8::shape(xi),
    }
}

fn ref_grad(cell: CellType, xi: &[f64]) -> Vec<Vec<f64>> {
    match cell {
        CellType::Line => Line2::grad(xi),
        CellType::Tri => Tri3::grad(xi),
        CellType::Quad => Quad4::grad(xi),
        CellType::Tet => Tet4::grad(xi),
        CellType::Hex => Hex8::grad(xi),
    }
}

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
    }
}

fn strain_dim(dim: usize) -> usize {
    match dim {
        1 => 1,
        2 => 3,
        3 => 6,
        _ => panic!("strain_dim: unsupported dimension {dim}"),
    }
}

/// Constitutive matrix `D` (`n_strain × n_strain`) for the given model.
fn constitutive(model: ElasticModel, e: f64, nu: f64, dim: usize) -> Vec<Vec<f64>> {
    match (model, dim) {
        (ElasticModel::BarAxial, 1) => vec![vec![e]],
        (ElasticModel::PlaneStress, 2) => {
            let c = e / (1.0 - nu * nu);
            vec![
                vec![c, c * nu, 0.0],
                vec![c * nu, c, 0.0],
                vec![0.0, 0.0, c * (1.0 - nu) / 2.0],
            ]
        }
        (ElasticModel::PlaneStrain, 2) => {
            let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            vec![
                vec![c * (1.0 - nu), c * nu, 0.0],
                vec![c * nu, c * (1.0 - nu), 0.0],
                vec![0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0],
            ]
        }
        (ElasticModel::Continuum3D, 3) => {
            let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            let m = (1.0 - nu) * c;
            let d = (1.0 - 2.0 * nu) / 2.0 * c;
            vec![
                vec![m, c * nu, c * nu, 0.0, 0.0, 0.0],
                vec![c * nu, m, c * nu, 0.0, 0.0, 0.0],
                vec![c * nu, c * nu, m, 0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, d, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, d, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, 0.0, d],
            ]
        }
        _ => panic!("constitutive: model {model:?} incompatible with dim {dim}"),
    }
}

/// Strain–displacement sub-matrix `B_i` (`n_strain × dim`) for node `i` with
/// physical gradients `g = [∂N/∂x, ∂N/∂y, ∂N/∂z]`.
fn b_matrix(g: &[f64], dim: usize) -> Vec<Vec<f64>> {
    match dim {
        1 => vec![vec![g[0]]],
        2 => vec![vec![g[0], 0.0], vec![0.0, g[1]], vec![g[1], g[0]]],
        3 => vec![
            vec![g[0], 0.0, 0.0],
            vec![0.0, g[1], 0.0],
            vec![0.0, 0.0, g[2]],
            vec![g[1], g[0], 0.0],
            vec![0.0, g[2], g[1]],
            vec![g[2], 0.0, g[0]],
        ],
        _ => panic!("b_matrix: unsupported dim {dim}"),
    }
}

/// Element stiffness matrix for the given elasticity model, returned in
/// `(node × dim)` local DOF order.
pub fn elasticity_element_matrix(
    mesh: &Mesh,
    eid: usize,
    model: ElasticModel,
    young: f64,
    poisson: f64,
    quad_order: usize,
) -> Vec<Vec<f64>> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let dim = ref_dim(cell);
    let n = elem.nodes.len();
    let nstr = strain_dim(dim);
    let d = constitutive(model, young, poisson, dim);
    let (qpts, qw) = cell_quad(cell, quad_order);

    let mut k = vec![vec![0.0; n * dim]; n * dim];
    for (qp, w) in qpts.iter().zip(&qw) {
        let local = ref_grad(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        let g: Vec<Vec<f64>> = local.iter().map(|g| map.physical_grad(g)).collect();
        let b: Vec<Vec<Vec<f64>>> = (0..n).map(|i| b_matrix(&g[i], dim)).collect();
        for i in 0..n {
            for j in 0..n {
                for a in 0..dim {
                    for bb in 0..dim {
                        let mut s = 0.0;
                        for s_ in 0..nstr {
                            for t in 0..nstr {
                                s += b[i][s_][a] * d[s_][t] * b[j][t][bb];
                            }
                        }
                        k[i * dim + a][j * dim + bb] += w * det * s;
                    }
                }
            }
        }
    }
    k
}

/// Element body-force vector (per node `dim` components) from `b(x)`.
pub fn elasticity_body_vector(
    mesh: &Mesh,
    eid: usize,
    body_force: impl Fn(&[f64]) -> Vec<f64>,
    quad_order: usize,
) -> Vec<f64> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let dim = ref_dim(cell);
    let n = elem.nodes.len();
    let (qpts, qw) = cell_quad(cell, quad_order);
    let mut f = vec![0.0; n * dim];
    for (qp, w) in qpts.iter().zip(&qw) {
        let ns = ref_shape(cell, qp);
        let local = ref_grad(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        let mut x = vec![0.0; phys[0].len()];
        for k in 0..n {
            for dd in 0..x.len() {
                x[dd] += ns[k] * phys[k][dd];
            }
        }
        let b = body_force(&x);
        for i in 0..n {
            for a in 0..dim {
                f[i * dim + a] += w * ns[i] * b[a] * det;
            }
        }
    }
    f
}

/// Solve a linear-elasticity problem.
///
/// `dirichlet` is a list of `(global_dof, value)` essential conditions on the
/// `node_id * dim + component` DOFs. Returns the displacement vector (one entry
/// per global DOF).
pub fn solve_elasticity(
    mesh: &Mesh,
    model: ElasticModel,
    young: f64,
    poisson: f64,
    quad_order: usize,
    body_force: impl Fn(&[f64]) -> Vec<f64>,
    dirichlet: &[(usize, f64)],
) -> Result<Vec<f64>, SparseError> {
    let dim = ref_dim(mesh.elements[0].cell_type);
    let ndof = mesh.node_count() * dim;
    let coo = assemble(mesh, dim, |eid, m| {
        elasticity_element_matrix(m, eid, model, young, poisson, quad_order)
    });
    let mut rhs = vec![0.0; ndof];
    for eid in 0..mesh.elements.len() {
        let f = elasticity_body_vector(mesh, eid, &body_force, quad_order);
        let elem = &mesh.elements[eid];
        for (i, &node) in elem.nodes.iter().enumerate() {
            for a in 0..dim {
                rhs[node * dim + a] += f[i * dim + a];
            }
        }
    }
    solve_with_dirichlet(&coo, &rhs, dirichlet)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    #[test]
    fn bar_axial_closed_form() {
        // Two bar elements on [0,1], EA=1, u(0)=0, tip load F=1 at x=1.
        // Analytic u(x)=x.
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0]);
        let n1 = b.add_node(vec![0.5]);
        let n2 = b.add_node(vec![1.0]);
        b.add_element(CellType::Line, vec![n0, n1]);
        b.add_element(CellType::Line, vec![n1, n2]);
        let mesh = b.build();
        let u = solve_elasticity(
            &mesh,
            ElasticModel::BarAxial,
            1.0,
            0.0,
            1,
            |_| vec![0.0],
            &[(n0, 0.0), (n2, 1.0)],
        )
        .unwrap();
        assert!((u[n0] - 0.0).abs() < 1e-12);
        assert!((u[n1] - 0.5).abs() < 1e-12);
        assert!((u[n2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn plane_stress_patch_test() {
        // Unit square split into four triangles around a centre node. Impose the
        // linear displacement field u = (0.1 x, 0.2 y) on the four corners; the
        // free centre node must reproduce u = (0.05, 0.1) exactly.
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

        let mut bcs = Vec::new();
        for (node, x, y) in [
            (c00, 0.0, 0.0),
            (c10, 0.1, 0.0),
            (c01, 0.0, 0.2),
            (c11, 0.1, 0.2),
        ] {
            bcs.push((node * 2, x));
            bcs.push((node * 2 + 1, y));
        }
        let u = solve_elasticity(
            &mesh,
            ElasticModel::PlaneStress,
            1.0,
            0.3,
            2,
            |_| vec![0.0, 0.0],
            &bcs,
        )
        .unwrap();
        assert!((u[mid * 2] - 0.05).abs() < 1e-9, "got {}", u[mid * 2]);
        assert!(
            (u[mid * 2 + 1] - 0.1).abs() < 1e-9,
            "got {}",
            u[mid * 2 + 1]
        );
    }
}
