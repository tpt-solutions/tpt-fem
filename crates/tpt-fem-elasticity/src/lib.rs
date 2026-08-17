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

use tpt_fem_assembly::{assemble, reduce_system, solve_with_dirichlet};
use tpt_fem_eigen::generalized_lanczos_eigs;
use tpt_fem_element::{
    Hex20, Hex27, Hex8, Line2, Map, Quad4, Quad8, Quad9, ReferenceElement, Tet10, Tet4, Tri3, Tri6,
};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};
use tpt_fem_sparse::{Coo, SparseError};

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
        CellType::Tri6 => Tri6::DIM,
        CellType::Quad8 => Quad8::DIM,
        CellType::Quad9 => Quad9::DIM,
        CellType::Tet10 => Tet10::DIM,
        CellType::Hex20 => Hex20::DIM,
        CellType::Hex27 => Hex27::DIM,
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
        // P2 elements need higher-order rules: the mass matrix involves
        // NᵀN (degree 4), so bump the tensor-product/1-D rules by one and use
        // the highest fixed simplex rules.
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

/// Consistent mass matrix `M = ∫ ρ Nᵀ N dΩ` for an elasticity model, returned as
/// a [`Coo`] with `dim` DOFs per node. `density` is the mass per reference
/// volume (or per unit length for [`ElasticModel::BarAxial`]).
pub fn elasticity_mass_matrix(
    mesh: &Mesh,
    _model: ElasticModel,
    density: f64,
    quad_order: usize,
) -> Coo {
    let dim = ref_dim(mesh.elements[0].cell_type);
    let mut coo = Coo::new();
    for elem in &mesh.elements {
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        let cell = elem.cell_type;
        let n = elem.nodes.len();
        let (qpts, qw) = cell_quad(cell, quad_order);
        // Scalar element mass m_ij = ρ ∫ N_i N_j dΩ.
        let mut m = vec![vec![0.0; n]; n];
        for (qp, w) in qpts.iter().zip(&qw) {
            let local = ref_grad(cell, qp);
            let map = Map::from_nodes_and_grad(&phys, &local);
            let det = map.determinant.abs();
            let ns = ref_shape(cell, qp);
            for i in 0..n {
                for j in 0..n {
                    m[i][j] += w * det * density * ns[i] * ns[j];
                }
            }
        }
        for i in 0..n {
            for j in 0..n {
                let v = m[i][j];
                if v != 0.0 {
                    for a in 0..dim {
                        coo.push(elem.nodes[i] * dim + a, elem.nodes[j] * dim + a, v);
                    }
                }
            }
        }
    }
    coo
}

/// Row-sum (lumped) mass matrix, with each node's mass spread equally across its
/// `dim` DOFs. Convenient for explicit dynamics where a diagonal mass is needed.
pub fn elasticity_lumped_mass(
    mesh: &Mesh,
    model: ElasticModel,
    density: f64,
    quad_order: usize,
) -> Coo {
    let dim = ref_dim(mesh.elements[0].cell_type);
    let consistent = elasticity_mass_matrix(mesh, model, density, quad_order);
    let csr = consistent.to_csr();
    let n = csr.nrows;
    let mut coo = Coo::new();
    for node in 0..n / dim {
        let mut total = 0.0;
        for a in 0..dim {
            let row = node * dim + a;
            for c in csr.row_ptrs[row]..csr.row_ptrs[row + 1] {
                total += csr.values[c];
            }
        }
        let lump = total / dim as f64;
        for a in 0..dim {
            coo.push(node * dim + a, node * dim + a, lump);
        }
    }
    coo
}

/// Solve the generalized eigenproblem `K Φ = ω² M Φ` (natural vibration modes).
///
/// Stiffness and mass are assembled, the Dirichlet conditions are condensed out
/// of *both* via [`reduce_system`] (so the reduced bases match), and the
/// resulting pencil is solved with [`generalized_lanczos_eigs`]. The returned
/// eigenpairs `(ω², φ)` are the squared natural frequencies and their mode
/// shapes, scattered back into the full (zero on fixed DOFs) displacement space.
#[allow(clippy::too_many_arguments)]
pub fn solve_modal(
    mesh: &Mesh,
    model: ElasticModel,
    young: f64,
    poisson: f64,
    density: f64,
    quad_order: usize,
    num_modes: usize,
    dirichlet: &[(usize, f64)],
) -> Result<Vec<(f64, Vec<f64>)>, SparseError> {
    let dim = ref_dim(mesh.elements[0].cell_type);
    let ndof = mesh.node_count() * dim;
    let k = assemble(mesh, dim, |eid, m| {
        elasticity_element_matrix(m, eid, model, young, poisson, quad_order)
    });
    let m = elasticity_mass_matrix(mesh, model, density, quad_order);

    let red_k = reduce_system(&k, &vec![0.0; ndof], dirichlet);
    let red_m = reduce_system(&m, &vec![0.0; ndof], dirichlet);
    let nred = red_k.free.len();
    let lanczos_dim = (num_modes + 8).min(nred.max(1));

    let modes = generalized_lanczos_eigs(&red_k.kred, &red_m.kred, 0.0, num_modes, lanczos_dim)?;
    Ok(modes
        .into_iter()
        .map(|(lam, v)| {
            let mut phi = vec![0.0; ndof];
            for (i, &dof) in red_k.free.iter().enumerate() {
                phi[dof] = v[i];
            }
            (lam, phi)
        })
        .collect())
}

// ---------------------------------------------------------------------------
// 2-D Euler–Bernoulli frame element (Hermite-cubic)
#[derive(Clone, Copy, Debug)]
pub struct BeamSection2D {
    /// Axial stiffness `E A`.
    pub ea: f64,
    /// Bending stiffness `E I`.
    pub ei: f64,
    /// Mass per unit length `ρ A` (used by the consistent-mass matrix).
    pub mass_per_length: f64,
}

impl BeamSection2D {
    /// Build from Young's modulus `e`, area `a`, second moment `i`, and density
    /// `rho` (mass density; `mass_per_length = rho * a`).
    pub fn from_eai(e: f64, a: f64, i: f64, rho: f64) -> Self {
        BeamSection2D {
            ea: e * a,
            ei: e * i,
            mass_per_length: rho * a,
        }
    }

    /// Build directly from `ea`, `ei`, and `mass_per_length`.
    pub fn from_ea_ei(ea: f64, ei: f64, mass_per_length: f64) -> Self {
        BeamSection2D {
            ea,
            ei,
            mass_per_length,
        }
    }
}

/// Local 6×6 Euler–Bernoulli frame stiffness in node-DOF order
/// `[u1, v1, θ1, u2, v2, θ2]` (axial `u`, transverse `v`, rotation `θ`).
pub fn beam2d_element_matrix(s: &BeamSection2D, length: f64) -> Vec<Vec<f64>> {
    let l = length;
    let mut k = vec![vec![0.0; 6]; 6];
    // Axial (DOFs 0 and 3).
    let c = s.ea / l;
    k[0][0] = c;
    k[0][3] = -c;
    k[3][0] = -c;
    k[3][3] = c;
    // Bending (DOFs 1, 2, 4, 5).
    let d = s.ei / (l * l * l);
    let bi: [[f64; 4]; 4] = [
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
    ];
    let bend = [1usize, 2, 4, 5];
    for i in 0..4 {
        for j in 0..4 {
            k[bend[i]][bend[j]] = d * bi[i][j];
        }
    }
    k
}

/// Local 6×6 consistent-mass matrix in node-DOF order `[u1, v1, θ1, u2, v2, θ2]`.
pub fn beam2d_consistent_mass(s: &BeamSection2D, length: f64) -> Vec<Vec<f64>> {
    let l = length;
    let mu = s.mass_per_length;
    let mut m = vec![vec![0.0; 6]; 6];
    // Axial consistent mass (DOFs 0 and 3): ρA L / 6 * [[2,1],[1,2]].
    let am = mu * l / 6.0;
    m[0][0] = 2.0 * am;
    m[0][3] = am;
    m[3][0] = am;
    m[3][3] = 2.0 * am;
    // Bending consistent mass (DOFs 1, 2, 4, 5): ρA L / 420 * [...].
    let bm = mu * l / 420.0;
    let bi: [[f64; 4]; 4] = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l * l, 13.0 * l, -3.0 * l * l],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l * l, -22.0 * l, 4.0 * l * l],
    ];
    let bend = [1usize, 2, 4, 5];
    for i in 0..4 {
        for j in 0..4 {
            m[bend[i]][bend[j]] = bm * bi[i][j];
        }
    }
    m
}

/// 6×6 local→global rotation for a frame element whose local axis makes angle
/// `phi` with the global X axis. Each node's `[u, v, θ]` block rotates by the
/// in-plane rotation `Q = [[cos, sin, 0], [-sin, cos, 0], [0, 0, 1]]`.
fn beam2d_rotation(phi: f64) -> Vec<Vec<f64>> {
    let c = phi.cos();
    let s = phi.sin();
    let mut r = vec![vec![0.0; 6]; 6];
    let q = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
    for b in 0..2 {
        for i in 0..3 {
            for j in 0..3 {
                r[b * 3 + i][b * 3 + j] = q[i][j];
            }
        }
    }
    r
}

/// Solve a 2-D frame (Euler–Bernoulli beam + axial) problem.
///
/// Each `Line2` element carries 3 DOFs per node `[u, v, θ]`; elements are
/// rotated into global coordinates by their orientation. `loads` returns the
/// *total* `[Fx, Fy, Mz]` applied at a node and is evaluated exactly once per
/// (unique) node, so a load at a node shared by several elements is **not**
/// double-counted. `dirichlet` is a list of `(global_dof, value)` essential
/// conditions (global dof `node * 3 + comp`). Returns the displacement vector
/// (one entry per global DOF).
pub fn solve_frame2d(
    mesh: &Mesh,
    section: &BeamSection2D,
    loads: impl Fn(usize, &[f64]) -> [f64; 3],
    dirichlet: &[(usize, f64)],
) -> Result<Vec<f64>, SparseError> {
    let ndof = mesh.node_count() * 3;
    let mut coo = Coo::new();
    let mut rhs = vec![0.0; ndof];

    for elem in &mesh.elements {
        assert_eq!(
            elem.cell_type,
            CellType::Line,
            "solve_frame2d needs Line2 elements"
        );
        let p0 = mesh.node_coords(elem.nodes[0]);
        let p1 = mesh.node_coords(elem.nodes[1]);
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let length = (dx * dx + dy * dy).sqrt();
        let phi = dy.atan2(dx);

        let kl = beam2d_element_matrix(section, length);
        let r = beam2d_rotation(phi);
        // K_global = R K_local Rᵀ  (rigid rotation of the element).
        let rt = transpose(&r);
        let temp = matmul(&r, &kl);
        let kg = matmul(&temp, &rt);

        let mut gdof = [0usize; 6];
        for i in 0..2 {
            for c in 0..3 {
                gdof[i * 3 + c] = elem.nodes[i] * 3 + c;
            }
        }
        for i in 0..6 {
            for j in 0..6 {
                let v = kg[i][j];
                if v != 0.0 {
                    coo.push(gdof[i], gdof[j], v);
                }
            }
        }
    }

    // Assemble the load vector once per unique node. (A node shared by two
    // elements would otherwise have its concentrated load counted twice.)
    let mut seen = std::collections::HashSet::new();
    for elem in &mesh.elements {
        for &node in &elem.nodes {
            if seen.insert(node) {
                let f = loads(node, mesh.node_coords(node));
                for c in 0..3 {
                    rhs[node * 3 + c] += f[c];
                }
            }
        }
    }

    solve_with_dirichlet(&coo, &rhs, dirichlet)
}

fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut t = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            t[i][j] = a[j][i];
        }
    }
    t
}

fn matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i][j] += aik * b[k][j];
            }
        }
    }
    c
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

    #[test]
    fn frame_cantilever_tip_load() {
        // Single horizontal cantilever, length 1, EI=1, tip transverse load
        // P=1 downward. Euler-Bernoulli gives delta = P L^3/(3EI) = 1/3 and
        // tip slope = P L^2/(2EI) = 1/2 (both negative for downward load).
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        b.add_element(CellType::Line, vec![n0, n1]);
        let mesh = b.build();
        let sec = BeamSection2D::from_ea_ei(1.0e9, 1.0, 1.0);
        let u = solve_frame2d(
            &mesh,
            &sec,
            |node, _| {
                if node == n1 {
                    [0.0, -1.0, 0.0]
                } else {
                    [0.0, 0.0, 0.0]
                }
            },
            &[(n0 * 3, 0.0), (n0 * 3 + 1, 0.0), (n0 * 3 + 2, 0.0)],
        )
        .unwrap();
        assert!(
            (u[n1 * 3 + 1] + (1.0 / 3.0)).abs() < 1e-9,
            "got {}",
            u[n1 * 3 + 1]
        );
        assert!((u[n1 * 3 + 2] + 0.5).abs() < 1e-9, "got {}", u[n1 * 3 + 2]);
    }

    #[test]
    fn frame_cantilever_rotated() {
        // Same cantilever rotated 90 degrees (vertical). A global horizontal tip
        // load P=1 should produce the same bending, now along global X.
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Line, vec![n0, n1]);
        let mesh = b.build();
        let sec = BeamSection2D::from_ea_ei(1.0e9, 1.0, 1.0);
        let u = solve_frame2d(
            &mesh,
            &sec,
            |node, _| {
                if node == n1 {
                    [-1.0, 0.0, 0.0]
                } else {
                    [0.0, 0.0, 0.0]
                }
            },
            &[(n0 * 3, 0.0), (n0 * 3 + 1, 0.0), (n0 * 3 + 2, 0.0)],
        )
        .unwrap();
        assert!((u[n1 * 3] + (1.0 / 3.0)).abs() < 1e-9, "got {}", u[n1 * 3]);
        assert!((u[n1 * 3 + 2] + 0.5).abs() < 1e-9, "got {}", u[n1 * 3 + 2]);
    }

    #[test]
    fn frame_simply_supported_central_load() {
        // Simply-supported beam (supports at both ends), length 1, EI=1, central
        // point load P=1. Midspan deflection delta = P L^3/(48EI) = 1/48. A
        // central point load produces a slope discontinuity that the C1-continuous
        // Hermite element smooths, so a fine mesh is used to converge.
        let n_elem = 20;
        let mut b = MeshBuilder::new();
        let mut nodes = Vec::new();
        for i in 0..=n_elem {
            nodes.push(b.add_node(vec![(i as f64) / (n_elem as f64), 0.0]));
        }
        for i in 0..n_elem {
            b.add_element(CellType::Line, vec![nodes[i], nodes[i + 1]]);
        }
        let mesh = b.build();
        let mid = nodes[n_elem / 2];
        let sec = BeamSection2D::from_ea_ei(1.0e9, 1.0, 1.0);
        let u = solve_frame2d(
            &mesh,
            &sec,
            |node, _| {
                if node == mid {
                    [0.0, -1.0, 0.0]
                } else {
                    [0.0, 0.0, 0.0]
                }
            },
            // Pinned at both ends (v=0) plus one axial roller (u=0 at node 0).
            &[
                (nodes[0] * 3, 0.0),
                (nodes[0] * 3 + 1, 0.0),
                (nodes[n_elem] * 3 + 1, 0.0),
            ],
        )
        .unwrap();
        assert!(
            (u[mid * 3 + 1] + (1.0 / 48.0)).abs() < 1e-4,
            "got {}",
            u[mid * 3 + 1]
        );
    }

    #[test]
    fn modal_bar_axial_fundamental() {
        // Uniform axial bar [0,1], fixed at x=0. Continuous fundamental
        // frequency w1 = (pi/2) sqrt(EA/(rhoA L^2)) = (pi/2) sqrt(EA/(rhoA)).
        // A P1 bar overestimates eigenfrequencies, so a fine mesh is used to
        // converge within tolerance.
        let n_elem = 30;
        let mut b = MeshBuilder::new();
        let mut nodes = Vec::new();
        for i in 0..=n_elem {
            nodes.push(b.add_node(vec![(i as f64) / (n_elem as f64)]));
        }
        for i in 0..n_elem {
            b.add_element(CellType::Line, vec![nodes[i], nodes[i + 1]]);
        }
        let mesh = b.build();
        let ea = 1.0;
        let rho_a = 1.0;
        let modes = solve_modal(
            &mesh,
            ElasticModel::BarAxial,
            ea,
            0.0,
            rho_a,
            2,
            1,
            &[(nodes[0], 0.0)],
        )
        .unwrap();
        let w1 = modes[0].0.sqrt();
        let expected = std::f64::consts::PI / 2.0 * (ea / rho_a).sqrt();
        assert!((w1 - expected).abs() < 1e-3, "got {w1} expected {expected}");
    }
}
