//! Element-to-global assembly and boundary-condition application for `tpt-fem`.
//!
//! This crate ties the core crates together into the standard FEM workflow:
//!
//! 1. [`assemble`] scatters per-element matrices (computed by a physics crate
//!    such as `tpt-fem-thermal`) into a global [`Coo`] triplet matrix.
//! 2. Natural boundary conditions are added through [`apply_neumann`] (load
//!    vector) and [`apply_robin`] (stiffness + load contributions).
//! 3. Essential conditions are enforced by [`solve_with_dirichlet`], which
//!    statically condenses the fixed degrees of freedom, solves the reduced
//!    system with `tpt-fem-sparse`, and scatters the solution back.
//!
//! All routines are dimension-agnostic across the five linear element types
//! `Line2`, `Tri3`, `Quad4`, `Tet4`, and `Hex8`, and support an arbitrary
//! (uniform) number of degrees of freedom per node.

use std::collections::{HashMap, HashSet};

use tpt_fem_element::{
    Hex20, Hex27, Hex8, Line2, Quad4, Quad8, Quad9, ReferenceElement, Tet10, Tet4, Tri3, Tri6,
};
use tpt_fem_mesh::{CellType, ElementId, Mesh};
use tpt_fem_quadrature::{gauss_legendre, tensor_square, triangle, TriangleRule};
use tpt_fem_sparse::{solve, Coo, SparseError};

/// Default quadrature order used for boundary-face integration.
const FACE_QUAD_ORDER: usize = 2;

/// A face of a reference element: the local node indices that bound it and the
/// reference element that parametrises the face (`Line` for an edge, `Tri`/`Quad`
/// for a surface).
struct FaceDef {
    nodes: &'static [usize],
    face: CellType,
}

const LINE_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[1],
        face: CellType::Line,
    },
];
const TRI_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[1, 2],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[2, 0],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[0, 1],
        face: CellType::Line,
    },
];
const QUAD_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 1],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[1, 2],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[2, 3],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[3, 0],
        face: CellType::Line,
    },
];
const TET_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[1, 2, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 2, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 1, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 1, 2],
        face: CellType::Tri,
    },
];
const HEX_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 1, 2, 3],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[4, 5, 6, 7],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 1, 5, 4],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[3, 2, 6, 7],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 4, 7, 3],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[1, 5, 6, 2],
        face: CellType::Quad,
    },
];

// Quadratic (P2) faces are keyed by their *corner* nodes only — the corner set
// uniquely identifies an edge/face, which is all `boundary_faces` needs to
// detect unshared (boundary) faces. Mid/center nodes are omitted. The corner
// local indices follow `tpt-fem-element`'s reference ordering: `Quad9` and
// `Hex27` interleave their interior/mid nodes, so their corners are not the
// first four/eight indices.
const TRI6_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[1, 2],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[2, 0],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[0, 1],
        face: CellType::Line,
    },
];
const QUAD8_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 1],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[1, 2],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[2, 3],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[3, 0],
        face: CellType::Line,
    },
];
const QUAD9_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 6],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[6, 8],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[8, 2],
        face: CellType::Line,
    },
    FaceDef {
        nodes: &[2, 0],
        face: CellType::Line,
    },
];
const TET10_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[1, 2, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 2, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 1, 3],
        face: CellType::Tri,
    },
    FaceDef {
        nodes: &[0, 1, 2],
        face: CellType::Tri,
    },
];
const HEX20_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 1, 2, 3],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[4, 5, 6, 7],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 1, 5, 4],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[3, 2, 6, 7],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 4, 7, 3],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[1, 5, 6, 2],
        face: CellType::Quad,
    },
];
const HEX27_FACES: &[FaceDef] = &[
    FaceDef {
        nodes: &[0, 18, 24, 6],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[2, 20, 26, 8],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 18, 20, 2],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[6, 24, 26, 8],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[0, 2, 8, 6],
        face: CellType::Quad,
    },
    FaceDef {
        nodes: &[18, 20, 26, 24],
        face: CellType::Quad,
    },
];

fn faces_of(cell: CellType) -> &'static [FaceDef] {
    match cell {
        CellType::Line => LINE_FACES,
        CellType::Tri => TRI_FACES,
        CellType::Quad => QUAD_FACES,
        CellType::Tet => TET_FACES,
        CellType::Hex => HEX_FACES,
        CellType::Tri6 => TRI6_FACES,
        CellType::Quad8 => QUAD8_FACES,
        CellType::Quad9 => QUAD9_FACES,
        CellType::Tet10 => TET10_FACES,
        CellType::Hex20 => HEX20_FACES,
        CellType::Hex27 => HEX27_FACES,
    }
}

/// Reference-element node coordinates for a cell type.
fn ref_nodes(cell: CellType) -> Vec<Vec<f64>> {
    match cell {
        CellType::Line => Line2::nodes(),
        CellType::Tri => Tri3::nodes(),
        CellType::Quad => Quad4::nodes(),
        CellType::Tet => Tet4::nodes(),
        CellType::Hex => Hex8::nodes(),
        CellType::Tri6 => Tri6::nodes(),
        CellType::Quad8 => Quad8::nodes(),
        CellType::Quad9 => Quad9::nodes(),
        CellType::Tet10 => Tet10::nodes(),
        CellType::Hex20 => Hex20::nodes(),
        CellType::Hex27 => Hex27::nodes(),
    }
}

/// Reference-element shape-function values at local coordinates `xi`.
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

/// Reference-element shape-function gradients at local coordinates `xi`.
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

/// Reference-element spatial dimension.
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

/// Quadrature points (in face-reference coordinates) and weights for a face.
fn face_quad(cell: CellType, order: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
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
        other => panic!("face_quad: unsupported face cell {other:?}"),
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn cross(a: &[f64], b: &[f64]) -> Vec<f64> {
    vec![
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mean(pts: &[Vec<f64>]) -> Vec<f64> {
    let dim = pts[0].len();
    let mut m = vec![0.0; dim];
    for p in pts {
        for (d, v) in p.iter().enumerate() {
            m[d] += v;
        }
    }
    for v in m.iter_mut() {
        *v /= pts.len() as f64;
    }
    m
}

/// Surface measure and outward unit normal from face tangent vectors.
///
/// * `face_dim == 1` (an edge): the normal is the in-plane perpendicular.
/// * `face_dim == 2` (a surface): the normal is the cross product.
fn surface(tangents: &[Vec<f64>], dim: usize) -> (f64, Vec<f64>) {
    match tangents.len() {
        1 => {
            let t = &tangents[0];
            let m = norm(t);
            let mut n = vec![0.0; dim];
            n[0] = -t[1];
            n[1] = t[0];
            let ln = norm(&n);
            for v in n.iter_mut() {
                *v /= ln;
            }
            (m, n)
        }
        2 => {
            let c = cross(&tangents[0], &tangents[1]);
            let m = norm(&c);
            let mut n = c;
            let ln = norm(&n);
            for v in n.iter_mut() {
                *v /= ln;
            }
            (m, n)
        }
        _ => panic!("surface: unsupported face dimension {}", tangents.len()),
    }
}

/// Fallible variant of [`assemble`].
///
/// Identical to [`assemble`], but `elem_matrix` may fail (returning an error of
/// type `E`), in which case the first failure is propagated rather than
/// panicking. Use this when the per-element operator can reject an invalid
/// model/dimension combination.
pub fn try_assemble<E>(
    mesh: &Mesh,
    dofs_per_node: usize,
    elem_matrix: impl Fn(usize, &Mesh) -> Result<Vec<Vec<f64>>, E>,
) -> Result<Coo, E> {
    let mut coo = Coo::new();
    for (eid, elem) in mesh.elements.iter().enumerate() {
        let k = elem_matrix(eid, mesh)?;
        let ndof = k.len();
        for i in 0..ndof {
            for j in 0..ndof {
                let kij = k[i][j];
                if kij == 0.0 {
                    continue;
                }
                let gi = elem.nodes[i / dofs_per_node] * dofs_per_node + i % dofs_per_node;
                let gj = elem.nodes[j / dofs_per_node] * dofs_per_node + j % dofs_per_node;
                coo.push(gi, gj, kij);
            }
        }
    }
    Ok(coo)
}

/// Scatter per-element matrices into a global [`Coo`] matrix.
///
/// `elem_matrix(eid, mesh)` returns the element matrix in *local DOF order*,
/// where local DOF `a` belongs to node `a / dofs_per_node` and component
/// `a % dofs_per_node`. Local DOF `a` maps to global DOF
/// `element.nodes[a / dofs_per_node] * dofs_per_node + a % dofs_per_node`.
pub fn assemble(
    mesh: &Mesh,
    dofs_per_node: usize,
    elem_matrix: impl Fn(usize, &Mesh) -> Vec<Vec<f64>>,
) -> Coo {
    let mut coo = Coo::new();
    for (eid, elem) in mesh.elements.iter().enumerate() {
        let k = elem_matrix(eid, mesh);
        let ndof = k.len();
        for i in 0..ndof {
            for j in 0..ndof {
                let kij = k[i][j];
                if kij == 0.0 {
                    continue;
                }
                let gi = elem.nodes[i / dofs_per_node] * dofs_per_node + i % dofs_per_node;
                let gj = elem.nodes[j / dofs_per_node] * dofs_per_node + j % dofs_per_node;
                coo.push(gi, gj, kij);
            }
        }
    }
    coo
}

/// The result of condensing the fixed (Dirichlet) DOFs out of a system.
///
/// `kred` and `fred` live in *reduced* DOF numbering `[0, free.len())`, where
/// reduced DOF `i` corresponds to the global DOF `free[i]`. This is the shared
/// low-level primitive behind [`solve_with_dirichlet`], the generalized
/// eigenproblem, and arc-length continuation: each of those just needs the
/// reduced matrix/vector and the `free` map back to global DOFs.
pub struct ReducedSystem {
    /// Reduced stiffness (or tangent) matrix, in reduced DOF numbering.
    pub kred: Coo,
    /// Reduced right-hand side, in reduced DOF numbering.
    pub fred: Vec<f64>,
    /// Global DOF index corresponding to each reduced DOF `i`.
    pub free: Vec<usize>,
}

/// Condense the fixed (Dirichlet) DOFs out of `A u = f`.
///
/// The prescribed values are moved to the right-hand side (the standard
/// "lifting" of the essential boundary condition), leaving a square reduced
/// system in the free DOFs. The returned [`ReducedSystem`] carries the `free`
/// map so the caller can scatter a reduced solution back to global DOFs.
pub fn reduce_system(coo: &Coo, rhs: &[f64], bcs: &[(usize, f64)]) -> ReducedSystem {
    let n = rhs.len();
    let csr = coo.to_csr();
    let fixed: HashSet<usize> = bcs.iter().map(|(i, _)| *i).collect();
    let free: Vec<usize> = (0..n).filter(|i| !fixed.contains(i)).collect();
    let free_idx: HashMap<usize, usize> = free.iter().enumerate().map(|(k, &v)| (v, k)).collect();

    let mut kred = Coo::new();
    let mut fred = vec![0.0; free.len()];
    for &fdof in &free {
        let mut r = rhs[fdof];
        for c in csr.row_ptrs[fdof]..csr.row_ptrs[fdof + 1] {
            let col = csr.col_ind[c];
            let v = csr.values[c];
            if fixed.contains(&col) {
                // `col` was collected directly from `bcs` into `fixed`, so this
                // DOF is guaranteed to appear in `bcs`. Contracted precondition
                // (mirrors the `mat_det`/`mat_inv` treatment from Phase 9b/10a):
                // converting the lookup to a `Result` was deferred because the
                // invariant is local and provably true. See `todo.md` (11b).
                debug_assert!(bcs.iter().any(|(i, _)| *i == col));
                let val = bcs.iter().find(|(i, _)| *i == col).unwrap().1;
                r -= v * val;
            } else {
                kred.push(free_idx[&fdof], free_idx[&col], v);
            }
        }
        fred[free_idx[&fdof]] = r;
    }

    ReducedSystem { kred, fred, free }
}

/// Solve `A u = f` with essential (Dirichlet) conditions.
///
/// Fixed DOFs are condensed out of the system (see [`reduce_system`]), the
/// reduced system is solved by `tpt-fem-sparse`, and the prescribed values
/// are scattered back into the solution vector (which has one entry per global
/// DOF).
pub fn solve_with_dirichlet(
    coo: &Coo,
    rhs: &[f64],
    bcs: &[(usize, f64)],
) -> Result<Vec<f64>, SparseError> {
    let n = rhs.len();
    let reduced = reduce_system(coo, rhs, bcs);
    let x = solve(&reduced.kred, &reduced.fred)?;
    let mut u = vec![0.0; n];
    for (i, &dof) in reduced.free.iter().enumerate() {
        u[dof] = x[i];
    }
    for (i, val) in bcs {
        u[*i] = *val;
    }
    Ok(u)
}

/// All `(element_id, face_index)` pairs lying on the mesh boundary (i.e. faces
/// not shared with any other element).
pub fn boundary_faces(mesh: &Mesh) -> Vec<(ElementId, usize)> {
    let mut counts: HashMap<Vec<usize>, usize> = HashMap::new();
    for elem in &mesh.elements {
        for f in faces_of(elem.cell_type) {
            let mut key: Vec<usize> = f.nodes.iter().map(|&i| elem.nodes[i]).collect();
            key.sort_unstable();
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    let mut out = Vec::new();
    for (eid, elem) in mesh.elements.iter().enumerate() {
        for (fi, f) in faces_of(elem.cell_type).iter().enumerate() {
            let mut key: Vec<usize> = f.nodes.iter().map(|&i| elem.nodes[i]).collect();
            key.sort_unstable();
            if counts[&key] == 1 {
                out.push((eid, fi));
            }
        }
    }
    out
}

/// Add a Neumann (natural) boundary load to `rhs`.
///
/// For each boundary face, `flux(x, n)` returns the outward flux `g(x)`; the
/// contribution `∫_Γ g v dΓ` is accumulated into the right-hand side, where `v`
/// is the test function. `dofs_per_node` tiles the scalar load across
/// components.
pub fn apply_neumann(
    mesh: &Mesh,
    dofs_per_node: usize,
    flux: impl Fn(&[f64], &[f64]) -> f64,
    rhs: &mut [f64],
) {
    apply_neumann_order(mesh, dofs_per_node, flux, rhs, FACE_QUAD_ORDER);
}

/// Like [`apply_neumann`] with an explicit quadrature order.
pub fn apply_neumann_order(
    mesh: &Mesh,
    dofs_per_node: usize,
    flux: impl Fn(&[f64], &[f64]) -> f64,
    rhs: &mut [f64],
    order: usize,
) {
    for (eid, fi) in boundary_faces(mesh) {
        let elem = &mesh.elements[eid];
        let f = &faces_of(elem.cell_type)[fi];
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        let (qpts, qw) = face_quad(f.face, order);
        let refp = ref_nodes(elem.cell_type);
        let dim = phys[0].len();
        let elem_centroid = mean(&phys);
        let face_coords: Vec<Vec<f64>> = f.nodes.iter().map(|&i| phys[i].clone()).collect();
        let face_centroid = mean(&face_coords);

        for (qp, w) in qpts.iter().zip(&qw) {
            let ns = ref_shape(f.face, qp);
            let ng = ref_grad(f.face, qp);
            let mut xr = vec![0.0; ref_dim(elem.cell_type)];
            for (k, &idx) in f.nodes.iter().enumerate() {
                for d in 0..xr.len() {
                    xr[d] += ns[k] * refp[idx][d];
                }
            }
            let mut xp = vec![0.0; dim];
            for (k, &idx) in f.nodes.iter().enumerate() {
                for d in 0..dim {
                    xp[d] += ns[k] * phys[idx][d];
                }
            }
            let fd = ng[0].len();
            let mut tangents = vec![vec![0.0; dim]; fd];
            for m in 0..fd {
                for (k, &idx) in f.nodes.iter().enumerate() {
                    for d in 0..dim {
                        tangents[m][d] += ng[k][m] * phys[idx][d];
                    }
                }
            }
            let (measure, mut normal) = surface(&tangents, dim);
            let mut dir = vec![0.0; dim];
            for d in 0..dim {
                dir[d] = face_centroid[d] - elem_centroid[d];
            }
            if dot(&normal, &dir) < 0.0 {
                for v in normal.iter_mut() {
                    *v = -*v;
                }
            }
            let g = flux(&xp, &normal);
            for (k, &idx) in f.nodes.iter().enumerate() {
                let base = elem.nodes[idx] * dofs_per_node;
                for c in 0..dofs_per_node {
                    rhs[base + c] += w * measure * g * ns[k];
                }
            }
            let _ = xr;
        }
    }
}

/// Add a Robin (mixed) boundary contribution to the global matrix `coo`.
///
/// The Robin condition `h u = g` contributes `∫_Γ h u v dΓ` to the bilinear
/// form, i.e. `coo += h N_i N_j` summed over the boundary-face quadrature. The
/// known flux `g` is handled by the caller through [`apply_neumann`].
pub fn apply_robin(
    mesh: &Mesh,
    dofs_per_node: usize,
    coeff: impl Fn(&[f64], &[f64]) -> f64,
    coo: &mut Coo,
) {
    apply_robin_order(mesh, dofs_per_node, coeff, coo, FACE_QUAD_ORDER);
}

/// Like [`apply_robin`] with an explicit quadrature order.
pub fn apply_robin_order(
    mesh: &Mesh,
    dofs_per_node: usize,
    coeff: impl Fn(&[f64], &[f64]) -> f64,
    coo: &mut Coo,
    order: usize,
) {
    for (eid, fi) in boundary_faces(mesh) {
        let elem = &mesh.elements[eid];
        let f = &faces_of(elem.cell_type)[fi];
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        let (qpts, qw) = face_quad(f.face, order);
        let refp = ref_nodes(elem.cell_type);
        let dim = phys[0].len();
        let elem_centroid = mean(&phys);
        let face_coords: Vec<Vec<f64>> = f.nodes.iter().map(|&i| phys[i].clone()).collect();
        let face_centroid = mean(&face_coords);

        for (qp, w) in qpts.iter().zip(&qw) {
            let ns = ref_shape(f.face, qp);
            let ng = ref_grad(f.face, qp);
            let mut xr = vec![0.0; ref_dim(elem.cell_type)];
            for (k, &idx) in f.nodes.iter().enumerate() {
                for d in 0..xr.len() {
                    xr[d] += ns[k] * refp[idx][d];
                }
            }
            let mut xp = vec![0.0; dim];
            for (k, &idx) in f.nodes.iter().enumerate() {
                for d in 0..dim {
                    xp[d] += ns[k] * phys[idx][d];
                }
            }
            let fd = ng[0].len();
            let mut tangents = vec![vec![0.0; dim]; fd];
            for m in 0..fd {
                for (k, &idx) in f.nodes.iter().enumerate() {
                    for d in 0..dim {
                        tangents[m][d] += ng[k][m] * phys[idx][d];
                    }
                }
            }
            let (measure, mut normal) = surface(&tangents, dim);
            let mut dir = vec![0.0; dim];
            for d in 0..dim {
                dir[d] = face_centroid[d] - elem_centroid[d];
            }
            if dot(&normal, &dir) < 0.0 {
                for v in normal.iter_mut() {
                    *v = -*v;
                }
            }
            let h = coeff(&xp, &normal);
            for (k, &ik) in f.nodes.iter().enumerate() {
                for (l, &il) in f.nodes.iter().enumerate() {
                    let base_i = elem.nodes[ik] * dofs_per_node;
                    let base_j = elem.nodes[il] * dofs_per_node;
                    for c in 0..dofs_per_node {
                        coo.push(base_i + c, base_j + c, w * measure * h * ns[k] * ns[l]);
                    }
                }
            }
            let _ = xr;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    fn square_tri_mesh() -> Mesh {
        // Unit square split by the diagonal (0,0)-(1,1) into two triangles,
        // plus a centre node at (0.5,0.5) split into four triangles.
        let mut b = MeshBuilder::new();
        let n00 = b.add_node(vec![0.0, 0.0]);
        let n10 = b.add_node(vec![1.0, 0.0]);
        let n01 = b.add_node(vec![0.0, 1.0]);
        let n11 = b.add_node(vec![1.0, 1.0]);
        // diagonal from n00 to n11
        b.add_element(CellType::Tri, vec![n00, n10, n11]);
        b.add_element(CellType::Tri, vec![n00, n11, n01]);
        b.build()
    }

    #[test]
    fn assemble_sizes_and_symmetry() {
        let mesh = square_tri_mesh();
        // Element matrix = identity-per-node (node block), 1 dof/node.
        let coo = assemble(&mesh, 1, |eid, m| {
            let n = m.elements[eid].nodes.len();
            let mut k = vec![vec![0.0; n]; n];
            for i in 0..n {
                k[i][i] = 1.0;
            }
            k
        });
        let csr = coo.to_csr();
        assert_eq!(csr.nrows, 4);
        assert_eq!(csr.ncols, 4);
        // Two triangles share nodes 0 and 3 (the diagonal), so the per-element
        // diagonal blocks at those nodes are summed; 4 distinct entries remain.
        assert_eq!(csr.nnz(), 4);
    }

    #[test]
    fn dirichlet_reduces_to_known_solution() {
        // 2-element line on [0,1] with unit conductivity, fix u(0)=0, u(1)=1.
        // The P1 solution is linear: u(0.5)=0.5.
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0]);
        let n1 = b.add_node(vec![0.5]);
        let n2 = b.add_node(vec![1.0]);
        b.add_element(CellType::Line, vec![n0, n1]);
        b.add_element(CellType::Line, vec![n1, n2]);
        let mesh = b.build();
        // Element stiffness for unit conductivity, length 0.5: 2*[[1,-1],[-1,1]].
        let coo = assemble(&mesh, 1, |_, _| vec![vec![2.0, -2.0], vec![-2.0, 2.0]]);
        let u = solve_with_dirichlet(&coo, &[0.0, 0.0, 0.0], &[(n0, 0.0), (n2, 1.0)]).unwrap();
        assert!((u[n0] - 0.0).abs() < 1e-12);
        assert!((u[n1] - 0.5).abs() < 1e-12);
        assert!((u[n2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn neumann_conserves_total_load() {
        // Single triangle: perimeter = 2 + sqrt(2). Constant unit flux over the
        // whole boundary. Sum of nodal loads = total measure (since sum of
        // shape functions on a face equals 1).
        let mut b = MeshBuilder::new();
        let a = b.add_node(vec![0.0, 0.0]);
        let c = b.add_node(vec![1.0, 0.0]);
        let d = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![a, c, d]);
        let mesh = b.build();
        let mut rhs = vec![0.0; 3];
        apply_neumann(&mesh, 1, |_, _| 1.0, &mut rhs);
        let total: f64 = rhs.iter().sum();
        let expected = 2.0 + 2.0_f64.sqrt();
        assert!(
            (total - expected).abs() < 1e-9,
            "got {total} expected {expected}"
        );
    }

    #[test]
    fn robin_conserves_total_measure() {
        // Robin coefficient h=1 over the boundary. Sum of all coo entries on a
        // single face = measure (sum_i sum_j N_i N_j = sum_i N_i = 1) times h.
        let mut b = MeshBuilder::new();
        let a = b.add_node(vec![0.0, 0.0]);
        let c = b.add_node(vec![1.0, 0.0]);
        let d = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![a, c, d]);
        let mesh = b.build();
        let mut coo = Coo::new();
        apply_robin(&mesh, 1, |_, _| 1.0, &mut coo);
        let total: f64 = coo.vals.iter().sum();
        let expected = 2.0 + 2.0_f64.sqrt();
        assert!(
            (total - expected).abs() < 1e-9,
            "got {total} expected {expected}"
        );
    }

    #[test]
    fn boundary_faces_of_square() {
        let mesh = square_tri_mesh();
        let bf = boundary_faces(&mesh);
        // Square has 4 edges, each an edge of exactly one triangle.
        assert_eq!(bf.len(), 4);
    }
}
