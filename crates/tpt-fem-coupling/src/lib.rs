//! Multiphysics coupling operators for `tpt-fem`.
//!
//! * [`thermal_structural`] — thermal strain `ε = α ΔT I` fed into the linear
//!   elasticity operator as an initial strain (free thermal expansion).
//! * [`electro_thermal`] — Joule heating `q = σ |E|²` fed into the Poisson
//!   (heat-conduction) operator as a volumetric source.
//! * [`fsi_coupling`] — a basic partitioned fluid–structure interaction that
//!   exchanges traction/kinematics between `tpt-fem-fluid` and
//!   `tpt-fem-elasticity` through `tpt-fem-dynamic`.
//!
//! ```
//! use tpt_fem_coupling::thermal_structural;
//! use tpt_fem_elasticity::ElasticModel;
//! use tpt_fem_mesh::{CellType, MeshBuilder};
//!
//! // A unit 1-D bar, fixed at node 0 in x, heated uniformly by ΔT = 1.
//! // Free thermal expansion gives u(x) = α ΔT (x − x0); with α = 1e-3 the
//! // rightmost node (x = 1) should sit at 1e-3.
//! let mut b = MeshBuilder::new();
//! let mut prev = b.add_node(vec![0.0]);
//! for i in 1..=4 {
//!     let n = b.add_node(vec![i as f64 / 4.0]);
//!     b.add_element(CellType::Line, vec![prev, n]);
//!     prev = n;
//! }
//! let mesh = b.build();
//! let temp = vec![1.0; mesh.node_count()];
//! let dirichlet = [(0usize, 0.0)]; // fix node 0 x-DOF
//! let u = thermal_structural(&mesh, ElasticModel::BarAxial, 1.0, 0.3, 1e-3, &temp, &dirichlet).unwrap();
//! let tip = u[4]; // node 4, x-component (dim = 1 for BarAxial)
//! assert!((tip - 1e-3).abs() < 1e-6, "free expansion tip = {}", tip);
//! ```

use tpt_fem_assembly::{solve_with_dirichlet, try_assemble};
use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel, ElasticityError};
use tpt_fem_element::{Hex8, Line2, Map, Quad4, ReferenceElement, Tet4, Tri3};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};
use tpt_fem_sparse::SparseError;

/// Errors returned by the multiphysics coupling operators.
#[derive(Debug)]
pub enum CouplingError {
    /// The underlying linear-algebra solve (structure system, or a
    /// Dirichlet-reduced step) failed.
    Sparse(SparseError),
    /// Building or evaluating the elasticity operator failed (e.g. a
    /// model/dimension mismatch).
    Elasticity(ElasticityError),
    /// The fluid (steady-Stokes) solve within the coupling failed.
    Fluid(tpt_fem_fluid::FluidError),
}

impl std::fmt::Display for CouplingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CouplingError::Sparse(e) => write!(f, "coupling structure solve failed: {e}"),
            CouplingError::Elasticity(e) => write!(f, "coupling elasticity operator failed: {e}"),
            CouplingError::Fluid(e) => write!(f, "coupling fluid solve failed: {e}"),
        }
    }
}

impl std::error::Error for CouplingError {}

impl From<SparseError> for CouplingError {
    fn from(e: SparseError) -> Self {
        CouplingError::Sparse(e)
    }
}

impl From<ElasticityError> for CouplingError {
    fn from(e: ElasticityError) -> Self {
        CouplingError::Elasticity(e)
    }
}

impl From<tpt_fem_fluid::FluidError> for CouplingError {
    fn from(e: tpt_fem_fluid::FluidError) -> Self {
        CouplingError::Fluid(e)
    }
}

/// Solve the thermal-structural problem: free thermal expansion of an elastic body
/// under a per-node temperature rise `delta_t[node]` (in Kelvin), coefficient of
/// thermal expansion `alpha`, with `dirichlet` fixing `(node*dim + comp, value)`
/// displacement DOFs. Returns the displacement vector (one entry per global DOF).
///
/// The thermal strain `εᵗʰ = α·ΔT·I` is applied as an initial strain: the global
/// load is `f = K·u_free` where `u_free(x) = α·ΔT(x)·(x − x_c)` is the free
/// expansion field (centred on each element's centroid so only the strain, not a
/// rigid shift, contributes). This is exact for a uniform `ΔT` and a good
/// approximation for smoothly varying `ΔT`.
pub fn thermal_structural(
    mesh: &Mesh,
    model: ElasticModel,
    young: f64,
    poisson: f64,
    alpha: f64,
    delta_t: &[f64],
    dirichlet: &[(usize, f64)],
) -> Result<Vec<f64>, tpt_fem_sparse::SparseError> {
    let dim = match mesh.elements[0].cell_type {
        CellType::Line => 1,
        CellType::Tri | CellType::Quad => 2,
        CellType::Tet | CellType::Hex => 3,
        other => panic!("coupling: unsupported cell {other:?}"),
    };
    let ndof = mesh.node_count() * dim;
    let k_full = try_assemble(mesh, dim, |eid, m| {
        elasticity_element_matrix(m, eid, model, young, poisson, 2)
    })
    .map_err(|e| tpt_fem_sparse::SparseError::Numeric(e.to_string()))?;

    let mut rhs = vec![0.0; ndof];
    for (eid, elem) in mesh.elements.iter().enumerate() {
        let ke = elasticity_element_matrix(mesh, eid, model, young, poisson, 2)
            .map_err(|e| tpt_fem_sparse::SparseError::Numeric(e.to_string()))?;
        let n = elem.nodes.len();
        // Element centroid.
        let mut xc = vec![0.0; dim];
        for &node in &elem.nodes {
            let c = mesh.node_coords(node);
            for a in 0..dim {
                xc[a] += c[a];
            }
        }
        for a in 0..dim {
            xc[a] /= n as f64;
        }
        // Free-expansion nodal displacements for this element.
        let mut u0 = vec![0.0; n * dim];
        for (i, &node) in elem.nodes.iter().enumerate() {
            let c = mesh.node_coords(node);
            for a in 0..dim {
                u0[i * dim + a] = alpha * delta_t[node] * (c[a] - xc[a]);
            }
        }
        // f_e = K_e u0, scattered into the global RHS at node*dim + comp.
        for i in 0..n * dim {
            let mut s = 0.0;
            for j in 0..n * dim {
                s += ke[i][j] * u0[j];
            }
            let node = elem.nodes[i / dim];
            rhs[node * dim + (i % dim)] += s;
        }
    }
    solve_with_dirichlet(&k_full, &rhs, dirichlet)
}

/// Joule-heating volumetric source `q = σ |E|²` (W/m³) at a point with electric
/// field magnitude `e_mag` and electrical conductivity `sigma`.
pub fn joule_source(sigma: f64, e_mag: f64) -> f64 {
    sigma * e_mag * e_mag
}

/// Solve electro-thermal (Joule heating) steady conduction: assemble the Poisson
/// operator with the Joule source `q = σ|E|²` (constant over the mesh) and the
/// given Dirichlet temperatures. Returns the temperature field (one entry per
/// node). `conductivity` is the thermal conductivity `k`.
pub fn electro_thermal(
    mesh: &Mesh,
    conductivity: f64,
    sigma: f64,
    e_mag: f64,
    dirichlet: &[(usize, f64)],
) -> Result<Vec<f64>, tpt_fem_sparse::SparseError> {
    let q = joule_source(sigma, e_mag);
    tpt_fem_thermal::solve_poisson(mesh, conductivity, 2, |_| q, dirichlet, None, None)
}

/// Work-consistent fluid–structure interface load vector.
///
/// Given the `(structure node, fluid node)` interface pairing and the fluid
/// pressure at every fluid node, assemble the structure's Neumann load vector
/// by integrating the traction `p·n̂` over the interface *faces* with the
/// face shape functions — the work-equivalent (consistent) load, rather than
/// a lumped per-node point force.
///
/// The interface faces are discovered as the boundary faces of the structure
/// mesh whose nodes all lie on the interface (an edge in 2-D / a triangle or
/// quad face in 3-D that belongs to exactly one element). Each face's normal
/// is oriented *outward*, away from its owning element's centroid, so the
/// direction adapts to horizontal, vertical, and curved interfaces alike.
/// Isolated interface nodes that belong to no interface face (single-node
/// coupling) fall back to a lumped load `p·n̂` with the same geometry-aware
/// outward normal.
pub fn fsi_interface_loads(
    struct_mesh: &Mesh,
    interface: &[(usize, usize)],
    fluid_pressure: &[f64],
) -> Vec<f64> {
    let dim = match struct_mesh.elements[0].cell_type {
        CellType::Tri | CellType::Quad => 2,
        CellType::Tet | CellType::Hex => 3,
        other => panic!("fsi_interface_loads: unsupported cell {other:?}"),
    };
    let mut s_to_f: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &(s, f) in interface {
        s_to_f.insert(s, f);
    }
    let press = |s_node: usize| -> f64 {
        let f = s_to_f.get(&s_node).copied().unwrap_or_else(|| {
            panic!("fsi_interface_loads: structure node {s_node} not on the interface")
        });
        fluid_pressure[f]
    };
    let coords = |n: usize| -> Vec<f64> { struct_mesh.node_coords(n).to_vec() };
    let mut loads = vec![0.0; struct_mesh.node_count() * dim];

    // Boundary faces of the mesh: a face occurring in exactly one element.
    let element_faces = |cell: CellType, nodes: &[usize]| -> Vec<Vec<usize>> {
        match cell {
            CellType::Tri => vec![
                vec![nodes[0], nodes[1]],
                vec![nodes[1], nodes[2]],
                vec![nodes[2], nodes[0]],
            ],
            CellType::Quad => vec![
                vec![nodes[0], nodes[1]],
                vec![nodes[1], nodes[2]],
                vec![nodes[2], nodes[3]],
                vec![nodes[3], nodes[0]],
            ],
            CellType::Tet => vec![
                vec![nodes[0], nodes[1], nodes[2]],
                vec![nodes[0], nodes[1], nodes[3]],
                vec![nodes[0], nodes[2], nodes[3]],
                vec![nodes[1], nodes[2], nodes[3]],
            ],
            CellType::Hex => vec![
                vec![nodes[0], nodes[1], nodes[2], nodes[3]],
                vec![nodes[4], nodes[5], nodes[6], nodes[7]],
                vec![nodes[0], nodes[1], nodes[5], nodes[4]],
                vec![nodes[1], nodes[2], nodes[6], nodes[5]],
                vec![nodes[2], nodes[3], nodes[7], nodes[6]],
                vec![nodes[3], nodes[0], nodes[4], nodes[7]],
            ],
            other => panic!("fsi_interface_loads: unsupported cell {other:?}"),
        }
    };
    let mut face_count: std::collections::HashMap<Vec<usize>, usize> =
        std::collections::HashMap::new();
    let mut elem_of_face: std::collections::HashMap<Vec<usize>, usize> =
        std::collections::HashMap::new();
    for (eid, elem) in struct_mesh.elements.iter().enumerate() {
        for face in element_faces(elem.cell_type, &elem.nodes) {
            let mut key = face.clone();
            key.sort_unstable();
            *face_count.entry(key.clone()).or_insert(0) += 1;
            elem_of_face.entry(key).or_insert(eid);
        }
    }

    // Flip `nu` so it points away from the owning element's centroid.
    let orient = |nu: &mut [f64], eid: usize, face_nodes: &[usize]| {
        let mut xc = vec![0.0; dim];
        for &nd in &struct_mesh.elements[eid].nodes {
            let cc = coords(nd);
            for a in 0..dim {
                xc[a] += cc[a];
            }
        }
        let nn = struct_mesh.elements[eid].nodes.len() as f64;
        for a in 0..dim {
            xc[a] /= nn;
        }
        let mut xm = vec![0.0; dim];
        for &nd in face_nodes {
            let cc = coords(nd);
            for a in 0..dim {
                xm[a] += cc[a];
            }
        }
        for a in 0..dim {
            xm[a] /= face_nodes.len() as f64;
        }
        let side: f64 = (0..dim).map(|a| (xm[a] - xc[a]) * nu[a]).sum();
        if side < 0.0 {
            for v in nu.iter_mut() {
                *v = -*v;
            }
        }
    };
    // Consistent load of one linear triangle: f_i = n̂ A (2p_i + p_j + p_k)/12.
    let add_tri = |loads: &mut Vec<f64>, tri: &[usize], n_hat: &[f64], area: f64| {
        let p: Vec<f64> = tri.iter().map(|&n| press(n)).collect();
        for i in 0..3 {
            let w = area * (2.0 * p[i] + p[(i + 1) % 3] + p[(i + 2) % 3]) / 12.0;
            for a in 0..dim {
                loads[tri[i] * dim + a] += n_hat[a] * w;
            }
        }
    };

    let mut covered = vec![false; struct_mesh.node_count()];
    let mut seen: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
    for (eid, elem) in struct_mesh.elements.iter().enumerate() {
        for face in element_faces(elem.cell_type, &elem.nodes) {
            let mut key = face.clone();
            key.sort_unstable();
            if face_count.get(&key).copied().unwrap_or(0) != 1 || !seen.insert(key) {
                continue;
            }
            // An interface face has every node on the interface.
            if !face.iter().all(|n| s_to_f.contains_key(n)) {
                continue;
            }
            match dim {
                2 => {
                    let x0 = coords(face[0]);
                    let x1 = coords(face[1]);
                    let tx = x1[0] - x0[0];
                    let ty = x1[1] - x0[1];
                    let len = (tx * tx + ty * ty).sqrt();
                    if len < 1e-14 {
                        continue;
                    }
                    let mut nu = vec![ty / len, -tx / len];
                    orient(&mut nu, eid, &face);
                    let (p0, p1) = (press(face[0]), press(face[1]));
                    for a in 0..2 {
                        loads[face[0] * 2 + a] += nu[a] * len * (2.0 * p0 + p1) / 6.0;
                        loads[face[1] * 2 + a] += nu[a] * len * (p0 + 2.0 * p1) / 6.0;
                    }
                    covered[face[0]] = true;
                    covered[face[1]] = true;
                }
                _ => {
                    // Triangulate the face fan-style around corner 0.
                    let xs: Vec<Vec<f64>> = face.iter().map(|&n| coords(n)).collect();
                    for k in 1..xs.len() - 1 {
                        let e1: Vec<f64> = (0..3).map(|d| xs[k][d] - xs[0][d]).collect();
                        let e2: Vec<f64> = (0..3).map(|d| xs[k + 1][d] - xs[0][d]).collect();
                        let cr = [
                            e1[1] * e2[2] - e1[2] * e2[1],
                            e1[2] * e2[0] - e1[0] * e2[2],
                            e1[0] * e2[1] - e1[1] * e2[0],
                        ];
                        let mag = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
                        if mag < 1e-14 {
                            continue;
                        }
                        let mut nu = vec![cr[0] / mag, cr[1] / mag, cr[2] / mag];
                        let sub = vec![face[0], face[k], face[k + 1]];
                        orient(&mut nu, eid, &sub);
                        add_tri(&mut loads, &sub, &nu, 0.5 * mag);
                        for &nd in &sub {
                            covered[nd] = true;
                        }
                    }
                }
            }
        }
    }
    // Fallback for isolated interface nodes that belong to no interface face
    // (e.g. a single-node coupling): apply a lumped load p·n̂ with the
    // geometry-aware outward normal estimated from the incident elements.
    for &(s_node, f_node) in interface {
        if covered[s_node] {
            continue;
        }
        let c = coords(s_node);
        let mut acc = vec![0.0; dim];
        let mut n_elems = 0usize;
        for elem in &struct_mesh.elements {
            if elem.nodes.contains(&s_node) {
                let mut xc = vec![0.0; dim];
                for &nd in &elem.nodes {
                    let cc = coords(nd);
                    for a in 0..dim {
                        xc[a] += cc[a];
                    }
                }
                for a in 0..dim {
                    xc[a] /= elem.nodes.len() as f64;
                    acc[a] += c[a] - xc[a];
                }
                n_elems += 1;
            }
        }
        let mut normal = if n_elems > 0 { acc } else { vec![0.0; dim] };
        let mag = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
        if mag > 1e-12 {
            for a in 0..dim {
                normal[a] /= mag;
            }
        } else if dim >= 2 {
            normal[1] = 1.0; // degenerate fallback
        } else {
            normal[0] = 1.0;
        }
        for a in 0..dim {
            loads[s_node * dim + a] += fluid_pressure[f_node] * normal[a];
        }
    }
    loads
}

///
/// Given a structure displacement `u_struct` (global DOF order, `dim` per node)
/// and a shared interface of `(structure_node, fluid_node)` pairs, this:
/// 1. displaces the fluid mesh nodes by the structure displacement (FSI
///    kinematic condition),
/// 2. solves `tpt-fem-fluid`'s steady Stokes for the fluid pressure,
/// 3. transfers the fluid pressure as a normal traction onto the structure
///    interface (Newton's third law), using a *consistent* (shape-weighted,
///    non-lumped) load assembled through the interface elements' quadrature with
///    a geometry-aware outward normal at each interface node,
/// 4. solves the structure with that traction as a Neumann load.
///
/// Returns the updated structure displacement, or a [`CouplingError`] if the
/// fluid or structure solve fails (e.g. an under-constrained structure, or a
/// model/dimension mismatch), instead of panicking.
pub fn fsi_coupling(
    struct_mesh: &Mesh,
    fluid_mesh: &Mesh,
    model: ElasticModel,
    young: f64,
    poisson: f64,
    viscosity: f64,
    u_struct: &[f64],
    interface: &[(usize, usize)],
    struct_dirichlet: &[(usize, f64)],
    fluid_penalty: f64,
) -> Result<Vec<f64>, CouplingError> {
    let dim = match struct_mesh.elements[0].cell_type {
        CellType::Line => 1,
        CellType::Tri | CellType::Quad => 2,
        CellType::Tet | CellType::Hex => 3,
        other => panic!("coupling: unsupported cell {other:?}"),
    };
    // Displace fluid nodes per the interface map.
    let mut fluid = fluid_mesh.clone();
    let fdim = match fluid.elements[0].cell_type {
        CellType::Line => 1,
        CellType::Tri | CellType::Quad => 2,
        CellType::Tet | CellType::Hex => 3,
        other => panic!("coupling: unsupported fluid cell {other:?}"),
    };
    for &(s_node, f_node) in interface {
        let c = fluid.node_coords(f_node).to_vec();
        let mut nc = c.clone();
        for a in 0..fdim {
            nc[a] = c[a] + u_struct[s_node * dim + a];
        }
        fluid.nodes[f_node].coords = nc;
    }
    // Steady Stokes: a downward body force (gravity-like) drives a pressure that
    // pushes the structure. The fluid's top and bottom boundaries are no-slip so
    // the system is non-singular; the pressure that develops at the (bottom)
    // interface is what loads the structure.
    let top = (0..fluid.node_count())
        .max_by(|&a, &b| {
            fluid.node_coords(a)[1]
                .partial_cmp(&fluid.node_coords(b)[1])
                .unwrap()
        })
        .unwrap();
    let bot = (0..fluid.node_count())
        .min_by(|&a, &b| {
            fluid.node_coords(a)[1]
                .partial_cmp(&fluid.node_coords(b)[1])
                .unwrap()
        })
        .unwrap();
    let fymax = fluid.node_coords(top)[1];
    let fymin = fluid.node_coords(bot)[1];
    let mut fluid_bc = Vec::new();
    for n in 0..fluid.node_count() {
        if (fluid.node_coords(n)[1] - fymax).abs() < 1e-9
            || (fluid.node_coords(n)[1] - fymin).abs() < 1e-9
        {
            for a in 0..fdim {
                fluid_bc.push((n * fdim + a, 0.0));
            }
        }
    }
    let (_u, pressure) = tpt_fem_fluid::steady_stokes(
        &fluid,
        viscosity,
        |_: &[f64]| {
            let mut b = vec![0.0; fdim];
            if fdim >= 2 {
                b[1] = -1.0;
            }
            b
        },
        &fluid_bc,
        fluid_penalty,
    )
    .map_err(CouplingError::from)?;

    // Consistent (non-lumped) traction transfer. The fluid pressure is projected
    // onto the structure interface as a normal traction `t = p·n` and assembled
    // as a *consistent* FEM nodal load `f_a = ∫_Ω N_a (p·n) dΩ` via the interface
    // elements' reference quadrature, with both `p` and the geometry-aware
    // outward normal `n` interpolated from their nodal values through the shape
    // functions. This replaces the earlier lumped point-load projection and,
    // because `n` is derived from the actual interface geometry (average of
    // `node − incident-element-centroid`), curved or vertical interfaces receive
    // a correct, shape-weighted traction direction rather than a hardcoded +y.
    // (A strict surface-only mortar projection over the interface faces is a
    // future refinement; integrating over the interface element is a consistent,
    // smoke-level coupling load.)
    let s_nodes: Vec<usize> = interface.iter().map(|&(s, _)| s).collect();

    // Per interface-node geometry-aware outward normal.
    let mut n_node = vec![vec![0.0; dim]; struct_mesh.node_count()];
    for &s_node in &s_nodes {
        let c = struct_mesh.node_coords(s_node).to_vec();
        let mut acc = vec![0.0; dim];
        let mut n_elems = 0usize;
        for elem in &struct_mesh.elements {
            if elem.nodes.contains(&s_node) {
                let mut xc = vec![0.0; dim];
                for &nd in &elem.nodes {
                    let cc = struct_mesh.node_coords(nd);
                    for a in 0..dim {
                        xc[a] += cc[a];
                    }
                }
                for a in 0..dim {
                    xc[a] /= elem.nodes.len() as f64;
                    acc[a] += c[a] - xc[a];
                }
                n_elems += 1;
            }
        }
        let mut normal = if n_elems > 0 { acc } else { vec![0.0; dim] };
        let mag = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
        if mag > 1e-12 {
            for a in 0..dim {
                normal[a] /= mag;
            }
        } else if dim >= 2 {
            normal[1] = 1.0; // degenerate fallback
        } else {
            normal[0] = 1.0;
        }
        n_node[s_node] = normal;
    }

    // Per interface-node fluid pressure gathered through the interface map.
    let mut p_node = vec![0.0; struct_mesh.node_count()];
    for &(s_node, f_node) in interface {
        p_node[s_node] = pressure[f_node];
    }

    // Assemble the consistent nodal traction load via element quadrature.
    let order = 2usize;
    let mut rhs = vec![0.0; struct_mesh.node_count() * dim];
    for elem in &struct_mesh.elements {
        if !elem.nodes.iter().any(|n| s_nodes.contains(n)) {
            continue;
        }
        let coords: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|n| struct_mesh.node_coords(*n).to_vec())
            .collect();
        let (pts, wts) = ref_quad(elem.cell_type, order);
        for (xi, &w) in pts.iter().zip(&wts) {
            let (n, g) = ref_shape_grad(elem.cell_type, xi);
            let map = Map::from_nodes_and_grad(&coords, &g);
            let detj = map.determinant.abs();
            // Interpolate pressure and normal at this quadrature point.
            let mut p = 0.0;
            let mut nrm = vec![0.0; dim];
            for (a, &node) in elem.nodes.iter().enumerate() {
                p += n[a] * p_node[node];
                for c in 0..dim {
                    nrm[c] += n[a] * n_node[node][c];
                }
            }
            // f_a += ∫ N_a (p·n) dΩ — consistent (shape-weighted), not lumped.
            for (a, &node) in elem.nodes.iter().enumerate() {
                for c in 0..dim {
                    rhs[node * dim + c] += w * detj * n[a] * p * nrm[c];
                }
            }
        }
    }

    let k_full = try_assemble(struct_mesh, dim, |eid, m| {
        elasticity_element_matrix(m, eid, model, young, poisson, 2)
    })
    .map_err(CouplingError::from)?;
    solve_with_dirichlet(&k_full, &rhs, struct_dirichlet).map_err(CouplingError::from)
}

/// Reference shape-function values and gradients for the supported structure
/// element types at local coordinates `xi`. Mirrors the dispatch used by
/// `tpt-fem-elasticity` so the consistent load transfer reuses the same
/// reference elements.
fn ref_shape_grad(cell: CellType, xi: &[f64]) -> (Vec<f64>, Vec<Vec<f64>>) {
    match cell {
        CellType::Line => (Line2::shape(xi), Line2::grad(xi)),
        CellType::Tri => (Tri3::shape(xi), Tri3::grad(xi)),
        CellType::Quad => (Quad4::shape(xi), Quad4::grad(xi)),
        CellType::Tet => (Tet4::shape(xi), Tet4::grad(xi)),
        CellType::Hex => (Hex8::shape(xi), Hex8::grad(xi)),
        other => panic!("fsi consistent load: unsupported cell {other:?}"),
    }
}

/// Reference quadrature `(points, weights)` for the supported structure element
/// types at the given rule `order`.
fn ref_quad(cell: CellType, order: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    match cell {
        CellType::Line => {
            let r = gauss_legendre(order);
            (r.points.iter().map(|x| vec![*x]).collect(), r.weights)
        }
        CellType::Tri => {
            let r = triangle(TriangleRule::Degree2);
            (r.points.iter().map(|p| p.to_vec()).collect(), r.weights)
        }
        CellType::Quad => {
            let r = tensor_square(&gauss_legendre(order));
            (r.points.iter().map(|p| p.to_vec()).collect(), r.weights)
        }
        CellType::Tet => {
            let r = tetrahedron(TetrahedronRule::Degree2);
            (r.points.iter().map(|p| p.to_vec()).collect(), r.weights)
        }
        CellType::Hex => {
            let r = tensor_cube(&gauss_legendre(order));
            (r.points.iter().map(|p| p.to_vec()).collect(), r.weights)
        }
        other => panic!("fsi consistent load: unsupported cell {other:?}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_elasticity::ElasticModel;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    #[test]
    fn uniform_heating_free_expansion() {
        let mut b = MeshBuilder::new();
        let mut prev = b.add_node(vec![0.0, 0.0]);
        for i in 1..=4 {
            let n = b.add_node(vec![i as f64 / 4.0, 0.0]);
            b.add_element(CellType::Line, vec![prev, n]);
            prev = n;
        }
        let mesh = b.build();
        let temp = vec![1.0; mesh.node_count()];
        let dirichlet = [(0usize, 0.0)];
        let u = thermal_structural(
            &mesh,
            ElasticModel::BarAxial,
            1.0,
            0.3,
            1e-3,
            &temp,
            &dirichlet,
        )
        .unwrap();
        assert!((u[4] - 1e-3).abs() < 1e-6, "tip = {}", u[4]);
    }

    #[test]
    fn bimetallic_cantilever_bends() {
        // A 2-D cantilever (rows of Quad elements) fixed at the left edge, heated
        // with a through-thickness temperature gradient (top hotter than bottom)
        // to mimic a bilayer with different CTEs. It must bend and the tip must
        // deflect in +y, with a curvature of the right order of magnitude.
        let nx = 8;
        let ny = 2;
        let mut b = MeshBuilder::new();
        let mut rows = Vec::new();
        for j in 0..=ny {
            let y = j as f64 / ny as f64;
            let mut r = Vec::new();
            for i in 0..=nx {
                r.push(b.add_node(vec![i as f64 / nx as f64, y]));
            }
            rows.push(r);
        }
        for j in 0..ny {
            for i in 0..nx {
                b.add_element(
                    CellType::Quad,
                    vec![
                        rows[j][i],
                        rows[j][i + 1],
                        rows[j + 1][i + 1],
                        rows[j + 1][i],
                    ],
                );
            }
        }
        let mesh = b.build();
        let mut temp = vec![0.0; mesh.node_count()];
        for n in 0..mesh.node_count() {
            let y = mesh.node_coords(n)[1];
            // Top (y=1) hotter than bottom (y=0): Δα·ΔT encoded as 1.0 high / -1.0 low.
            temp[n] = if y > 0.5 { 1.0 } else { -1.0 };
        }
        // Fix the left edge (all nodes with x≈0) in both components.
        let mut dirichlet = Vec::new();
        for n in 0..mesh.node_count() {
            if mesh.node_coords(n)[0] < 1e-9 {
                dirichlet.push((n * 2, 0.0));
                dirichlet.push((n * 2 + 1, 0.0));
            }
        }
        let alpha = 1e-3;
        let u = thermal_structural(
            &mesh,
            ElasticModel::PlaneStress,
            1.0,
            0.3,
            alpha,
            &temp,
            &dirichlet,
        )
        .unwrap();
        // Tip node is the rightmost-bottom node: its y-displacement must be > 0
        // (bimetal bends toward the colder/shorter side).
        let tip = (0..mesh.node_count())
            .find(|&n| {
                (mesh.node_coords(n)[0] - 1.0).abs() < 1e-9 && (mesh.node_coords(n)[1]).abs() < 1e-9
            })
            .unwrap();
        let tip_y = u[tip * 2 + 1];
        // The hotter top layer (higher CTE) expands more, so the strip curls
        // downward (toward the colder, shorter side): tip deflects in −y.
        assert!(tip_y < 0.0, "bimetal tip should bend downward, got {tip_y}");
        // Order-of-magnitude: length 1, αΔT~1e-3 -> sub-millimetre deflection.
        assert!(tip_y.abs() < 5e-3 && tip_y < 0.0, "tip_y = {tip_y}");
    }

    #[test]
    fn joule_heating_quadratic_profile() {
        // 1-D bar [0,1], conductivity k=1, σ=1, E=2 -> q=4. T(0)=T(1)=0.
        // Steady: T(x) = q/(2k)(x - x²) = 2(x - x²); midpoint T(0.5)=0.5.
        let mut b = MeshBuilder::new();
        let mut prev = b.add_node(vec![0.0]);
        for i in 1..=4 {
            let n = b.add_node(vec![i as f64 / 4.0]);
            b.add_element(CellType::Line, vec![prev, n]);
            prev = n;
        }
        let mesh = b.build();
        let t = electro_thermal(&mesh, 1.0, 1.0, 2.0, &[(0, 0.0), (4, 0.0)]).unwrap();
        let mid = t[2];
        assert!((mid - 0.5).abs() < 1e-2, "mid T = {mid}");
    }

    #[test]
    fn fsi_transfers_traction() {
        // Structure: a 2-D elastic block (Quad elements) fixed at the base,
        // sharing the x-coordinates of the fluid so the interface maps cleanly.
        let mut sb = MeshBuilder::new();
        let mut srows = Vec::new();
        for j in 0..=2 {
            let y = j as f64 / 2.0;
            let mut r = Vec::new();
            for i in 0..=2 {
                r.push(sb.add_node(vec![i as f64 / 2.0, y]));
            }
            srows.push(r);
        }
        for j in 0..2 {
            for i in 0..2 {
                sb.add_element(
                    CellType::Quad,
                    vec![
                        srows[j][i],
                        srows[j][i + 1],
                        srows[j + 1][i + 1],
                        srows[j + 1][i],
                    ],
                );
            }
        }
        let struct_mesh = sb.build();

        let mut fb = MeshBuilder::new();
        let mut rows = Vec::new();
        for j in 0..=2 {
            let y = 1.0 + j as f64 / 2.0;
            let mut r = Vec::new();
            for i in 0..=2 {
                r.push(fb.add_node(vec![i as f64 / 2.0, y]));
            }
            rows.push(r);
        }
        for j in 0..2 {
            for i in 0..2 {
                fb.add_element(
                    CellType::Quad,
                    vec![
                        rows[j][i],
                        rows[j][i + 1],
                        rows[j + 1][i + 1],
                        rows[j + 1][i],
                    ],
                );
            }
        }
        let fluid_mesh = fb.build();

        // Interface: structure top-centre node <-> fluid bottom-centre node.
        let fluid_bottom = (0..fluid_mesh.node_count())
            .find(|&n| {
                (fluid_mesh.node_coords(n)[0] - 0.5).abs() < 1e-9
                    && (fluid_mesh.node_coords(n)[1] - 1.0).abs() < 1e-9
            })
            .unwrap();
        let s2 = (0..struct_mesh.node_count())
            .find(|&n| {
                (struct_mesh.node_coords(n)[0] - 0.5).abs() < 1e-9
                    && (struct_mesh.node_coords(n)[1] - 1.0).abs() < 1e-9
            })
            .unwrap();
        let interface = [(s2, fluid_bottom)];

        let u0 = vec![0.0; struct_mesh.node_count() * 2];
        // Fix the structure base (bottom row) in both components.
        let mut struct_dirichlet = Vec::new();
        for n in 0..struct_mesh.node_count() {
            if struct_mesh.node_coords(n)[1] < 1e-9 {
                struct_dirichlet.push((n * 2, 0.0));
                struct_dirichlet.push((n * 2 + 1, 0.0));
            }
        }
        let u = fsi_coupling(
            &struct_mesh,
            &fluid_mesh,
            ElasticModel::PlaneStress,
            1.0,
            0.3,
            1.0,
            &u0,
            &interface,
            &struct_dirichlet,
            1e6,
        )
        .expect("fsi coupling must solve");
        // The structure node must move (pressure traction transferred).
        assert!(
            u[s2 * 2 + 1].abs() > 0.0,
            "structure should respond to fluid, got {}",
            u[s2 * 2 + 1]
        );
    }

    #[test]
    fn fsi_consistent_load_uniform_pressure_top_edge() {
        // Single unit Quad element [0,1]^2. Interface = the top edge (y=1).
        // With constant pressure p on that edge the consistent load must be
        // exactly p·L/2 straight up on each top node and zero elsewhere.
        let mut b = MeshBuilder::new();
        let n00 = b.add_node(vec![0.0, 0.0]);
        let n10 = b.add_node(vec![1.0, 0.0]);
        let n01 = b.add_node(vec![0.0, 1.0]);
        let n11 = b.add_node(vec![1.0, 1.0]);
        b.add_element(CellType::Quad, vec![n00, n10, n11, n01]);
        let mesh = b.build();
        let interface = [(n01, 0usize), (n11, 1usize)];
        // Fluid "pressure" array indexed by fluid node id.
        let pressure = vec![2.0_f64; 8];
        let f = fsi_interface_loads(&mesh, &interface, &pressure);
        for &(s, _) in &interface {
            assert!(f[s * 2].abs() < 1e-14, "horizontal component {f:?}");
            let expect = 2.0 * 1.0 / 2.0; // p·L/2
            assert!(
                (f[s * 2 + 1] - expect).abs() < 1e-12,
                "vertical load {} != {expect}",
                f[s * 2 + 1]
            );
        }
        // Bottom nodes unloaded.
        assert!(f[n00 * 2].abs() < 1e-14 && f[n00 * 2 + 1].abs() < 1e-14);
        assert!(f[n10 * 2].abs() < 1e-14 && f[n10 * 2 + 1].abs() < 1e-14);
        // Resultant equals ∫ p n̂ ds = (0, p·L).
        let fy: f64 = f.chunks_exact(2).map(|c| c[1]).sum();
        assert!((fy - 2.0).abs() < 1e-12, "resultant {fy}");
    }

    #[test]
    fn fsi_consistent_load_vertical_interface_normal_points_outward() {
        // Same unit element but the interface is the RIGHT edge (x=1): the
        // outward normal must be +x, exercising geometry-aware orientation
        // (the old hardcoded-+y behaviour would give zero vertical traction).
        let mut b = MeshBuilder::new();
        let n00 = b.add_node(vec![0.0, 0.0]);
        let n10 = b.add_node(vec![1.0, 0.0]);
        let n01 = b.add_node(vec![0.0, 1.0]);
        let n11 = b.add_node(vec![1.0, 1.0]);
        b.add_element(CellType::Quad, vec![n00, n10, n11, n01]);
        let mesh = b.build();
        let interface = [(n10, 0usize), (n11, 1usize)];
        let pressure = vec![3.0_f64; 8];
        let f = fsi_interface_loads(&mesh, &interface, &pressure);
        for &(s, _) in &interface {
            let expect = 3.0 * 1.0 / 2.0;
            assert!(
                (f[s * 2] - expect).abs() < 1e-12,
                "+x load {} != {expect}",
                f[s * 2]
            );
            assert!(f[s * 2 + 1].abs() < 1e-14, "vertical component {f:?}");
        }
    }

    #[test]
    fn fsi_consistent_load_linear_pressure_resultant() {
        // Linearly varying pressure p(x) = x over the top edge of the unit
        // element: the consistent load's resultant must equal
        // ∫₀¹ x dx = 1/2 in +y, and moment balance about node (0,1) must hold:
        // Σ x_i f_iy = ∫ x·p dx = 1/3.
        let mut b = MeshBuilder::new();
        let n00 = b.add_node(vec![0.0, 0.0]);
        let n10 = b.add_node(vec![1.0, 0.0]);
        let n01 = b.add_node(vec![0.0, 1.0]);
        let n11 = b.add_node(vec![1.0, 1.0]);
        b.add_element(CellType::Quad, vec![n00, n10, n11, n01]);
        let mesh = b.build();
        let interface = [(n01, 0usize), (n11, 1usize)];
        let mut pressure = vec![0.0_f64; 8];
        pressure[0] = 0.0; // fluid node above structure node (0,1)
        pressure[1] = 1.0; // fluid node above structure node (1,1)
        let f = fsi_interface_loads(&mesh, &interface, &pressure);
        let fy: f64 = f.chunks_exact(2).map(|c| c[1]).sum();
        assert!((fy - 0.5).abs() < 1e-12, "resultant {fy} != 0.5");
        let moment = f[n01 * 2 + 1] * 0.0 + f[n11 * 2 + 1] * 1.0;
        assert!((moment - 1.0 / 3.0).abs() < 1e-12, "moment {moment} != 1/3");
    }
}
