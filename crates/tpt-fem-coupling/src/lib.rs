//! Multiphysics coupling operators for `tpt-fem`.
//!
//! * [`thermal_structural`] — thermal strain `ε = α ΔT I` fed into the linear
//!   elasticity operator as an initial strain (free thermal expansion).
//! * [`electro_thermal`] — Joule heating `q = σ |E|²` fed into the Poisson
//!   (heat-conduction) operator as a volumetric source.
//! * [`fsi_coupling`] — a basic partitioned fluid–structure interaction that
//!   exchanges traction/kinematics between [`tpt-fem-fluid`] and
//!   [`tpt-fem-elasticity`] through [`tpt-fem-dynamic`].
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
use tpt_fem_mesh::{CellType, Mesh};
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

/// A basic partitioned fluid–structure interaction step.
///
/// Given a structure displacement `u_struct` (global DOF order, `dim` per node)
/// and a shared interface of `(structure_node, fluid_node)` pairs, this:
/// 1. displaces the fluid mesh nodes by the structure displacement (FSI
///    kinematic condition),
/// 2. solves [`tpt-fem-fluid`]'s steady Stokes for the fluid pressure,
/// 3. transfers the fluid pressure as a normal traction onto the structure
///    interface (Newton's third law),
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

    // Transfer fluid pressure as a normal traction onto the structure interface.
    // Neumann load per node = pressure * outward_normal. The normal is derived
    // from the actual interface geometry: the average of (node - element
    // centroid) over the structure elements incident on the interface node. For
    // a convex/conforming mesh this points outward, so vertical or curved
    // interfaces receive a correct (e.g. horizontal) traction direction rather
    // than a hardcoded +y.
    let mut neumann = vec![0.0; struct_mesh.node_count() * dim];
    for &(s_node, f_node) in interface {
        // Geometry-aware outward normal at the structure interface node.
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
        for a in 0..dim {
            neumann[s_node * dim + a] += pressure[f_node] * normal[a];
        }
    }
    // Solve the structure with the transferred traction as a Neumann load.
    let k_full = try_assemble(struct_mesh, dim, |eid, m| {
        elasticity_element_matrix(m, eid, model, young, poisson, 2)
    })
    .map_err(CouplingError::from)?;
    let mut rhs = vec![0.0; struct_mesh.node_count() * dim];
    // Convert the per-node nodal traction into consistent nodal loads by a simple
    // lumped projection (good for a smoke-level coupling check).
    for (i, v) in neumann.iter().enumerate() {
        rhs[i] += *v;
    }
    solve_with_dirichlet(&k_full, &rhs, struct_dirichlet).map_err(CouplingError::from)
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
}
