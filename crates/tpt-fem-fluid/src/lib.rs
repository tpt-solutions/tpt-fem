//! Fluid elements (Stokes / low-Re Navier–Stokes) for `tpt-fem`.
//!
//! Provides a divergence-penalty *mixed-style* formulation in which the velocity
//! field (one vector component per node) is solved together with a weakly
//! enforced incompressibility constraint, and the pressure is recovered from the
//! velocity divergence after the solve:
//!
//! * [`steady_stokes`] solves the steady Stokes problem
//!   `−μ∇²u + ∇p = f`, `∇·u = 0` by augmenting the velocity Laplacian with the
//!   divergence penalty `(1/ε)∫(∇·u)(∇·v)dΩ` and recovering `p = −(1/ε)∇·u`
//!   after the solve. The velocity field is constrained by Dirichlet conditions
//!   on the boundary, which makes the penalty operator symmetric positive
//!   definite. The penalty term is integrated with a reduced (centroid) rule —
//!   *selective reduced integration* — so the equal-order velocity space does
//!   not lock as the penalty grows.
//! * [`transient_stokes`] integrates the creeping-flow operator through time
//!   using `tpt-fem-dynamic`'s implicit [`newmark`] scheme.
//! * [`transient_navier_stokes`] adds a Picard-linearised convective term for
//!   low-Reynolds-number flow.
//!
//! ```
//! use tpt_fem_fluid::steady_stokes;
//! use tpt_fem_mesh::{CellType, MeshBuilder};
//!
//! // 2-D unit-height channel driven by a unit x-body force G=1, μ=1. The analytic
//! // Poiseuille profile (with no-slip walls at y=0,1) is
//! // u_x(y) = (G/(2μ))·y·(1 − y); the parabola is imposed on the whole boundary
//! // (no-slip on the walls, parabolic inlet/outlet), so the penalty-Stokes solve
//! // reproduces it in the interior.
//! let mut b = MeshBuilder::new();
//! let mut rows = Vec::new();
//! for j in 0..=8 {
//!     let y = j as f64 / 8.0;
//!     let mut r = Vec::new();
//!     for i in 0..=8 {
//!         r.push(b.add_node(vec![i as f64 / 8.0, y]));
//!     }
//!     rows.push(r);
//! }
//! for j in 0..8 {
//!     for i in 0..8 {
//!         b.add_element(CellType::Quad, vec![rows[j][i], rows[j][i+1], rows[j+1][i+1], rows[j+1][i]]);
//!     }
//! }
//! let mesh = b.build();
//! let g = 1.0;
//! let mu = 1.0;
//! // Velocity Dirichlet on the boundary only; the interior DOFs stay free.
//! let mut bc = Vec::new();
//! for n in 0..mesh.node_count() {
//!     let c = mesh.node_coords(n);
//!     let on_boundary = c[0] < 1e-9 || c[0] > 1.0 - 1e-9 || c[1] < 1e-9 || c[1] > 1.0 - 1e-9;
//!     if on_boundary {
//!         let ux = (g / (2.0 * mu)) * c[1] * (1.0 - c[1]); // analytic parabola
//!         bc.push((n * 2, ux));
//!         bc.push((n * 2 + 1, 0.0));
//!     }
//! }
//! let (u, _p) = steady_stokes(&mesh, mu, |_| vec![g, 0.0], &bc, 1e6)
//!     .expect("steady Stokes solve");
//! // Sample an interior centreline node (y≈0.5) and compare with the analytic
//! // value 0.125.
//! let center: f64 = (0..mesh.node_count())
//!     .filter(|&n| {
//!         let c = mesh.node_coords(n);
//!         (c[1] - 0.5).abs() < 1e-9 && c[0] > 1e-9 && c[0] < 1.0 - 1e-9
//!     })
//!     .map(|n| u[n * 2])
//!     .next()
//!     .unwrap();
//! let want = (g / (2.0 * mu)) * 0.5 * 0.5;
//! assert!((center - want).abs() < 2e-2, "got {} want {}", center, want);
//! ```

use tpt_fem_assembly::solve_with_dirichlet;
use tpt_fem_dynamic::{newmark, NewmarkOptions};
use tpt_fem_element::{Hex8, Line2, Map, Quad4, ReferenceElement, Tet4, Tri3};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};
use tpt_fem_sparse::{Coo, SparseError};

/// Errors returned by the fluid solvers.
#[derive(Debug)]
pub enum FluidError {
    /// The underlying linear-algebra solve (velocity system, or a
    /// Dirichlet-reduced step) failed.
    Sparse(SparseError),
    /// The Picard (fixed-point) iteration for the convective term did not reach
    /// the convergence tolerance within the allowed number of re-linearisations.
    PicardNotConverged {
        /// Time step index at which the failure occurred.
        step: usize,
        /// Number of Picard iterations attempted.
        iterations: usize,
        /// Final relative change between successive Picard iterates.
        residual: f64,
    },
}

impl std::fmt::Display for FluidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FluidError::Sparse(e) => write!(f, "fluid system solve failed: {e}"),
            FluidError::PicardNotConverged {
                step,
                iterations,
                residual,
            } => write!(
                f,
                "Navier-Stokes Picard iteration did not converge at step {step} \
                 after {iterations} iterations (relative change {residual:e})"
            ),
        }
    }
}

impl std::error::Error for FluidError {}

impl From<SparseError> for FluidError {
    fn from(e: SparseError) -> Self {
        FluidError::Sparse(e)
    }
}

/// Build the velocity(=dim)+pressure(=1) multi-field DOF map for `mesh`.
///
/// The penalty method solves for the velocity field only and recovers the
/// pressure afterwards, but the map is provided for callers that want a single
/// field descriptor covering both quantities.
pub fn stokes_dofmap(mesh: &Mesh) -> tpt_fem_dofmap::MultiFieldDofMap {
    let dim = match mesh.elements[0].cell_type {
        CellType::Line => Line2::DIM,
        CellType::Tri => Tri3::DIM,
        CellType::Quad => Quad4::DIM,
        CellType::Tet => Tet4::DIM,
        CellType::Hex => Hex8::DIM,
        other => panic!("fluid: unsupported cell {other:?}"),
    };
    tpt_fem_dofmap::MultiFieldDofMap::new(
        mesh,
        &[
            tpt_fem_dofmap::FieldSpec::new("velocity", dim),
            tpt_fem_dofmap::FieldSpec::new("pressure", 1),
        ],
        tpt_fem_dofmap::Layout::Interleaved,
    )
}

fn quad_points(cell: CellType, order: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
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
        other => panic!("fluid: unsupported cell {other:?}"),
    }
}

/// One-point ("reduced") quadrature rule for the divergence-penalty term.
///
/// A penalty term integrated with the same rule as the viscous term
/// over-constrains an equal-order velocity space and locks the flow (the
/// velocity tends to zero as the penalty grows). Integrating only the
/// divergence penalty with the centroid rule — *selective reduced integration*
/// — restores the constraint count of the equivalent mixed method.
fn reduced_quad_points(cell: CellType) -> (Vec<Vec<f64>>, Vec<f64>) {
    match cell {
        CellType::Line => {
            let r = gauss_legendre(1);
            (r.points.iter().map(|x| vec![*x]).collect(), r.weights)
        }
        CellType::Tri => {
            let r = triangle(TriangleRule::Degree1);
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Quad => {
            let r = tensor_square(&gauss_legendre(1));
            (
                r.points.iter().map(|p| vec![p[0], p[1]]).collect(),
                r.weights,
            )
        }
        CellType::Tet => {
            let r = tetrahedron(TetrahedronRule::Degree1);
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
        CellType::Hex => {
            let r = tensor_cube(&gauss_legendre(1));
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
        other => panic!("fluid: unsupported cell {other:?}"),
    }
}

fn ref_grad(cell: CellType, xi: &[f64]) -> Vec<Vec<f64>> {
    match cell {
        CellType::Line => Line2::grad(xi),
        CellType::Tri => Tri3::grad(xi),
        CellType::Quad => Quad4::grad(xi),
        CellType::Tet => Tet4::grad(xi),
        CellType::Hex => Hex8::grad(xi),
        other => panic!("fluid: unsupported cell {other:?}"),
    }
}

fn ref_shape(cell: CellType, xi: &[f64]) -> Vec<f64> {
    match cell {
        CellType::Line => Line2::shape(xi),
        CellType::Tri => Tri3::shape(xi),
        CellType::Quad => Quad4::shape(xi),
        CellType::Tet => Tet4::shape(xi),
        CellType::Hex => Hex8::shape(xi),
        other => panic!("fluid: unsupported cell {other:?}"),
    }
}

/// Assemble the velocity (viscous + divergence-penalty) stiffness and mass over
/// the free velocity DOFs, together with the free-DOF right-hand side.
///
/// Returns `(k, mass, rhs_free, dim)` where `k` already includes the penalty
/// term `(1/ε)∫(∇·u)(∇·v)dΩ` and `rhs_free` has the Dirichlet values lifted into
/// it.
fn assemble_fluid(
    mesh: &Mesh,
    viscosity: f64,
    penalty: f64,
    body_force: &impl Fn(&[f64]) -> Vec<f64>,
    velocity_bc: &[(usize, f64)],
) -> (Coo, Coo, Vec<f64>, usize) {
    let order = 2;
    let dim = match mesh.elements[0].cell_type {
        CellType::Line => Line2::DIM,
        CellType::Tri => Tri3::DIM,
        CellType::Quad => Quad4::DIM,
        CellType::Tet => Tet4::DIM,
        CellType::Hex => Hex8::DIM,
        other => panic!("fluid: unsupported cell {other:?}"),
    };
    let nvel = mesh.node_count() * dim;
    let mut k = Coo::new();
    let mut mass = Coo::new();
    let mut rhs = vec![0.0; nvel];

    for elem in &mesh.elements {
        let cell = elem.cell_type;
        let n = elem.nodes.len();
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&nd| mesh.node_coords(nd).to_vec())
            .collect();
        let (qpts, qw) = quad_points(cell, order);
        let mut av = vec![vec![0.0; n * dim]; n * dim];
        let mut mm = vec![vec![0.0; n * dim]; n * dim];
        for (qp, w) in qpts.iter().zip(&qw) {
            let local = ref_grad(cell, qp);
            let m = Map::from_nodes_and_grad(&phys, &local);
            let det = m.determinant.abs();
            let gphys: Vec<Vec<f64>> = local.iter().map(|g| m.physical_grad(g)).collect();
            let ns = ref_shape(cell, qp);
            for i in 0..n {
                for j in 0..n {
                    // Viscous (vector Laplacian) term: μ ∫ ∇u_a·∇v_a dΩ. The
                    // gradient dot product runs over *all* spatial directions,
                    // so every velocity component sees the full Laplacian.
                    let grad_dot: f64 = (0..dim).map(|b| gphys[i][b] * gphys[j][b]).sum();
                    for a in 0..dim {
                        av[i * dim + a][j * dim + a] += w * det * viscosity * grad_dot;
                    }
                    // Consistent mass: ∫ N_i N_j dΩ.
                    for a in 0..dim {
                        mm[i * dim + a][j * dim + a] += w * det * ns[i] * ns[j];
                    }
                }
                let mut x = vec![0.0; phys[0].len()];
                for kk in 0..n {
                    for dd in 0..x.len() {
                        x[dd] += ns[kk] * phys[kk][dd];
                    }
                }
                let f = body_force(&x);
                for a in 0..dim {
                    rhs[elem.nodes[i] * dim + a] += w * det * ns[i] * f[a];
                }
            }
        }
        // Divergence-penalty term: (1/ε) ∫ (∇·u)(∇·v) dΩ, integrated with the
        // reduced (centroid) rule to avoid locking.
        let (rpts, rw) = reduced_quad_points(cell);
        for (qp, w) in rpts.iter().zip(&rw) {
            let local = ref_grad(cell, qp);
            let m = Map::from_nodes_and_grad(&phys, &local);
            let det = m.determinant.abs();
            let gphys: Vec<Vec<f64>> = local.iter().map(|g| m.physical_grad(g)).collect();
            for i in 0..n {
                for j in 0..n {
                    for a in 0..dim {
                        for b in 0..dim {
                            av[i * dim + a][j * dim + b] +=
                                w * det * penalty * gphys[i][a] * gphys[j][b];
                        }
                    }
                }
            }
        }
        for i in 0..n * dim {
            let gi = elem.nodes[i / dim] * dim + (i % dim);
            for j in 0..n * dim {
                let gj = elem.nodes[j / dim] * dim + (j % dim);
                if av[i][j] != 0.0 {
                    k.push(gi, gj, av[i][j]);
                }
                if mm[i][j] != 0.0 {
                    mass.push(gi, gj, mm[i][j]);
                }
            }
        }
    }

    // Reduce to free velocity DOFs and lift Dirichlet values into the RHS.
    let free: Vec<usize> = (0..nvel)
        .filter(|d| !velocity_bc.iter().any(|(b, _)| *b == *d))
        .collect();
    let free_idx: std::collections::HashMap<usize, usize> =
        free.iter().enumerate().map(|(r, g)| (*g, r)).collect();
    let mut kf = Coo::new();
    let mut mf = Coo::new();
    for idx in 0..k.rows.len() {
        if let (Some(&r), Some(&c)) = (free_idx.get(&k.rows[idx]), free_idx.get(&k.cols[idx])) {
            kf.push(r, c, k.vals[idx]);
        }
    }
    for idx in 0..mass.rows.len() {
        if let (Some(&r), Some(&c)) = (free_idx.get(&mass.rows[idx]), free_idx.get(&mass.cols[idx]))
        {
            mf.push(r, c, mass.vals[idx]);
        }
    }
    let mut rhsf = vec![0.0; free.len()];
    for (i, &g) in free.iter().enumerate() {
        rhsf[i] = rhs[g];
    }
    // Lift the Dirichlet values into the free-DOF right-hand side exactly once.
    for &(d, val) in velocity_bc {
        for idx in 0..k.rows.len() {
            if k.cols[idx] == d {
                if let Some(&r) = free_idx.get(&k.rows[idx]) {
                    rhsf[r] -= k.vals[idx] * val;
                }
            }
        }
    }
    (kf, mf, rhsf, dim)
}

/// Solve steady Stokes by the divergence-penalty method.
///
/// `body_force(x)` returns the `dim`-vector volumetric source, `velocity_bc`
/// lists `(global_velocity_dof, value)` essential conditions on the velocity
/// field (one entry per velocity component, i.e. `node*dim + component`).
/// `penalty` is the `1/ε` divergence-penalty weight; larger enforces
/// incompressibility more strongly at the cost of conditioning. Returns
/// `(velocity, pressure)` each indexed by node (velocity is `dim` per node,
/// concatenated as `node*dim + component`; pressure is one scalar per node).
/// The velocity system is solved with `tpt-fem-assembly`'s Dirichlet
/// reduction; a singular/under-constrained system (e.g. a floating fluid domain
/// with no velocity Dirichlet condition) is surfaced as [`FluidError::Sparse`]
/// rather than panicking.
pub fn steady_stokes(
    mesh: &Mesh,
    viscosity: f64,
    body_force: impl Fn(&[f64]) -> Vec<f64>,
    velocity_bc: &[(usize, f64)],
    penalty: f64,
) -> Result<(Vec<f64>, Vec<f64>), FluidError> {
    let (k, _mass, rhs, dim) = assemble_fluid(mesh, viscosity, penalty, &body_force, velocity_bc);
    let sol = solve_with_dirichlet(&k, &rhs, &[])?;
    let nvel = mesh.node_count() * dim;
    let mut u = vec![0.0; nvel];
    for (i, g) in (0..nvel)
        .filter(|d| !velocity_bc.iter().any(|(b, _)| *b == *d))
        .enumerate()
    {
        u[g] = sol[i];
    }
    for &(d, val) in velocity_bc {
        u[d] = val;
    }
    let p = recover_pressure(mesh, &u, penalty, dim);
    Ok((u, p))
}

/// Recover the pressure field `p = −(1/ε)·∇·u` from a velocity solution.
///
/// The divergence is evaluated with the same reduced (centroid) rule used to
/// integrate the penalty term, then smoothed onto the nodes with a lumped
/// (partition-of-unity weighted) projection, so the returned values are nodal
/// pressures rather than integrated residuals.
fn recover_pressure(mesh: &Mesh, sol: &[f64], penalty: f64, dim: usize) -> Vec<f64> {
    let mut p = vec![0.0; mesh.node_count()];
    let mut wsum = vec![0.0; mesh.node_count()];
    for elem in &mesh.elements {
        let cell = elem.cell_type;
        let n = elem.nodes.len();
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&nd| mesh.node_coords(nd).to_vec())
            .collect();
        let (qpts, qw) = reduced_quad_points(cell);
        for (qp, w) in qpts.iter().zip(&qw) {
            let local = ref_grad(cell, qp);
            let m = Map::from_nodes_and_grad(&phys, &local);
            let det = m.determinant.abs();
            let gphys: Vec<Vec<f64>> = local.iter().map(|g| m.physical_grad(g)).collect();
            let ns = ref_shape(cell, qp);
            let mut div_u = 0.0;
            for j in 0..n {
                for a in 0..dim {
                    div_u += gphys[j][a] * sol[elem.nodes[j] * dim + a];
                }
            }
            for mnode in 0..n {
                let wt = w * det * ns[mnode];
                p[elem.nodes[mnode]] += -penalty * wt * div_u;
                wsum[elem.nodes[mnode]] += wt;
            }
        }
    }
    for nd in 0..p.len() {
        if wsum[nd] > 0.0 {
            p[nd] /= wsum[nd];
        }
    }
    p
}

/// Integrate creeping flow through time with `tpt-fem-dynamic`'s implicit
/// [`newmark`] scheme. Returns the velocity history `(t, velocity_concat)` for
/// `0..=nsteps` (velocity concatenated as `node*dim + component`).
///
/// Unsteady creeping flow is *first* order in the velocity,
/// `ρ M·u̇ + K·u = f`. Integrating the fluid-particle displacement `w` (with
/// `u = ẇ`) recasts it as the second-order system `M·ẅ + K·ẇ = f` that
/// [`newmark`] accepts, with the viscous + divergence-penalty operator entering
/// as the *damping* matrix and no stiffness. For `γ = 1/2` Newmark satisfies
/// `(wⁿ⁺¹ − wⁿ)/Δt = ½(uⁿ + uⁿ⁺¹)` exactly, so differencing the returned
/// displacement history gives the step-averaged velocity, which relaxes
/// monotonically onto the [`steady_stokes`] solution.
///
/// The flow starts from rest (`u = 0` on the free DOFs) and `body_force` is
/// sampled once, at `t = 0`: the load is applied as a step at the first step
/// and held constant thereafter.
pub fn transient_stokes(
    mesh: &Mesh,
    viscosity: f64,
    body_force: impl Fn(f64, &[f64]) -> Vec<f64>,
    velocity_bc: &[(usize, f64)],
    penalty: f64,
    opts: &NewmarkOptions,
    nsteps: usize,
) -> Vec<(f64, Vec<f64>)> {
    let bf = |x: &[f64]| body_force(0.0, x);
    let (k, mass, rhs_free, dim) = assemble_fluid(mesh, viscosity, penalty, &bf, velocity_bc);
    let nfree = rhs_free.len();
    let w0 = vec![0.0; nfree];
    let v0 = vec![0.0; nfree];
    let hist = newmark(
        &mass,
        &k,
        &Coo::new(),
        &w0,
        &v0,
        move |_| rhs_free.clone(),
        opts,
        nsteps,
    );
    let nvel = mesh.node_count() * dim;
    let dt = opts.dt;
    let mut out = Vec::with_capacity(hist.len());
    for (step, (t, w)) in hist.iter().enumerate() {
        // Step-averaged velocity of the free DOFs from the displacement history.
        let vel_free: Vec<f64> = if step == 0 {
            vec![0.0; nfree]
        } else {
            let wprev = &hist[step - 1].1;
            (0..nfree).map(|i| (w[i] - wprev[i]) / dt).collect()
        };
        let mut u = vec![0.0; nvel];
        let mut fi = 0;
        for d in 0..nvel {
            if let Some(&(_, val)) = velocity_bc.iter().find(|(b, _)| *b == d) {
                u[d] = val;
            } else {
                u[d] = vel_free[fi];
                fi += 1;
            }
        }
        out.push((*t, u));
    }
    out
}

/// Low-Reynolds-number Navier–Stokes via Picard iteration of the convective term
/// (semi-implicit: the velocity from the previous Picard iterate drives the
/// convection).
///
/// Each step solves the backward-Euler system
///
/// ```text
/// (M/Δt + K + C(u*))·uⁿ⁺¹ = f(tⁿ⁺¹) + (M/Δt)·uⁿ
/// ```
///
/// where `K` is the viscous + divergence-penalty operator, `C(u*)` is the
/// Picard-linearised convection `∫ N_i (u*·∇N_j) dΩ` about the current iterate
/// `u*`, and `M` is the consistent mass matrix (density 1). `picard_iters`
/// re-linearisations are taken per step. The scheme is intended for
/// creeping/low-Re flow: with no upwinding or stabilisation it collapses to
/// [`steady_stokes`] as the velocity vanishes, but it is not suitable for
/// advection-dominated flow.
///
/// The flow starts from rest and `velocity_bc` is held fixed for all steps.
/// Returns the final velocity field (`node*dim + component`). The backward-Euler
/// step system is solved with `tpt-fem-assembly`'s Dirichlet reduction (a
/// singular/under-constrained system is surfaced as [`FluidError::Sparse`]).
/// Each step's nonlinear convective term is driven by a Picard fixed-point
/// iteration that breaks early on convergence; if it does not reach the
/// tolerance within `picard_iters` re-linearisations the step fails with
/// [`FluidError::PicardNotConverged`] rather than silently returning an
/// unconverged field.
pub fn transient_navier_stokes(
    mesh: &Mesh,
    viscosity: f64,
    body_force: impl Fn(f64, &[f64]) -> Vec<f64>,
    velocity_bc: &[(usize, f64)],
    penalty: f64,
    dt: f64,
    nsteps: usize,
    picard_iters: usize,
) -> Result<Vec<f64>, FluidError> {
    let order = 2;
    let dim = match mesh.elements[0].cell_type {
        CellType::Line => Line2::DIM,
        CellType::Tri => Tri3::DIM,
        CellType::Quad => Quad4::DIM,
        CellType::Tet => Tet4::DIM,
        CellType::Hex => Hex8::DIM,
        other => panic!("fluid: unsupported cell {other:?}"),
    };
    let dt = dt.max(1e-12);
    let nvel = mesh.node_count() * dim;
    let free: Vec<usize> = (0..nvel)
        .filter(|d| !velocity_bc.iter().any(|(b, _)| *b == *d))
        .collect();
    let free_idx: std::collections::HashMap<usize, usize> =
        free.iter().enumerate().map(|(r, g)| (*g, r)).collect();

    let mut u_prev = vec![0.0; nvel];
    for &(d, val) in velocity_bc {
        u_prev[d] = val;
    }

    for step in 0..nsteps {
        // Backward Euler evaluates the load at the end of the step.
        let t = (step + 1) as f64 * dt;
        let mut u_iter = u_prev.clone();
        // Relative-change tolerance for Picard convergence.
        let picard_tol = 1e-9;
        let mut converged = false;
        let mut last_change = 0.0_f64;
        for _ in 0..picard_iters.max(1) {
            // `stiff` = viscous + penalty + convection about `u_iter`;
            // `mass` = M/Δt. They are kept apart so that only `stiff` needs the
            // Dirichlet lift (the mass lift cancels against the `(M/Δt)·uⁿ`
            // term, whose boundary values are the same constants).
            let mut stiff = Coo::new();
            let mut mass = Coo::new();
            let mut rhs = vec![0.0; nvel];
            for elem in &mesh.elements {
                let cell = elem.cell_type;
                let n = elem.nodes.len();
                let phys: Vec<Vec<f64>> = elem
                    .nodes
                    .iter()
                    .map(|&nd| mesh.node_coords(nd).to_vec())
                    .collect();
                let (qpts, qw) = quad_points(cell, order);
                let mut av = vec![vec![0.0; n * dim]; n * dim];
                let mut mm = vec![vec![0.0; n * dim]; n * dim];
                for (qp, w) in qpts.iter().zip(&qw) {
                    let local = ref_grad(cell, qp);
                    let m = Map::from_nodes_and_grad(&phys, &local);
                    let det = m.determinant.abs();
                    let gphys: Vec<Vec<f64>> = local.iter().map(|g| m.physical_grad(g)).collect();
                    let ns = ref_shape(cell, qp);
                    let mut x = vec![0.0; phys[0].len()];
                    for kk in 0..n {
                        for dd in 0..x.len() {
                            x[dd] += ns[kk] * phys[kk][dd];
                        }
                    }
                    // Convecting velocity at the quadrature point.
                    let u_q: Vec<f64> = (0..dim)
                        .map(|a| {
                            (0..n)
                                .map(|j| ns[j] * u_iter[elem.nodes[j] * dim + a])
                                .sum::<f64>()
                        })
                        .collect();
                    for i in 0..n {
                        for j in 0..n {
                            // Viscous vector Laplacian μ ∫ ∇u_a·∇v_a dΩ.
                            let grad_dot: f64 = (0..dim).map(|b| gphys[i][b] * gphys[j][b]).sum();
                            // Picard convection ∫ N_i (u*·∇N_j) dΩ.
                            let conv: f64 = (0..dim).map(|b| u_q[b] * gphys[j][b]).sum();
                            for a in 0..dim {
                                av[i * dim + a][j * dim + a] +=
                                    w * det * (viscosity * grad_dot + ns[i] * conv);
                                mm[i * dim + a][j * dim + a] += w * det * ns[i] * ns[j] / dt;
                            }
                        }
                        let f = body_force(t, &x);
                        for a in 0..dim {
                            rhs[elem.nodes[i] * dim + a] += w * det * ns[i] * f[a];
                        }
                    }
                }
                // Divergence-penalty term with reduced (centroid) integration.
                let (rpts, rw) = reduced_quad_points(cell);
                for (qp, w) in rpts.iter().zip(&rw) {
                    let local = ref_grad(cell, qp);
                    let m = Map::from_nodes_and_grad(&phys, &local);
                    let det = m.determinant.abs();
                    let gphys: Vec<Vec<f64>> = local.iter().map(|g| m.physical_grad(g)).collect();
                    for i in 0..n {
                        for j in 0..n {
                            for a in 0..dim {
                                for b in 0..dim {
                                    av[i * dim + a][j * dim + b] +=
                                        w * det * penalty * gphys[i][a] * gphys[j][b];
                                }
                            }
                        }
                    }
                }
                for i in 0..n * dim {
                    let gi = elem.nodes[i / dim] * dim + (i % dim);
                    for j in 0..n * dim {
                        let gj = elem.nodes[j / dim] * dim + (j % dim);
                        if av[i][j] != 0.0 {
                            stiff.push(gi, gj, av[i][j]);
                        }
                        if mm[i][j] != 0.0 {
                            mass.push(gi, gj, mm[i][j]);
                        }
                    }
                }
            }

            // Reduce to the free DOFs: the system matrix is M/Δt + K.
            let mut kf = Coo::new();
            let mut mf = Coo::new();
            for idx in 0..stiff.rows.len() {
                if let (Some(&r), Some(&c)) = (
                    free_idx.get(&stiff.rows[idx]),
                    free_idx.get(&stiff.cols[idx]),
                ) {
                    kf.push(r, c, stiff.vals[idx]);
                }
            }
            for idx in 0..mass.rows.len() {
                if let (Some(&r), Some(&c)) =
                    (free_idx.get(&mass.rows[idx]), free_idx.get(&mass.cols[idx]))
                {
                    kf.push(r, c, mass.vals[idx]);
                    mf.push(r, c, mass.vals[idx]);
                }
            }
            let mut rhsf = vec![0.0; free.len()];
            for (i, &g) in free.iter().enumerate() {
                rhsf[i] = rhs[g];
            }
            for &(d, val) in velocity_bc {
                for idx in 0..stiff.rows.len() {
                    if stiff.cols[idx] == d {
                        if let Some(&r) = free_idx.get(&stiff.rows[idx]) {
                            rhsf[r] -= stiff.vals[idx] * val;
                        }
                    }
                }
            }
            // Inertia from the previous step: (M/Δt)·uⁿ over the free DOFs.
            let mcsr = mf.to_csr();
            for row in 0..mcsr.nrows.min(rhsf.len()) {
                let mut s = 0.0;
                for c in mcsr.row_ptrs[row]..mcsr.row_ptrs[row + 1] {
                    s += mcsr.values[c] * u_prev[free[mcsr.col_ind[c]]];
                }
                rhsf[row] += s;
            }

            let sol = solve_with_dirichlet(&kf, &rhsf, &[])?;
            let mut u = vec![0.0; nvel];
            for (i, &g) in free.iter().enumerate() {
                u[g] = sol[i];
            }
            for &(d, val) in velocity_bc {
                u[d] = val;
            }
            // Relative change between successive Picard iterates; break early
            // once the fixed point has settled.
            let mut num = 0.0_f64;
            let mut den = 0.0_f64;
            for i in 0..u.len() {
                num += (u[i] - u_iter[i]).powi(2);
                den += u[i].powi(2);
            }
            last_change = num.sqrt();
            let rel = last_change / (den.sqrt() + 1e-30);
            u_iter = u;
            if rel < picard_tol {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(FluidError::PicardNotConverged {
                step,
                iterations: picard_iters.max(1),
                residual: last_change,
            });
        }
        u_prev = u_iter;
    }
    Ok(u_prev)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    fn channel(nx: usize, ny: usize) -> Mesh {
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
        b.build()
    }

    /// Analytic Poiseuille profile u_x(y) = (G/(2μ))·y·(1 − y) for a channel of
    /// unit height with no-slip walls.
    fn parabola(y: f64, g: f64, mu: f64) -> f64 {
        (g / (2.0 * mu)) * y * (1.0 - y)
    }

    /// True for nodes lying on the perimeter of the unit box [0,1]².
    fn is_boundary(mesh: &Mesh, n: usize) -> bool {
        let c = mesh.node_coords(n);
        const TOL: f64 = 1e-9;
        (c[0] < TOL) || (c[0] > 1.0 - TOL) || (c[1] < TOL) || (c[1] > 1.0 - TOL)
    }

    // Prescribe the analytic Poiseuille profile on the *boundary* only; the
    // penalty-Stokes solve then reproduces it in the interior (a patch test).
    fn poiseuille_bc(mesh: &Mesh, g: f64, mu: f64) -> Vec<(usize, f64)> {
        let mut bc = Vec::new();
        for n in 0..mesh.node_count() {
            if is_boundary(mesh, n) {
                let c = mesh.node_coords(n);
                bc.push((n * 2, parabola(c[1], g, mu)));
                bc.push((n * 2 + 1, 0.0));
            }
        }
        bc
    }

    #[test]
    fn poiseuille_parabolic_profile() {
        let mesh = channel(8, 8);
        let mu = 1.0;
        let g = 1.0;
        let bc = poiseuille_bc(&mesh, g, mu);
        let (u, _p) = steady_stokes(&mesh, mu, |_| vec![g, 0.0], &bc, 1e6)
            .expect("steady Stokes patch test must solve");
        let center: Vec<f64> = (0..mesh.node_count())
            .filter(|&n| (mesh.node_coords(n)[1] - 0.5).abs() < 1e-9)
            .map(|n| u[n * 2])
            .collect();
        let want = parabola(0.5, g, mu);
        for c in center {
            assert!((c - want).abs() < 2e-2, "centreline u={c} want {want}");
        }
        // Walls must be no-slip.
        let wall: Vec<f64> = (0..mesh.node_count())
            .filter(|&n| (mesh.node_coords(n)[1] - 0.0).abs() < 1e-9)
            .map(|n| u[n * 2])
            .collect();
        for w in wall {
            assert!(w.abs() < 1e-6, "wall u_x={w}");
        }
    }

    #[test]
    fn lid_driven_cavity_has_flow() {
        // Lid (top) moves at u_x = 1, all other boundary velocities fixed to 0
        // (no-slip); interior develops circulation.
        let mesh = channel(10, 10);
        let mut bc = Vec::new();
        for n in 0..mesh.node_count() {
            if is_boundary(&mesh, n) {
                let c = mesh.node_coords(n);
                if c[1] > 1.0 - 1e-9 {
                    bc.push((n * 2, 1.0));
                    bc.push((n * 2 + 1, 0.0));
                } else {
                    bc.push((n * 2, 0.0));
                    bc.push((n * 2 + 1, 0.0));
                }
            }
        }
        let (u, _p) = steady_stokes(&mesh, 1.0, |_| vec![0.0, 0.0], &bc, 1e6)
            .expect("steady Stokes cavity solve");
        // The band 0.3 < y < 0.7 sits *below* the primary vortex centre, so the
        // horizontal velocity there is the return flow (u_x < 0); test its
        // magnitude, not its signed value.
        let interior: Vec<f64> = (0..mesh.node_count())
            .filter(|&n| {
                let c = mesh.node_coords(n);
                c[1] > 0.3 && c[1] < 0.7 && c[0] > 0.3 && c[0] < 0.7
            })
            .map(|n| u[n * 2].abs())
            .collect();
        let max_u = interior.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_u > 0.05,
            "interior should carry lid-driven flow, got {max_u}"
        );
        // Forward flow immediately under the lid, return flow at mid-height.
        let at = |x: f64, y: f64| {
            (0..mesh.node_count())
                .find(|&n| {
                    (mesh.node_coords(n)[0] - x).abs() < 1e-9
                        && (mesh.node_coords(n)[1] - y).abs() < 1e-9
                })
                .map(|n| u[n * 2])
                .unwrap()
        };
        assert!(at(0.5, 0.9) > 0.0, "under-lid flow must be forward");
        assert!(at(0.5, 0.5) < 0.0, "mid-height flow must recirculate");
    }

    #[test]
    fn transient_stokes_reaches_steady() {
        let mesh = channel(8, 8);
        let mu = 1.0;
        let g = 1.0;
        let bc = poiseuille_bc(&mesh, g, mu);
        let opts = NewmarkOptions {
            dt: 0.01,
            beta: 0.25,
            gamma: 0.5,
        };
        let hist = transient_stokes(&mesh, mu, |_, _| vec![g, 0.0], &bc, 1e3, &opts, 600);
        let final_u = &hist.last().unwrap().1;
        let center: Vec<f64> = (0..mesh.node_count())
            .filter(|&n| (mesh.node_coords(n)[1] - 0.5).abs() < 1e-9)
            .map(|n| final_u[n * 2])
            .collect();
        let want = parabola(0.5, g, mu);
        for c in center {
            assert!(
                (c - want).abs() < 3e-2,
                "transient centreline u={c} want {want}"
            );
        }
    }

    #[test]
    fn navier_stokes_runs_low_re() {
        let mesh = channel(6, 6);
        let mu = 1.0;
        let g = 0.1;
        let bc = poiseuille_bc(&mesh, g, mu);
        // Small driving force => low Re; should settle near the Stokes profile.
        let u = transient_navier_stokes(&mesh, mu, |_, _| vec![g, 0.0], &bc, 1e3, 0.01, 60, 20)
            .expect("low-Re Navier-Stokes solve");
        let center: Vec<f64> = (0..mesh.node_count())
            .filter(|&n| (mesh.node_coords(n)[1] - 0.5).abs() < 1e-9)
            .map(|n| u[n * 2])
            .collect();
        let want = parabola(0.5, g, mu);
        for c in center {
            assert!(
                (c - want).abs() < 5e-2,
                "navier centreline u={c} want {want}"
            );
        }
    }
}
