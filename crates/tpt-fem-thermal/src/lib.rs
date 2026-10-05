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

use tpt_fem_assembly::{apply_neumann, apply_robin, solve_with_dirichlet, try_assemble};
use tpt_fem_element::{
    Hex20, Hex27, Hex8, Line2, Map, Quad4, Quad8, Quad9, ReferenceElement, Tet10, Tet4, Tri3, Tri6,
};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    tensor_cube, tensor_square, tetrahedron, triangle, try_gauss_legendre, QuadratureError,
    TetrahedronRule, TriangleRule,
};
use tpt_fem_sparse::SparseError;

/// Errors returned by this crate's quadrature-dependent element operators.
#[derive(Debug)]
pub enum ThermalError {
    /// The requested Gauss–Legendre `quad_order` is out of range.
    Quadrature(QuadratureError),
    /// The sparse linear solve of a time step failed.
    Sparse(SparseError),
    /// A caller-supplied parameter is invalid (non-positive `dt`, wrong-length
    /// initial field, `theta` outside `[0, 1]`, ...).
    InvalidInput(String),
}

impl std::fmt::Display for ThermalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThermalError::Quadrature(e) => write!(f, "{e}"),
            ThermalError::Sparse(e) => write!(f, "thermal solve failed: {e}"),
            ThermalError::InvalidInput(m) => write!(f, "thermal: invalid input: {m}"),
        }
    }
}

impl std::error::Error for ThermalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ThermalError::Quadrature(e) => Some(e),
            ThermalError::Sparse(e) => Some(e),
            ThermalError::InvalidInput(_) => None,
        }
    }
}

impl From<SparseError> for ThermalError {
    fn from(e: SparseError) -> Self {
        ThermalError::Sparse(e)
    }
}

impl From<QuadratureError> for ThermalError {
    fn from(e: QuadratureError) -> Self {
        ThermalError::Quadrature(e)
    }
}

impl From<ThermalError> for SparseError {
    fn from(e: ThermalError) -> Self {
        match e {
            ThermalError::Sparse(inner) => inner,
            other => SparseError::Numeric(other.to_string()),
        }
    }
}

/// Quadrature points (reference coordinates) and weights for a cell type.
///
/// Returns [`ThermalError`] if `order` (or, for P2 element types, `order + 1`)
/// is outside the supported `1..=5` range, rather than panicking.
fn cell_quad(cell: CellType, order: usize) -> Result<(Vec<Vec<f64>>, Vec<f64>), ThermalError> {
    Ok(match cell {
        CellType::Line => {
            let r = try_gauss_legendre(order)?;
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
            let r = tensor_square(&try_gauss_legendre(order)?);
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
            let r = tensor_cube(&try_gauss_legendre(order)?);
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
            let r = tensor_square(&try_gauss_legendre(order + 1)?);
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
            let r = tensor_cube(&try_gauss_legendre(order + 1)?);
            (
                r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect(),
                r.weights,
            )
        }
    })
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
///
/// Returns [`ThermalError`] if `quad_order` is out of range, rather than
/// panicking.
pub fn poisson_element_matrix(
    mesh: &Mesh,
    eid: usize,
    conductivity: f64,
    quad_order: usize,
) -> Result<Vec<Vec<f64>>, ThermalError> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let (qpts, qw) = cell_quad(cell, quad_order)?;
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
    Ok(k)
}

/// Element load vector from a source term `f(x)`, returned in node order.
///
/// Returns [`ThermalError`] if `quad_order` is out of range, rather than
/// panicking.
pub fn poisson_source_vector(
    mesh: &Mesh,
    eid: usize,
    source: impl Fn(&[f64]) -> f64,
    quad_order: usize,
) -> Result<Vec<f64>, ThermalError> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let (qpts, qw) = cell_quad(cell, quad_order)?;
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
    Ok(f)
}

/// Element heat-capacity (consistent mass) matrix `C_e = rho_c * int N^T N dV`,
/// returned in node order (1 DOF/node).
///
/// Returns [`ThermalError`] if `quad_order` is out of range.
pub fn heat_capacity_element_matrix(
    mesh: &Mesh,
    eid: usize,
    rho_c: f64,
    quad_order: usize,
) -> Result<Vec<Vec<f64>>, ThermalError> {
    let elem = &mesh.elements[eid];
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let cell = elem.cell_type;
    let (qpts, qw) = cell_quad(cell, quad_order)?;
    let n = elem.nodes.len();
    let mut c = vec![vec![0.0; n]; n];
    for (qp, w) in qpts.iter().zip(&qw) {
        let ns = ref_shape(cell, qp);
        let local = ref_grad(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        for i in 0..n {
            for j in 0..n {
                c[i][j] += w * rho_c * ns[i] * ns[j] * det;
            }
        }
    }
    Ok(c)
}

/// Parameters of the transient heat-conduction solver
/// [`solve_transient_heat`].
#[derive(Clone, Copy, Debug)]
pub struct TransientHeatOptions {
    /// Constant conductivity `k`.
    pub conductivity: f64,
    /// Volumetric heat capacity `rho * c`.
    pub rho_c: f64,
    /// Gauss-Legendre order for the element integrals.
    pub quad_order: usize,
    /// Time step (must be positive).
    pub dt: f64,
    /// Number of steps to take.
    pub nsteps: usize,
    /// Time-integration parameter `theta` in `[0, 1]`: `1` = backward Euler
    /// (unconditionally stable, first order), `0.5` = Crank-Nicolson (second
    /// order), `0` = forward Euler (conditionally stable).
    pub theta: f64,
}

impl Default for TransientHeatOptions {
    fn default() -> Self {
        TransientHeatOptions {
            conductivity: 1.0,
            rho_c: 1.0,
            quad_order: 2,
            dt: 0.01,
            nsteps: 10,
            theta: 1.0,
        }
    }
}

/// Transient heat conduction `rho_c * dT/dt - div(k grad T) = f(x, t)` by the
/// theta-method.
///
/// Each step solves
/// `(C/dt + theta K) T_{n+1} = (C/dt - (1-theta) K) T_n + theta f_{n+1} + (1-theta) f_n`
/// with the `dirichlet` DOFs held at their prescribed `(dof, value)`. `source`
/// is `f(x, t)`; `initial` is the nodal temperature at `t = 0`. Returns the
/// nodal temperature history `(t, T)` for steps `0..=nsteps` (step 0 is
/// `initial`, with the Dirichlet values imposed).
pub fn solve_transient_heat(
    mesh: &Mesh,
    opts: &TransientHeatOptions,
    initial: &[f64],
    source: impl Fn(&[f64], f64) -> f64,
    dirichlet: &[(usize, f64)],
) -> Result<Vec<(f64, Vec<f64>)>, ThermalError> {
    let n = mesh.node_count();
    if initial.len() != n {
        return Err(ThermalError::InvalidInput(format!(
            "initial field has {} entries, mesh has {n} nodes",
            initial.len()
        )));
    }
    if !(opts.dt.is_finite() && opts.dt > 0.0) {
        return Err(ThermalError::InvalidInput(format!(
            "dt must be finite and positive (got {})",
            opts.dt
        )));
    }
    if !(0.0..=1.0).contains(&opts.theta) {
        return Err(ThermalError::InvalidInput(format!(
            "theta must be in [0, 1] (got {})",
            opts.theta
        )));
    }
    if mesh.elements.is_empty() {
        return Err(ThermalError::InvalidInput("mesh has no elements".into()));
    }
    if let Some(&(d, _)) = dirichlet.iter().find(|(d, _)| *d >= n) {
        return Err(ThermalError::InvalidInput(format!(
            "Dirichlet DOF {d} is out of range for {n} nodes"
        )));
    }

    let k_coo = try_assemble(mesh, 1, |eid, m| {
        poisson_element_matrix(m, eid, opts.conductivity, opts.quad_order)
    })?;
    let c_coo = try_assemble(mesh, 1, |eid, m| {
        heat_capacity_element_matrix(m, eid, opts.rho_c, opts.quad_order)
    })?;
    let (dt, th) = (opts.dt, opts.theta);

    // lhs = C/dt + th*K
    let mut lhs = tpt_fem_sparse::Coo::new();
    for i in 0..c_coo.rows.len() {
        lhs.push(c_coo.rows[i], c_coo.cols[i], c_coo.vals[i] / dt);
    }
    for i in 0..k_coo.rows.len() {
        lhs.push(k_coo.rows[i], k_coo.cols[i], th * k_coo.vals[i]);
    }
    // rhs operator = C/dt - (1-th)*K, applied by matvec.
    let mut rop = tpt_fem_sparse::Coo::new();
    for i in 0..c_coo.rows.len() {
        rop.push(c_coo.rows[i], c_coo.cols[i], c_coo.vals[i] / dt);
    }
    for i in 0..k_coo.rows.len() {
        rop.push(k_coo.rows[i], k_coo.cols[i], -(1.0 - th) * k_coo.vals[i]);
    }
    let rop = rop.to_csr();

    let load = |t: f64| -> Result<Vec<f64>, ThermalError> {
        let mut f = vec![0.0; n];
        for eid in 0..mesh.elements.len() {
            let fe = poisson_source_vector(mesh, eid, |x| source(x, t), opts.quad_order)?;
            for (i, &node) in mesh.elements[eid].nodes.iter().enumerate() {
                f[node] += fe[i];
            }
        }
        Ok(f)
    };

    let mut temp = initial.to_vec();
    for &(d, v) in dirichlet {
        temp[d] = v;
    }
    let mut history = Vec::with_capacity(opts.nsteps + 1);
    history.push((0.0, temp.clone()));
    let mut f_prev = load(0.0)?;
    for step in 1..=opts.nsteps {
        let t = step as f64 * dt;
        let f_next = load(t)?;
        let mut rhs = rop.matvec(&temp);
        rhs.truncate(n);
        for i in 0..n {
            rhs[i] += th * f_next[i] + (1.0 - th) * f_prev[i];
        }
        temp = solve_with_dirichlet(&lhs, &rhs, dirichlet)?;
        history.push((t, temp.clone()));
        f_prev = f_next;
    }
    Ok(history)
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
    let mut coo = try_assemble(mesh, 1, |eid, m| {
        poisson_element_matrix(m, eid, conductivity, quad_order)
    })?;

    let mut rhs = vec![0.0; ndof];
    for eid in 0..mesh.elements.len() {
        let f = poisson_source_vector(mesh, eid, &source, quad_order)?;
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

    fn rod(nx: usize) -> Mesh {
        let mut b = MeshBuilder::new();
        let nodes: Vec<usize> = (0..=nx)
            .map(|i| b.add_node(vec![i as f64 / nx as f64]))
            .collect();
        for w in nodes.windows(2) {
            b.add_element(CellType::Line, vec![w[0], w[1]]);
        }
        b.build()
    }

    #[test]
    fn transient_heat_matches_analytic_decay() {
        // u_t = u_xx, u(0)=u(1)=0, u0 = sin(pi x)  =>  u = exp(-pi^2 t) sin(pi x).
        let nx = 40;
        let mesh = rod(nx);
        let pi = std::f64::consts::PI;
        let u0: Vec<f64> = (0..=nx)
            .map(|i| (pi * i as f64 / nx as f64).sin())
            .collect();
        let opts = TransientHeatOptions {
            dt: 0.002,
            nsteps: 50,
            theta: 0.5,
            ..TransientHeatOptions::default()
        };
        let hist =
            solve_transient_heat(&mesh, &opts, &u0, |_, _| 0.0, &[(0, 0.0), (nx, 0.0)]).unwrap();
        let (t, u) = hist.last().unwrap();
        let decay = (-pi * pi * t).exp();
        let mid = nx / 2;
        assert!((u[mid] - decay).abs() < 2e-3, "{} vs {}", u[mid], decay);
        // Monotone decay in time at the midpoint.
        assert!(hist.windows(2).all(|w| w[1].1[mid] <= w[0].1[mid] + 1e-12));
    }

    #[test]
    fn transient_heat_reaches_steady_state_with_source() {
        // u_t = u_xx + 1, u(0)=u(1)=0 relaxes to 0.5 x (1 - x).
        let nx = 20;
        let mesh = rod(nx);
        let opts = TransientHeatOptions {
            dt: 0.05,
            nsteps: 100,
            theta: 1.0,
            ..TransientHeatOptions::default()
        };
        let u0 = vec![0.0; nx + 1];
        let hist =
            solve_transient_heat(&mesh, &opts, &u0, |_, _| 1.0, &[(0, 0.0), (nx, 0.0)]).unwrap();
        let u = &hist.last().unwrap().1;
        let x = 0.5;
        assert!((u[nx / 2] - 0.5 * x * (1.0 - x)).abs() < 1e-3);
    }

    #[test]
    fn transient_heat_rejects_bad_input() {
        let mesh = rod(4);
        let ok = TransientHeatOptions::default();
        let u0 = vec![0.0; 5];
        let none = |_: &[f64], _: f64| 0.0;
        assert!(solve_transient_heat(&mesh, &ok, &[0.0; 3], none, &[]).is_err());
        let bad_dt = TransientHeatOptions { dt: 0.0, ..ok };
        assert!(solve_transient_heat(&mesh, &bad_dt, &u0, none, &[]).is_err());
        let bad_theta = TransientHeatOptions { theta: 1.5, ..ok };
        assert!(solve_transient_heat(&mesh, &bad_theta, &u0, none, &[]).is_err());
        assert!(solve_transient_heat(&mesh, &ok, &u0, none, &[(99, 0.0)]).is_err());
    }

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
        // Four Tri6 elements around an interior centre node tile the unit
        // square. `u = x` is harmonic and linear, so it lies in the P2 space;
        // with Dirichlet imposed on every boundary node (the four corners and
        // the four edge midpoints) the Galerkin solution must equal `u = x`
        // everywhere — in particular at the interior centre and the four
        // "radius" midpoints. This exercises the full P2 assembly + solve path
        // end to end through a genuine interior node (a single Tri6, or two
        // Tri6 forming a square, has no interior node, so corner-only Dirichlet
        // does not pin the edge-mid values to the linear field).
        let mut b = MeshBuilder::new();
        // Boundary corners.
        let a = b.add_node(vec![0.0, 0.0]); // u = 0
        let bp = b.add_node(vec![1.0, 0.0]); // u = 1
        let c = b.add_node(vec![1.0, 1.0]); // u = 1
        let d = b.add_node(vec![0.0, 1.0]); // u = 0
                                            // Boundary edge midpoints.
        let eab = b.add_node(vec![0.5, 0.0]); // u = 0.5
        let ebc = b.add_node(vec![1.0, 0.5]); // u = 1.0
        let ecd = b.add_node(vec![0.5, 1.0]); // u = 0.5
        let eda = b.add_node(vec![0.0, 0.5]); // u = 0.0
                                              // Interior centre node + the four "radius" midpoints.
        let o = b.add_node(vec![0.5, 0.5]); // u = 0.5
        let moa = b.add_node(vec![0.25, 0.25]); // u = 0.25
        let mob = b.add_node(vec![0.75, 0.25]); // u = 0.75
        let moc = b.add_node(vec![0.75, 0.75]); // u = 0.75
        let mod_ = b.add_node(vec![0.25, 0.75]); // u = 0.25
                                                 // Four Tri6 elements, each a corner-corner-centre triangle with the
                                                 // appropriate edge midpoints in reference (corner0, corner1, corner2,
                                                 // mid(c0-c1), mid(c1-c2), mid(c2-c0)) order.
        b.add_element(CellType::Tri6, vec![a, bp, o, eab, mob, moa]);
        b.add_element(CellType::Tri6, vec![bp, c, o, ebc, moc, mob]);
        b.add_element(CellType::Tri6, vec![c, d, o, ecd, mod_, moc]);
        b.add_element(CellType::Tri6, vec![d, a, o, eda, moa, mod_]);
        let mesh = b.build();
        let u = solve_poisson(
            &mesh,
            1.0,
            4,
            |_| 0.0,
            &[
                (a, 0.0),
                (bp, 1.0),
                (c, 1.0),
                (d, 0.0),
                (eab, 0.5),
                (ebc, 1.0),
                (ecd, 0.5),
                (eda, 0.0),
            ],
            None,
            None,
        )
        .unwrap();
        // Interior nodes must equal the exact linear field `u = x`.
        assert!((u[o] - 0.5).abs() < 1e-9, "centre got {}", u[o]);
        assert!((u[moa] - 0.25).abs() < 1e-9, "moa got {}", u[moa]);
        assert!((u[mob] - 0.75).abs() < 1e-9, "mob got {}", u[mob]);
        assert!((u[moc] - 0.75).abs() < 1e-9, "moc got {}", u[moc]);
        assert!((u[mod_] - 0.25).abs() < 1e-9, "mod got {}", u[mod_]);
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

    #[test]
    fn quad_order_out_of_range_errors_instead_of_panicking() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0]);
        let n1 = b.add_node(vec![1.0]);
        b.add_element(CellType::Line, vec![n0, n1]);
        let mesh = b.build();

        assert!(poisson_element_matrix(&mesh, 0, 1.0, 6).is_err());
        assert!(poisson_source_vector(&mesh, 0, |_| 1.0, 6).is_err());
        assert!(solve_poisson(&mesh, 1.0, 6, |_| 1.0, &[(n0, 0.0)], None, None).is_err());
    }
}
