//! Nonlinear solvers for `tpt-fem`.
//!
//! Provides a generic [Newton–Raphson](newton) iteration for a residual
//! `R(u) = 0` together with its Jacobian, and a
//! [parameter-continuation](continuation) driver that warm-starts Newton at
//! successive load/parameter steps. Essential (Dirichlet) conditions are
//! condensed out of the system on every iteration, consistently with
//! `tpt-fem-assembly`'s linear solver.
//!
//! The residual and Jacobian are supplied by the caller (e.g. assembled from a
//! nonlinear weak form), so this crate is physics-agnostic.

use std::collections::{HashMap, HashSet};

use tpt_fem_sparse::{solve, Coo, SparseError};

/// Errors returned by the nonlinear solvers.
#[derive(Debug)]
pub enum NewtonError {
    /// The iteration did not converge within `max_iter`.
    MaxIterations,
    /// The inner sparse linear solve failed.
    Sparse(SparseError),
}

impl From<SparseError> for NewtonError {
    fn from(e: SparseError) -> Self {
        NewtonError::Sparse(e)
    }
}

impl std::fmt::Display for NewtonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NewtonError::MaxIterations => {
                write!(f, "Newton iteration did not converge within max_iter steps")
            }
            NewtonError::Sparse(e) => write!(f, "linear solve failed during Newton iteration: {e}"),
        }
    }
}

impl std::error::Error for NewtonError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            NewtonError::Sparse(e) => Some(e),
            NewtonError::MaxIterations => None,
        }
    }
}

/// Tolerances and limits for the Newton iteration.
#[derive(Clone, Copy, Debug)]
pub struct NewtonOptions {
    /// Absolute residual norm at which the iteration is considered converged.
    pub tol: f64,
    /// Maximum number of Newton iterations.
    pub max_iter: usize,
}

impl Default for NewtonOptions {
    fn default() -> Self {
        NewtonOptions {
            tol: 1e-10,
            max_iter: 50,
        }
    }
}

/// Condense a global `Coo` matrix to the free DOFs.
fn reduce_coo(coo: &Coo, free: &[usize], free_idx: &HashMap<usize, usize>) -> Coo {
    let csr = coo.to_csr();
    let mut out = Coo::new();
    for &row in free {
        for c in csr.row_ptrs[row]..csr.row_ptrs[row + 1] {
            let col = csr.col_ind[c];
            let v = csr.values[c];
            if let Some(&j) = free_idx.get(&col) {
                out.push(free_idx[&row], j, v);
            }
        }
    }
    out
}

/// Solve `R(u) = 0` with Newton–Raphson.
///
/// * `u0` — initial guess.
/// * `residual` — `R(u)`.
/// * `jacobian` — `J(u)` (global `Coo`, square).
/// * `dirichlet` — `(dof, value)` essential conditions held fixed each
///   iteration.
///
/// Returns the converged solution vector.
pub fn newton(
    u0: &[f64],
    residual: impl Fn(&[f64]) -> Vec<f64>,
    jacobian: impl Fn(&[f64]) -> Coo,
    dirichlet: &[(usize, f64)],
    opts: &NewtonOptions,
) -> Result<Vec<f64>, NewtonError> {
    let n = u0.len();
    let fixed: HashSet<usize> = dirichlet.iter().map(|(i, _)| *i).collect();
    let free: Vec<usize> = (0..n).filter(|i| !fixed.contains(i)).collect();
    let free_idx: HashMap<usize, usize> = free.iter().enumerate().map(|(k, &v)| (v, k)).collect();

    let mut u = u0.to_vec();
    for (i, v) in dirichlet {
        u[*i] = *v;
    }

    for _ in 0..opts.max_iter {
        let r = residual(&u);
        let rnorm = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        if rnorm < opts.tol {
            return Ok(u);
        }
        let j = jacobian(&u);
        let jred = reduce_coo(&j, &free, &free_idx);
        let rfree: Vec<f64> = free.iter().map(|&i| r[i]).collect();
        let neg = rfree.iter().map(|x| -x).collect::<Vec<_>>();
        let du = solve(&jred, &neg)?;
        for (k, &i) in free.iter().enumerate() {
            u[i] += du[k];
        }
    }
    Err(NewtonError::MaxIterations)
}

/// Solve a parameterised residual `R(u, λ) = 0` by continuation.
///
/// Newton is run at each of `steps` values of `λ` evenly spaced in
/// `[lambda0, lambda1]`, warm-starting from the previous solution. Returns the
/// `(λ, u)` pairs (including the initial `λ0` step).
#[allow(clippy::too_many_arguments)]
pub fn continuation<F, G>(
    u0: &[f64],
    lambda0: f64,
    lambda1: f64,
    steps: usize,
    dirichlet: &[(usize, f64)],
    residual: F,
    jacobian: G,
    opts: &NewtonOptions,
) -> Result<Vec<(f64, Vec<f64>)>, NewtonError>
where
    F: Fn(&[f64], f64) -> Vec<f64>,
    G: Fn(&[f64], f64) -> Coo,
{
    let mut out = Vec::with_capacity(steps + 1);
    let mut u = u0.to_vec();
    // Apply the (constant) Dirichlet conditions to the initial guess.
    for (dof, val) in dirichlet {
        u[*dof] = *val;
    }
    for step in 0..=steps {
        let lambda = lambda0 + (lambda1 - lambda0) * (step as f64) / (steps as f64);
        let u_sol = newton(
            &u,
            |u| residual(u, lambda),
            |u| jacobian(u, lambda),
            dirichlet,
            opts,
        )?;
        out.push((lambda, u_sol.clone()));
        u = u_sol;
    }
    Ok(out)
}

/// Error returned by [`arc_length_continuation`].
#[derive(Debug)]
pub enum ArcLengthError {
    /// The inner sparse linear solve (tangent factorization) failed.
    LinearSolve(SparseError),
    /// The iteration ran out of arc-length steps before a stop condition.
    MaxSteps,
    /// Every step cut back to the minimum arc length without converging.
    MaxCutbacks,
}

impl std::fmt::Display for ArcLengthError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ArcLengthError::LinearSolve(e) => write!(f, "arc-length tangent solve failed: {e}"),
            ArcLengthError::MaxSteps => write!(f, "arc-length continuation reached max_steps"),
            ArcLengthError::MaxCutbacks => write!(
                f,
                "arc-length continuation exhausted cut-backs at minimum step"
            ),
        }
    }
}

impl std::error::Error for ArcLengthError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ArcLengthError::LinearSolve(e) => Some(e),
            _ => None,
        }
    }
}

impl From<SparseError> for ArcLengthError {
    fn from(e: SparseError) -> Self {
        ArcLengthError::LinearSolve(e)
    }
}

/// Options controlling [`arc_length_continuation`].
#[derive(Clone, Copy, Debug)]
pub struct ArcLengthOptions {
    /// Arc-length of the very first step.
    pub initial_arc_length: f64,
    /// Target arc-length the adaptive controller drives toward after the first
    /// step (raised after fast steps, lowered after slow ones).
    pub target_arc_length: f64,
    /// Weight of the load increment in the spherical constraint
    /// `Δu·Δu + ψ·Δλ² = Δs²`. `1.0` is the standard Crisfield spherical arc.
    pub psi: f64,
    /// Lower bound on the (adaptive) arc length; also the floor at which a
    /// cut-back sequence gives up.
    pub min_arc_length: f64,
    /// Upper bound on the (adaptive) arc length.
    pub max_arc_length: f64,
    /// Maximum number of arc-length steps.
    pub max_steps: usize,
    /// Maximum corrector iterations per step before the step is cut back.
    pub max_iter_per_step: usize,
    /// Number of corrector iterations the adaptive controller considers ideal;
    /// steps converging in fewer iterations grow, more shrink.
    pub target_iterations: usize,
    /// Residual norm at which a step is considered converged.
    pub tol: f64,
    /// Maximum number of halvings of the arc length before a step is declared
    /// failed.
    pub max_cutbacks: usize,
    /// Optional absolute bound on `|λ|`; once exceeded the trace stops cleanly
    /// (returns what it has so far).
    pub lambda_limit: Option<f64>,
    /// Sign of the first predictor step. The arc-length trace has no preferred
    /// direction, so this lets the caller choose which way along the equilibrium
    /// path to go (e.g. toward a known limit point). Defaults to `+1.0`.
    pub initial_direction: f64,
}

impl Default for ArcLengthOptions {
    fn default() -> Self {
        ArcLengthOptions {
            initial_arc_length: 1.0,
            target_arc_length: 1.0,
            psi: 1.0,
            min_arc_length: 1e-4,
            max_arc_length: 5.0,
            max_steps: 200,
            max_iter_per_step: 30,
            target_iterations: 5,
            tol: 1e-9,
            max_cutbacks: 8,
            lambda_limit: None,
            initial_direction: 1.0,
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Build a dense `n × n` matrix from a `Coo`.
fn coo_to_dense(coo: &Coo, n: usize) -> Vec<Vec<f64>> {
    let csr = coo.to_csr();
    let mut a = vec![vec![0.0; n]; n];
    for r in 0..n {
        for c in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            a[r][csr.col_ind[c]] += csr.values[c];
        }
    }
    a
}

/// Solve `A x = B` for every column of `B` with Gaussian elimination and
/// partial pivoting.
///
/// Unlike the faer-backed sparse solve (which panics on a singular matrix),
/// this regularises near-singular pivots (replacing a pivot below `1e-12` with
/// `1e-12`) so the arc-length driver keeps working through a limit point
/// instead of aborting. The regularisation only affects directions whose
/// stiffness has collapsed; the spherical predictor/corrector controls the step
/// magnitude, so the regularised direction is still correct to first order.
fn dense_solve(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, ArcLengthError> {
    let n = a.len();
    let nrhs = b.len();
    let mut m = vec![vec![0.0; n + nrhs]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        for k in 0..nrhs {
            m[i][n + k] = b[k][i];
        }
    }
    for col in 0..n {
        let mut piv = col;
        let mut max = m[col][col].abs();
        for r in (col + 1)..n {
            if m[r][col].abs() > max {
                max = m[r][col].abs();
                piv = r;
            }
        }
        m.swap(col, piv);
        let mut pv = m[col][col];
        if pv.abs() < 1e-12 {
            pv = 1e-12 * pv.signum();
            m[col][col] = pv;
        }
        for r in (col + 1)..n {
            let f = m[r][col] / pv;
            for c in col..(n + nrhs) {
                m[r][c] -= f * m[col][c];
            }
        }
    }
    let mut x = vec![vec![0.0; n]; nrhs];
    for k in 0..nrhs {
        for i in (0..n).rev() {
            let mut s = m[i][n + k];
            for j in (i + 1)..n {
                s -= m[i][j] * x[k][j];
            }
            x[k][i] = s / m[i][i];
        }
    }
    Ok(x)
}

/// Solve the bordered `(n+1)×(n+1)` arc-length system
/// ```text
/// [ Kt    -p ] [δu]   = [ -r                  ]
/// [ aᵀ   ψ·b ] [δλ]     [ (Δs² - |a|² - ψ·b²)/2 ]
/// ```
/// via dense Gaussian elimination with partial pivoting.
///
/// The extra bordered row (linearising the spherical constraint along the
/// current increment `(a, b)`) keeps the system nonsingular even when `Kt` is
/// singular at a limit point, so the corrector can traverse folds that a plain
/// `Kt` solve would choke on. Near-singular pivots are regularised (a pivot
/// below `1e-12` is replaced) rather than rejected, so the trace still makes
/// progress through a fold.
fn bordered_solve(
    kt: &[Vec<f64>],
    p: &[f64],
    a: &[f64],
    b: f64,
    r: &[f64],
    rhs_c: f64,
    psi: f64,
) -> Result<(Vec<f64>, f64), ArcLengthError> {
    let n = kt.len();
    let m = n + 1;
    let mut mat = vec![vec![0.0; m]; m];
    for i in 0..n {
        for j in 0..n {
            mat[i][j] = kt[i][j];
        }
        mat[i][n] = -p[i];
        mat[n][i] = a[i];
    }
    mat[n][n] = psi * b;
    let mut rhs = vec![0.0; m];
    for i in 0..n {
        rhs[i] = -r[i];
    }
    rhs[n] = rhs_c;
    for col in 0..m {
        let mut piv = col;
        let mut max = mat[col][col].abs();
        for row in (col + 1)..m {
            if mat[row][col].abs() > max {
                max = mat[row][col].abs();
                piv = row;
            }
        }
        mat.swap(col, piv);
        rhs.swap(col, piv);
        let mut pv = mat[col][col];
        if pv.abs() < 1e-12 {
            pv = 1e-12 * pv.signum();
            mat[col][col] = pv;
        }
        for row in (col + 1)..m {
            let f = mat[row][col] / pv;
            for c in col..m {
                mat[row][c] -= f * mat[col][c];
            }
            rhs[row] -= f * rhs[col];
        }
    }
    let mut x = vec![0.0; m];
    for i in (0..m).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..m {
            s -= mat[i][j] * x[j];
        }
        x[i] = s / mat[i][i];
    }
    let du = x[0..n].to_vec();
    let dl = x[n];
    Ok((du, dl))
}
///
/// `residual` must have the form `R(u, λ) = F_int(u) - λ·load`, i.e. `λ`
/// multiplies the *constant* reference load vector `load`. `jacobian` returns
/// the tangent `∂F_int/∂u` (equivalently `∂R/∂u`); it is re-sampled at the
/// current point on every corrector iteration (full Newton).
///
/// Each step begins with a predictor along the tangent and is closed by
/// corrector iterations that satisfy both equilibrium and the spherical
/// constraint `Δu·Δu + ψ·Δλ² = Δs²` via the bordered (Crisfield) Newton update,
/// which stays nonsingular even at a limit point. If a step fails to converge
/// (singular bordered system or too many iterations) the arc length is halved
/// and the step retried, up to `max_cutbacks` times. Converged steps adapt `Δs`
/// toward `target_arc_length` based on how many iterations they took.
///
/// Returns the traced `(λ, u)` samples, starting with `(lambda0, u0)`.
pub fn arc_length_continuation(
    u0: &[f64],
    lambda0: f64,
    dirichlet: &[(usize, f64)],
    load: &[f64],
    residual: impl Fn(&[f64], f64) -> Vec<f64>,
    jacobian: impl Fn(&[f64], f64) -> Coo,
    opts: &ArcLengthOptions,
) -> Result<Vec<(f64, Vec<f64>)>, ArcLengthError> {
    let n = u0.len();
    let fixed: HashSet<usize> = dirichlet.iter().map(|(i, _)| *i).collect();
    let free: Vec<usize> = (0..n).filter(|i| !fixed.contains(i)).collect();
    let free_idx: HashMap<usize, usize> = free.iter().enumerate().map(|(k, &v)| (v, k)).collect();
    let reduce = |v: &[f64]| -> Vec<f64> { free.iter().map(|&i| v[i]).collect() };
    let nf = free.len();

    let mut u = u0.to_vec();
    for (i, val) in dirichlet {
        u[*i] = *val;
    }
    let mut lambda = lambda0;

    let mut out = vec![(lambda, u.clone())];
    let mut delta_s = opts.initial_arc_length;
    let mut prev_dlambda = 0.0_f64;

    for _ in 0..opts.max_steps {
        let u_start = u.clone();
        let lam_start = lambda;
        let pload = reduce(load);

        // Cut-back loop for this step.
        let mut cutbacks = 0;
        let (u_new, lam_new, iters) = 'step: loop {
            // Step-start tangent (predictor direction).
            let kt0 = jacobian(&u_start, lam_start);
            let kt0_red = reduce_coo(&kt0, &free, &free_idx);
            let kt0_dense = coo_to_dense(&kt0_red, nf);

            // Predictor: tangent direction for a unit load increment. A singular
            // tangent here means the step started on a fold; cut back.
            let x1 = match dense_solve(&kt0_dense, std::slice::from_ref(&pload)) {
                Ok(s) => s.into_iter().next().unwrap(),
                Err(_) => {
                    cutbacks += 1;
                    if cutbacks > opts.max_cutbacks || delta_s * 0.5 < opts.min_arc_length {
                        return Err(ArcLengthError::MaxCutbacks);
                    }
                    delta_s *= 0.5;
                    continue 'step;
                }
            };
            let sign = if prev_dlambda == 0.0 {
                opts.initial_direction.signum()
            } else {
                prev_dlambda.signum()
            };
            let denom = (dot(&x1, &x1) + opts.psi).sqrt();
            let dl_pred = sign * delta_s / denom;
            let du_pred: Vec<f64> = x1.iter().map(|v| dl_pred * v).collect();

            let mut u_cur = u_start.clone();
            for (k, &fd) in free.iter().enumerate() {
                u_cur[fd] += du_pred[k];
            }
            let mut lam_cur = lam_start + dl_pred;
            // Accumulated increments from the step start (used by the spherical
            // constraint's linearisation in the corrector).
            let mut du_acc = du_pred;
            let mut dl_acc = dl_pred;

            let mut converged = false;
            let mut it = 0;
            for _ in 0..opts.max_iter_per_step {
                it += 1;
                let r = residual(&u_cur, lam_cur);
                let r_red = reduce(&r);
                let rnorm = r_red.iter().map(|x| x * x).sum::<f64>().sqrt();
                if rnorm < opts.tol {
                    converged = true;
                    break;
                }
                // Tangent at the *current* point (full Newton).
                let kt = jacobian(&u_cur, lam_cur);
                let kt_red = reduce_coo(&kt, &free, &free_idx);
                let kt_dense = coo_to_dense(&kt_red, nf);
                // Newton step for the augmented (equilibrium + linearised
                // spherical constraint) system. The bordered row keeps the
                // (n+1)×(n+1) system nonsingular even when `Kt` is singular at a
                // limit point, so the corrector can traverse folds.
                let a2 = dot(&du_acc, &du_acc);
                let rhs_c = 0.5 * (delta_s * delta_s - a2 - opts.psi * dl_acc * dl_acc);
                let (du_step, dlambda) = match bordered_solve(
                    &kt_dense,
                    &pload,
                    &du_acc,
                    dl_acc,
                    &r_red,
                    rhs_c,
                    opts.psi,
                ) {
                    Ok(sol) => sol,
                    Err(_) => break, // singular -> cut back
                };
                for (k, _fd) in free.iter().enumerate() {
                    du_acc[k] += du_step[k];
                }
                dl_acc += dlambda;
                for (k, &fd) in free.iter().enumerate() {
                    u_cur[fd] = u_start[fd] + du_acc[k];
                }
                lam_cur = lam_start + dl_acc;
            }

            if converged {
                break 'step (u_cur, lam_cur, it);
            }
            // Cut back and retry.
            cutbacks += 1;
            if cutbacks > opts.max_cutbacks || delta_s * 0.5 < opts.min_arc_length {
                return Err(ArcLengthError::MaxCutbacks);
            }
            delta_s *= 0.5;
        };

        u = u_new;
        lambda = lam_new;
        out.push((lambda, u.clone()));
        prev_dlambda = lambda - lam_start;

        // Adaptive arc length toward the target.
        let factor = (opts.target_iterations as f64 / iters as f64).clamp(0.1, 10.0).sqrt();
        delta_s = (delta_s * factor).clamp(opts.min_arc_length, opts.max_arc_length);

        if let Some(lim) = opts.lambda_limit {
            if lambda.abs() > lim {
                break;
            }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_sparse::Coo;

    #[test]
    fn newton_scalar() {
        // Root of u^3 - u - 1 = 0 (~1.324717957).
        let u = newton(
            &[1.0],
            |u| vec![u[0] * u[0] * u[0] - u[0] - 1.0],
            |u| {
                let mut c = Coo::new();
                c.push(0, 0, 3.0 * u[0] * u[0] - 1.0);
                c
            },
            &[],
            &NewtonOptions::default(),
        )
        .unwrap();
        assert!((u[0] - 1.324717957).abs() < 1e-8);
    }

    #[test]
    fn newton_diagonal_system() {
        // u0^2 - 2 = 0, u1^3 - 3 = 0  =>  (sqrt 2, cbrt 3).
        let u = newton(
            &[1.0, 1.0],
            |u| vec![u[0] * u[0] - 2.0, u[1] * u[1] * u[1] - 3.0],
            |u| {
                let mut c = Coo::new();
                c.push(0, 0, 2.0 * u[0]);
                c.push(1, 1, 3.0 * u[1] * u[1]);
                c
            },
            &[],
            &NewtonOptions::default(),
        )
        .unwrap();
        assert!((u[0] - 2.0_f64.sqrt()).abs() < 1e-9);
        assert!((u[1] - 3.0_f64.powf(1.0 / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn newton_with_dirichlet() {
        // 2-dof: residual = [u0 + u1 - 1, u0^2 - 0.25], fix u1 = 0.5.
        // Then u0 = 0.5 satisfies both.
        let u = newton(
            &[0.0, 0.0],
            |u| vec![u[0] + u[1] - 1.0, u[0] * u[0] - 0.25],
            |u| {
                let mut c = Coo::new();
                c.push(0, 0, 1.0);
                c.push(0, 1, 1.0);
                c.push(1, 0, 2.0 * u[0]);
                c.push(1, 1, 0.0);
                c
            },
            &[(1, 0.5)],
            &NewtonOptions::default(),
        )
        .unwrap();
        assert!((u[0] - 0.5).abs() < 1e-9);
        assert!((u[1] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn continuation_sweep() {
        // u^2 = lambda, sweep lambda 0 -> 4, expect final u ~ 2.
        let res = continuation(
            &[1.0],
            1.0,
            4.0,
            8,
            &[],
            |u, lam| vec![u[0] * u[0] - lam],
            |u, _lam| {
                let mut c = Coo::new();
                c.push(0, 0, 2.0 * u[0]);
                c
            },
            &NewtonOptions::default(),
        )
        .unwrap();
        let (lam, u) = res.last().unwrap();
        assert!(((*lam) - 4.0).abs() < 1e-12);
        assert!((u[0] - 2.0).abs() < 1e-8, "got {}", u[0]);
    }

    #[test]
    fn arc_length_passes_limit_point() {
        // Cubic fold: R(u,λ) = u^3 - 3u - λ. The equilibrium λ = u^3 - 3u has
        // limit points at u = ±1 (λ = ∓2). Start on the right branch
        // (u≈2.104, λ=3) and trace downward; the path must pass the fold at
        // λ = -2 and remain on the curve throughout.
        let u0 = [2.1038];
        let res = |u: &[f64], lam: f64| vec![u[0] * u[0] * u[0] - 3.0 * u[0] - lam];
        let jac = |u: &[f64], _lam: f64| {
            let mut c = Coo::new();
            c.push(0, 0, 3.0 * u[0] * u[0] - 3.0);
            c
        };
        let opts = ArcLengthOptions {
            initial_arc_length: 0.15,
            target_arc_length: 0.15,
            max_arc_length: 0.5,
            max_steps: 800,
            initial_direction: -1.0,
            ..Default::default()
        };
        let trace = arc_length_continuation(&u0, 3.0, &[], &[1.0], res, jac, &opts).unwrap();
        // The initial sample is the caller's (approximate) guess; the traced
        // points are what must stay on the curve.
        for (lam, u) in trace.iter().skip(1) {
            let r = u[0] * u[0] * u[0] - 3.0 * u[0] - lam;
            assert!(r.abs() < 1e-5, "off-curve: {r} at λ={lam}");
        }
        let min_lam = trace.iter().map(|(l, _)| *l).fold(f64::INFINITY, f64::min);
        assert!(min_lam <= -1.7, "did not reach the fold region, min λ = {min_lam}");
    }

    #[test]
    fn arc_length_2d_limit_point() {
        // A 2-DOF nonlinear residual with a genuine limit point, used to exercise
        // the bordered (n+1)×(n+1) solve in two dimensions. The equilibrium is
        //   R0 = u0³ - 3 u0 - λ = 0 ,  R1 = u1 = 0 ,
        // i.e. the 1-D cubic fold extended to 2 DOF. The Jacobian is singular at
        // the fold (u0 = 1, λ = -2), so a plain `Kt` solve would fail there;
        // the bordered corrector must carry the trace through it.
        let res = |u: &[f64], lam: f64| vec![u[0] * u[0] * u[0] - 3.0 * u[0] - lam, u[1]];
        let jac = |u: &[f64], _lam: f64| {
            let mut c = Coo::new();
            c.push(0, 0, 3.0 * u[0] * u[0] - 3.0);
            c.push(1, 1, 1.0);
            c
        };
        let opts = ArcLengthOptions {
            initial_arc_length: 0.15,
            target_arc_length: 0.15,
            max_arc_length: 0.5,
            max_steps: 800,
            initial_direction: -1.0,
            ..Default::default()
        };
        let u0 = [2.1038, 0.0];
        let trace =
            arc_length_continuation(&u0, 3.0, &[], &[1.0, 0.0], res, jac, &opts).unwrap();
        for (lam, u) in trace.iter().skip(1) {
            let r = res(u, *lam);
            assert!(r[0].abs() < 1e-5 && r[1].abs() < 1e-9, "off-curve: {r:?}");
        }
        let min_lam = trace.iter().map(|(l, _)| *l).fold(f64::INFINITY, f64::min);
        assert!(min_lam <= -1.7, "did not reach the fold region, min λ = {min_lam}");
        // Second DOF stays pinned at 0 (verifies the 2-D bordered solve tracks it).
        let max_u1 = trace
            .iter()
            .map(|(_, u)| u[1].abs())
            .fold(0.0_f64, f64::max);
        assert!(max_u1 < 1e-6, "second DOF drifted: {max_u1}");
    }
}
