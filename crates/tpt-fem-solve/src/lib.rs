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
}
