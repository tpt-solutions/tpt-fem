//! Sparse eigenvalue solvers for `tpt-fem`.
//!
//! Provides:
//!
//! * [`power_iteration`] — the dominant (largest-magnitude) eigenpair,
//! * [`inverse_iteration`] — the eigenpair nearest a target shift (shift-invert),
//! * [`lanczos_eigs`] — a few extreme eigenpairs via Lanczos tridiagonalisation
//!   with a dense Jacobi eigensolve of the projected matrix.
//!
//! All routines operate on a [`Coo`] matrix from `tpt-fem-sparse` and are
//! physics-agnostic.

use tpt_fem_sparse::{solve, Coo, SparseError};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Matrix–vector product `y = A x` for a [`Coo`] matrix.
pub fn matvec(coo: &Coo, x: &[f64]) -> Vec<f64> {
    let csr = coo.to_csr();
    let mut y = vec![0.0; csr.nrows];
    for r in 0..csr.nrows {
        let mut s = 0.0;
        for c in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            s += csr.values[c] * x[csr.col_ind[c]];
        }
        y[r] = s;
    }
    y
}

/// Rayleigh quotient `xᵀ A x / xᵀ x`.
pub fn rayleigh(coo: &Coo, x: &[f64]) -> f64 {
    dot(x, &matvec(coo, x)) / dot(x, x)
}

/// Return `A - shift·I` as a new [`Coo`].
fn shifted(coo: &Coo, shift: f64) -> Coo {
    let csr = coo.to_csr();
    let n = csr.nrows;
    let mut out = Coo::new();
    for r in 0..n {
        let mut diag = 0.0;
        for c in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            let col = csr.col_ind[c];
            let v = csr.values[c];
            if col == r {
                diag = v - shift;
            } else {
                out.push(r, col, v);
            }
        }
        out.push(r, r, diag);
    }
    out
}

/// Compute the dominant (largest-magnitude) eigenpair by the power method.
///
/// Returns `(λ, x)` with `||x|| = 1`. `max_iter` and `tol` bound the iteration
/// (converged when the relative change in `λ` drops below `tol`).
pub fn power_iteration(coo: &Coo, max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = coo.to_csr().nrows;
    let mut x: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let s = norm(&x);
    x.iter_mut().for_each(|v| *v /= s);
    let mut lambda = 0.0;
    for _ in 0..max_iter {
        let y = matvec(coo, &x);
        let ny = norm(&y);
        for i in 0..n {
            x[i] = y[i] / ny;
        }
        let new = rayleigh(coo, &x);
        if (new - lambda).abs() < tol * new.abs().max(1.0) {
            lambda = new;
            break;
        }
        lambda = new;
    }
    (lambda, x)
}

/// Compute the eigenpair nearest `shift` by inverse iteration (shift-invert).
///
/// Each step solves `(A - shift·I) x = x_prev` with `tpt-fem-sparse`, so the
/// iteration converges to the eigenvalue closest to `shift`.
pub fn inverse_iteration(
    coo: &Coo,
    shift: f64,
    max_iter: usize,
    tol: f64,
) -> Result<(f64, Vec<f64>), SparseError> {
    let n = coo.to_csr().nrows;
    let a_shift = shifted(coo, shift);
    let mut x: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let s = norm(&x);
    x.iter_mut().for_each(|v| *v /= s);
    let mut lambda = rayleigh(coo, &x);
    for _ in 0..max_iter {
        let y = solve(&a_shift, &x)?;
        let ny = norm(&y);
        if ny < 1e-14 {
            break;
        }
        for i in 0..n {
            x[i] = y[i] / ny;
        }
        let new = rayleigh(coo, &x);
        if (new - lambda).abs() < tol * new.abs().max(1.0) {
            lambda = new;
            break;
        }
        lambda = new;
    }
    Ok((lambda, x))
}

/// Which end of the spectrum to extract.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EigWhich {
    /// Smallest eigenvalues.
    Smallest,
    /// Largest eigenvalues.
    Largest,
}

/// Lanczos tridiagonalisation of `A`, returning the orthonormal basis
/// `V` (`n × m`), the diagonal `α` (`m`) and sub-/super-diagonal `β` (`m-1`).
fn lanczos(coo: &Coo, m: usize) -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    let n = coo.to_csr().nrows;
    let mut v0 = vec![0.0; n];
    v0[0] = 1.0;
    let nv = norm(&v0);
    v0.iter_mut().for_each(|x| *x /= nv);
    let mut basis = vec![v0];
    let mut alpha = Vec::new();
    let mut beta = Vec::new();
    for j in 0..m {
        let mut w = matvec(coo, &basis[j]);
        let a = dot(&w, &basis[j]);
        // Full re-orthogonalisation against all previous Lanczos vectors.
        for k in 0..=j {
            let p = dot(&w, &basis[k]);
            for d in 0..n {
                w[d] -= p * basis[k][d];
            }
        }
        alpha.push(a);
        if j + 1 == m {
            break;
        }
        let b = norm(&w);
        if b < 1e-12 {
            break;
        }
        beta.push(b);
        basis.push(w.iter().map(|x| x / b).collect());
    }
    (basis, alpha, beta)
}

/// Eigen-decomposition of a dense symmetric matrix by the cyclic Jacobi method.
/// Returns the eigenvalues and the eigenvectors as columns.
fn jacobi(mut a: Vec<Vec<f64>>, max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    for _ in 0..max_sweeps {
        let mut changed = false;
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() > 1e-12 {
                    changed = true;
                    let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                    let t = if theta >= 0.0 {
                        1.0 / (theta + (1.0 + theta * theta).sqrt())
                    } else {
                        1.0 / (theta - (1.0 + theta * theta).sqrt())
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    for i in 0..n {
                        let aip = a[i][p];
                        let aiq = a[i][q];
                        a[i][p] = c * aip - s * aiq;
                        a[i][q] = s * aip + c * aiq;
                    }
                    for i in 0..n {
                        let api = a[p][i];
                        let aqi = a[q][i];
                        a[p][i] = c * api - s * aqi;
                        a[q][i] = s * api + c * aqi;
                    }
                    for i in 0..n {
                        let vip = v[i][p];
                        let viq = v[i][q];
                        v[i][p] = c * vip - s * viq;
                        v[i][q] = s * vip + c * viq;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    let eig = (0..n).map(|i| a[i][i]).collect();
    (eig, v)
}

/// Compute `num` extreme eigenpairs of a symmetric `Coo` by Lanczos.
///
/// Builds a Lanczos basis of dimension `lanczos_dim`, projects `A` onto the
/// tridiagonal `T`, and returns the Ritz pairs from a Jacobi eigensolve of `T`.
pub fn lanczos_eigs(
    coo: &Coo,
    num: usize,
    which: EigWhich,
    lanczos_dim: usize,
) -> Vec<(f64, Vec<f64>)> {
    let (basis, alpha, beta) = lanczos(coo, lanczos_dim);
    let m = basis.len();
    let mut t = vec![vec![0.0; m]; m];
    for i in 0..m {
        t[i][i] = alpha[i];
    }
    for i in 0..(m - 1) {
        t[i][i + 1] = beta[i];
        t[i + 1][i] = beta[i];
    }
    let (eig, q) = jacobi(t, 200);
    let mut idx: Vec<usize> = (0..m).collect();
    idx.sort_by(|&a, &b| {
        eig[a]
            .partial_cmp(&eig[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let take = idx.iter().take(num).cloned().collect::<Vec<usize>>();
    let order: Vec<usize> = match which {
        EigWhich::Smallest => take,
        EigWhich::Largest => idx.iter().rev().take(num).cloned().collect(),
    };
    order
        .into_iter()
        .map(|i| {
            let mut vec = vec![0.0; basis[0].len()];
            for (r, vr) in basis.iter().enumerate() {
                let c = q[r][i];
                for (d, v) in vr.iter().enumerate() {
                    vec[d] += c * v;
                }
            }
            let nv = norm(&vec);
            for v in vec.iter_mut() {
                *v /= nv;
            }
            (eig[i], vec)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sym2(a: f64, b: f64) -> Coo {
        // [[a, b], [b, a]]
        let mut c = Coo::new();
        c.push(0, 0, a);
        c.push(0, 1, b);
        c.push(1, 0, b);
        c.push(1, 1, a);
        c
    }

    #[test]
    fn power_iteration_dominant() {
        // [[2,1],[1,2]] has eigenvalues 1 and 3; dominant is 3.
        let coo = sym2(2.0, 1.0);
        let (lam, v) = power_iteration(&coo, 200, 1e-12);
        assert!((lam - 3.0).abs() < 1e-8, "got {lam}");
        let av = matvec(&coo, &v);
        let res = av
            .iter()
            .zip(&v)
            .map(|(a, b)| a - lam * b)
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(res < 1e-6, "eigenvector residual {res}");
    }

    #[test]
    fn inverse_iteration_smallest() {
        // Shift 0 -> eigenpair nearest 0, i.e. smallest eigenvalue 1.
        let coo = sym2(2.0, 1.0);
        let (lam, v) = inverse_iteration(&coo, 0.0, 200, 1e-12).unwrap();
        assert!((lam - 1.0).abs() < 1e-8, "got {lam}");
        let av = matvec(&coo, &v);
        let res = av
            .iter()
            .zip(&v)
            .map(|(a, b)| a - lam * b)
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(res < 1e-6, "eigenvector residual {res}");
    }

    #[test]
    fn lanczos_1d_laplacian() {
        // Discrete 1-D Laplacian (n=10): eigenvalues 2 - 2 cos(kπ/(n+1)).
        let n = 10;
        let mut c = Coo::new();
        for i in 0..n {
            c.push(i, i, 2.0);
            if i + 1 < n {
                c.push(i, i + 1, -1.0);
                c.push(i + 1, i, -1.0);
            }
        }
        let smallest = lanczos_eigs(&c, 1, EigWhich::Smallest, n);
        let expected = 2.0 - 2.0 * (std::f64::consts::PI / (n as f64 + 1.0)).cos();
        assert!(
            (smallest[0].0 - expected).abs() < 1e-6,
            "got {}",
            smallest[0].0
        );

        let largest = lanczos_eigs(&c, 1, EigWhich::Largest, n);
        let expected_max = 2.0 - 2.0 * (std::f64::consts::PI * n as f64 / (n as f64 + 1.0)).cos();
        assert!(
            (largest[0].0 - expected_max).abs() < 1e-6,
            "got {}",
            largest[0].0
        );
    }
}
