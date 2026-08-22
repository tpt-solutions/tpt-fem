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

/// Cholesky factor `M = L Lᵀ` (lower `L`). Returns `None` if `M` is not
/// positive-definite.
fn cholesky(m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = m[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i][j] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(l)
}

fn solve_lower(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i][k] * x[k];
        }
        x[i] = s / l[i][i];
    }
    x
}

fn solve_upper(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in (i + 1)..n {
            s -= l[k][i] * x[k];
        }
        x[i] = s / l[i][i];
    }
    x
}

/// Solve the dense symmetric system `(A - σ·I) x = b` by Gaussian elimination
/// with partial pivoting.
///
/// When the shift σ sits at or between clustered eigenvalues, `A - σ·I` is
/// near-singular; without pivoting the elimination can divide by a pivot that
/// has collapsed to (or through) zero and produce inf/NaN, poisoning the whole
/// Lanczos recurrence. Pivoting keeps the solve finite, and a pivot that is
/// still tiny after row exchange is regularised to `±1e-12` (mirroring
/// `tpt-fem-solve`'s `dense_solve`) so the shift-invert operator stays usable:
/// the huge but finite inverse maps the iteration onto the eigenpairs nearest
/// σ, which is exactly the desired shift-invert behaviour.
fn solve_shifted(a: &[Vec<f64>], sigma: f64, b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            mat[i][j] = a[i][j] - if i == j { sigma } else { 0.0 };
        }
    }
    let mut x = b.to_vec();
    for col in 0..n {
        // Partial pivoting: pick the largest remaining entry in this column.
        let piv_row = (col..n)
            .reduce(|best, r| {
                if mat[r][col].abs() > mat[best][col].abs() {
                    r
                } else {
                    best
                }
            })
            .unwrap_or(col);
        if piv_row != col {
            mat.swap(piv_row, col);
            x.swap(piv_row, col);
        }
        // Regularise a collapsed pivot instead of dividing by ~0.
        let piv = if mat[col][col].abs() < 1e-12 {
            if mat[col][col] < 0.0 {
                -1e-12
            } else {
                1e-12
            }
        } else {
            mat[col][col]
        };
        for r in (col + 1)..n {
            let f = mat[r][col] / piv;
            for c in col..n {
                mat[r][c] -= f * mat[col][c];
            }
            x[r] -= f * x[col];
        }
    }
    for col in (0..n).rev() {
        let piv = if mat[col][col].abs() < 1e-12 {
            if mat[col][col] < 0.0 {
                -1e-12
            } else {
                1e-12
            }
        } else {
            mat[col][col]
        };
        let mut s = x[col];
        for c in (col + 1)..n {
            s -= mat[col][c] * x[c];
        }
        x[col] = s / piv;
    }
    x
}

/// Standard symmetric Lanczos with the shift-invert operator `(A - σ·I)⁻¹`,
/// returning the orthonormal basis, diagonal `α` and off-diagonal `β`.
fn shifted_lanczos(a: &[Vec<f64>], sigma: f64, m: usize) -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    let n = a.len();
    let mut v0 = vec![0.0; n];
    v0[0] = 1.0;
    let nv = dot(&v0, &v0).sqrt();
    for v in v0.iter_mut() {
        *v /= nv;
    }
    let mut basis = vec![v0];
    let mut alpha = Vec::new();
    let mut beta = Vec::new();
    for j in 0..m {
        let w = solve_shifted(a, sigma, &basis[j]);
        let an = dot(&w, &basis[j]);
        let mut w = w;
        // Two full Gram-Schmidt passes ("twice is enough"): one pass loses
        // orthogonality whenever the shifted-inverse nearly cancels the
        // current basis vector (which happens for clustered spectra, where
        // the new direction can be orders of magnitude smaller than the
        // rounding noise in the cancelled components).
        for _pass in 0..2 {
            for k in 0..=j {
                let p = dot(&w, &basis[k]);
                for d in 0..n {
                    w[d] -= p * basis[k][d];
                }
            }
        }
        alpha.push(an);
        if j + 1 == m {
            break;
        }
        let b = dot(&w, &w).sqrt();
        if b < 1e-12 {
            break;
        }
        beta.push(b);
        basis.push(w.iter().map(|x| x / b).collect());
    }
    (basis, alpha, beta)
}

/// Solve the symmetric generalized eigenproblem `K x = λ M x` by shift-invert
/// Lanczos.
///
/// The pencil is first reduced to the standard symmetric problem `K' y = λ y`
/// with `K' = L⁻¹ K L⁻ᵀ` (where `M = L Lᵀ` is the Cholesky factor of the mass
/// matrix), so the existing symmetric Lanczos machinery applies directly. The
/// shift `σ` is applied through `(K' - σ)⁻¹`, which targets the eigenpairs
/// nearest `σ` (a shift near `0` therefore extracts the low structural
/// frequencies). Returned eigenpairs `(λ, x)` are sorted ascending by `λ` and
/// the eigenvectors `x` live in the input `K`/`M` DOF space.
pub fn generalized_lanczos_eigs(
    k: &Coo,
    m: &Coo,
    sigma: f64,
    num: usize,
    lanczos_dim: usize,
) -> Result<Vec<(f64, Vec<f64>)>, SparseError> {
    let n = k.to_csr().nrows;
    let kd = coo_to_dense(k, n);
    let md = coo_to_dense(m, n);
    let l = cholesky(&md).ok_or_else(|| {
        SparseError::Numeric(
            "generalized eigenproblem: mass matrix is not positive-definite".into(),
        )
    })?;

    // K' = L⁻¹ K L⁻ᵀ:  X = L⁻¹ K (column-wise), then K' = X L⁻ᵀ.
    let mut xmat = vec![vec![0.0; n]; n];
    for c in 0..n {
        let b: Vec<f64> = (0..n).map(|r| kd[r][c]).collect();
        let sol = solve_lower(&l, &b);
        for r in 0..n {
            xmat[r][c] = sol[r];
        }
    }
    let mut linv_t = vec![vec![0.0; n]; n];
    for c in 0..n {
        let b: Vec<f64> = (0..n).map(|r| if r == c { 1.0 } else { 0.0 }).collect();
        let sol = solve_upper(&l, &b);
        for r in 0..n {
            linv_t[r][c] = sol[r];
        }
    }
    let kprime = matmul(&xmat, &linv_t);

    let (basis, alpha, beta) = shifted_lanczos(&kprime, sigma, lanczos_dim);
    let mdim = basis.len();
    let mut t = vec![vec![0.0; mdim]; mdim];
    for i in 0..mdim {
        t[i][i] = alpha[i];
    }
    for i in 0..(mdim - 1) {
        t[i][i + 1] = beta[i];
        t[i + 1][i] = beta[i];
    }
    let (eig, qmat) = jacobi(t, 200);
    // Ritz values μ of (K'-σ)⁻¹ satisfy λ = σ + 1/μ.
    let mut pairs: Vec<(f64, usize)> = (0..mdim)
        .map(|i| {
            let mu = eig[i];
            let lam = sigma + 1.0 / mu;
            (lam, i)
        })
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let take = pairs.iter().take(num).cloned().collect::<Vec<_>>();
    Ok(take
        .into_iter()
        .map(|(lam, i)| {
            // Eigenvector in the transformed space, then back to x = L⁻ᵀ y.
            let mut y = vec![0.0; n];
            for (r, vr) in basis.iter().enumerate() {
                let c = qmat[r][i];
                for (d, v) in vr.iter().enumerate() {
                    y[d] += c * v;
                }
            }
            let x = solve_upper(&l, &y);
            (lam, x)
        })
        .collect())
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

    #[test]
    fn generalized_with_identity_mass() {
        // K = [[2,1],[1,2]], M = I  =>  λ = 1, 3.
        let mut k = Coo::new();
        k.push(0, 0, 2.0);
        k.push(0, 1, 1.0);
        k.push(1, 0, 1.0);
        k.push(1, 1, 2.0);
        let mut m = Coo::new();
        m.push(0, 0, 1.0);
        m.push(1, 1, 1.0);
        let eigs = generalized_lanczos_eigs(&k, &m, 0.0, 2, 8).unwrap();
        let mut vals: Vec<f64> = eigs.iter().map(|(l, _)| *l).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((vals[0] - 1.0).abs() < 1e-6, "got {}", vals[0]);
        assert!((vals[1] - 3.0).abs() < 1e-6, "got {}", vals[1]);
        for (_, v) in &eigs {
            let mv = matvec(&m, v);
            assert!((dot(v, &mv) - 1.0).abs() < 1e-4);
        }
    }

    #[test]
    fn generalized_with_diagonal_mass() {
        // K = [[2,1],[1,2]], M = diag(2,2)  =>  λ = 0.5, 1.5.
        let mut k = Coo::new();
        k.push(0, 0, 2.0);
        k.push(0, 1, 1.0);
        k.push(1, 0, 1.0);
        k.push(1, 1, 2.0);
        let mut m = Coo::new();
        m.push(0, 0, 2.0);
        m.push(1, 1, 2.0);
        let eigs = generalized_lanczos_eigs(&k, &m, 0.0, 2, 8).unwrap();
        let mut vals: Vec<f64> = eigs.iter().map(|(l, _)| *l).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((vals[0] - 0.5).abs() < 1e-5, "got {}", vals[0]);
        assert!((vals[1] - 1.5).abs() < 1e-5, "got {}", vals[1]);
        assert!((vals[0] - 0.5).abs() < 1e-5, "got {}", vals[0]);
        assert!((vals[1] - 1.5).abs() < 1e-5, "got {}", vals[1]);
    }

    /// Stiffness `I + delta·L` where `L` is the n-point free-free discrete
    /// 1-D Laplacian. Its eigenvalues are `1 + delta·(2 - 2cos(kπ/(n+1)))`,
    /// i.e. a tight cluster of width `~delta` around 1.
    fn clustered_stiffness(n: usize, delta: f64) -> Coo {
        let mut c = Coo::new();
        for i in 0..n {
            c.push(i, i, 1.0 + delta * 2.0);
            if i + 1 < n {
                c.push(i, i + 1, -delta);
                c.push(i + 1, i, -delta);
            }
        }
        c
    }

    fn identity_mass(n: usize) -> Coo {
        let mut c = Coo::new();
        for i in 0..n {
            c.push(i, i, 1.0);
        }
        c
    }

    #[test]
    fn generalized_closely_clustered_eigenvalues() {
        // Stress case flagged in todo.md (11c): six eigenvalues packed into a
        // band of width ~4e-6. The shift-invert operator maps them to very
        // large, very close Ritz targets, stressing both the (K' - σ)⁻¹ solve
        // and the projected Jacobi eigensolve. A full-dimension Lanczos basis
        // must still resolve every member of the cluster to high accuracy.
        let n = 6;
        let delta = 1e-6;
        let k = clustered_stiffness(n, delta);
        let m = identity_mass(n);
        let eigs = generalized_lanczos_eigs(&k, &m, 0.0, n, n).unwrap();
        assert_eq!(eigs.len(), n, "expected the full spectrum back");
        let mut vals: Vec<f64> = eigs.iter().map(|(l, _)| *l).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, lam) in vals.iter().enumerate() {
            let expected = 1.0
                + delta
                    * (2.0
                        - 2.0 * (std::f64::consts::PI * (i + 1) as f64 / (n as f64 + 1.0)).cos());
            assert!(
                (lam - expected).abs() < 1e-9,
                "cluster member {i}: got {lam}, want {expected}"
            );
        }
        // Eigenpair residuals ||K x - lambda M x|| and M-orthonormality.
        for (lam, x) in &eigs {
            let mx = matvec(&m, x);
            let kx = matvec(&k, x);
            let res: f64 = kx
                .iter()
                .zip(&mx)
                .map(|(a, b)| (a - lam * b) * (a - lam * b))
                .sum::<f64>()
                .sqrt();
            assert!(res < 1e-6, "eigenpair residual {res} at λ={lam}");
            assert!(
                (dot(x, &mx) - 1.0).abs() < 1e-6,
                "not M-normalised at λ={lam}"
            );
        }
    }

    #[test]
    fn generalized_shift_inside_cluster_is_accurate() {
        // Same clustered pencil, but with the shift σ placed *inside* the
        // cluster: the shifted matrix K - σI is then near-singular (its
        // smallest eigenvalue magnitude is ~1e-6), which is the hardest input
        // for a shift-invert solve short of an exactly singular pencil. The
        // recovered spectrum must still match the closed form, to a looser
        // (but still strict) tolerance reflecting the ~1e6 condition number.
        let n = 6;
        let delta = 1e-6;
        let k = clustered_stiffness(n, delta);
        let m = identity_mass(n);
        // Mid-cluster shift: halfway between members k=3 and k=4.
        let s3 = 2.0 - 2.0 * (std::f64::consts::PI * 3.0 / (n as f64 + 1.0)).cos();
        let s4 = 2.0 - 2.0 * (std::f64::consts::PI * 4.0 / (n as f64 + 1.0)).cos();
        let sigma = 1.0 + delta * 0.5 * (s3 + s4);
        let eigs = generalized_lanczos_eigs(&k, &m, sigma, 3, n).unwrap();
        assert_eq!(eigs.len(), 3);
        let mut vals: Vec<f64> = eigs.iter().map(|(l, _)| *l).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, lam) in vals.iter().enumerate() {
            let expected = 1.0
                + delta
                    * (2.0
                        - 2.0 * (std::f64::consts::PI * (i + 1) as f64 / (n as f64 + 1.0)).cos());
            assert!(
                (lam - expected).abs() < 1e-4,
                "cluster member {i}: got {lam}, want {expected}"
            );
        }
    }

    #[test]
    fn generalized_lanczos_clustered_eigenvalues() {
        // Two eigenvalues tightly clustered near 1.0 (1 ± 1e-3) in a *coupled*
        // block, plus two well-separated ones. The coupling keeps the Krylov
        // start vector `[1,0,0,0]` off the eigenspace so the subspace does not
        // collapse to one dimension. Shift-invert Lanczos must recover BOTH
        // clustered values (not merge them into a single mode) and return
        // accurate Ritz vectors. Regression for todo.md:716.
        let b = 1e-3;
        let mut k = Coo::new();
        // Block [[1, b],[b, 1]] → eigenvalues 1±b, coupled so v0 is not an eigenvector.
        k.push(0, 0, 1.0);
        k.push(0, 1, b);
        k.push(1, 0, b);
        k.push(1, 1, 1.0);
        k.push(2, 2, 3.0);
        k.push(3, 3, 7.0);
        let mut m = Coo::new();
        for i in 0..4 {
            m.push(i, i, 1.0);
        }
        let eigs = generalized_lanczos_eigs(&k, &m, 0.0, 2, 8).unwrap();
        assert_eq!(eigs.len(), 2);
        let mut lam: Vec<f64> = eigs.iter().map(|(l, _)| *l).collect();
        lam.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((lam[0] - (1.0 - b)).abs() < 1e-4, "got {}", lam[0]);
        assert!((lam[1] - (1.0 + b)).abs() < 1e-4, "got {}", lam[1]);
        // Each Ritz vector must satisfy K x ≈ λ M x.
        let kd = coo_to_dense(&k, 4);
        let md = coo_to_dense(&m, 4);
        for (l, x) in &eigs {
            let kx: Vec<f64> = (0..4)
                .map(|r| (0..4).map(|c| kd[r][c] * x[c]).sum())
                .collect();
            let mx: Vec<f64> = (0..4)
                .map(|r| (0..4).map(|c| md[r][c] * x[c]).sum())
                .collect();
            let res = (0..4)
                .map(|r| kx[r] - l * mx[r])
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt();
            assert!(res < 1e-6, "eigenvector residual {res}");
        }
    }
}
