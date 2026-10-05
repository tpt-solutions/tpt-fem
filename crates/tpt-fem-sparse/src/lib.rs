//! FEM-specific sparse-matrix assembly adapter with a `tpt-math-linalg-dense`-backed solve.
#![allow(clippy::needless_range_loop)]
//!
//! Finite-element assembly naturally produces a sparse global matrix as a bag of
//! `(row, col, value)` triplets, where the same entry may be written many times
//! (once per element that touches it). [`Coo`] is a growable coordinate-list
//! accumulator that supports duplicate-summing assembly; [`Coo::to_csr`]
//! collapses it into a canonical compressed-sparse-row [`Csr`] matrix.
//!
//! [`solve`] assembles the matrix into a dense [`tpt_math_linalg_dense::DMatrix`]
//! and solves `A x = b` with the in-house partial-pivot LU decomposition, so the
//! workspace carries no Apache-2.0-only linear-algebra dependency.
//!
//! # Example
//!
//! ```
//! use tpt_fem_sparse::Coo;
//!
//! // Assemble [[2, 1], [1, 3]] by writing each entry twice then summing.
//! let mut c = Coo::new();
//! c.push(0, 0, 1.0);
//! c.push(0, 0, 1.0);
//! c.push(0, 1, 1.0);
//! c.push(1, 0, 1.0);
//! c.push(1, 1, 1.5);
//! c.push(1, 1, 1.5);
//! let csr = c.to_csr();
//! assert_eq!(csr.nnz(), 4);
//! assert_eq!(csr.row_ptrs, vec![0, 2, 4]);
//! assert_eq!(csr.values, vec![2.0, 1.0, 1.0, 3.0]);
//! ```

use tpt_math_linalg_dense::{DMatrix, DVector};

/// A coordinate-list (triplet) accumulator for sparse matrix assembly.
///
/// Entries written to the same `(row, col)` are summed when the list is
/// collapsed with [`Coo::to_csr`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Coo {
    /// Row indices.
    pub rows: Vec<usize>,
    /// Column indices.
    pub cols: Vec<usize>,
    /// Values.
    pub vals: Vec<f64>,
}

/// A compressed-sparse-row matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Csr {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row pointers (length `nrows + 1`); row `r` occupies
    /// `col_ind[row_ptrs[r]..row_ptrs[r+1]]`.
    pub row_ptrs: Vec<usize>,
    /// Column indices, grouped by row.
    pub col_ind: Vec<usize>,
    /// Non-zero values, parallel to `col_ind`.
    pub values: Vec<f64>,
}

/// Errors produced while building or solving a sparse system.
#[derive(Debug)]
pub enum SparseError {
    /// The matrix could not be constructed from triplets.
    Creation(String),
    /// The symbolic factorization failed.
    Symbolic(String),
    /// The numeric factorization or solve failed.
    Numeric(String),
}

impl std::fmt::Display for SparseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SparseError::Creation(m) => write!(f, "failed to build sparse matrix: {m}"),
            SparseError::Symbolic(m) => write!(f, "sparse symbolic factorization failed: {m}"),
            SparseError::Numeric(m) => write!(f, "sparse numeric factorization/solve failed: {m}"),
        }
    }
}

impl std::error::Error for SparseError {}

impl Coo {
    /// Create an empty accumulator.
    pub fn new() -> Self {
        Coo::default()
    }

    /// Create an empty accumulator with space reserved for `capacity` entries.
    pub fn with_capacity(capacity: usize) -> Self {
        Coo {
            rows: Vec::with_capacity(capacity),
            cols: Vec::with_capacity(capacity),
            vals: Vec::with_capacity(capacity),
        }
    }

    /// Number of stored entries (before duplicate summing).
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// True if no entries have been stored.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Append a `(row, col, value)` entry.
    pub fn push(&mut self, row: usize, col: usize, value: f64) {
        self.rows.push(row);
        self.cols.push(col);
        self.vals.push(value);
    }

    /// Collapse into a canonical CSR matrix, summing duplicate `(row, col)`
    /// entries. Rows are stored in ascending column order.
    pub fn to_csr(&self) -> Csr {
        let n = self.rows.len();
        let nrows = self
            .rows
            .iter()
            .cloned()
            .chain(std::iter::once(0))
            .max()
            .map(|m| m + 1)
            .unwrap_or(0);
        let ncols = self
            .cols
            .iter()
            .cloned()
            .chain(std::iter::once(0))
            .max()
            .map(|m| m + 1)
            .unwrap_or(0);

        let mut entries: Vec<(usize, usize, f64)> = (0..n)
            .map(|i| (self.rows[i], self.cols[i], self.vals[i]))
            .collect();
        entries.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut col_ind = Vec::new();
        let mut values = Vec::new();
        let mut row_ptrs = vec![0usize; nrows + 1];

        let mut k = 0;
        while k < entries.len() {
            let (r, c, v) = entries[k];
            let mut sum = v;
            let mut kk = k + 1;
            while kk < entries.len() && entries[kk].0 == r && entries[kk].1 == c {
                sum += entries[kk].2;
                kk += 1;
            }
            col_ind.push(c);
            values.push(sum);
            row_ptrs[r + 1] += 1;
            k = kk;
        }
        for r in 0..nrows {
            row_ptrs[r + 1] += row_ptrs[r];
        }

        Csr {
            nrows,
            ncols,
            row_ptrs,
            col_ind,
            values,
        }
    }
}

impl Csr {
    /// Number of stored non-zeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Matrix–vector product `y = A x` using the compressed rows.
    ///
    /// The row-contiguous storage lets the scalar accumulation loop be
    /// auto-vectorised, and — unlike [`Coo::to_csr`] — the conversion cost is
    /// paid once by the caller, not on every call. Prefer this in time-stepping
    /// or iterative loops: convert with [`Coo::to_csr`] once, then call
    /// [`Csr::matvec`] per iteration.
    ///
    /// `y` has `max(nrows, x.len())` entries (so a placeholder/empty matrix
    /// still yields a correctly-sized zero vector); `x` must have at least
    /// `ncols` entries (extra entries are ignored, a short `x` panics on index).
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.nrows.max(x.len())];
        for r in 0..self.nrows {
            let mut s = 0.0;
            for c in self.row_ptrs[r]..self.row_ptrs[r + 1] {
                s += self.values[c] * x[self.col_ind[c]];
            }
            y[r] = s;
        }
        y
    }
}

/// Solve the square linear system `A x = b`, where `A` is supplied as a
/// [`Coo] accumulator, returning the solution vector `x`.
///
/// Duplicate `(row, col)` entries in `coo` are summed (via [`Coo::to_csr`]),
/// the matrix is assembled into a dense
/// [`tpt_math_linalg_dense::DMatrix`] and factored with the in-house
/// partial-pivot LU decomposition, and the system is solved.
///
/// # Cost note
///
/// Systems with at least 200 unknowns that are symmetric positive-definite are
/// automatically routed to the sparse envelope Cholesky ([`solve_skyline`]);
/// large unsymmetric or indefinite systems whose reordered band is narrow use
/// the banded LU ([`solve_banded`]); everything else uses the dense path described
/// below. The default dense backend assembles a **dense** `n×n` matrix regardless of `A`'s
/// sparsity and factors it in `O(n³)` time with `O(n²)` storage. This is fine
/// for the small, hand-built meshes used to validate the crate, but it does
/// **not** scale to large sparse problems. For those, enable the optional
/// `russell` feature (which dispatches to the `russell_sparse`
/// UMFPACK/MUMPS direct solvers) — see `solve_russell`. As a rule of thumb,
/// prefer `russell` once `A` has more than a few thousand rows, or [`solve_cg`]
/// (pure Rust, `O(nnz)`) when `A` is symmetric positive-definite.
pub fn solve(coo: &Coo, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
    let sols = solve_multi(coo, std::slice::from_ref(&rhs.to_vec()))?;
    sols.into_iter()
        .next()
        .ok_or_else(|| SparseError::Numeric("solve received an empty right-hand side".into()))
}

/// System size from which [`solve`]/[`solve_multi`] try the sparse envelope
/// Cholesky ([`solve_skyline`]) before the dense LU.
const SKYLINE_THRESHOLD: usize = 200;

/// Options for the iterative [`solve_cg`] solver.
#[derive(Clone, Copy, Debug)]
pub struct CgOptions {
    /// Relative residual tolerance: stop when `‖r‖ ≤ tol · ‖b‖`.
    pub tol: f64,
    /// Maximum number of iterations (`0` means `10 · n`).
    pub max_iter: usize,
}

impl Default for CgOptions {
    fn default() -> Self {
        CgOptions {
            tol: 1e-10,
            max_iter: 0,
        }
    }
}

/// Result of a successful [`solve_cg`].
#[derive(Clone, Debug)]
pub struct CgSolution {
    /// The solution vector.
    pub x: Vec<f64>,
    /// Iterations taken.
    pub iterations: usize,
    /// Final relative residual `‖r‖ / ‖b‖` (absolute `‖r‖` if `b = 0`).
    pub relative_residual: f64,
}

/// Solve the **symmetric positive-definite** system `A x = b` with
/// Jacobi-preconditioned conjugate gradients.
///
/// Unlike [`solve`] this works on the CSR form directly — `O(nnz)` memory and
/// `O(nnz)` work per iteration — so it scales to systems far larger than the
/// dense-LU backend allows, without any external toolchain. `A` must be SPD
/// (e.g. a Dirichlet-reduced stiffness matrix); for indefinite or
/// unsymmetric systems use [`solve`] or the `russell` backend. Returns
/// [`SparseError::Numeric`] if the matrix is not square, a diagonal entry is
/// non-positive, the iteration breaks down (non-SPD input), or it does not
/// reach `opts.tol` within `opts.max_iter` iterations.
pub fn solve_cg(coo: &Coo, rhs: &[f64], opts: &CgOptions) -> Result<CgSolution, SparseError> {
    let a = coo.to_csr();
    let n = a.nrows;
    if a.ncols != n {
        return Err(SparseError::Numeric(format!(
            "solve_cg requires a square matrix, got {n} x {}",
            a.ncols
        )));
    }
    if rhs.len() != n {
        return Err(SparseError::Numeric(format!(
            "rhs length {} does not match matrix dimension {n}",
            rhs.len()
        )));
    }
    if !(opts.tol.is_finite() && opts.tol > 0.0) {
        return Err(SparseError::Numeric(format!(
            "solve_cg tolerance must be finite and positive (got {})",
            opts.tol
        )));
    }

    // Jacobi preconditioner M⁻¹ = diag(A)⁻¹.
    let mut inv_diag = vec![0.0; n];
    for r in 0..n {
        let mut d = 0.0;
        for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            if a.col_ind[c] == r {
                d += a.values[c];
            }
        }
        if !(d.is_finite() && d > 0.0) {
            return Err(SparseError::Numeric(format!(
                "solve_cg requires a positive diagonal (row {r} has {d})"
            )));
        }
        inv_diag[r] = 1.0 / d;
    }

    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(a, b)| a * b).sum::<f64>();
    let bnorm = dot(rhs, rhs).sqrt();
    let scale = if bnorm > 0.0 { bnorm } else { 1.0 };
    let max_iter = if opts.max_iter == 0 {
        10 * n.max(1)
    } else {
        opts.max_iter
    };

    let mut x = vec![0.0; n];
    let mut r = rhs.to_vec();
    let mut res = dot(&r, &r).sqrt() / scale;
    if res <= opts.tol {
        return Ok(CgSolution {
            x,
            iterations: 0,
            relative_residual: res,
        });
    }
    let mut z: Vec<f64> = r.iter().zip(&inv_diag).map(|(r, m)| r * m).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    for it in 1..=max_iter {
        let ap = a.matvec(&p);
        let pap = dot(&p, &ap);
        if !(pap.is_finite() && pap > 0.0) {
            return Err(SparseError::Numeric(
                "solve_cg broke down: matrix is not positive-definite".into(),
            ));
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        res = dot(&r, &r).sqrt() / scale;
        if res <= opts.tol {
            return Ok(CgSolution {
                x,
                iterations: it,
                relative_residual: res,
            });
        }
        for i in 0..n {
            z[i] = r[i] * inv_diag[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    Err(SparseError::Numeric(format!(
        "solve_cg did not converge in {max_iter} iterations (relative residual {res:e})"
    )))
}

/// Reverse Cuthill-McKee ordering of the symmetrised sparsity pattern of `a`.
///
/// Returns `perm` with `perm[new] = old`. Each connected component is started
/// from a minimum-degree node and visited breadth-first (neighbours by
/// ascending degree); the final order is reversed. This clusters the non-zeros
/// near the diagonal, which is what bounds the fill of an envelope Cholesky.
fn reverse_cuthill_mckee(a: &Csr) -> Vec<usize> {
    let n = a.nrows;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for r in 0..n {
        for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            let col = a.col_ind[c];
            if col != r && col < n {
                adj[r].push(col);
                adj[col].push(r);
            }
        }
    }
    for list in &mut adj {
        list.sort_unstable();
        list.dedup();
    }
    let degree: Vec<usize> = adj.iter().map(Vec::len).collect();
    let mut by_degree: Vec<usize> = (0..n).collect();
    by_degree.sort_by_key(|&i| degree[i]);

    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for &start in &by_degree {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let head = order.len();
        order.push(start);
        let mut qi = head;
        while qi < order.len() {
            let v = order[qi];
            qi += 1;
            let mut next: Vec<usize> = adj[v].iter().copied().filter(|&w| !visited[w]).collect();
            next.sort_by_key(|&w| degree[w]);
            for w in next {
                visited[w] = true;
                order.push(w);
            }
        }
    }
    order.reverse();
    order
}

/// Solve the **symmetric positive-definite** system `A x = b` with a sparse
/// direct method: reverse Cuthill-McKee reordering followed by an envelope
/// (skyline) Cholesky factorisation `P A P^T = L L^T`.
///
/// Memory is `O(n * bandwidth)` and time `O(n * bandwidth^2)` after
/// reordering, which for mesh-based FEM matrices is far below the dense
/// `O(n^2)` / `O(n^3)` of [`solve`] -- e.g. a 2-D Laplacian with 10^4 unknowns
/// factors in milliseconds instead of needing ~800 MB. Pure Rust, no external
/// toolchain. Only the lower triangle of `A` is read (duplicates summed); `A`
/// must be square and SPD, otherwise [`SparseError::Numeric`] reports a
/// non-positive pivot. Unlike [`solve_cg`] the result is exact to round-off and
/// the cost does not depend on conditioning.
pub fn solve_skyline(coo: &Coo, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
    if coo.is_empty() && rhs.is_empty() {
        return Ok(Vec::new());
    }
    let a = coo.to_csr();
    if a.ncols != a.nrows {
        return Err(SparseError::Numeric(format!(
            "solve_skyline requires a square matrix, got {} x {}",
            a.nrows, a.ncols
        )));
    }
    if rhs.len() != a.nrows {
        return Err(SparseError::Numeric(format!(
            "rhs length {} does not match matrix dimension {}",
            rhs.len(),
            a.nrows
        )));
    }
    SkylineFactor::new(&a)?.solve(rhs)
}

/// Envelope Cholesky factor `P A P^T = L L^T` of an SPD matrix.
struct SkylineFactor {
    perm: Vec<usize>,
    first: Vec<usize>,
    start: Vec<usize>,
    l: Vec<f64>,
}

impl SkylineFactor {
    /// Factor the square CSR matrix `a` (only its lower triangle is read).
    fn new(a: &Csr) -> Result<Self, SparseError> {
        let n = a.nrows;
        let perm = reverse_cuthill_mckee(a);
        let mut inv = vec![0usize; n];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }

        // Envelope of the permuted lower triangle: first[i] = leftmost column of row i.
        let mut first: Vec<usize> = (0..n).collect();
        for r in 0..n {
            for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
                let (i, j) = (inv[r], inv[a.col_ind[c]]);
                let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
                first[hi] = first[hi].min(lo);
            }
        }
        let mut start = vec![0usize; n + 1];
        for i in 0..n {
            start[i + 1] = start[i] + (i - first[i] + 1);
        }
        let mut l = vec![0.0f64; start[n]];
        let at = |i: usize, j: usize| start[i] + (j - first[i]);
        for r in 0..n {
            for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
                // Read each stored entry of the original lower triangle once and
                // mirror it into the permuted lower triangle.
                if a.col_ind[c] <= r {
                    let (i, j) = (inv[r], inv[a.col_ind[c]]);
                    let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
                    l[at(hi, lo)] += a.values[c];
                }
            }
        }

        // In-place envelope Cholesky.
        for i in 0..n {
            for j in first[i]..=i {
                let k0 = first[i].max(first[j]);
                let mut sum = l[at(i, j)];
                for k in k0..j {
                    sum -= l[at(i, k)] * l[at(j, k)];
                }
                if j < i {
                    l[at(i, j)] = sum / l[at(j, j)];
                } else {
                    if !(sum.is_finite() && sum > 0.0) {
                        return Err(SparseError::Numeric(format!(
                            "solve_skyline: non-positive pivot at row {i}; matrix is not positive-definite"
                        )));
                    }
                    l[at(i, i)] = sum.sqrt();
                }
            }
        }
        Ok(SkylineFactor {
            perm,
            first,
            start,
            l,
        })
    }

    /// Solve `A x = rhs` with the stored factor.
    fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
        let n = self.perm.len();
        let at = |i: usize, j: usize| self.start[i] + (j - self.first[i]);
        // L y = P b, then L^T z = y, then x = P^T z.
        let mut y: Vec<f64> = self.perm.iter().map(|&old| rhs[old]).collect();
        for i in 0..n {
            let mut sum = y[i];
            for k in self.first[i]..i {
                sum -= self.l[at(i, k)] * y[k];
            }
            y[i] = sum / self.l[at(i, i)];
        }
        for i in (0..n).rev() {
            y[i] /= self.l[at(i, i)];
            let yi = y[i];
            for k in self.first[i]..i {
                y[k] -= self.l[at(i, k)] * yi;
            }
        }
        let mut x = vec![0.0; n];
        for (new, &old) in self.perm.iter().enumerate() {
            x[old] = y[new];
        }
        Ok(x)
    }
}

impl Csr {
    /// `true` if the matrix equals its transpose to within `rel_tol` of the
    /// largest entry magnitude (`O(nnz log nnz)`).
    pub fn is_symmetric(&self, rel_tol: f64) -> bool {
        if self.nrows != self.ncols {
            return false;
        }
        let scale = self.values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let tol = rel_tol * scale;
        for r in 0..self.nrows {
            for c in self.row_ptrs[r]..self.row_ptrs[r + 1] {
                let col = self.col_ind[c];
                let lo = self.row_ptrs[col];
                let hi = self.row_ptrs[col + 1];
                // Columns within a row are sorted (see `Coo::to_csr`).
                let mirror = match self.col_ind[lo..hi].binary_search(&r) {
                    Ok(k) => self.values[lo + k],
                    Err(_) => 0.0,
                };
                if (self.values[c] - mirror).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
}

/// Banded LU factorisation with partial pivoting (LAPACK `gbtf2` layout) of a
/// reverse-Cuthill-McKee-permuted matrix. Handles unsymmetric and indefinite
/// systems.
struct BandedLu {
    perm: Vec<usize>,
    n: usize,
    kl: usize,
    kv: usize,
    ldab: usize,
    ab: Vec<f64>,
    ipiv: Vec<usize>,
}

impl BandedLu {
    /// Half-bandwidths `(kl, ku)` of `a` after the RCM permutation, plus the
    /// permutation itself.
    fn bandwidths(a: &Csr) -> (Vec<usize>, Vec<usize>, usize, usize) {
        let perm = reverse_cuthill_mckee(a);
        let mut inv = vec![0usize; a.nrows];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }
        let (mut kl, mut ku) = (0usize, 0usize);
        for r in 0..a.nrows {
            for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
                let (i, j) = (inv[r], inv[a.col_ind[c]]);
                if i > j {
                    kl = kl.max(i - j);
                } else {
                    ku = ku.max(j - i);
                }
            }
        }
        (perm, inv, kl, ku)
    }

    fn new(
        a: &Csr,
        perm: Vec<usize>,
        inv: &[usize],
        kl: usize,
        ku: usize,
    ) -> Result<Self, SparseError> {
        let n = a.nrows;
        let kv = kl + ku;
        let ldab = 2 * kl + ku + 1;
        let mut ab = vec![0.0f64; ldab * n];
        for r in 0..n {
            for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
                let (i, j) = (inv[r], inv[a.col_ind[c]]);
                ab[(kv + i - j) + j * ldab] += a.values[c];
            }
        }
        let mut ipiv = vec![0usize; n];
        let mut ju = 0usize;
        for j in 0..n {
            let km = kl.min(n - 1 - j);
            // Pivot search in column j, rows j..=j+km.
            let mut jp = 0usize;
            let mut best = ab[kv + j * ldab].abs();
            for i in 1..=km {
                let v = ab[kv + i + j * ldab].abs();
                if v > best {
                    best = v;
                    jp = i;
                }
            }
            ipiv[j] = jp;
            if !(best.is_finite() && best > 0.0) {
                return Err(SparseError::Numeric(format!(
                    "banded LU: singular matrix (zero pivot in column {j})"
                )));
            }
            ju = ju.max((j + ku + jp).min(n - 1));
            if jp != 0 {
                for c in 0..=(ju - j) {
                    let a_idx = (kv + jp - c) + (j + c) * ldab;
                    let b_idx = (kv - c) + (j + c) * ldab;
                    ab.swap(a_idx, b_idx);
                }
            }
            if km > 0 {
                let piv = ab[kv + j * ldab];
                for i in 1..=km {
                    ab[kv + i + j * ldab] /= piv;
                }
                for c in 1..=(ju - j) {
                    let u = ab[(kv - c) + (j + c) * ldab];
                    if u != 0.0 {
                        for i in 1..=km {
                            let l = ab[kv + i + j * ldab];
                            ab[(kv + i - c) + (j + c) * ldab] -= l * u;
                        }
                    }
                }
            }
        }
        Ok(BandedLu {
            perm,
            n,
            kl,
            kv,
            ldab,
            ab,
            ipiv,
        })
    }

    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let (n, kl, kv, ldab) = (self.n, self.kl, self.kv, self.ldab);
        let mut b: Vec<f64> = self.perm.iter().map(|&old| rhs[old]).collect();
        for j in 0..n.saturating_sub(1) {
            let lm = kl.min(n - 1 - j);
            let l = j + self.ipiv[j];
            if l != j {
                b.swap(l, j);
            }
            let bj = b[j];
            for i in 1..=lm {
                b[j + i] -= self.ab[kv + i + j * ldab] * bj;
            }
        }
        for j in (0..n).rev() {
            b[j] /= self.ab[kv + j * ldab];
            let bj = b[j];
            let i1 = j.saturating_sub(kv);
            for i in i1..j {
                b[i] -= self.ab[(kv + i - j) + j * ldab] * bj;
            }
        }
        let mut x = vec![0.0; n];
        for (new, &old) in self.perm.iter().enumerate() {
            x[old] = b[new];
        }
        x
    }
}

/// Solve a general (possibly unsymmetric or indefinite) sparse system `A x = b`
/// with a banded LU after reverse Cuthill-McKee reordering and partial
/// pivoting.
///
/// Memory is `O(n * (2 kl + ku + 1))` for the half-bandwidths `kl`, `ku` of the
/// reordered matrix, so it scales like [`solve_skyline`] on mesh-based
/// matrices but needs no symmetry or definiteness: it is the pure-Rust sparse
/// direct path for e.g. convection-dominated or saddle-point systems. Returns
/// [`SparseError::Numeric`] for a non-square matrix, mismatched `rhs`, or a
/// singular matrix.
pub fn solve_banded(coo: &Coo, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
    if coo.is_empty() && rhs.is_empty() {
        return Ok(Vec::new());
    }
    let a = coo.to_csr();
    if a.ncols != a.nrows {
        return Err(SparseError::Numeric(format!(
            "solve_banded requires a square matrix, got {} x {}",
            a.nrows, a.ncols
        )));
    }
    if rhs.len() != a.nrows {
        return Err(SparseError::Numeric(format!(
            "rhs length {} does not match matrix dimension {}",
            rhs.len(),
            a.nrows
        )));
    }
    let (perm, inv, kl, ku) = BandedLu::bandwidths(&a);
    Ok(BandedLu::new(&a, perm, &inv, kl, ku)?.solve(rhs))
}

/// Solve `A x_k = rhs[k]` for every right-hand side in `rhs` against the
/// *same* matrix `A`.
///
/// Equivalent to calling [`solve`] once per right-hand side. Callers with
/// multiple RHS vectors against an unchanged `A` (e.g. an arc-length
/// continuation corrector, which needs both a tangent and a
/// residual-correction solve per iteration) should prefer this over repeated
/// `solve` calls.
///
/// Like [`solve`], the default backend assembles a **dense** `n×n` matrix and
/// factors it in `O(n³)` / `O(n²)` time/storage; enable the `russell` feature
/// for true sparse-direct solves on large systems.
pub fn solve_multi(coo: &Coo, rhs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, SparseError> {
    let csr = coo.to_csr();
    let n = csr.nrows;
    if csr.ncols != n {
        return Err(SparseError::Numeric(format!(
            "solve requires a square matrix, got {n} x {}",
            csr.ncols
        )));
    }
    for r in rhs {
        if r.len() != n {
            return Err(SparseError::Numeric(format!(
                "rhs length {} does not match matrix dimension {n}",
                r.len()
            )));
        }
    }

    // Large symmetric systems: try the sparse envelope Cholesky first. It only
    // succeeds on an SPD matrix (any non-positive pivot is an error), in which
    // case it is both exact and far cheaper than the dense LU; otherwise fall
    // through to the dense path below.
    if n >= SKYLINE_THRESHOLD && csr.is_symmetric(1e-12) {
        if let Ok(factor) = SkylineFactor::new(&csr) {
            return rhs.iter().map(|r| factor.solve(r)).collect();
        }
    }

    // Large unsymmetric / indefinite systems: banded LU when the reordered band
    // is narrow enough to be a clear win over the dense `n x n` matrix.
    if n >= SKYLINE_THRESHOLD {
        let (perm, inv, kl, ku) = BandedLu::bandwidths(&csr);
        if (2 * kl + ku + 1) * 4 <= n {
            if let Ok(lu) = BandedLu::new(&csr, perm, &inv, kl, ku) {
                return Ok(rhs.iter().map(|r| lu.solve(r)).collect());
            }
        }
    }

    let mut data = vec![0.0_f64; n * n];
    for r in 0..n {
        for idx in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            data[r + csr.col_ind[idx] * n] += csr.values[idx];
        }
    }
    let mat = DMatrix::from_vec(n, n, data);

    let mut out = Vec::with_capacity(rhs.len());
    for b in rhs {
        let x = mat
            .solve(&DVector::from_vec(b.clone()))
            .map_err(|e| SparseError::Numeric(format!("{e}")))?;
        out.push(x.iter().copied().collect());
    }
    Ok(out)
}

/// Optional `russell_sparse` (MUMPS/UMFPACK) backend for large-scale problems.
///
/// Enabled by the `russell` feature. `russell_sparse` wraps external
/// Fortran/C solvers and its build requires a SuiteSparse/MUMPS toolchain
/// (see <https://github.com/cpmech/russell>), so this module is compiled only
/// when that feature is active.
#[cfg(feature = "russell")]
mod russell {
    use super::*;
    use russell_lab::Vector;
    use russell_sparse::prelude::*;

    /// Number of rows (unknowns) described by `coo`.
    fn dim(coo: &Coo) -> usize {
        coo.rows
            .iter()
            .cloned()
            .chain(std::iter::once(0))
            .max()
            .map(|m| m + 1)
            .unwrap_or(0)
    }

    /// Build a `russell_sparse` COO matrix from our duplicate-summing `Coo`.
    fn to_russell_coo(coo: &Coo) -> Result<CooMatrix, SparseError> {
        let n = dim(coo);
        let mut rc = CooMatrix::new(n, n, coo.len(), Sym::No)
            .map_err(|e| SparseError::Creation(e.to_string()))?;
        for i in 0..coo.len() {
            rc.put(coo.rows[i], coo.cols[i], coo.vals[i])
                .map_err(|e| SparseError::Creation(e.to_string()))?;
        }
        Ok(rc)
    }

    /// Solve `A x = b` with the `russell_sparse` direct solver.
    ///
    /// `genie` selects the underlying library (`Genie::Umfpack` is the usual
    /// default; `Genie::Mumps` helps for very large or symmetric systems).
    pub fn solve_russell(coo: &Coo, rhs: &[f64], genie: Genie) -> Result<Vec<f64>, SparseError> {
        let n = dim(coo);
        if rhs.len() != n {
            return Err(SparseError::Numeric(format!(
                "rhs length {} does not match matrix dimension {n}",
                rhs.len()
            )));
        }
        let rc = to_russell_coo(coo)?;
        let b = Vector::from(&rhs);
        let mut x = Vector::new(n);
        LinSolver::compute(genie, &mut x, &rc, &b, None)
            .map_err(|e| SparseError::Numeric(e.to_string()))?;
        Ok(x.as_data().to_vec())
    }

    /// Solve `A x_k = rhs[k]` for several right-hand sides against the same `A`,
    /// performing a single factorization (one factorization, multiple solves) —
    /// the sparse analogue of [`solve_multi`].
    pub fn solve_russell_multi(
        coo: &Coo,
        rhs: &[Vec<f64>],
        genie: Genie,
    ) -> Result<Vec<Vec<f64>>, SparseError> {
        let n = dim(coo);
        let rc = to_russell_coo(coo)?;
        let mut solver = LinSolver::new(genie).map_err(|e| SparseError::Numeric(e.to_string()))?;
        solver
            .actual
            .factorize(&rc, None)
            .map_err(|e| SparseError::Numeric(e.to_string()))?;
        let mut out = Vec::with_capacity(rhs.len());
        for b in rhs {
            if b.len() != n {
                return Err(SparseError::Numeric(format!(
                    "rhs length {} does not match matrix dimension {n}",
                    b.len()
                )));
            }
            let rhs_vec = Vector::from(b);
            let mut x = Vector::new(n);
            solver
                .actual
                .solve(&mut x, &rhs_vec, false)
                .map_err(|e| SparseError::Numeric(e.to_string()))?;
            out.push(x.as_data().to_vec());
        }
        Ok(out)
    }
}

#[cfg(feature = "russell")]
pub use russell::{solve_russell, solve_russell_multi};

#[cfg(test)]
mod tests {
    use super::*;

    /// 1-D Laplacian (tridiagonal 2/-1) with `n` unknowns.
    fn laplacian_1d(n: usize) -> Coo {
        let mut c = Coo::new();
        for i in 0..n {
            c.push(i, i, 2.0);
            if i + 1 < n {
                c.push(i, i + 1, -1.0);
                c.push(i + 1, i, -1.0);
            }
        }
        c
    }

    #[test]
    fn cg_matches_dense_lu() {
        let n = 40;
        let a = laplacian_1d(n);
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin() + 1.0).collect();
        let direct = solve(&a, &b).unwrap();
        let cg = solve_cg(&a, &b, &CgOptions::default()).unwrap();
        for (x, y) in direct.iter().zip(&cg.x) {
            assert!((x - y).abs() < 1e-7, "{x} vs {y}");
        }
        assert!(cg.relative_residual <= 1e-10);
    }

    /// 2-D 5-point Laplacian on an `m x m` grid.
    fn laplacian_2d(m: usize) -> Coo {
        let mut c = Coo::new();
        let id = |i: usize, j: usize| i * m + j;
        for i in 0..m {
            for j in 0..m {
                c.push(id(i, j), id(i, j), 4.0);
                if i + 1 < m {
                    c.push(id(i, j), id(i + 1, j), -1.0);
                    c.push(id(i + 1, j), id(i, j), -1.0);
                }
                if j + 1 < m {
                    c.push(id(i, j), id(i, j + 1), -1.0);
                    c.push(id(i, j + 1), id(i, j), -1.0);
                }
            }
        }
        c
    }

    #[test]
    fn skyline_matches_dense_lu() {
        let a = laplacian_2d(9);
        let b: Vec<f64> = (0..81).map(|i| (i as f64 * 0.21).cos() + 2.0).collect();
        let dense = solve(&a, &b).unwrap();
        let sky = solve_skyline(&a, &b).unwrap();
        for (x, y) in dense.iter().zip(&sky) {
            assert!((x - y).abs() < 1e-10, "{x} vs {y}");
        }
    }

    #[test]
    fn skyline_scales_beyond_dense_lu() {
        // 120 x 120 = 14 400 unknowns: ~1.6 GB dense, modest for the envelope.
        let m = 120;
        let a = laplacian_2d(m);
        let b = vec![1.0; m * m];
        let x = solve_skyline(&a, &b).unwrap();
        let r = a.to_csr().matvec(&x);
        let err = r
            .iter()
            .zip(&b)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-8, "max residual {err}");
    }

    #[test]
    fn skyline_handles_scrambled_ordering_and_rejects_bad_input() {
        // A permuted tridiagonal matrix: RCM must recover the band.
        let n = 30;
        let p: Vec<usize> = (0..n).map(|i| (i * 7) % n).collect();
        let mut c = Coo::new();
        for i in 0..n {
            c.push(p[i], p[i], 2.0);
            if i + 1 < n {
                c.push(p[i], p[i + 1], -1.0);
                c.push(p[i + 1], p[i], -1.0);
            }
        }
        let b: Vec<f64> = (0..n).map(|i| i as f64 + 1.0).collect();
        let x = solve_skyline(&c, &b).unwrap();
        let r = c.to_csr().matvec(&x);
        assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-9));

        let mut indef = Coo::new();
        indef.push(0, 0, 1.0);
        indef.push(1, 1, -1.0);
        assert!(solve_skyline(&indef, &[1.0, 1.0]).is_err());
        assert!(solve_skyline(&laplacian_2d(2), &[1.0]).is_err());
        assert!(solve_skyline(&Coo::new(), &[]).unwrap().is_empty());
    }

    #[test]
    fn solve_auto_dispatches_large_spd_and_falls_back_otherwise() {
        // Large SPD: `solve` should handle 3 600 unknowns quickly (the envelope
        // path; the dense path would assemble a 3 600 x 3 600 matrix).
        let m = 60;
        let a = laplacian_2d(m);
        let b = vec![1.0; m * m];
        let x = solve(&a, &b).unwrap();
        let r = a.to_csr().matvec(&x);
        assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-8));
        assert!(a.to_csr().is_symmetric(1e-12));

        // Large *unsymmetric* tridiagonal system must still be solved correctly
        // (the skyline path must not silently read only a triangle).
        let n = 250;
        let mut c = Coo::new();
        for i in 0..n {
            c.push(i, i, 4.0);
            if i + 1 < n {
                c.push(i, i + 1, -1.0);
                c.push(i + 1, i, -2.0);
            }
        }
        assert!(!c.to_csr().is_symmetric(1e-12));
        let b: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let x = solve(&c, &b).unwrap();
        let r = c.to_csr().matvec(&x);
        assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-8));

        // Large symmetric *indefinite* system falls back to the dense LU.
        let mut d = Coo::new();
        for i in 0..n {
            d.push(i, i, if i % 2 == 0 { 1.0 } else { -1.0 });
        }
        let x = solve(&d, &b).unwrap();
        assert!((x[1] + 1.0).abs() < 1e-12 && (x[2] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn banded_lu_handles_unsymmetric_indefinite_and_scrambled() {
        // Unsymmetric, diagonally weak (needs pivoting), scrambled ordering.
        let n = 120;
        let p: Vec<usize> = (0..n).map(|i| (i * 37) % n).collect();
        let mut c = Coo::new();
        for i in 0..n {
            c.push(p[i], p[i], if i % 3 == 0 { 0.1 } else { 3.0 });
            if i + 1 < n {
                c.push(p[i], p[i + 1], 2.0);
                c.push(p[i + 1], p[i], -1.5);
            }
            if i + 2 < n {
                c.push(p[i + 2], p[i], 0.7);
            }
        }
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.3).sin() + 1.0).collect();
        let x = solve_banded(&c, &b).unwrap();
        let r = c.to_csr().matvec(&x);
        assert!(r.iter().zip(&b).all(|(u, v)| (u - v).abs() < 1e-8));
        let dense = solve(&c, &b).unwrap();
        assert!(x.iter().zip(&dense).all(|(u, v)| (u - v).abs() < 1e-7));

        // Singular and malformed inputs are errors, not panics.
        let mut sing = Coo::new();
        sing.push(0, 0, 1.0);
        sing.push(1, 0, 1.0);
        assert!(solve_banded(&sing, &[1.0, 1.0]).is_err());
        assert!(solve_banded(&c, &[1.0]).is_err());
        assert!(solve_banded(&Coo::new(), &[]).unwrap().is_empty());
    }

    #[test]
    fn banded_lu_scales_for_unsymmetric_2d_convection() {
        // 100 x 100 convection-diffusion stencil (unsymmetric, diagonally
        // dominant so the system is well conditioned): 10 000 unknowns.
        let m = 100;
        let id = |i: usize, j: usize| i * m + j;
        let mut c = Coo::new();
        for i in 0..m {
            for j in 0..m {
                c.push(id(i, j), id(i, j), 9.0);
                if i + 1 < m {
                    c.push(id(i, j), id(i + 1, j), -1.0);
                    c.push(id(i + 1, j), id(i, j), -1.6);
                }
                if j + 1 < m {
                    c.push(id(i, j), id(i, j + 1), -1.0);
                    c.push(id(i, j + 1), id(i, j), -1.0);
                }
            }
        }
        let b = vec![1.0; m * m];
        // `solve` should take the banded path (dense would need 800 MB).
        let x = solve(&c, &b).unwrap();
        let r = c.to_csr().matvec(&x);
        let err = r
            .iter()
            .zip(&b)
            .map(|(u, v)| (u - v).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-8, "max residual {err}");
    }

    #[test]
    fn cg_scales_beyond_dense_lu() {
        // 4 000 unknowns is already ~128 MB as a dense matrix and O(n³) to factor.
        let n = 4_000;
        let a = laplacian_1d(n);
        let b = vec![1.0; n];
        let opts = CgOptions {
            tol: 1e-8,
            max_iter: 0,
        };
        let cg = solve_cg(&a, &b, &opts).unwrap();
        // Residual check against the original system.
        let r = a.to_csr().matvec(&cg.x);
        let err = r
            .iter()
            .zip(&b)
            .map(|(p, q)| (p - q).abs())
            .fold(0.0, f64::max);
        assert!(err < 1e-3, "max residual {err}");
    }

    #[test]
    fn cg_rejects_bad_input() {
        let a = laplacian_1d(3);
        assert!(solve_cg(&a, &[1.0, 2.0], &CgOptions::default()).is_err());
        let mut neg = Coo::new();
        neg.push(0, 0, -1.0);
        assert!(solve_cg(&neg, &[1.0], &CgOptions::default()).is_err());
        let tiny = CgOptions {
            tol: 1e-14,
            max_iter: 2,
        };
        assert!(solve_cg(&laplacian_1d(50), &vec![1.0; 50], &tiny).is_err());
        let zero = solve_cg(&a, &[0.0; 3], &CgOptions::default()).unwrap();
        assert_eq!(zero.iterations, 0);
    }

    #[test]
    fn coo_sums_duplicates_and_sorts() {
        let mut c = Coo::new();
        c.push(1, 1, 1.0);
        c.push(0, 1, 1.0);
        c.push(0, 0, 1.0);
        c.push(0, 0, 1.0);
        c.push(1, 0, 1.0);
        let csr = c.to_csr();
        assert_eq!(csr.nrows, 2);
        assert_eq!(csr.ncols, 2);
        assert_eq!(csr.row_ptrs, vec![0, 2, 4]);
        assert_eq!(csr.col_ind, vec![0, 1, 0, 1]);
        assert_eq!(csr.values, vec![2.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn solve_2x2() {
        // [[2, 1], [1, 3]] x = [3, 5]  =>  x = [0.8, 1.4]
        let mut c = Coo::new();
        c.push(0, 0, 2.0);
        c.push(0, 1, 1.0);
        c.push(1, 0, 1.0);
        c.push(1, 1, 3.0);
        let x = solve(&c, &[3.0, 5.0]).expect("solve");
        assert!((x[0] - 0.8).abs() < 1e-10);
        assert!((x[1] - 1.4).abs() < 1e-10);
    }

    #[test]
    fn solve_3x3_diagonally_dominant() {
        // [[4, -1, 0], [-1, 4, -1], [0, -1, 4]] x = [3, 2, 3]
        let mut c = Coo::new();
        c.push(0, 0, 4.0);
        c.push(0, 1, -1.0);
        c.push(1, 0, -1.0);
        c.push(1, 1, 4.0);
        c.push(1, 2, -1.0);
        c.push(2, 1, -1.0);
        c.push(2, 2, 4.0);
        let x = solve(&c, &[3.0, 2.0, 3.0]).expect("solve");
        // Hand-checked: x0 = 1, x1 = 1, x2 = 1.
        for v in x {
            assert!((v - 1.0).abs() < 1e-9, "got {v}");
        }
    }

    #[test]
    fn solve_multi_matches_repeated_solve() {
        let mut c = Coo::new();
        c.push(0, 0, 2.0);
        c.push(0, 1, 1.0);
        c.push(1, 0, 1.0);
        c.push(1, 1, 3.0);
        let x1 = solve(&c, &[3.0, 5.0]).unwrap();
        let x2 = solve(&c, &[1.0, 1.0]).unwrap();
        let both = solve_multi(&c, &[vec![3.0, 5.0], vec![1.0, 1.0]]).unwrap();
        assert_eq!(both.len(), 2);
        for (a, b) in x1.iter().zip(&both[0]) {
            assert!((a - b).abs() < 1e-10);
        }
        for (a, b) in x2.iter().zip(&both[1]) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn solve_with_duplicates() {
        // Same 2x2 system written in duplicate pieces; summing must recover it.
        let mut c = Coo::new();
        c.push(0, 0, 1.0);
        c.push(0, 0, 1.0);
        c.push(0, 1, 0.5);
        c.push(0, 1, 0.5);
        c.push(1, 0, 0.5);
        c.push(1, 0, 0.5);
        c.push(1, 1, 1.5);
        c.push(1, 1, 1.5);
        let x = solve(&c, &[3.0, 5.0]).expect("solve");
        assert!((x[0] - 0.8).abs() < 1e-10);
        assert!((x[1] - 1.4).abs() < 1e-10);
    }
}

#[cfg(all(test, feature = "russell"))]
mod russell_tests {
    use super::*;
    use russell_sparse::Genie;

    #[test]
    fn russell_solve_2x2() {
        // [[2, 1], [1, 3]] x = [3, 5]  =>  x = [0.8, 1.4]
        let mut c = Coo::new();
        c.push(0, 0, 2.0);
        c.push(0, 1, 1.0);
        c.push(1, 0, 1.0);
        c.push(1, 1, 3.0);
        let x = solve_russell(&c, &[3.0, 5.0], Genie::Umfpack).expect("solve");
        assert!((x[0] - 0.8).abs() < 1e-10);
        assert!((x[1] - 1.4).abs() < 1e-10);
    }

    #[test]
    fn russell_solve_multi_matches() {
        let mut c = Coo::new();
        c.push(0, 0, 2.0);
        c.push(0, 1, 1.0);
        c.push(1, 0, 1.0);
        c.push(1, 1, 3.0);
        let x1 = solve_russell(&c, &[3.0, 5.0], Genie::Umfpack).unwrap();
        let x2 = solve_russell(&c, &[1.0, 1.0], Genie::Umfpack).unwrap();
        let both =
            solve_russell_multi(&c, &[vec![3.0, 5.0], vec![1.0, 1.0]], Genie::Umfpack).unwrap();
        assert_eq!(both.len(), 2);
        for (a, b) in x1.iter().zip(&both[0]) {
            assert!((a - b).abs() < 1e-10);
        }
        for (a, b) in x2.iter().zip(&both[1]) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn csr_matvec_matches_dense() {
        // A = [[2, 1], [0, 3]] (row 1 has an implicit zero), x = [1, -2]:
        // A·x = [0, -6].
        let mut c = Coo::new();
        c.push(0, 0, 2.0);
        c.push(0, 1, 1.0);
        c.push(1, 1, 3.0);
        let csr = c.to_csr();
        assert_eq!(csr.nnz(), 3);
        let y = csr.matvec(&[1.0, -2.0]);
        assert_eq!(y.len(), 2);
        assert!(y[0].abs() < 1e-14);
        assert!((y[1] + 6.0).abs() < 1e-14);
    }

    #[test]
    fn csr_matvec_empty_matrix_is_zero_vector() {
        // A placeholder/empty matrix must return a zero vector of the input's
        // length — the contract `coo_matvec` relies on for optional operators
        // like damping.
        let csr = Coo::new().to_csr();
        let y = csr.matvec(&[1.0, 2.0, 3.0]);
        assert_eq!(y.len(), 3);
        assert!(y.iter().all(|&v| v == 0.0));
    }
}
