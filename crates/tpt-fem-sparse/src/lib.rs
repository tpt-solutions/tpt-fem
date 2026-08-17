//! FEM-specific sparse-matrix assembly adapter with a [`tpt-math-linalg-dense`]-backed solve.
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
/// The default backend assembles a **dense** `n×n` matrix regardless of `A`'s
/// sparsity and factors it in `O(n³)` time with `O(n²)` storage. This is fine
/// for the small, hand-built meshes used to validate the crate, but it does
/// **not** scale to large sparse problems. For those, enable the optional
/// `russell` feature (which dispatches to the `russell_sparse`
/// UMFPACK/MUMPS direct solvers) — see [`solve_russell`]. As a rule of thumb,
/// prefer `russell` once `A` has more than a few thousand rows.
pub fn solve(coo: &Coo, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
    let sols = solve_multi(coo, std::slice::from_ref(&rhs.to_vec()))?;
    sols.into_iter()
        .next()
        .ok_or_else(|| SparseError::Numeric("solve received an empty right-hand side".into()))
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
        let b = Vector::from(rhs);
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
}
