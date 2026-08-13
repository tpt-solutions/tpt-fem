//! FEM-specific sparse-matrix assembly adapter with a [`faer`]-backed solve.
#![allow(clippy::needless_range_loop)]
//!
//! Finite-element assembly naturally produces a sparse global matrix as a bag of
//! `(row, col, value)` triplets, where the same entry may be written many times
//! (once per element that touches it). [`Coo`] is a growable coordinate-list
//! accumulator that supports duplicate-summing assembly; [`Coo::to_csr`]
//! collapses it into a canonical compressed-sparse-row [`Csr`] matrix.
//!
//! [`solve`] factors the assembled matrix with [`faer`]'s sparse LU
//! decomposition and solves `A x = b`.
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

use faer::sparse::linalg::solvers::{Lu, SpSolver, SymbolicLu};
use faer::sparse::SparseColMat;

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
/// [`Coo`] accumulator, returning the solution vector `x`.
///
/// Duplicate `(row, col)` entries in `coo` are summed (via [`Coo::to_csr`]),
/// the matrix is factored with [`faer`]'s sparse LU decomposition, and the
/// system is solved.
pub fn solve(coo: &Coo, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
    let mut sols = solve_multi(coo, std::slice::from_ref(&rhs.to_vec()))?;
    Ok(sols
        .pop()
        .expect("solve_multi returns one solution per rhs"))
}

/// Solve `A x_k = rhs[k]` for every right-hand side in `rhs` against the
/// *same* matrix `A`, factoring it only once.
///
/// Equivalent to calling [`solve`] once per right-hand side, but far
/// cheaper when there is more than one: callers with multiple RHS vectors
/// against an unchanged `A` (e.g. an arc-length continuation corrector,
/// which needs both a tangent and a residual-correction solve per
/// iteration) should prefer this over repeated `solve` calls.
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

    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(csr.nnz());
    for r in 0..n {
        for idx in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            triplets.push((r, csr.col_ind[idx], csr.values[idx]));
        }
    }

    let mat = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
        .map_err(|e| SparseError::Creation(format!("{e:?}")))?;
    let symbolic =
        SymbolicLu::try_new(mat.symbolic()).map_err(|e| SparseError::Symbolic(format!("{e:?}")))?;
    let lu = Lu::try_new_with_symbolic(symbolic, mat.as_ref())
        .map_err(|e| SparseError::Numeric(format!("{e:?}")))?;

    let ncols = rhs.len();
    let mut b = faer::Mat::<f64>::from_fn(n, ncols, |i, j| rhs[j][i]);
    lu.solve_in_place(b.as_mut());
    Ok((0..ncols)
        .map(|j| (0..n).map(|i| b.read(i, j)).collect())
        .collect())
}

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
