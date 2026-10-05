# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `solve` / `solve_multi` now route systems with ≥ 200 unknowns that are symmetric positive-definite to the envelope Cholesky automatically (every FEM solver in the workspace benefits, no API change); unsymmetric or indefinite systems, and any Cholesky failure, still use the dense LU. Results agree with the dense path to round-off.

### Added

- `solve_cg` / `CgOptions` / `CgSolution` — Jacobi-preconditioned conjugate-gradient solver for symmetric positive-definite systems. Pure Rust, `O(nnz)` memory and work per iteration, so SPD systems scale far beyond the dense-LU `solve` without needing the `russell` SuiteSparse/MUMPS toolchain.
- `solve_skyline` — pure-Rust sparse *direct* solver for SPD systems: reverse Cuthill-McKee reordering + envelope (skyline) Cholesky. `O(n·bandwidth)` memory, exact to round-off, independent of conditioning; a 14 400-unknown 2-D Laplacian solves in well under a second where the dense LU would need ~1.6 GB.
- `Csr::is_symmetric(rel_tol)`.
- `solve_banded` — pure-Rust sparse direct solver for general (unsymmetric / indefinite) systems: reverse Cuthill-McKee reordering + banded LU with partial pivoting. `solve` / `solve_multi` use it automatically for large non-SPD systems whose reordered band is narrow enough (`(2·kl+ku+1)·4 ≤ n`).

## [0.1.0] - 2026-08-23

### Added

- `Coo` coordinate-list accumulator with duplicate-summing `push`.
- `Csr` compressed-sparse-row matrix produced by `Coo::to_csr`.
- `solve` and `solve_multi` backed by an in-house `tpt-math-linalg-dense`
  partial-pivot LU factorisation.
- `SparseError` error type for singular / non-finite systems.
- Optional `russell` feature exposing `solve_russell`/`solve_russell_multi`,
  a `russell_sparse` (SuiteSparse/MUMPS) sparse-direct backend for
  large-scale problems.

### Changed

- Default `solve`/`solve_multi` backend swapped from `faer` sparse LU to the
  in-house dense LU, dropping the Apache-2.0-only `faer` dependency.
- `tpt-math-linalg-dense` dependency now points at the published `0.1.0`
  crates.io release instead of the git/vendored copy.

### Documented after release

_These items were recorded under `[Unreleased]` but are already present in the published 0.1.0 tarball._

#### Added

- `Csr::matvec` — matrix–vector product on the compressed rows: the
  conversion cost is paid once by the caller instead of per call, and the
  row-contiguous accumulation loop is auto-vectoriser friendly. Intended for
  time-stepping / iterative loops that previously re-ran `to_csr()` on every
  matvec.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-sparse-0.1.0
