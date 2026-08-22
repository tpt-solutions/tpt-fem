# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `Csr::matvec` — matrix–vector product on the compressed rows: the
  conversion cost is paid once by the caller instead of per call, and the
  row-contiguous accumulation loop is auto-vectoriser friendly. Intended for
  time-stepping / iterative loops that previously re-ran `to_csr()` on every
  matvec.

## [0.1.0] - 2026-08-13

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

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-sparse-0.1.0
