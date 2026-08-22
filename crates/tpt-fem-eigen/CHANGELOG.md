# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `generalized_lanczos_eigs`: the shift-invert Lanczos now runs two full
  reorthogonalization passes per step ("twice is enough"). With a single pass,
  closely-clustered eigenvalues (Ritz targets separated by ~1e-6) lost basis
  orthogonality and polluted the projected tridiagonal's extreme Ritz values
  by ~1e-2; the clustered spectrum is now resolved to ~1e-12.

### Fixed

- `solve_shifted` (internal): Gaussian elimination now uses partial pivoting
  and regularises collapsed pivots to `±1e-12`. Previously a shift placed at
  or between clustered eigenvalues (near-singular `K - σM`) could produce
  inf/NaN that poisoned the whole Lanczos recurrence.

### Added

- Regression tests `generalized_closely_clustered_eigenvalues` (six
  eigenvalues packed into a ~4e-6 band, checked against the closed form) and
  `generalized_shift_inside_cluster_is_accurate` (shift mid-cluster, stressing
  the near-singular shift-invert solve).

## [0.1.0] - 2026-08-13

### Added

- `matvec` / `rayleigh` sparse matvec and Rayleigh-quotient helpers.
- `power_iteration` dominant-eigenpair solver.
- `inverse_iteration` shift-invert eigenpair solver.
- `EigWhich` eigenvalue selection enum.
- `lanczos_eigs` Lanczos tridiagonalisation eigensovle.
- `generalized_lanczos_eigs` for the `K x = λ M x` eigenproblem.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-eigen-0.1.0
