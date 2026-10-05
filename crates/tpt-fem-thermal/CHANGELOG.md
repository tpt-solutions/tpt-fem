# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-10-06

### Changed

- Element assembly (`solve_poisson`, `solve_transient_heat`) now runs in parallel via `try_assemble_parallel`; results are unchanged.

### Added

- `solve_transient_heat`, `TransientHeatOptions` and `heat_capacity_element_matrix`: transient heat conduction `ρc ∂T/∂t − ∇·(k∇T) = f(x, t)` by the θ-method (backward Euler / Crank–Nicolson / forward Euler), with Dirichlet conditions and a time-dependent source; verified against the analytic `exp(−π²t)·sin(πx)` decay and the steady-state limit.
- `ThermalError::Sparse` and `ThermalError::InvalidInput` variants (with `From<SparseError>`).

### Changed

- **Breaking:** `poisson_element_matrix` and `poisson_source_vector` now
  return `Result<_, ThermalError>` instead of panicking when `quad_order`
  (or, for `Quad8`/`Quad9`/`Hex20`/`Hex27` elements, `quad_order + 1`) is
  outside `tpt-fem-quadrature`'s supported `1..=5` range. `solve_poisson`'s
  signature is unchanged (still `Result<_, SparseError>`); `ThermalError`
  converts into `SparseError` automatically.

## [0.1.0] - 2026-08-13

### Added

- `poisson_element_matrix` element stiffness integration.
- `poisson_source_vector` element load integration.
- `solve_poisson` end-to-end steady Poisson / heat-conduction solve.
- Support for Dirichlet, Neumann, and Robin boundary conditions.
- Method-of-Manufactured-Solutions (MMS) convergence tests asserting P1
  `L2` (order 2) and `H1` (order 1) rates.

[0.2.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-thermal-0.2.0
[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-thermal-0.1.0
