# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `poisson_element_matrix` element stiffness integration.
- `poisson_source_vector` element load integration.
- `solve_poisson` end-to-end steady Poisson / heat-conduction solve.
- Support for Dirichlet, Neumann, and Robin boundary conditions.
- Method-of-Manufactured-Solutions (MMS) convergence tests asserting P1
  `L2` (order 2) and `H1` (order 1) rates.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-thermal-0.1.0
