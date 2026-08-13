# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `matvec` / `rayleigh` sparse matvec and Rayleigh-quotient helpers.
- `power_iteration` dominant-eigenpair solver.
- `inverse_iteration` shift-invert eigenpair solver.
- `EigWhich` eigenvalue selection enum.
- `lanczos_eigs` Lanczos tridiagonalisation eigensovle.
- `generalized_lanczos_eigs` for the `K x = λ M x` eigenproblem.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-eigen-0.1.0
