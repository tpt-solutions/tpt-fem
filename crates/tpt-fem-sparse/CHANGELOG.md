# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `Coo` coordinate-list accumulator with duplicate-summing `push`.
- `Csr` compressed-sparse-row matrix produced by `Coo::to_csr`.
- `solve` and `solve_multi` backed by `faer` sparse LU factorisation.
- `SparseError` error type for singular / non-finite systems.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-sparse-0.1.0
