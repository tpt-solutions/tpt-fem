# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `solve_darcy` — steady Darcy flow: assembles `K = ∫ k ∇Nᵀ∇N dΩ` for a scalar pressure field over `Line`/`Tri`/`Quad`/`Tet`/`Hex` cells and solves `−∇·(k∇p) = f` with Dirichlet (fixed-head) nodes; element sources supplied as node-order vectors, returns `Result<Vec<f64>, SparseError>`.
- `terzaghi_consolidation` — 1-D Terzaghi consolidation of a saturated column (`Line2` mesh spanning `[0, H]`), integrating `∂u/∂t = cᵥ ∂²u/∂z²` with a backward-Euler step; returns the `(t, settlement, max_u)` history and asserts the explicit-stability limit `dt ≤ Δz²/(2·cᵥ)`.
- Conductivity / mass element matrix helper supporting the scalar-field cell types used by the solvers.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-porous-0.1.0
