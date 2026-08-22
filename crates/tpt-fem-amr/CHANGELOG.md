# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned

- 3-D octree refinement with hanging-node constraints.
- Coarsening (de-refinement) rounds between solve cycles.

## [0.1.0] - 2026-08-13

### Added

- `QuadTree` / `CellKey`: leaf set of a quadtree over `[0,1]²` with
  `refine` and `balance` (1-irregularity enforcement).
- `build_mesh`: conforming `HangingMesh` extraction with hanging-node
  constraints `u_h = (u_a + u_b)/2`.
- `solve_poisson`: Q1 finite-element solve of `−Δu = f` with hanging-node
  elimination (2×2 Gauss quadrature).
- `zz_estimates`: Zienkiewicz–Zhu gradient-recovery per-element error
  estimates.
- `solve_adaptive` / `AmrOptions` / `AdaptiveSolution`: full
  identify–mark–refine loop with Dörfler bulk marking and an element budget.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-amr-0.1.0
