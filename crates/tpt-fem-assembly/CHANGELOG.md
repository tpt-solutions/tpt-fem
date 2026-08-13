# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `assemble` to scatter per-element matrices into a global `Coo`.
- `reduce_system` / `ReducedSystem` static condensation of fixed DOFs.
- `solve_with_dirichlet` essential-condition enforcement + sparse solve.
- `boundary_faces` enumeration of mesh boundary faces.
- `apply_neumann` / `apply_neumann_order` natural (flux) boundary loads.
- `apply_robin` / `apply_robin_order` convective boundary contributions.
- Dimension-agnostic support for all five linear element types and arbitrary
  DOFs-per-node.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-assembly-0.1.0
