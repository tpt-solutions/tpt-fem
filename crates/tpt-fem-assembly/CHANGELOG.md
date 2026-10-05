# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `try_assemble_parallel` — multi-threaded element assembly (`std::thread::scope`, no new dependency). The triplet lists are concatenated in element order, so the result is bit-for-bit identical to `try_assemble` for any thread count; `threads = 0` uses all cores and small meshes stay single-threaded.

### Fixed

- `solve_with_dirichlet` no longer fails with an opaque dimension error when every DOF is prescribed (it returns the prescribed values), and reports an out-of-range Dirichlet DOF as `SparseError` instead of panicking on an index.

### Changed

- **Breaking:** `apply_neumann_order` / `apply_robin_order` now return
  `Result<(), SparseError>` instead of panicking when `order` is outside
  `tpt-fem-quadrature`'s supported `1..=5` range. `apply_neumann` /
  `apply_robin` (the fixed-order convenience wrappers) are unchanged.

### Fixed

- `reduce_system` looks up prescribed DOF values through a `HashMap` instead of a linear scan of the BC list per matrix entry (removes an `unwrap` and an `O(nnz · n_bc)` cost on large Dirichlet sets); a DOF listed twice now uses the last value.

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

### Documented after release

_These items were recorded under `[Unreleased]` but are already present in the published 0.1.0 tarball._

#### Added

- `try_assemble`, a fallible variant of `assemble` whose element-matrix closure may
  fail (e.g. an invalid model/dimension combination), propagating the first error
  instead of panicking.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-assembly-0.1.0
