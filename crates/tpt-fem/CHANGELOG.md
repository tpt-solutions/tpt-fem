# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Breaking (via re-exports):** the umbrella re-exports the breaking `Result`
  signature changes in `tpt-fem-element`, `-assembly`, `-thermal`,
  `-elasticity`, `-eigen`, `-contact`, `-dynamic`, `-fluid`, `-porous`,
  `-modal` and `-coupling`. See `docs/MIGRATING-0.2.md`.
- Minimum supported Rust version is now 1.85 (the previously declared 1.75
  could not build the dependency tree).

### Added

- `tpt_fem::Error` gains `From` conversions for the new per-crate error enums.

## [0.1.0] - 2026-08-23

### Added

- Umbrella crate re-exporting all `tpt-fem-*` core crates behind Cargo features.
- Feature flags for each constituent crate (all enabled by default).
- `prelude` module for convenient glob imports.
- `thermal_solve` example (mesh → solve → VTK).
- Integration tests `tests/patch_test.rs` and `tests/end_to_end.rs`.
- `amr` feature wiring `tpt-fem-amr`: `build_mesh`, `solve_adaptive`,
  `zz_estimates`, `AdaptiveSolution`, `AmrOptions`, `CellKey`, `HangingMesh`,
  `QuadTree` re-exported from the crate root and `prelude`; the full
  namespace (including `tpt_fem_amr::solve_poisson`, which clashes with
  `tpt-fem-thermal`'s and is therefore excluded from the flat re-export) is
  reachable via the `amr` module.
- `tests/manifest_drift.rs` guarding workspace membership and
  `[workspace.dependencies]` consistency.
- Three end-to-end example showcases chaining multiple constituent crates:
  `thermoelastic_plate_showcase` (conduction → thermal expansion → VTK),
  `contact_column_showcase` (one-sided contact of a bar against a rigid
  stop), and `cantilever_dynamics_showcase` (modal analysis vs. beam theory
  plus Newmark/modal-superposition transient response).

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-0.1.0
