# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `ElasticityError::Quadrature` variant, and `elasticity_body_vector`,
  `elasticity_mass_matrix`, `elasticity_lumped_mass` now return
  `Result<_, ElasticityError>` instead of panicking when `quad_order` (or,
  for P2 element types, `quad_order + 1`) is outside `tpt-fem-quadrature`'s
  supported `1..=5` range. `solve_elasticity` / `solve_modal` are unchanged
  (still `Result<_, SparseError>`); `ElasticityError` converts into
  `SparseError` automatically.

### Changed

- **Breaking:** the three points above change the return type of
  `elasticity_body_vector`, `elasticity_mass_matrix`, and
  `elasticity_lumped_mass` from a bare value to `Result`.

## [0.1.0] - 2026-08-13

### Added

- `ElasticModel` for axial bar, plane-stress, plane-strain, and 3-D continua.
- `elasticity_element_matrix`, `elasticity_body_vector` element operators.
- `elasticity_mass_matrix` / `elasticity_lumped_mass` consistent and lumped mass.
- `solve_elasticity` end-to-end static solve.
- `solve_modal` modal analysis via `tpt-fem-eigen`.
- 2-D Euler–Bernoulli beam support: `BeamSection2D`, `beam2d_element_matrix`,
  `beam2d_consistent_mass`, `solve_frame2d`.

### Documented after release

_These items were recorded under `[Unreleased]` but are already present in the published 0.1.0 tarball._

#### Added

- `ElasticityError` surfaced when a model/dimension combination has no
  constitutive definition (e.g. `PlaneStress` on a 3-D mesh), replacing the
  previous `panic!` in `strain_dim` / `constitutive` / `b_matrix`.
- `elasticity_element_matrix` now returns `Result<_, ElasticityError>`; the
  assembly primitive gained a fallible `try_assemble` so element-matrix failures
  propagate instead of panicking.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-elasticity-0.1.0
