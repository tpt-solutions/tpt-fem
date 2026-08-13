# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `ElasticModel` for axial bar, plane-stress, plane-strain, and 3-D continua.
- `elasticity_element_matrix`, `elasticity_body_vector` element operators.
- `elasticity_mass_matrix` / `elasticity_lumped_mass` consistent and lumped mass.
- `solve_elasticity` end-to-end static solve.
- `solve_modal` modal analysis via `tpt-fem-eigen`.
- 2-D Euler–Bernoulli beam support: `BeamSection2D`, `beam2d_element_matrix`,
  `beam2d_consistent_mass`, `solve_frame2d`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-elasticity-0.1.0
