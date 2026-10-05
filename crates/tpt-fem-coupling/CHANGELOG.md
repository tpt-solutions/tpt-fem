# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `fsi_interface_loads` validates that every interface structure node and fluid node is in range and returns `CouplingError::Interface` instead of panicking on an out-of-bounds index.

### Fixed

- `thermal_structural` and the FSI coupling operator return an error for an empty structure/fluid mesh instead of panicking.

### Changed

- **Breaking:** `fsi_interface_loads` now returns `Result<Vec<f64>, CouplingError>`: an empty structure mesh, an unsupported cell type, or a malformed interface pairing is reported as `CouplingError::Interface` instead of panicking.
- The consistent-load quadrature uses `try_gauss_legendre`, so an out-of-range rule order surfaces as the new `CouplingError::Quadrature` variant (`From<QuadratureError>`) rather than aborting.

## [0.1.0] - 2026-08-20

### Added

- `thermal_structural` thermal expansion / stress from a temperature field.
- `joule_source` Ohmic volumetric heat `σ·|E|²`.
- `electro_thermal` steady heat conduction driven by a Joule source.
- `fsi_coupling` explicit fluid→structure substep with pressure-traction transfer.
- Examples: `thermal_bar_expansion`, `thermal_bimetal_strip`, `joule_heating`,
  `fsi_coupling`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-coupling-0.1.0
