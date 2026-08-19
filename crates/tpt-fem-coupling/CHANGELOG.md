# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `thermal_structural` thermal expansion / stress from a temperature field.
- `joule_source` Ohmic volumetric heat `σ·|E|²`.
- `electro_thermal` steady heat conduction driven by a Joule source.
- `fsi_coupling` explicit fluid→structure substep with pressure-traction transfer.
- Examples: `thermal_bar_expansion`, `thermal_bimetal_strip`, `joule_heating`,
  `fsi_coupling`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-coupling-0.1.0
