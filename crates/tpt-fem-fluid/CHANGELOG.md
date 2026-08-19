# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `stokes_dofmap` mixed velocity/pressure DOF-map builder.
- `steady_stokes` Stokes solver with divergence-penalty incompressibility.
- `transient_stokes` time-dependent Stokes relaxation via Newmark.
- `transient_navier_stokes` low-Re Navier–Stokes via semi-implicit Picard.
- Selective reduced integration of the penalty term and lumped nodal pressure recovery.
- Examples: `stokes_dofmap`, `steady_stokes_poiseuille`, `lid_driven_cavity`,
  `transient_stokes_startup`, `navier_stokes_low_re`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-fluid-0.1.0
