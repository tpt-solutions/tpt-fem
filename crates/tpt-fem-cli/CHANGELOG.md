# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `tpt-fem` binary with `solve` subcommand (TOML-configured steady
  Poisson/heat-conduction run).
- `mesh info` subcommand for mesh summary statistics.
- `mesh convert` subcommand for Gmsh `.msh` → ParaView `.vtk` conversion.
- Human-readable error reporting via core-crate `Display` impls.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-cli-0.1.0
