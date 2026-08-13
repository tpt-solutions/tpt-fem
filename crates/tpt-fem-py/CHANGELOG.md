# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `Mesh` Python class: `load`, `box_mesh`, `coords`, `nodes_on_plane`,
  `nodes_in_box`, `write_vtk`.
- `solve_poisson` accepting a constant source or a Python callable `f(x, y, z)`.
- `pyo3` / `maturin` bindings excluded from the Cargo workspace (dev-only this
  pass).
- Core-crate errors surfaced as Python exceptions via `Display` impls.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-py-0.1.0
