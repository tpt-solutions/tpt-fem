# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `PointData` named per-node scalar field.
- `mesh_to_vtk` to build a `vtkio` dataset from a mesh + point data.
- `write_vtk` / `write_vtk_ascii` legacy VTK export.
- `write_vtk_with_data` one-call export with point data.
- `VtkError` error type.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-vtk-0.1.0
