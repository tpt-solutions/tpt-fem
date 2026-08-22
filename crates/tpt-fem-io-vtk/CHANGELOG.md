# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-23

### Added

- `PointData` named per-node scalar field.
- `mesh_to_vtk` to build a `vtkio` dataset from a mesh + point data.
- `write_vtk` / `write_vtk_ascii` legacy VTK export.
- `write_vtk_with_data` one-call export with point data.
- `VtkError` error type.
- `read_vtk` / `mesh_from_vtk` to import a `Mesh` from a legacy `.vtk`/`.vtu`
  file, promoted from the `tpt-fem-cli`-only reader so VTK round-tripping is
  reachable from the umbrella crate and Python bindings. Supports the linear
  and quadratic (`P2`) cell types written by `mesh_to_vtk`.
- `VtkError::Parse` variant for malformed input (non-unstructured grids,
  unsupported cell types, non-`f64` coordinate buffers).

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-vtk-0.1.0
