# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-10-06

### Added

- `Mesh::to_stl_string` — ASCII STL of the boundary surface with outward normals (3-D meshes export boundary faces, 2-D meshes the elements themselves; quads split into triangles).
- `Mesh::to_xdmf_string` — XDMF 3 document with embedded XML data and nodal attributes (single linear cell type), openable directly in ParaView without HDF5.
- `ExportError` variants `UnsupportedCell`, `MixedCells`, `Empty`.

### Added

- `Mesh::to_msh_string` — ASCII Gmsh MSH 4.1 writer (inverse of `from_msh_bytes`, including P2 node reordering); regions are not written.
- `Mesh::nodal_csv` — CSV export of node coordinates plus scalar/vector nodal fields.
- `ExportError` for the two export helpers.

### Changed

- Depend on `tpt-fem-element` 0.2.0 (internal use only — no public API change).

## [0.1.0] - 2026-08-13

### Added

- `Node`, `Element`, `Mesh`, and `DofMap` mesh model types.
- `CellType` enumeration for the five supported linear elements.
- `MeshBuilder` API for constructing meshes in code.
- `Mesh::number_dofs` for configurable DOF-per-node numbering.
- `Mesh::from_msh_bytes` Gmsh `.msh` v4.1 ASCII import via `mshio`.
- `MeshError` including `UnsupportedElementType` for non-linear elements.

[0.1.1]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-mesh-0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-mesh-0.1.0
