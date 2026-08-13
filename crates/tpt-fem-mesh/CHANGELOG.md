# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `Node`, `Element`, `Mesh`, and `DofMap` mesh model types.
- `CellType` enumeration for the five supported linear elements.
- `MeshBuilder` API for constructing meshes in code.
- `Mesh::number_dofs` for configurable DOF-per-node numbering.
- `Mesh::from_msh_bytes` Gmsh `.msh` v4.1 ASCII import via `mshio`.
- `MeshError` including `UnsupportedElementType` for non-linear elements.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-mesh-0.1.0
