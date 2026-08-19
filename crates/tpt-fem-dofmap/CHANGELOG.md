# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `FieldSpec` — a named physical field with a per-node component count.
- `FieldSpec::new` — construct a field spec (`components` must be non-zero).
- `Layout` — `Block` (field-major) and `Interleaved` (node-major) numbering.
- `MultiFieldDofMap` — multi-field DOF map over a `Mesh`.
- `MultiFieldDofMap::new` — build the map for a mesh, field list, and layout.
- `MultiFieldDofMap::node_field_dof` — global DOF for `(node, field, component)`.
- `MultiFieldDofMap::components` — per-node component count of a field.
- `MultiFieldDofMap::field_range` — `(start, count)` of a field's contiguous block.
- `MultiFieldDofMap::dofs_of` — all global DOFs owned by a node.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-dofmap-0.1.0
