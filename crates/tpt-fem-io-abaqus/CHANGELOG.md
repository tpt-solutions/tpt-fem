# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `read_inp` parser for the `*NODE` / `*ELEMENT` subset of Abaqus decks.
- `write_inp` serialiser for `Mesh` → `.inp`.
- Mapping of common linear Abaqus types (T2D2, CPS3/CPS4, C3D4, C3D8) to
  `tpt-fem-mesh` `CellType`s.
- `InpError` error type.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-abaqus-0.1.0
