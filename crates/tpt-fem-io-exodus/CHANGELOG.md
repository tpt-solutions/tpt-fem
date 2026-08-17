# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `mesh_to_exodus_bytes` now returns `Result<Vec<u8>, ExodusError>`; the NetCDF-3
  encoder (`encode_nc3` / `build_header`) propagates dimension-product overflow as
  an error instead of `.expect()`-ing on the writer path, so encoding can no longer
  panic.

## [0.1.0] - 2026-08-13

### Added

- Minimal NetCDF-3 classic (CDF-1, big-endian) codec.
- `read_exodus` / `bytes_to_mesh` Exodus II readers.
- `write_exodus` / `mesh_to_exodus_bytes` Exodus II writers.
- Round-trip of linear meshes (`coords`, `connectN`, element-block metadata,
  numbering maps, `time_whole`).
- `ExodusError` error type.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-exodus-0.1.0
