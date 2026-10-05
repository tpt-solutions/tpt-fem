# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-10-06

### Fixed

- A `connect*` variable whose suffix is not a number is now rejected with `ExodusError::Parse` instead of being silently attributed to block 1 (which could mis-assign cell types and corrupt connectivity).

## [0.1.0] - 2026-08-13

### Added

- Minimal NetCDF-3 classic (CDF-1, big-endian) codec.
- `read_exodus` / `bytes_to_mesh` Exodus II readers.
- `write_exodus` / `mesh_to_exodus_bytes` Exodus II writers.
- Round-trip of linear meshes (`coords`, `connectN`, element-block metadata,
  numbering maps, `time_whole`).
- `ExodusError` error type.

### Documented after release

_These items were recorded under `[Unreleased]` but are already present in the published 0.1.0 tarball._

#### Changed

- `mesh_to_exodus_bytes` now returns `Result<Vec<u8>, ExodusError>`; the NetCDF-3
  encoder (`encode_nc3` / `build_header`) propagates dimension-product overflow as
  an error instead of `.expect()`-ing on the writer path, so encoding can no longer
  panic.

[0.1.1]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-exodus-0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-io-exodus-0.1.0
