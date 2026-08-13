# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- Umbrella crate re-exporting all `tpt-fem-*` core crates behind Cargo features.
- Feature flags for each constituent crate (all enabled by default).
- `prelude` module for convenient glob imports.
- `thermal_solve` example (mesh → solve → VTK).
- Integration tests `tests/patch_test.rs` and `tests/end_to_end.rs`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-0.1.0
