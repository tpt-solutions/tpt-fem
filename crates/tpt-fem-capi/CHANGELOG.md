# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Initial C ABI (`cdylib` + `staticlib`) with a cbindgen-generated
  `include/tpt_fem.h`: mesh (`tpt_mesh_load` / `tpt_mesh_box` / coords / node
  queries / `tpt_mesh_write_vtk`), `tpt_solve_poisson` (constant or callback
  source), `tpt_solve_elasticity`, `tpt_solve_modal`, `tpt_topopt_cantilever`.
- Status codes + thread-local `tpt_last_error_message()`; panics are caught at
  the FFI boundary.
- `examples/poisson.c` and FFI integration tests.
