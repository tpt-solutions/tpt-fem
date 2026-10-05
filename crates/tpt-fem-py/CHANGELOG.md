# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-23

### Added

- `topopt_cantilever` + `TopOptSolution`: SIMP topology optimization of a 2-D
  cantilever with `to_numpy()` density grid and compliance history.
- `Mesh` Python class: `load`, `box_mesh`, `coords`, `nodes_on_plane`,
  `nodes_in_box`, `write_vtk`.
- `solve_poisson` accepting a constant source or a Python callable `f(x, y, z)`.
- `pyo3` / `maturin` bindings excluded from the Cargo workspace (dev-only this
  pass).
- Core-crate errors surfaced as Python exceptions via `Display` impls.
- `solve_poisson` / `solve_elasticity` / `solve_modal` now return Jupyter-friendly
  result objects (`PoissonSolution`, `ElasticitySolution`, `ModalSolution` of
  `ModeShape`) instead of bare lists, with rich `__repr__` / `_repr_html_`
  display and `to_numpy()` / `to_pyvista()` interop accessors (the `viz` extra
  pulls in `numpy` and `pyvista`).
- `tpt_fem.pyi` type stub covering `Mesh` and the result/solver API, shipped
  alongside the compiled extension via `[tool.maturin] include`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-py-0.1.0
