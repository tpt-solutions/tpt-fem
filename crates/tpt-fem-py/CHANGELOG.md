# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `solve_transient_heat` and the `HeatHistory` result object (`times`, indexing to a per-step field, `to_numpy()` of shape `(nsteps + 1, n_nodes)`): transient heat conduction by the theta-method. Type stubs and a pytest case included.
- `j2_uniaxial_response` (J2 elastic-plastic stress along a monotonic strain path, isotropic + kinematic hardening) and `neo_hookean_uniaxial` (incompressible neo-Hookean nominal stress): material-point drivers over `tpt-fem-plasticity` / `tpt-fem-hyperelastic`, for calibrating against test data. Stubs and pytest cases included.
- `solve_darcy` (steady Darcy pressure field, over `tpt-fem-porous`) and `solve_stokes` (steady Stokes flow by the penalty method, returns `(velocity, pressure)`, over `tpt-fem-fluid`).
- `laminate_abd` (classical lamination theory ABD matrix over `tpt-fem-composite`) and `newmark` (implicit Newmark-beta integration of a small dense `M u'' + C u' + K u = f(t)` system, constant or callable load, over `tpt-fem-dynamic`).
- `solve_thermal_structural` (free thermal expansion, over `tpt-fem-coupling`) and `contact_pairs` (octree nearest-point pairing of two point sets, over `tpt-fem-contact`; returns the index into `b`).
- `contact_augmented_lagrangian`: unilateral (non-penetration) contact on a small dense system by the augmented-Lagrangian method, returning displacements and contact forces.
- `fsi_interface_loads` (work-consistent fluid-structure interface load vector, over `tpt-fem-coupling`).
- Optional `gpu` Cargo feature (`maturin develop --features gpu`): `gpu_enabled()`, `gpu_adapter()` and `gpu_solve_cg(triplets, rhs, tol)` — GPU-resident Jacobi-PCG with f64 refinement over the new `tpt-fem-gpu` crate (3.6x-7.7x faster than CPU PCG from 250k to 1.4M unknowns on an RTX 3050).

### Added

- `solve_transient_heat` and the `HeatHistory` result object (`times`, indexing to a per-step field, `to_numpy()` of shape `(nsteps + 1, n_nodes)`): transient heat conduction by the θ-method. Type stubs and a pytest case included.

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
