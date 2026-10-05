# Migrating from 0.1 to 0.2

0.2 removes the remaining *reachable* panics from the public solver APIs.
Functions that used to abort on bad input (an out-of-range quadrature order, a
singular matrix, an empty or unsupported mesh) now return a `Result` with a
typed error. The numerical behaviour for valid input is unchanged.

Mechanical fix for almost every site: add `?` (or `.unwrap()` in examples and
tests). All of the new error enums implement `std::error::Error`, and the
`tpt-fem` umbrella's `tpt_fem::Error` converts from every one of them, so a
single `-> Result<_, tpt_fem::Error>` covers a whole pipeline.

## Crates that need a version bump in your `Cargo.toml`

Breaking (0.1 → 0.2): `tpt-fem`, `-element`, `-assembly`, `-thermal`,
`-elasticity`, `-eigen`, `-contact`, `-dynamic`, `-fluid`, `-porous`, `-modal`,
`-coupling`, `-cli`.
Compatible patch releases (0.1.0 → 0.1.1): `-quadrature`, `-sparse`, `-mesh`,
`-mesh-gen`, `-io-abaqus`, `-io-exodus`, `-composite`, `-topopt`.
Unchanged: `-amr`, `-dofmap`, `-hyperelastic`, `-io-vtk`, `-plasticity`, `-solve`.

## Signature changes

| Crate | Item | 0.1 | 0.2 |
|-------|------|-----|-----|
| element | `line_rule`, `quad_rule`, `hex_rule` | `-> Quad1D/2D/3D` | `-> Result<_, ElementError>` |
| assembly | `apply_neumann_order`, `apply_robin_order` | `()` | `Result<(), SparseError>` (`apply_neumann`/`apply_robin` unchanged) |
| thermal | `poisson_element_matrix`, `poisson_source_vector` | bare value | `Result<_, ThermalError>` |
| elasticity | `elasticity_body_vector`, `elasticity_mass_matrix`, `elasticity_lumped_mass` | bare value | `Result<_, ElasticityError>` |
| eigen | `lanczos_eigs` | `Vec<(f64, Vec<f64>)>` | `Result<Vec<(f64, Vec<f64>)>, SparseError>` |
| contact | `augmented_lagrangian` | `(Vec<f64>, Vec<f64>)` | `Result<_, ContactError>` |
| dynamic | `newmark` | `Vec<(f64, Vec<f64>)>` | `Result<_, DynamicError>` |
| modal | `ModalData::modal_superposition` | `Vec<(f64, Vec<f64>)>` | `Result<_, DynamicError>` |
| fluid | `stokes_dofmap` | `MultiFieldDofMap` | `Result<_, FluidError>` |
| fluid | `transient_stokes` | `Vec<(f64, Vec<f64>)>` | `Result<_, FluidError>` |
| porous | `solve_darcy` | `Result<_, SparseError>` | `Result<_, PorousError>` |
| porous | `terzaghi_consolidation` | `Vec<(f64, f64, f64)>` | `Result<_, PorousError>` |
| coupling | `fsi_interface_loads` | `Vec<f64>` | `Result<_, CouplingError>` |

New and extended error variants: `FluidError::{EmptyMesh, UnsupportedCell, Dynamic}`,
`CouplingError::Quadrature`, `ElasticityError::Quadrature`, `ContactError`,
`PorousError`.

## Before / after

```rust
// 0.1
let rule = quad_rule(3);
let hist = newmark(&m, &c, &k, &u0, &v0, load, &opts, nsteps);

// 0.2
let rule = quad_rule(3)?;
let hist = newmark(&m, &c, &k, &u0, &v0, load, &opts, nsteps)?;
```

## Behaviour fixes worth knowing about

- `tpt-fem-quadrature`: the `TetrahedronRule::Keast4` 11-point table had
  points outside the reference tetrahedron; `Tet10` results computed with it in
  0.1.0 were wrong (the table is now regression-tested).
- `tpt-fem-eigen`: shift-invert Lanczos now reorthogonalises twice per step, so
  closely clustered eigenvalues resolve correctly.
- Minimum supported Rust version is **1.85** (0.1 declared 1.75 but could not
  build).
- New (additive): `tpt_fem_sparse::solve_cg` — a pure-Rust Jacobi-preconditioned
  conjugate-gradient solver for SPD systems that scales to far larger problems
  than the dense-LU `solve`.
