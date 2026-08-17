# tpt-fem

Umbrella crate re-exporting the [tpt-fem](https://github.com/tpt-solutions/tpt-fem)
core crates behind Cargo features — the mesh-based finite element method (FEM)
library from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate re-exports the `tpt-fem-*` crates behind Cargo features, so
downstream users can depend on a single crate and opt into only the pieces they
need:

| Feature       | Re-exports                       | Backing crate          |
|---------------|---------------------------------|-----------------------|
| `quadrature`  | `tpt_fem_quadrature::*`         | `tpt-fem-quadrature`  |
| `element`     | `tpt_fem_element::*`            | `tpt-fem-element`     |
| `mesh`        | `tpt_fem_mesh::*`               | `tpt-fem-mesh`        |
| `sparse`      | `tpt_fem_sparse::*`             | `tpt-fem-sparse`      |
| `assembly`    | `tpt_fem_assembly::*`           | `tpt-fem-assembly`    |
| `thermal`     | `tpt_fem_thermal::*`            | `tpt-fem-thermal`     |
| `io-vtk`      | `tpt_fem_io_vtk::*`             | `tpt-fem-io-vtk`      |
| `elasticity`  | `tpt_fem_elasticity::*`         | `tpt-fem-elasticity`  |
| `solve`       | `tpt_fem_solve::*`              | `tpt-fem-solve`       |
| `eigen`       | `tpt_fem_eigen::*`              | `tpt-fem-eigen`       |
| `io-abaqus`   | `tpt_fem_io_abaqus::*`          | `tpt-fem-io-abaqus`   |
| `io-exodus`   | `tpt_fem_io_exodus::*`          | `tpt-fem-io-exodus`   |
| `mesh-gen`    | `tpt_fem_mesh_gen::*`           | `tpt-fem-mesh-gen`    |

All features are enabled by default. A `prelude` module re-exports the most
commonly used items for convenient glob imports.

The end-to-end pipeline — reference-element shape functions and gradients,
quadrature, triplet assembly, and a sparse solve — is exercised by the
integration test `tests/patch_test.rs`.

## Installation

```toml
[dependencies]
tpt-fem = "0.1"
# or opt into only what you need:
# default-features = false, features = ["mesh", "thermal", "io-vtk"]
```

## Usage

```rust
use tpt_fem::prelude::*;

// Build a box mesh, solve a 3D Poisson problem, and export to VTK — all
// through the umbrella prelude.
let mesh = box_mesh(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0), [8, 8, 8]);
let solution = solve_poisson::<_>(mesh.clone(), 1.0, 1.0, &[(0, 0.0)]);
let vtk = mesh_to_vtk(&mesh, &[PointData::new("u", solution.clone())]);
write_vtk_with_data("thermal_solve.vtk", &mesh, &[PointData::new("u", solution)]).unwrap();
```

A runnable version of this example lives at
`crates/tpt-fem/examples/thermal_solve.rs` (`cargo run -p tpt-fem --example
thermal_solve`).

## Constituent crates

- `tpt-fem-quadrature` — Gauss quadrature rules.
- `tpt-fem-element` — reference elements & shape functions.
- `tpt-fem-mesh` — mesh model & Gmsh import.
- `tpt-fem-sparse` — COO/CSR assembly & in-house dense-LU solve (optional `russell` sparse-direct backend).
- `tpt-fem-assembly` — assembly & boundary conditions.
- `tpt-fem-thermal` — Poisson / heat conduction.
- `tpt-fem-elasticity` — linear elasticity & beams.
- `tpt-fem-solve` — nonlinear & continuation solvers.
- `tpt-fem-eigen` — sparse eigenvalue solvers.
- `tpt-fem-io-vtk` / `tpt-fem-io-abaqus` / `tpt-fem-io-exodus` — mesh I/O.
- `tpt-fem-mesh-gen` — tetrahedral mesh generation.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
