# tpt-fem

A modular finite-element method (FEM) workspace for Rust, part of the
`tpt-solutions` ecosystem. It provides a from-the-ground-up FEM core — reference
elements, quadrature, mesh I/O, sparse assembly, physics (thermal, elasticity),
eigen- and continuation-solvers, a native tet-mesh generator, and a CLI driver —
all under a permissive `MIT OR Apache-2.0` policy.

Crates.io publishing is intentionally **out of scope** for this repository; the
crates are consumed as path dependencies within the workspace.

## Crates

| Crate | Purpose | Depends on |
|-------|---------|------------|
| [`tpt-fem-quadrature`](crates/tpt-fem-quadrature) | Fixed-order Gauss quadrature on reference elements (line / tri / quad / tet / hex). | — |
| [`tpt-fem-element`](crates/tpt-fem-element) | Reference elements, Lagrange `P1` shape functions, isoparametric Jacobian map. | `tpt-fem-quadrature` |
| [`tpt-fem-mesh`](crates/tpt-fem-mesh) | Mesh data model, DOF numbering, Gmsh `.msh` v4.1 import (via `mshio`), region selectors. | `mshio` |
| [`tpt-fem-sparse`](crates/tpt-fem-sparse) | FEM-specific COO/CSR assembly adapter with an in-house dense LU solve (optional `russell` sparse-direct backend for large problems). | `tpt-math-linalg-dense` |
| [`tpt-fem-assembly`](crates/tpt-fem-assembly) | Element-to-global scatter, Dirichlet/Neumann/Robin boundary conditions, reduced-system extraction. | `tpt-fem-{sparse,element,mesh}` |
| [`tpt-fem-thermal`](crates/tpt-fem-thermal) | Heat-conduction / Poisson elements (`-∇·(k∇u) = f`). | `tpt-fem-assembly` |
| [`tpt-fem-io-vtk`](crates/tpt-fem-io-vtk) | ParaView `.vtk` / `.vtu` export (via `vtkio`). | `tpt-fem-mesh`, `vtkio` |
| [`tpt-fem-elasticity`](crates/tpt-fem-elasticity) | 2-D Euler–Bernoulli frames, and 2-D/3-D continuum linear elasticity. | `tpt-fem-assembly` |
| [`tpt-fem-solve`](crates/tpt-fem-solve) | Newton–Raphson and real arc-length (Crisfield) continuation. | `tpt-fem-assembly` |
| [`tpt-fem-eigen`](crates/tpt-fem-eigen) | Sparse shift-invert Lanczos/Arnoldi eigensolver, including the generalized `Kx = λMx` problem. | `tpt-fem-sparse` |
| [`tpt-fem-io-abaqus`](crates/tpt-fem-io-abaqus) | Abaqus `.inp` reader/writer. | `tpt-fem-mesh` |
| [`tpt-fem-io-exodus`](crates/tpt-fem-io-exodus) | Exodus II reader/writer (hand-rolled NetCDF-3 codec). | `tpt-fem-mesh` |
| [`tpt-fem-mesh-gen`](crates/tpt-fem-mesh-gen) | Native 3-D tetrahedral mesh generation (Delaunay + structured box), dependency-free. | `tpt-fem-mesh` |
| [`tpt-fem-dofmap`](crates/tpt-fem-dofmap) | Multi-field DOF map: per-field `dofs_per_node`, block/interleaved numbering. | `tpt-fem-mesh` |
| [`tpt-fem-dynamic`](crates/tpt-fem-dynamic) | Time integration: implicit Newmark-β/HHT-α and explicit central-difference over generic M/C/K. | `tpt-fem-{sparse,assembly}` |
| [`tpt-fem-plasticity`](crates/tpt-fem-plasticity) | J2 (von Mises) plasticity with isotropic/kinematic hardening (radial return). | `tpt-fem-{assembly,solve}` |
| [`tpt-fem-hyperelastic`](crates/tpt-fem-hyperelastic) | Neo-Hookean, Mooney-Rivlin, Ogden soft-tissue models. | `tpt-fem-{assembly,solve}` |
| [`tpt-fem-composite`](crates/tpt-fem-composite) | Classical lamination theory and cohesive-zone delamination. | `tpt-fem-elasticity` |
| [`tpt-fem-porous`](crates/tpt-fem-porous) | Biot consolidation and Darcy flow (steady + transient). | `tpt-fem-{assembly,dynamic}` |
| [`tpt-fem-contact`](crates/tpt-fem-contact) | Surface-to-surface contact (penalty + augmented Lagrangian), wear. | `tpt-fem-{assembly,dofmap}` |
| [`tpt-fem-fluid`](crates/tpt-fem-fluid) | Stokes and low-Re Navier-Stokes (mixed elements, transient). | `tpt-fem-{dofmap,dynamic,assembly,solve}` |
| [`tpt-fem-coupling`](crates/tpt-fem-coupling) | Multiphysics: thermal-structural, electro-thermal, fluid-structure. | `tpt-fem-{thermal,elasticity,fluid,dofmap,dynamic}` |
| [`tpt-fem-modal`](crates/tpt-fem-modal) | Modal analysis + frequency-response / modal superposition. | `tpt-fem-{eigen,sparse,dynamic}` |
| [`tpt-fem-topopt`](crates/tpt-fem-topopt) | SIMP topology optimization for 2-D linear elasticity. | `tpt-fem-{sparse,assembly,element,quadrature}` |
| [`tpt-fem-amr`](crates/tpt-fem-amr) | Adaptive h-refinement: 1-irregular quadtree Poisson with hanging-node elimination and ZZ error estimation. | `tpt-fem-sparse` |
| [`tpt-fem`](crates/tpt-fem) | Umbrella crate re-exporting all of the above behind Cargo features, plus `prelude` and end-to-end tests. | all of the above |
| [`tpt-fem-cli`](crates/tpt-fem-cli) | Command-line driver: `solve`, `elasticity`, `modal`, `amr`, `mesh info`, `mesh convert`. | `tpt-fem` |

Every crate is tracked as `git` in the sibling
[`tpt-rust-map/registry.toml`](https://github.com/tpt-solutions/tpt-rust-map).

## Umbrella feature flags

The `tpt-fem` umbrella re-exports each constituent behind a Cargo feature
(`quadrature`, `element`, `mesh`, `sparse`, `assembly`, `thermal`, `io-vtk`,
`elasticity`, `solve`, `eigen`, `io-abaqus`, `io-exodus`, `mesh-gen`). **All are
enabled by default.** The advanced Phase 12+ crates
(`dofmap`, `dynamic`, `plasticity`, `hyperelastic`, `composite`, `porous`,
`contact`, `fluid`, `coupling`, `modal`, `topopt`, `amr`) are opt-in features, not
enabled by default. `use tpt_fem::prelude::*;` pulls in the public API of every
enabled crate.

## Quick start

Add the umbrella (or individual crates) as a path dependency, then:

```rust
use tpt_fem::prelude::*;

// Build a structured 3-D tet mesh of the unit cube and solve -∇²u = 1.
let mesh = box_mesh([0.0; 3], [1.0; 3], [8, 8, 8]);
let mut bcs = Vec::new();
for axis in 0..3 {
    for &c in &[0.0, 1.0] {
        for n in mesh.nodes_on_plane(axis, c, 1e-9) {
            bcs.push((n, 0.0));
        }
    }
}
let u = solve_poisson(&mesh, 1.0, 2, |_| 1.0, &bcs, None, None).unwrap();
write_vtk_with_data(&mesh, &[PointData::new("u", u)], "out.vtk").unwrap();
```

A runnable version lives in
[`crates/tpt-fem/examples/thermal_solve.rs`](crates/tpt-fem/examples/thermal_solve.rs)
(`cargo run -p tpt-fem --example thermal_solve`).

See also:

* [`crates/tpt-fem/tests/patch_test.rs`](crates/tpt-fem/tests/patch_test.rs) —
  low-level element-loop patch test.
* [`crates/tpt-fem/tests/end_to_end.rs`](crates/tpt-fem/tests/end_to_end.rs) —
  real `mesh → solve → VTK` path.
* [`crates/tpt-fem-thermal/tests/mms_convergence.rs`](crates/tpt-fem-thermal/tests/mms_convergence.rs) —
  Method-of-Manufactured-Solutions convergence (L2/H1 rates).

## Getting started

Every adoption path begins with a clone — crates.io / PyPI publishing are out of
scope for this repository, so there is nothing to `cargo add`. From a checkout:

```bash
# 1. Clone the workspace.
git clone https://github.com/tpt-solutions/tpt-fem.git
cd tpt-fem

# 2. Run the end-to-end Poisson example (uses only the umbrella prelude).
cargo run -p tpt-fem --example thermal_solve
```

Expected output (values vary slightly with mesh count):

```text
Poisson solution on 864 tets: u in [0.000000, 0.013778]
Wrote thermal_solve.vtk
```

Open `thermal_solve.vtk` in ParaView to inspect the temperature field. From
here, explore the other examples under
[`crates/tpt-fem/examples/`](crates/tpt-fem/examples) (frame elasticity, modal
analysis, box mesh generation, Abaqus import), or drive a problem from a TOML
config with the CLI:

```bash
cargo run -p tpt-fem-cli -- solve examples/poisson.toml
```

## Command-line driver

```bash
# Scaffold a starter problem config for the chosen problem type.
cargo run -p tpt-fem-cli -- init poisson problem.toml

# Solve a Poisson problem from a TOML config (see crates/tpt-fem-cli for schema).
cargo run -p tpt-fem-cli -- solve problem.toml

# Solve an elasticity / modal problem from a TOML config.
cargo run -p tpt-fem-cli -- elasticity problem.toml
cargo run -p tpt-fem-cli -- modal problem.toml

# Adaptive h-refinement Poisson solve on [0,1]^2 (quadtree + ZZ estimator).
cargo run -p tpt-fem-cli -- amr --max-elements 512 --output amr.vtk

# Inspect a mesh (.msh / .vtk / .inp / .ex).
cargo run -p tpt-fem-cli -- mesh info mesh.msh

# Convert a Gmsh .msh mesh to a ParaView .vtk file.
cargo run -p tpt-fem-cli -- mesh convert mesh.msh mesh.vtk
```

A `problem.toml` looks like:

```toml
[problem]
type = "poisson"            # or "elasticity" / "modal" (see crates/tpt-fem-cli/examples)

[mesh]
dim = 2                    # 2 or 3; or set `file = "mesh.msh"` to import
min = [0.0, 0.0]
max = [1.0, 1.0]
n   = [20, 20]

[material]
conductivity = 1.0

[source]
constant = 1.0             # volumetric source f(x)

[[bc]]
value = 0.0
boundary = true            # selectors: nodes / plane / box / region / boundary

[output]
vtk = "solution.vtk"
```

## Building & testing

```bash
cargo build  --workspace --all-features
cargo test   --workspace --all-features
cargo fmt    --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo deny   check
```

A `Justfile` wraps these (and the `cargo run` examples) for convenience.

## Status

All phases are implemented: Phases 1–4 (core, assembly, first physics,
structural/nonlinear, ecosystem-gap crates), Phase 5 error ergonomics &
validation, Phase 6 physics completeness, and Phase 7–8 (CLI, MMS convergence
suite, fuzz targets, Python bindings).

## Examples

Run with the umbrella crate's re-exported prelude (`use tpt_fem::prelude::*;`):

| Example | Command | Description |
|---------|---------|-------------|
| `thermal_solve` | `cargo run -p tpt-fem --example thermal_solve` | Box-mesh Poisson solve and ParaView `.vtk` export. |
| `modal_analysis` | `cargo run -p tpt-fem --example modal_analysis` | Natural-vibration modes of a clamped plate. |
| `mesh_gen_box` | `cargo run -p tpt-fem --example mesh_gen_box` | Native tetrahedral box-mesh generation. |
| `elasticity_frame` | `cargo run -p tpt-fem --example elasticity_frame` | 2-D Euler–Bernoulli frame cantilever. |
| `abaqus_import` | `cargo run -p tpt-fem --example abaqus_import` | Round-trip an Abaqus `.inp` mesh import. |
| `quadrature_demo` | `cargo run -p tpt-fem --example quadrature_demo` | Quadrature weights = area and shape-function partition of unity. |
| `elasticity_bar_demo` | `cargo run -p tpt-fem --example elasticity_bar_demo` | Axial-bar element stiffness via the re-exported elasticity API. |

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at
your option.
