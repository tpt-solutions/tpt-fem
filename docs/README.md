# tpt-fem Documentation

This folder is the single browsable entry point for the `tpt-fem` workspace.
The workspace is a layered finite-element-method core; each crate has its own
`README.md` (linked below), and this page ties them together by capability.

> Crates.io / PyPI publishing are intentionally **out of scope** for this
> repository — every crate is consumed as a path dependency within the
> workspace. Start from a clone (see the root `README.md` "Getting started").

## Capability map

| Capability | Crate(s) | Example |
|------------|----------|---------|
| Quadrature rules (line / tri / quad / tet / hex) | `tpt-fem-quadrature` | — |
| Reference elements, P1/P2 shape functions, isoparametric Jacobian | `tpt-fem-element` | — |
| Mesh model, Gmsh import, region selectors, validation | `tpt-fem-mesh` | `abaqus_import` |
| COO/CSR assembly, in-house dense-LU solve (+ optional `russell`) | `tpt-fem-sparse` | — |
| Assembly, Dirichlet/Neumann/Robin BCs, reduced-system extraction | `tpt-fem-assembly` | — |
| Steady heat conduction / Poisson | `tpt-fem-thermal` | `thermal_solve` |
| Linear elasticity & 2-D frames | `tpt-fem-elasticity` | `elasticity_frame` |
| Newton–Raphson + arc-length continuation | `tpt-fem-solve` | — |
| Generalized eigenproblem `K φ = ω² M φ` | `tpt-fem-eigen` | `modal_analysis` |
| ParaView VTK export (+ VTK import) | `tpt-fem-io-vtk` | all examples |
| Abaqus `.inp` reader/writer | `tpt-fem-io-abaqus` | `abaqus_import` |
| Exodus II reader/writer (NetCDF-3) | `tpt-fem-io-exodus` | — |
| Native 3-D tet mesh generation | `tpt-fem-mesh-gen` | `mesh_gen_box` |
| DOF maps for multi-field problems | `tpt-fem-dofmap` | all examples |
| Transient dynamics (Newmark) | `tpt-fem-dynamic` | all examples |
| J2 (von Mises) plasticity, return mapping | `tpt-fem-plasticity` | `uniaxial_sweep` |
| Hyperelastic models (Neo-Hookean, etc.) | `tpt-fem-hyperelastic` | all examples |
| Composite / layered materials | `tpt-fem-composite` | all examples |
| Poroelasticity / coupled flow–mechanics | `tpt-fem-porous` | all examples |
| Contact & friction interface conditions | `tpt-fem-contact` | all examples |
| Stokes / Navier–Stokes (penalty method) | `tpt-fem-fluid` | — |
| Multi-physics coupling driver | `tpt-fem-coupling` | all examples |
| Umbrella re-exports + prelude | `tpt-fem` | all examples |
| Command-line driver | `tpt-fem-cli` | — |
| Python bindings (maturin/pyo3) | `tpt-fem-py` | `walkthrough.py` |

## Per-crate READMEs

- [tpt-fem-quadrature](../crates/tpt-fem-quadrature/README.md)
- [tpt-fem-element](../crates/tpt-fem-element/README.md)
- [tpt-fem-mesh](../crates/tpt-fem-mesh/README.md)
- [tpt-fem-sparse](../crates/tpt-fem-sparse/README.md)
- [tpt-fem-assembly](../crates/tpt-fem-assembly/README.md)
- [tpt-fem-thermal](../crates/tpt-fem-thermal/README.md)
- [tpt-fem-io-vtk](../crates/tpt-fem-io-vtk/README.md)
- [tpt-fem-elasticity](../crates/tpt-fem-elasticity/README.md)
- [tpt-fem-solve](../crates/tpt-fem-solve/README.md)
- [tpt-fem-eigen](../crates/tpt-fem-eigen/README.md)
- [tpt-fem-io-abaqus](../crates/tpt-fem-io-abaqus/README.md)
- [tpt-fem-io-exodus](../crates/tpt-fem-io-exodus/README.md)
- [tpt-fem-mesh-gen](../crates/tpt-fem-mesh-gen/README.md)
- [tpt-fem-dofmap](../crates/tpt-fem-dofmap/README.md)
- [tpt-fem-dynamic](../crates/tpt-fem-dynamic/README.md)
- [tpt-fem-plasticity](../crates/tpt-fem-plasticity/README.md)
- [tpt-fem-hyperelastic](../crates/tpt-fem-hyperelastic/README.md)
- [tpt-fem-composite](../crates/tpt-fem-composite/README.md)
- [tpt-fem-porous](../crates/tpt-fem-porous/README.md)
- [tpt-fem-contact](../crates/tpt-fem-contact/README.md)
- [tpt-fem-fluid](../crates/tpt-fem-fluid/README.md)
- [tpt-fem-coupling](../crates/tpt-fem-coupling/README.md)
- [tpt-fem (umbrella)](../crates/tpt-fem/README.md)
- [tpt-fem-cli](../crates/tpt-fem-cli/README.md)
- [tpt-fem-py](../crates/tpt-fem-py/README.md)

## End-to-end examples

Runnable Rust examples under [`crates/tpt-fem/examples/`](../crates/tpt-fem/examples):

- `thermal_solve.rs` — Poisson on a unit cube → VTK.
- `elasticity_frame.rs` — 2-D Euler–Bernoulli cantilever → VTK.
- `modal_analysis.rs` — natural frequencies of a 3-D cantilever block → VTK.
- `mesh_gen_box.rs` — structured box + Delaunay tetrahedralisation → VTK.
- `abaqus_import.rs` — parse an Abaqus `.inp` deck and report metadata.

Python walkthrough: [`crates/tpt-fem-py/examples/walkthrough.py`](../crates/tpt-fem-py/examples/walkthrough.py).

## Developer docs

- [CONTRIBUTING.md](../CONTRIBUTING.md) — `Justfile` recipes, the crate
  dependency DAG, and the `todo.md` phase workflow.
- [todo.md](../todo.md) — the phased build/audit backlog (source of truth for
  what is implemented and what is planned).
