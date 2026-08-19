# tpt-fem-mesh

Mesh data structures, degree-of-freedom (DOF) numbering, and Gmsh `.msh`
v4.1 import for [tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the
mesh-based finite element core from
[tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate provides the in-memory mesh containers used by the rest of the
workspace:

- `Node` / `Element` / `Mesh` — the in-memory mesh model.
- `MeshBuilder` — a manual mesh-builder API for constructing meshes in code.
- `Mesh::number_dofs` — configurable (DOFs-per-node) degree-of-freedom
  numbering.
- `Mesh::from_msh_bytes` — import of Gmsh `.msh` version 4.1 ASCII files via
  the [`mshio`](https://crates.io/crates/mshio) crate.

Only the five linear element types `Line2`, `Tri3`, `Quad4`, `Tet4`, and
`Hex8` are imported; other Gmsh element types produce
`MeshError::UnsupportedElementType`. Higher-order element types are a tracked
follow-up.

## Installation

```toml
[dependencies]
tpt-fem-mesh = "0.1"
```

## Usage

```rust
use tpt_fem_mesh::{MeshBuilder, CellType};

let mut b = MeshBuilder::new();
let n0 = b.add_node(vec![0.0, 0.0]);
let n1 = b.add_node(vec![1.0, 0.0]);
let n2 = b.add_node(vec![0.0, 1.0]);
b.add_element(CellType::Tri, vec![n0, n1, n2]);
let mesh = b.build();

assert_eq!(mesh.node_count(), 3);
assert_eq!(mesh.element_count(), 1);

// Number one DOF per node and query the global system size.
mesh.number_dofs(1);
let ndof = mesh.dof_count();
```

## API highlights

| Item | Description |
|------|-------------|
| `CellType` | Linear element kinds (`Line`, `Tri`, `Quad`, `Tet`, `Hex`). |
| `Node` / `Element` / `Mesh` | Core mesh model. |
| `DofMap` | DOF-per-node numbering and lookups. |
| `MeshBuilder` | Programmatic mesh construction. |
| `Mesh::from_msh_bytes` / `Mesh::number_dofs` | Gmsh import and DOF numbering. |
| `MeshError` | Error type, including `UnsupportedElementType`. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-assembly / tpt-fem-thermal / tpt-fem-elasticity
     │
     └──► tpt-fem-io-* (readers/writers) and tpt-fem-mesh-gen
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `build_tri` | `cargo run -p tpt-fem-mesh --example build_tri` | Build a single `Tri3` mesh and assert its node/element counts and coordinates. |
| `quad4_grid` | `cargo run -p tpt-fem-mesh --example quad4_grid` | Build a 2×2 `Quad4` grid and assert the derived node/element counts. |
| `dof_numbering` | `cargo run -p tpt-fem-mesh --example dof_numbering` | Number DOFs per node and check the global system size and per-node lookups. |
| `node_selectors` | `cargo run -p tpt-fem-mesh --example node_selectors` | Select nodes on a plane and inside a box and check the results. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
