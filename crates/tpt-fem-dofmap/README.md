# tpt-fem-dofmap

Multi-field degree-of-freedom numbering for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

The rest of the core (`tpt-fem-mesh::Mesh::number_dofs`) hard-codes a single
*uniform* `dofs_per_node` across the whole mesh — fine for a single physical
field (temperature, or one displacement vector) but insufficient for coupled
problems where different fields carry different component counts. Examples
include displacement with 3 components and pressure with 1, or a fluid
velocity with 2 components alongside a scalar temperature.

[`MultiFieldDofMap`] generalises numbering to N named fields, each with its own
per-node component count, laid out on a single [`tpt-fem-mesh::Mesh`]:

```text
node n, field f, component c  ──►  single global DOF index
```

The layout may be either field-major (`Layout::Block`, all DOFs of field 0 for
every node first) or node-major (`Layout::Interleaved`, every field's
component block contiguous per node). With a single field the map degenerates
exactly to `Mesh::number_dofs`.

> Out of scope (intentionally not implemented): the map is purely topological —
> it assigns indices but performs no assembly, element traversal, or constraint
> (e.g. tie/periodic) handling, and it does not reorder for bandwidth minimisation.

## Installation

```toml
[dependencies]
tpt-fem-dofmap = "0.1"
# A caller typically also needs the mesh the map is built over:
tpt-fem-mesh = "0.1"
```

## Usage

```rust
use tpt_fem_dofmap::{FieldSpec, Layout, MultiFieldDofMap};
use tpt_fem_mesh::{CellType, MeshBuilder};

let mut b = MeshBuilder::new();
let a = b.add_node(vec![0.0, 0.0]);
let c = b.add_node(vec![1.0, 0.0]);
b.add_element(CellType::Line, vec![a, c]);
let mesh = b.build();

// One displacement field (2 components) + one pressure field (1 component).
let map = MultiFieldDofMap::new(
    &mesh,
    &[FieldSpec::new("u", 2), FieldSpec::new("p", 1)],
    Layout::Interleaved,
);
// Node 0 owns dofs 0,1 (u) and 2 (p); node 1 owns 3,4 and 5.
assert_eq!(map.node_field_dof(0, 0, 1), 1);
assert_eq!(map.node_field_dof(0, 1, 0), 2);
assert_eq!(map.node_field_dof(1, 1, 0), 5);
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `single_field` | `cargo run -p tpt-fem-dofmap --example single_field` | Shows a single-field map degenerating to `Mesh::number_dofs`. |
| `two_field` | `cargo run -p tpt-fem-dofmap --example two_field` | Builds a 3-component velocity + 1-component pressure map and prints its global DOF layout. |
| `layout_compare` | `cargo run -p tpt-fem-dofmap --example layout_compare` | Compares `Layout::Block` vs `Layout::Interleaved` index tables side by side. |

## API highlights

| Item | Description |
|------|-------------|
| `FieldSpec` / `FieldSpec::new` | A named field and its per-node component count (`components` must be non-zero). |
| `Layout` | `Block` (field-major) or `Interleaved` (node-major) global numbering. |
| `MultiFieldDofMap` | Multi-field DOF map over a mesh. |
| `MultiFieldDofMap::new` | Build the map for a mesh, fields, and layout. |
| `node_field_dof` | Global DOF for `(node, field, component)`. |
| `components` | Per-node component count of a field. |
| `field_range` | `(start, count)` of a field's `Block`-style contiguous DOF block. |
| `dofs_of` | All global DOFs owned by a node. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-dofmap
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
