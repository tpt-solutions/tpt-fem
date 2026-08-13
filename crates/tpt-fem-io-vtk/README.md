# tpt-fem-io-vtk

ParaView-compatible result export for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

Wraps [`vtkio`](https://crates.io/crates/vtkio) to write the linear element
types of `tpt-fem-mesh` (`Line2`, `Tri3`, `Quad4`, `Tet4`, `Hex8`) as an
unstructured-grid `.vtk` (legacy) or `.vtu` (XML) file, with optional
per-node scalar fields such as a computed temperature or displacement
magnitude.

## Installation

```toml
[dependencies]
tpt-fem-io-vtk = "0.1"
```

## Usage

```rust
use tpt_fem_io_vtk::PointData;
use tpt_fem_mesh::{CellType, MeshBuilder};

let mut b = MeshBuilder::new();
let n0 = b.add_node(vec![0.0, 0.0]);
let n1 = b.add_node(vec![1.0, 0.0]);
let n2 = b.add_node(vec![0.0, 1.0]);
b.add_element(CellType::Tri, vec![n0, n1, n2]);
let mesh = b.build();

let vtk = tpt_fem_io_vtk::mesh_to_vtk(&mesh, &[PointData::new("u", vec![0.0, 1.0, 1.0])]);
assert!(matches!(vtk.data, vtkio::model::DataSet::UnstructuredGrid { .. }));

// Write to disk (binary `.vtk` or ASCII via `write_vtk_ascii`).
tpt_fem_io_vtk::write_vtk(&mesh, "result.vtk").unwrap();
tpt_fem_io_vtk::write_vtk_with_data("result.vtu", &mesh, &[PointData::new("u", vec![0.0, 1.0, 1.0])]).unwrap();
```

## API highlights

| Item | Description |
|------|-------------|
| `PointData` | Named per-node scalar field (`new(name, values)`). |
| `mesh_to_vtk` | Build a `vtkio` dataset from a mesh + point data. |
| `write_vtk` / `write_vtk_ascii` | Write legacy `.vtk` (binary / ASCII). |
| `write_vtk_with_data` | Write `.vtk`/`.vtu` with point data in one call. |
| `VtkError` | Error type for export failures. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-io-vtk ──► ParaView / VTK tooling
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
