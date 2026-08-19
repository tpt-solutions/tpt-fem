# tpt-fem-mesh-gen

Native (dependency-free) 3D tetrahedral mesh generation for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate provides dependency-free mesh generation, written from scratch so
it can ship under the repo's `MIT OR Apache-2.0` policy where wrapped
alternatives (e.g. TetGen-based tooling) are AGPL and therefore disallowed.

Two generators are provided:

- `delaunay_3d` — an incremental (Bowyer–Watson) Delaunay tetrahedralisation
  of an arbitrary point cloud. It meshes the convex hull of the input points
  and is useful for filling a volume bounded by a set of seed nodes.
- `box_mesh` — a structured mesher that splits every axis-aligned brick of a
  bounding box into six tetrahedra. This path is *guaranteed* to produce a
  valid, intersection-free, positively-oriented mesh with no external
  dependency and no robustness caveats, so it is the recommended route when a
  quality grid is needed.

Both return a `tpt_fem_mesh::Mesh` of `CellType::Tet` elements. Quality can be
inspected with `tet_quality` and improved with `laplacian_smooth`.

### Robustness note

The Delaunay predicates (`orient3`, `in_sphere`) are evaluated in `f64` with a
relative tolerance. They are exact for points in general position (no four
coplanar, no five cospherical) and degrade gracefully otherwise; coincident
input points are de-duplicated automatically. For highest quality on arbitrary
domains, the `box_mesh` + `laplacian_smooth` route is preferred.

## Installation

```toml
[dependencies]
tpt-fem-mesh-gen = "0.1"
```

## Usage

```rust
use tpt_fem_mesh_gen::{box_mesh, delaunay_3d, laplacian_smooth, tet_quality, Point3};

// Structured box mesh: split each brick of the [min, max] box into 6 tets.
let mesh = box_mesh(
    Point3::new(0.0, 0.0, 0.0),
    Point3::new(1.0, 1.0, 1.0),
    [4, 4, 4],
);

// Improve element quality with a few Laplacian-smoothing passes.
let worst = laplacian_smooth(&mut mesh, 10);
let quality = tet_quality(&mesh);
```

## API highlights

| Item | Description |
|------|-------------|
| `delaunay_3d` | Bowyer–Watson Delaunay tetrahedralisation of a point cloud. |
| `box_mesh` | Robust structured box → tetrahedra mesher. |
| `orient3` / `in_sphere` | Exact-in-`f64` Delaunay predicates. |
| `all_positively_oriented` | Orientation sanity check. |
| `tet_quality` / `TetQuality` | Minimum/mean aspect-ratio metrics. |
| `laplacian_smooth` | Iterative node smoothing; returns worst quality. |

## Position in the crate stack

```text
tpt-fem-mesh ◄── tpt-fem-mesh-gen (produces Tet meshes)
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `box_mesh_demo` | `cargo run -p tpt-fem-mesh-gen --example box_mesh_demo` | Structured box mesh; asserts node/element counts and positive orientation. |
| `delaunay_cloud` | `cargo run -p tpt-fem-mesh-gen --example delaunay_cloud` | Delaunay tetrahedralisation of a cube corner cloud; checks counts and orientation. |
| `tet_quality_demo` | `cargo run -p tpt-fem-mesh-gen --example tet_quality_demo` | Reports tetrahedral quality metrics and runs Laplacian smoothing. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
