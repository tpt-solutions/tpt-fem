# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `delaunay_3d` incremental (Bowyer–Watson) Delaunay tetrahedralisation.
- `box_mesh` robust structured box → tetrahedra mesher.
- `orient3` / `in_sphere` exact-in-`f64` Delaunay predicates with
  coincident-point de-duplication.
- `all_positively_oriented` orientation check.
- `tet_quality` / `TetQuality` aspect-ratio metrics.
- `laplacian_smooth` iterative node smoothing.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-mesh-gen-0.1.0
