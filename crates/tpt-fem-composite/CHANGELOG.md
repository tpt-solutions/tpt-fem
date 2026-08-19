# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `Ply` — orthotropic ply definition (moduli, Poisson ratio, thickness, fibre angle in degrees).
- `laminate_abd` — assembly of the `6×6` extensional–bending–coupling (`A`/`B`/`D`) matrix from a bottom-to-top ply stack (row-major ordering `[A B; B D]`).
- `CohesiveLaw` — bilinear traction–separation model for cohesive-zone delamination interfaces.
- `CohesiveLaw::from_toughness` — construction from fracture toughness `Gc`, peak traction and critical opening `δc`, deriving `δf` so the triangle area equals `Gc`.
- `CohesiveLaw::traction` — traction for an effective opening on the monotonic loading branch.
- `CohesiveLaw::toughness` — fracture toughness as the area under the traction–separation curve.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-composite-0.1.0
