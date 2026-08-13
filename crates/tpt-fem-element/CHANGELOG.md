# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `ReferenceElement` trait unifying the five linear (`P1`) Lagrange elements.
- Reference-element types `Line2`, `Tri3`, `Quad4`, `Tet4`, `Hex8` with shape
  functions and reference gradients.
- `Map` for isoparametric Jacobian assembly and local→physical gradient
  mapping.
- Element-to-quadrature helpers `line_rule`, `quad_rule`, `hex_rule`,
  `tri_rule`, `tet_rule`.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-element-0.1.0
