# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- Gauss–Legendre quadrature rules of orders 1–5 on `[-1, 1]` and `[0, 1]`
  (`gauss_legendre`, `gauss_legendre_unit`).
- Tensor-product quadrature on the reference square and cube
  (`tensor_square`, `tensor_cube`).
- Fixed low-order quadrature rules on the reference triangle
  (`triangle`) and tetrahedron (`tetrahedron`).
- `Quad1D` / `Quad2D` / `Quad3D` point-and-weight containers.
- Unit tests verifying polynomial-exactness against closed-form monomial
  integrals.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-quadrature-0.1.0
