# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-10-06

### Fixed

- **`TetrahedronRule::Keast4`'s 11-point table was wrong**: several of its
  points had barycentric coordinates summing to more than 1 (up to `1.2607`),
  placing them outside the reference tetrahedron. The rule still integrated
  monomials exactly (moment-matching is a purely algebraic property that
  doesn't require the points to lie anywhere in particular), so the bug was
  silent — it was reported by a downstream consumer after a `Tet10` solve
  using this rule came out 42x the closed-form answer. The table has been
  re-derived from Burkardt's public-domain reference implementation of
  Keast (1986) and is now covered by a regression test
  (`every_rule_stays_inside_its_reference_element`) that checks every rule's
  points against its reference element's domain, not just its moments.

### Added

- `QuadratureError` and fallible `try_gauss_legendre` / `try_gauss_legendre_unit`
  constructors, returning an error instead of panicking when `order` is
  outside the supported `1..=5` range. Prefer these over `gauss_legendre` /
  `gauss_legendre_unit` whenever `order` isn't a fixed, known-valid literal.

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

[0.1.1]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-quadrature-0.1.1
[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-quadrature-0.1.0
