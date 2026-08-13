# tpt-fem-quadrature

Gauss quadrature rules for finite-element reference elements, part of the
[`tpt-fem`](https://github.com/tpt-solutions/tpt-fem) workspace — a mesh-based
finite element method (FEM) core for [tpt-solutions](https://github.com/tpt-solutions).

This crate has **no dependencies** and is the lowest-level building block of the
stack: it provides the integration rules consumed by `tpt-fem-element`,
`tpt-fem-thermal`, `tpt-fem-elasticity`, and `tpt-fem-assembly`.

## Overview

Numerical integration of a reference-element weak form requires quadrature
points `xᵢ` and weights `wᵢ` such that

```text
∫_Ω f(x) dΩ  ≈  Σ_i wᵢ f(xᵢ)
```

with the rule exact for polynomials up to a stated degree. This crate provides
exact fixed-order rules on every reference domain used by `tpt-fem`:

- **1-D Gauss–Legendre** of orders 1–5 on `[-1, 1]` and on `[0, 1]`.
- **Tensor-product** rules on the square `[-1, 1]²` and cube `[-1, 1]³`.
- **Fixed low-order rules** on the reference triangle
  `(0,0),(1,0),(0,1)`.
- **Fixed low-order rules** on the reference tetrahedron
  `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`.

Every rule is verified by the unit tests, which integrate monomials against
their closed-form values.

## Installation

```toml
[dependencies]
tpt-fem-quadrature = "0.1"
```

## Usage

```rust
use tpt_fem_quadrature::gauss_legendre;

// Order-2 Gauss-Legendre is exact for cubics on [-1, 1].
let rule = gauss_legendre(2);
let approx: f64 = rule.weights.iter().zip(&rule.points)
    .map(|(w, x)| w * x * x * x)
    .sum();
assert!((approx - 0.0).abs() < 1e-12);
```

## API highlights

| Item | Description |
|------|-------------|
| `Quad1D` / `Quad2D` / `Quad3D` | Point/weight containers for 1-, 2-, and 3-D rules. |
| `gauss_legendre(order)` | Gauss–Legendre rule of order 1–5 on `[-1, 1]`. |
| `gauss_legendre_unit(order)` | Same, but on `[0, 1]`. |
| `tensor_square(rule)` / `tensor_cube(rule)` | Tensor-product 2-D / 3-D rules. |
| `TriangleRule` / `triangle(rule)` | Rule selection + rule on the reference triangle. |
| `TetrahedronRule` / `tetrahedron(rule)` | Rule selection + rule on the reference tetrahedron. |

## Position in the crate stack

```text
tpt-fem-quadrature  ──►  tpt-fem-element
       │
       └──►  tpt-fem-thermal / tpt-fem-elasticity  ──►  tpt-fem-assembly
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
