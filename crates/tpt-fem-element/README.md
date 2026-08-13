# tpt-fem-element

Reference elements, Lagrange shape functions, and the isoparametric Jacobian
mapping for [tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the
mesh-based finite element method core from
[tpt-solutions](https://github.com/tpt-solutions).

This crate depends only on `tpt-fem-quadrature` and sits directly above it in
the dependency stack.

## Overview

Five linear (`P1`) Lagrange reference elements are provided:

| Element | Spatial dim | Reference domain                              |
|---------|-------------|----------------------------------------------|
| `Line2` | 1           | `[-1, 1]`                                    |
| `Tri3`  | 2           | `(0,0),(1,0),(0,1)`                          |
| `Quad4` | 2           | `[-1, 1]²`                                   |
| `Tet4`  | 3           | `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`            |
| `Hex8`  | 3           | `[-1, 1]³`                                   |

Each element exposes its reference-node coordinates, its shape-function values
`Nᵢ(ξ)`, and the shape-function gradients `∂Nᵢ/∂ξⱼ` with respect to the
reference coordinates. The `Map` type assembles the isoparametric Jacobian
from physical node coordinates and the local gradients, and maps local
gradients to their physical-space counterparts.

> Quadratic (`P2`) elements are a tracked follow-up and are intentionally not
> implemented in this pass.

## Installation

```toml
[dependencies]
tpt-fem-element = "0.1"
```

## Usage

```rust
use tpt_fem_element::{Tri3, ReferenceElement, Map};

// Evaluate the shape functions and their reference gradients at a point.
let xi = [0.3, 0.3];
let n = Tri3::shape(&xi);
let dn = Tri3::grad(&xi);

// Build the isoparametric Jacobian for a physical triangle.
let nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
let map = Map::new::<Tri3>(&nodes);
let (jacobian, det, inv) = map.at(&xi);
```

Convenience constructors return the matching quadrature rule for an element:

| Function | Returns |
|----------|---------|
| `line_rule(order)` | `tpt_fem_quadrature::Quad1D` |
| `quad_rule(order)` | `tpt_fem_quadrature::Quad2D` |
| `hex_rule(order)` | `tpt_fem_quadrature::Quad3D` |
| `tri_rule(rule)` | `tpt_fem_quadrature::Quad2D` |
| `tet_rule(rule)` | `tpt_fem_quadrature::Quad3D` |

## API highlights

| Item | Description |
|------|-------------|
| `ReferenceElement` (trait) | Uniform shape/derivative interface for all elements. |
| `Line2` / `Tri3` / `Quad4` / `Tet4` / `Hex8` | Zero-sized reference-element types. |
| `Map` | Isoparametric Jacobian assembly and local→physical gradient mapping. |
| `line_rule` / `quad_rule` / `hex_rule` / `tri_rule` / `tet_rule` | Element-to-quadrature helpers. |

## Position in the crate stack

```text
tpt-fem-quadrature ──► tpt-fem-element ──► tpt-fem-thermal / tpt-fem-elasticity
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
