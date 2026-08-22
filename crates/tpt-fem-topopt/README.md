# tpt-fem-topopt

SIMP topology optimization for 2-D linear elasticity on the `tpt-fem` stack.

This crate provides a self-contained minimum-compliance (stiffness
maximization) optimizer using the **Solid Isotropic Material with
Penalization (SIMP)** model, a density (sensitivity) filter, and the
**optimality-criteria (OC)** update. It sits on top of the existing `tpt-fem`
crates:

- `tpt-fem-element` — `Quad4` bilinear shape functions and reference-coordinate
  derivatives.
- `tpt-fem-quadrature` — Gauss-Legendre quadrature.
- `tpt-fem-sparse` — `Coo` triplet accumulation.
- `tpt-fem-assembly` — `solve_with_dirichlet` for the per-iteration
  essential-boundary-condition solve.

The optimizer minimizes the structural compliance `c = uᵀ K(ρ) u` subject to a
volume fraction `vol_frac`, where the Young's modulus of each element is
interpolated as `E(ρ) = E_min + ρᵖ (E₀ − E_min)` (`ρ ∈ [0, 1]`, `p ≥ 1` the
SIMP penalty). This is the classic Bendsøe/Sigmund formulation.

## Example

```rust
use tpt_fem_topopt::{cantilever_load, topopt_simp, TopOptParams, Grid};

let grid = Grid::new(20, 10, 1.0);
let (f, bcs) = cantilever_load(&grid, 1.0);
let params = TopOptParams {
    grid: grid.clone(),
    e0: 1.0,
    nu: 0.3,
    vol_frac: 0.5,
    penal: 3.0,
    filter_radius: 2.0,
    max_iter: 40,
    move_limit: 0.2,
};
let res = topopt_simp(&params, &f, &bcs).unwrap();
// The optimized design is lighter (lower compliance) than the uniform start,
// while still consuming exactly `vol_frac` of the domain.
assert!(res.compliance.last().unwrap() < &res.compliance[0]);
```

License: MIT OR Apache-2.0.
