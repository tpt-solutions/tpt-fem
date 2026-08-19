# tpt-fem-porous

Steady Darcy flow and transient Biot consolidation for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate provides saturated porous-media solvers built on a scalar
(pressure / excess-pore-pressure) field, reusing the isoparametric gradient
machinery of the displacement elements.

- [`solve_darcy`] assembles the conductivity matrix
  `K = ∫ k ∇Nᵀ∇N dΩ` for the steady-state equation
  `−∇·(k∇p) = f` (Darcy's law `q = −(k/μ)∇p`, with `μ` absorbed into the
  scalar permeability `k` passed here) and solves it with `tpt-fem-sparse`.
- [`terzaghi_consolidation`] integrates Terzaghi's 1-D consolidation equation
  `∂u/∂t = cᵥ ∂²u/∂z²` (a diffusion problem) with a stable backward-Euler step
  and returns the settlement history together with the maximum excess pore
  pressure at each step. The final drained settlement equals the closed form
  `q0·H / E_v`.

The coupling is the classical (uncoupled) Terzaghi form: the pore pressure
diffuses independently and the solid settles by `1/E_v` times the relaxed
pore pressure. Fully coupled Biot poroelasticity (solid–fluid matrix with
`α` and the drained/skeleton split) is a tracked follow-up and is **not**
implemented here.

Sign / unit conventions:

- `permeability` `k` is the scalar conductivity (length²/time, e.g. m²/s when
  Darcy velocity is in m/s); the surrounding fluid viscosity is absorbed.
- Excess pore pressure `u` and head `p` are reported in the same units as the
  `source` and Dirichlet values supplied by the caller (typically Pa).
- In [`terzaghi_consolidation`] the column spans `[0, H]`; the top
  (`z = H`, largest coordinate) is drained (`u = 0`) and the base (`z = 0`)
  is impermeable (natural zero-flux boundary).
- Settlement is positive downward and equals `(1/E_v) ∫ (q0 − u) dΩ`;
  `E_v` is the 1-D drained modulus and `cᵥ = k·E_v / γ_w` the consolidation
  coefficient (length²/time).
- The time step `dt` must satisfy the explicit-stability limit
  `dt ≤ Δz² / (2·cᵥ)`; the routine asserts this.

## Installation

```toml
[dependencies]
tpt-fem-porous = "0.1"
# The solvers assemble/evaluate on a mesh and solve via the sparse backend:
tpt-fem-mesh = "0.1"
tpt-fem-sparse = "0.1"
```

## Usage

```rust
use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_porous::solve_darcy;

// 1-D column of unit length, permeability k=1, drained at both ends (p=0),
// zero source: the solution is the trivial p = 0.
let mut b = MeshBuilder::new();
let a = b.add_node(vec![0.0]);
let m = b.add_node(vec![0.5]);
let c = b.add_node(vec![1.0]);
b.add_element(CellType::Line, vec![a, m]);
b.add_element(CellType::Line, vec![m, c]);
let mesh = b.build();
let src: Vec<Vec<f64>> = Vec::new();
let p = solve_darcy(&mesh, 1.0, &src, &[(a, 0.0), (c, 0.0)]).unwrap();
assert!(p.iter().all(|&x| x.abs() < 1e-9));
```

## Examples

- `darcy_steady` — steady 1-D Darcy seepage with fixed heads at both ends,
  checked against the linear analytic pressure profile `p(x) = 1 − x`.

  ```text
  cargo run -p tpt-fem-porous --example darcy_steady
  ```

- `terzaghi_consolidation` — transient 1-D consolidation compared against the
  closed-form average degree of consolidation `U(T_v)` at several dimensionless
  time factors `T_v`.

  ```text
  cargo run -p tpt-fem-porous --example terzaghi_consolidation
  ```

- `permeability_study` — parameter study confirming the pressure field is
  independent of `k` while the Darcy flux `q = −k∇p` scales linearly with `k`.

  ```text
  cargo run -p tpt-fem-porous --example permeability_study
  ```

## API highlights

| Item | Description |
|------|-------------|
| `solve_darcy` | Assembles `K = ∫ k ∇Nᵀ∇N dΩ` and solves `−∇·(k∇p) = f` with Dirichlet heads. |
| `terzaghi_consolidation` | Backward-Euler 1-D consolidation; returns `(t, settlement, max_u)` history. |
| `Mesh` / `CellType` (re-exported dependency) | Geometry the solvers operate on. |
| `tpt_fem_sparse::SparseError` | Error type returned by `solve_darcy`. |

## Position in the crate stack

```text
tpt-fem-quadrature ──► tpt-fem-element ──► tpt-fem-assembly ──► tpt-fem-porous
tpt-fem-mesh        ──► tpt-fem-porous
tpt-fem-sparse      ──► tpt-fem-porous
```

## References

- K. Terzaghi, *Theoretical Soil Mechanics*, Wiley, 1943 — 1-D consolidation.
- M. A. Biot, "General Theory of Three-Dimensional Consolidation", *J. Appl.
  Phys.*, 12(2), 1941 — coupled poroelasticity (the full coupling is a
  follow-up to this crate's uncoupled Terzaghi solver).

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
