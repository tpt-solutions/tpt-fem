# tpt-fem-fluid

Stokes and low-Reynolds-number Navier–Stokes finite-element elements
(mixed velocity/pressure) for [tpt-fem](https://github.com/tpt-solutions/tpt-fem)
— the mesh-based finite element core from
[tpt-solutions](https://github.com/tpt-solutions).

## Overview

The crate assembles both steady and transient incompressible flow problems on an
unstructured mesh using a velocity/pressure mixed formulation with a
divergence-penalty (regularised incompressibility) stabilisation. That keeps the
element count to the velocity and pressure fields alone — no separate Lagrange
multiplier field is required.

- `steady_stokes` — Stokes (linearised, creeping) flow with a body force.
- `transient_stokes` — time-dependent Stokes relaxation via Newmark integration
  of the fluid-particle displacement.
- `transient_navier_stokes` — low-Re Navier–Stokes with a semi-implicit Picard
  linearisation of the convective term.
- `stokes_dofmap` — builds the mixed DOF layout (interleaved velocity then
  pressure) used by the solvers.

The penalty operator is integrated with a **selective reduced-integration**
(centroid) rule to avoid pressure locking, and the pressure is recovered with the
same reduced rule plus a lumped nodal projection, so the returned field is a
genuine nodal pressure `p = −(1/ε)·∇·u`.

## Installation

```toml
[dependencies]
tpt-fem-fluid = "0.1"
```

## Usage

```rust
use tpt_fem_fluid::steady_stokes;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

// A 2-D unit square triangulated into four cells, with no-slip walls.
let mut b = MeshBuilder::new();
let n = |x, y| b.add_node(vec![x, y]);
let a = n(0.0, 0.0);
let c = n(1.0, 0.0);
let d = n(1.0, 1.0);
let e = n(0.0, 1.0);
let m = n(0.5, 0.5);
b.add_element(CellType::Tri, vec![a, c, m]);
b.add_element(CellType::Tri, vec![c, d, m]);
b.add_element(CellType::Tri, vec![d, e, m]);
b.add_element(CellType::Tri, vec![e, a, m]);
let mesh = b.build();

// Drive the flow with a downward body force; pin every node to zero velocity.
let bc: Vec<(usize, f64)> = (0..mesh.node_count()).flat_map(|n| {
    let k = n * 2;
    [k, k + 1].into_iter().map(|d| (d, 0.0))
}).collect();
let body = |_x: &[f64]| vec![0.0, -1.0];

let (velocity, pressure) = steady_stokes(&mesh, 1.0, body, &bc, 1.0e6);
assert_eq!(velocity.len(), mesh.node_count() * 2);
assert_eq!(pressure.len(), mesh.node_count());
```

## Examples

Runnable examples live in `examples/` and are a good starting point:

- `stokes_dofmap` — inspects the mixed velocity/pressure DOF layout.
- `steady_stokes_poiseuille` — plane Poiseuille flow against the analytic
  parabolic profile.
- `lid_driven_cavity` — lid-driven cavity and a penalty-parameter sweep.
- `transient_stokes_startup` — transient Stokes relaxation to the steady state.
- `navier_stokes_low_re` — low-Re Poiseuille consistency and a cavity Re sweep.

Run one with:

```bash
cargo run -p tpt-fem-fluid --example steady_stokes_poiseuille
```

## API highlights

| Item | Description |
|------|-------------|
| `steady_stokes` | Solve steady Stokes; returns `(velocity, pressure)` with velocity indexed `node*dim + component`. |
| `transient_stokes` | Time-dependent Stokes relaxation; returns the `(time, velocity)` history. |
| `transient_navier_stokes` | Low-Re Navier–Stokes with Picard iteration; returns the final velocity. |
| `stokes_dofmap` | Mixed `MultiFieldDofMap` (interleaved velocity then pressure). |

## Position in the crate stack

```text
tpt-fem-mesh / -dofmap / -element / -quadrature
        │
        ▼
   tpt-fem-assembly ──► tpt-fem-fluid ──► tpt-fem-coupling
        │
        ▼
   tpt-fem-sparse / -dynamic / -solve
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
