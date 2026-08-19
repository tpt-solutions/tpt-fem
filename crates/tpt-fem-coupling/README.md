# tpt-fem-coupling

Multiphysics coupling operators (thermal-structural, electro-thermal, FSI) for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate wires the single-physics solvers of `tpt-fem` into coupled
multiphysics problems. Each operator takes the relevant sub-domain meshes and
material data and returns the response field of the driven sub-problem:

- `thermal_structural` — free thermal expansion / thermal stress. A temperature
  field becomes an initial strain `ε_th = α·ΔT·I`, solved with the supplied
  elastic model (bar-axial, plane-stress, or plane-strain).
- `joule_source` — the volumetric Ohmic dissipation `q = σ·|E|²` (W/m³) from a
  conductivity and an applied field magnitude.
- `electro_thermal` — steady heat conduction `−k·∇²T = q` driven by that
  constant source (a Poisson solve from `tpt-fem-thermal`).
- `fsi_coupling` — one explicit fluid-structure substep: displace the fluid
  interface by the current structure motion, solve steady Stokes, transfer the
  recovered nodal pressure as a normal traction onto the structure, and return
  the updated structural displacement.

## Installation

```toml
[dependencies]
tpt-fem-coupling = "0.1"
```

## Usage

```rust
use tpt_fem_coupling::{thermal_structural, electro_thermal, joule_source};
use tpt_fem_elasticity::ElasticModel;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

// A 1-D bar with a uniform temperature rise dT -- stress-free, so it just grows.
let mut b = MeshBuilder::new();
let mut prev = b.add_node(vec![0.0]);
for i in 1..=4 {
    let n = b.add_node(vec![i as f64 / 4.0]);
    b.add_element(CellType::Line, vec![prev, n]);
    prev = n;
}
let mesh = b.build();
let temp = vec![1.0; mesh.node_count()];
let dirichlet = [(0usize, 0.0)]; // pin the left end
let u = thermal_structural(
    &mesh, ElasticModel::BarAxial, 1.0, 0.3, 1e-3, &temp, &dirichlet,
).unwrap();
// tip extension == alpha * dT * L
assert!((u[4] - 1e-3).abs() < 1e-6);

// Joule heating of a bar held at T = 0 at both ends.
let sigma = 5.8e7; // copper, S/m
let e_mag = 0.05;  // V/m
let q = joule_source(sigma, e_mag); // == sigma * |E|^2
let bc = [(0usize, 0.0), (mesh.node_count() - 1, 0.0)];
let _t = electro_thermal(&mesh, 400.0, sigma, e_mag, &bc).unwrap();
```

## Examples

Runnable examples live in `examples/`:

- `thermal_bar_expansion` — stress-free growth, `ΔL = α·ΔT·L`, swept over `ΔT`
  and length.
- `thermal_bimetal_strip` — a through-thickness gradient bends a clamped strip,
  compared with beam-theory curvature.
- `joule_heating` — `q = σ·|E|²` fed into conduction; `ΔT` matches the
  `q·x·(L−x)/(2k)` parabola and scales as `|E|²`.
- `fsi_coupling` — an elastic block under a fluid layer, relaxed by repeated
  `fsi_coupling` substeps until the coupling converges.

Run one with:

```bash
cargo run -p tpt-fem-coupling --example thermal_bar_expansion
```

## API highlights

| Item | Description |
|------|-------------|
| `thermal_structural` | Thermal expansion / stress from a temperature field and CTE. |
| `joule_source` | Ohmic volumetric heat `σ·|E|²`. |
| `electro_thermal` | Steady conduction `−k∇²T = σ|E|²`. |
| `fsi_coupling` | One explicit fluid→structure substep, returning the new structural displacement. |

## Position in the crate stack

```text
tpt-fem-thermal / -elasticity / -fluid
        │
        ▼
   tpt-fem-coupling
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
