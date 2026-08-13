# tpt-fem-elasticity

Linear-elasticity element formulations for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

Supports:

- `Line2` — 1-D axial bar (`EA` stiffness),
- `Tri3` / `Quad4` — 2-D plane-stress and plane-strain continua,
- `Tet4` / `Hex8` — 3-D isotropic continua.

Each element stiffness is `K_e = ∫_Ω Bᵀ D B dΩ`, with `B` the
strain–displacement matrix and `D` the isotropic constitutive matrix. The
per-element matrices are scattered by `tpt-fem-assembly` and solved with
`tpt-fem-sparse`. Modal (eigenvalue) analysis reuses `tpt-fem-eigen`.

The crate also provides 2-D Euler–Bernoulli beam elements (`BeamSection2D`,
`beam2d_element_matrix`, `beam2d_consistent_mass`, `solve_frame2d`).

## Installation

```toml
[dependencies]
tpt-fem-elasticity = "0.1"
```

## Usage

```rust
use tpt_fem_elasticity::{solve_elasticity, ElasticModel};

// Assemble and solve a linear elasticity problem for an isotropic continuum.
let solution = solve_elasticity::<Tet4, _>(
    mesh,
    ElasticModel::ThreeD,   // or PlaneStress / PlaneStrain / AxialBar
    youngs_modulus,
    poissons_ratio,
    &body_force,
    &dirichlet,
);
```

Modal analysis:

```rust
use tpt_fem_elasticity::solve_modal;

// Recover a few eigenpairs of the stiffness/mass generalised eigenproblem.
let modes = solve_modal::<Tet4, _>(mesh, elastic_model, e, nu, n_modes);
```

## API highlights

| Item | Description |
|------|-------------|
| `ElasticModel` | `AxialBar`, `PlaneStress`, `PlaneStrain`, `ThreeD`. |
| `elasticity_element_matrix` | Element stiffness `K_e = ∫ Bᵀ D B dΩ`. |
| `elasticity_body_vector` | Element body-force load vector. |
| `elasticity_mass_matrix` / `elasticity_lumped_mass` | Consistent / lumped mass. |
| `solve_elasticity` | End-to-end static solve. |
| `solve_modal` | Modal eigenproblem via `tpt-fem-eigen`. |
| `BeamSection2D`, `beam2d_element_matrix`, `beam2d_consistent_mass`, `solve_frame2d` | 2-D frame/beam analysis. |

## Position in the crate stack

```text
tpt-fem-assembly + tpt-fem-sparse + tpt-fem-eigen ──► tpt-fem-elasticity
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
