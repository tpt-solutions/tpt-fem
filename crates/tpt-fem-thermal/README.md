# tpt-fem-thermal

Heat-conduction and Poisson finite-element formulations for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

The steady-state scalar problem solved here is

```text
-∇·(k ∇u) = f   in Ω,
```

equipped with Dirichlet (essential), Neumann (outward-flux), and Robin
(convective) boundary conditions. The element stiffness and load vectors are
integrated with the reference-element quadrature from `tpt-fem-quadrature` and
the isoparametric Jacobian from `tpt-fem-element`, then scattered into a global
system and solved by `tpt-fem-assembly` + `tpt-fem-sparse`.

Scalar fields (one degree of freedom per node) are assumed.

## Installation

```toml
[dependencies]
tpt-fem-thermal = "0.1"
```

## Usage

```rust
use tpt_fem_thermal::solve_poisson;

// Solve -∇·(k ∇u) = f on `mesh` with constant conductivity `k` and source `f`,
// applying Dirichlet conditions (dof, value). `S` selects the sparse backend.
let solution = solve_poisson::<S>(mesh, k, f, &[(boundary_dof, value)]);
```

Lower-level building blocks are also exposed:

```rust
use tpt_fem_thermal::{poisson_element_matrix, poisson_source_vector};

let ke = poisson_element_matrix::<Tri3>(&map, &quad, conductivity);
let fe = poisson_source_vector::<Tri3>(&map, &quad, source);
```

## API highlights

| Item | Description |
|------|-------------|
| `poisson_element_matrix` | Element stiffness `K_e = ∫ k ∇Nᵀ∇N dΩ`. |
| `poisson_source_vector` | Element load vector for source `f`. |
| `solve_poisson` | End-to-end Poisson/heat-conduction solve. |

Convergence of the P1 discretisation is verified by the method-of-manufactured-
solutions (MMS) integration tests (`tests/mms_convergence.rs`), asserting `L2`
order 2 and `H1` order 1 rates.

## Position in the crate stack

```text
tpt-fem-quadrature + tpt-fem-element ──► tpt-fem-thermal ──► tpt-fem-assembly ──► tpt-fem-sparse
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `tri3_conductivity` | `cargo run -p tpt-fem-thermal --example tri3_conductivity` | Single Tri3 conductivity matrix vs. the hand-computed `k·A·(∇Nᵢ·∇Nⱼ)`. |
| `steady_1d` | `cargo run -p tpt-fem-thermal --example steady_1d` | 1-D `-u'' = 1` Poisson solve converging to the analytic midpoint `0.125`. |
| `quad4_stiffness` | `cargo run -p tpt-fem-thermal --example quad4_stiffness` | Quad4 conductivity matrix symmetry and rigid-body (zero row-sum) check. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
