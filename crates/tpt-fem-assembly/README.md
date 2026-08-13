# tpt-fem-assembly

Element-to-global assembly and boundary-condition application for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate ties the core crates together into the standard FEM workflow:

1. `assemble` scatters per-element matrices (computed by a physics crate such
   as `tpt-fem-thermal`) into a global `Coo` triplet matrix.
2. Natural boundary conditions are added through `apply_neumann` (load vector)
   and `apply_robin` (stiffness + load contributions).
3. Essential conditions are enforced by `solve_with_dirichlet`, which
   statically condenses the fixed degrees of freedom, solves the reduced system
   with `tpt-fem-sparse`, and scatters the solution back.

All routines are dimension-agnostic across the five linear element types
`Line2`, `Tri3`, `Quad4`, `Tet4`, and `Hex8`, and support an arbitrary
(uniform) number of degrees of freedom per node.

## Installation

```toml
[dependencies]
tpt-fem-assembly = "0.1"
```

## Usage

```rust
use tpt_fem_assembly::{assemble, solve_with_dirichlet, ReducedSystem};

// `element_matrices` is a Vec<(element_id, local_matrix)> from a physics crate.
let coo = assemble(&mesh, &element_matrices);

// Apply Dirichlet conditions (dof, value) and solve the reduced system.
let ReducedSystem { solution, .. } =
    solve_with_dirichlet(&coo, &rhs, &[(0, 0.0), (last_dof, 1.0)]);
```

Natural boundary conditions:

```rust
use tpt_fem_assembly::{apply_neumann, apply_robin};

// Add an outward-flux Neumann load on a face.
apply_neumann(&mut coo, &mut rhs, mesh, element_id, face_index, flux);

// Add a convective Robin condition with coefficient `h` and ambient `u_inf`.
apply_robin(&mut coo, &mut rhs, mesh, element_id, face_index, h, u_inf);
```

## API highlights

| Item | Description |
|------|-------------|
| `assemble` | Scatter per-element matrices into a global `Coo`. |
| `reduce_system` / `ReducedSystem` | Static condensation of fixed DOFs. |
| `solve_with_dirichlet` | Dirichlet BCS + sparse solve + scatter-back. |
| `boundary_faces` | Enumerate `(element_id, face_index)` on the mesh boundary. |
| `apply_neumann` / `apply_neumann_order` | Natural (flux) boundary loads. |
| `apply_robin` / `apply_robin_order` | Convective boundary contributions. |

## Position in the crate stack

```text
tpt-fem-thermal / tpt-fem-elasticity ──► tpt-fem-assembly ──► tpt-fem-sparse
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
