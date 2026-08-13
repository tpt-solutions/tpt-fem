# tpt-fem-eigen

Sparse eigenvalue solvers (power iteration, shift-invert, Lanczos) for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

Provides:

- `power_iteration` — the dominant (largest-magnitude) eigenpair,
- `inverse_iteration` — the eigenpair nearest a target shift (shift-invert),
- `lanczos_eigs` — a few extreme eigenpairs via Lanczos tridiagonalisation
  with a dense Jacobi eigensolve of the projected matrix.

All routines operate on a `Coo` matrix from `tpt-fem-sparse` and are
physics-agnostic (used by `tpt-fem-elasticity` for modal analysis).

## Installation

```toml
[dependencies]
tpt-fem-eigen = "0.1"
```

## Usage

```rust
use tpt_fem_eigen::{power_iteration, inverse_iteration, lanczos_eigs, EigWhich};

// Largest-magnitude eigenpair.
let (lambda, v) = power_iteration(&coo, 1000, 1e-10);

// Eigenpair nearest a shift (shift-invert / inverse iteration).
let (lambda, v) = inverse_iteration(&coo, 1.0, 1000, 1e-10);

// A few extreme eigenpairs via Lanczos. `EigWhich` selects smallest/largest
// magnitude or algebraic value.
let pairs = lanczos_eigs(&coo, EigWhich::SmallestMagnitude, 6, 200, 1e-8);
```

A generalised eigenproblem `K x = λ M x` is also supported:

```rust
use tpt_fem_eigen::generalized_lanczos_eigs;

let pairs = generalized_lanczos_eigs(&stiffness, &mass, EigWhich::SmallestMagnitude, 6, 200, 1e-8);
```

## API highlights

| Item | Description |
|------|-------------|
| `matvec` / `rayleigh` | Sparse matrix-vector product and Rayleigh quotient. |
| `power_iteration` | Dominant eigenpair. |
| `inverse_iteration` | Shift-invert eigenpair nearest a target. |
| `EigWhich` | Eigenvalue selection (smallest/largest magnitude or value). |
| `lanczos_eigs` | Extreme eigenpairs via Lanczos. |
| `generalized_lanczos_eigs` | `K x = λ M x` generalised eigenproblem. |

## Position in the crate stack

```text
tpt-fem-sparse ──► tpt-fem-eigen ──► tpt-fem-elasticity (modal)
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
