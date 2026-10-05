# tpt-fem-sparse

A FEM-specific sparse-matrix assembly adapter with a `tpt-math`-backed solve,
part of [tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the
mesh-based finite element core from
[tpt-solutions](https://github.com/tpt-solutions).

## Overview

Finite-element assembly naturally produces a sparse global matrix as a bag of
`(row, col, value)` triplets, where the same entry may be written many times
(once per element that touches it). `Coo` is a growable coordinate-list
accumulator that supports duplicate-summing assembly; `Coo::to_csr` collapses
it into a canonical compressed-sparse-row `Csr` matrix.

`solve` assembles the matrix into a dense
[`tpt-math-linalg-dense`](https://github.com/tpt-solutions/tpt-math) matrix and
factors it with the in-house partial-pivot LU decomposition, so the workspace
carries no Apache-2.0-only linear-algebra dependency.

## Installation

```toml
[dependencies]
tpt-fem-sparse = "0.1"
```

## Usage

```rust
use tpt_fem_sparse::Coo;

// Assemble [[2, 1], [1, 3]] by writing each entry twice then summing.
let mut c = Coo::new();
c.push(0, 0, 1.0);
c.push(0, 0, 1.0);
c.push(0, 1, 1.0);
c.push(1, 0, 1.0);
c.push(1, 1, 1.5);
c.push(1, 1, 1.5);
let csr = c.to_csr();
assert_eq!(csr.nnz(), 4);
assert_eq!(csr.row_ptrs, vec![0, 2, 4]);
assert_eq!(csr.values, vec![2.0, 1.0, 1.0, 3.0]);

// Solve the system.
let x = tpt_fem_sparse::solve(&c, &[3.0, 5.0]).unwrap();
```

## API highlights

| Item | Description |
|------|-------------|
| `Coo` | Growable coordinate-list accumulator with duplicate summing. |
| `Coo::push` / `Coo::to_csr` | Add entries and collapse to CSR. |
| `Csr` | Compressed-sparse-row matrix (`row_ptrs`, `col_idxs`, `values`). |
| `solve` / `solve_multi` | Single/multiple RHS. Symmetric positive-definite systems with ≥ 200 unknowns are routed to the sparse envelope Cholesky automatically; large unsymmetric/indefinite systems with a narrow reordered band use the banded LU; everything else uses the dense LU from `tpt-math-linalg-dense`. |
| `solve_banded` | Pure-Rust sparse direct solver for general systems (RCM + banded LU with partial pivoting); used automatically for large unsymmetric/indefinite systems with a narrow band. |
| `solve_skyline` | Pure-Rust sparse direct solver for SPD systems (reverse Cuthill-McKee + envelope Cholesky). |
| `solve_cg` | Jacobi-preconditioned conjugate gradients for SPD systems (`O(nnz)` memory). |
| `SparseError` | Error type for singular / non-finite systems. |

## Position in the crate stack

```text
tpt-fem-sparse ◄── tpt-fem-assembly ◄── tpt-fem-thermal / tpt-fem-elasticity
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `coo_to_csr` | `cargo run -p tpt-fem-sparse --example coo_to_csr` | Builds a `Coo` with duplicate entries and checks the summed `Csr` layout. |
| `solve_linear` | `cargo run -p tpt-fem-sparse --example solve_linear` | Solves `[[2,1],[1,3]] x = [3,5]` and checks `x = [0.8, 1.4]`. |
| `solve_multi` | `cargo run -p tpt-fem-sparse --example solve_multi` | Solves two right-hand sides against the same matrix and checks both solutions. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
