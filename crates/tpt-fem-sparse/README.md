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
| `solve` / `solve_multi` | Dense LU factorisation and solve via `tpt-math-linalg-dense` (single/multiple RHS). |
| `SparseError` | Error type for singular / non-finite systems. |

## Position in the crate stack

```text
tpt-fem-sparse ◄── tpt-fem-assembly ◄── tpt-fem-thermal / tpt-fem-elasticity
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
