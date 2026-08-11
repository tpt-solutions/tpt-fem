# tpt-fem

A modular finite-element method (FEM) workspace for Rust, part of the
`tpt-solutions` ecosystem. This repository currently hosts **Phase 1** of the
planned `tpt-fem` stack: the five foundational, dependency-light core crates
used by every higher-level FEM crate.

## Crates (Phase 1)

| Crate | Purpose | Depends on |
|-------|---------|------------|
| [`tpt-fem-quadrature`](crates/tpt-fem-quadrature) | Fixed-order Gauss quadrature on reference elements (line / tri / quad / tet / hex). | — |
| [`tpt-fem-element`](crates/tpt-fem-element) | Reference elements, Lagrange `P1` shape functions, isoparametric Jacobian map. | `tpt-fem-quadrature` |
| [`tpt-fem-mesh`](crates/tpt-fem-mesh) | Mesh data model, DOF numbering, Gmsh `.msh` v4.1 import (via `mshio`). | `mshio` |
| [`tpt-fem-sparse`](crates/tpt-fem-sparse) | FEM-specific COO/CSR assembly adapter with a `faer`-backed sparse LU solve. | `faer` |
| [`tpt-fem`](crates/tpt-fem) | Umbrella crate re-exporting the above behind Cargo features, plus the end-to-end patch test. | all of the above |

The umbrella crate's feature flags are: `quadrature`, `element`, `mesh`, and
`sparse` (all enabled by default).

## Pipeline

A standard steady-diffusion solve flows through the crates as:

```
reference-element shape functions + gradients   (tpt-fem-element)
        + quadrature rules                        (tpt-fem-quadrature)
        -> per-element stiffness matrices
        -> triplet accumulation into a COO        (tpt-fem-sparse)
        -> CSR, sparse LU factorization + solve    (tpt-fem-sparse / faer)
```

The integration test `tpt-fem/tests/patch_test.rs` drives this entire pipeline
on a hand-built multi-element mesh and checks the result against an analytical
solution.

## Building & testing

```bash
cargo build --workspace --all-features
cargo test  --workspace --all-features
cargo fmt   --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo deny  check
```

## Status

Phase 1 (`tpt-fem-quadrature`, `tpt-fem-element`, `tpt-fem-mesh`,
`tpt-fem-sparse`, `tpt-fem`) is implemented and tracked as `git` in
[`registry.toml`](../../registry.toml). Follow-up phases (assembly, I/O,
solvers, elasticity, etc.) remain `planned`.

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at
your option.
