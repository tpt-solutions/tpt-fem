# tpt-fem-dynamic

Time integration for [tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the
mesh-based finite element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

Provides two standard schemes for the generic second-order system

```text
M·ü + C·v + K·u = f(t)
```

assembled from any physics crate (e.g. the mass and stiffness matrices produced
by `tpt-fem-elasticity`):

* `newmark` — implicit [Newmark-beta](https://en.wikipedia.org/wiki/Newmark-beta_method)
  integration (defaults to the unconditionally-stable average-acceleration rule,
  `β = 0.25`, `γ = 0.5`), and
* `central_difference` — explicit central-difference integration, which lumps
  the mass matrix to a diagonal so that `M⁻¹` is a cheap component-wise divide.

Both routines take the matrices as `Coo` accumulators from `tpt-fem-sparse` and
return the displacement history as `(t, u)` pairs for steps `0..=nsteps` (step 0
is the initial state). Small sparse-matrix helpers `coo_scale`, `coo_add`, and
`coo_matvec` are also exposed for building effective-operator right-hand sides.

> Out of scope (intentionally not implemented): only the Newmark-beta and
> explicit central-difference families are provided — no HHT-α, no generalized-α,
> no nonlinear (Newton) corrector, and no parallel/large-scale sparse solve
> (the `russell` backend is not used here).

## Installation

```toml
[dependencies]
tpt-fem-dynamic = "0.1"
# The matrices are Coo accumulators from tpt-fem-sparse:
tpt-fem-sparse = "0.1"
# and are usually produced by a physics crate such as:
tpt-fem-elasticity = "0.1"
```

## Usage

```rust
use tpt_fem_dynamic::{central_difference, newmark, NewmarkOptions};
use tpt_fem_sparse::Coo;

// SDOF oscillator: M=1, K=4 (ω=2). Free vibration from u0=1, v0=0.
let m = Coo { rows: vec![0], cols: vec![0], vals: vec![1.0] };
let k = Coo { rows: vec![0], cols: vec![0], vals: vec![4.0] };
let c = Coo::new();
let nsteps = 200;
let opts = NewmarkOptions { dt: 0.01, beta: 0.25, gamma: 0.5 };
let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);
let (t, u) = hist[nsteps].clone();
// Closed form u(t) = cos(ω t), ω = 2.
let want = (2.0 * t).cos();
assert!((u[0] - want).abs() < 1e-3, "got {} want {}", u[0], want);
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `newmark_sdof` | `cargo run -p tpt-fem-dynamic --example newmark_sdof` | SDOF oscillator integrated with `newmark`; compared to the closed form `cos(ωt)`. |
| `central_difference_sdof` | `cargo run -p tpt-fem-dynamic --example central_difference_sdof` | Same oscillator with explicit `central_difference`, noting the `Δt < 2/ω` stability limit. |
| `damped_sdof` | `cargo run -p tpt-fem-dynamic --example damped_sdof` | Damped SDOF response `e^{-ζωₙt}·(…)` compared to its analytical decay. |
| `two_dof_chain` | `cargo run -p tpt-fem-dynamic --example two_dof_chain` | A 2-DOF spring chain integrated and checked against the two-mode closed form. |
| `coo_helpers` | `cargo run -p tpt-fem-dynamic --example coo_helpers` | Standalone demo of `coo_scale`, `coo_add`, and `coo_matvec`. |

## API highlights

| Item | Description |
|------|-------------|
| `NewmarkOptions` | `dt`, `β`, `γ` for the implicit integrator (defaults: 0.01, 0.25, 0.5). |
| `newmark` | Implicit Newmark-beta integration of `M·ü + C·v + K·u = f(t)`. |
| `CentralOptions` | `dt` for the explicit integrator (default 0.01). |
| `central_difference` | Explicit central-difference integration; internally lumps the mass. |
| `coo_scale` | `s·A` as a new `Coo`. |
| `coo_add` | `A + B` as a new `Coo` (duplicates summed on collapse). |
| `coo_matvec` | `y = A x` for a `Coo`; empty matrix yields the zero vector. |

## Position in the crate stack

```text
tpt-fem-sparse ──► tpt-fem-dynamic
```

At the code level `tpt-fem-dynamic` only needs `tpt-fem-sparse` (for `Coo` and
`solve`); it is designed to consume the `M`, `C`, and `K` matrices produced by
`tpt-fem-elasticity` / `tpt-fem-assembly`.

## References

* N. M. Newmark, "A Method of Computation for Structural Dynamics", *Journal of
  the Engineering Mechanics Division*, ASCE, 85(3):67–94, 1959.
* T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite
  Element Analysis*, Dover, 2000.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
