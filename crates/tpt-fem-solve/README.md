# tpt-fem-solve

Nonlinear solvers for [tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the
mesh-based finite element core from
[tpt-solutions](https://github.com/tpt-solutions).

## Overview

Provides a generic Newton–Raphson iteration for a residual `R(u) = 0` together
with its Jacobian, and a parameter-continuation driver that warm-starts Newton
at successive load/parameter steps. Essential (Dirichlet) conditions are
condensed out of the system on every iteration, consistently with
`tpt-fem-assembly`'s linear solver.

The residual and Jacobian are supplied by the caller (e.g. assembled from a
nonlinear weak form), so this crate is physics-agnostic.

## Installation

```toml
[dependencies]
tpt-fem-solve = "0.1"
```

## Usage

```rust
use tpt_fem_solve::{newton, NewtonOptions};

// residual(u) -> R, jacobian(u) -> J
let options = NewtonOptions { tol: 1e-10, max_iter: 50, ..Default::default() };
let u = newton(initial_guess, &residual, &jacobian, options)?;
```

Parameter continuation:

```rust
use tpt_fem_solve::continuation;

// Walk a load parameter from `start` to `end` in `steps`, warm-starting each
// Newton solve with the previous solution.
let path = continuation(start, end, steps, &residual_at, &jacobian_at, options)?;
```

Arc-length continuation:

```rust
use tpt_fem_solve::{arc_length_continuation, ArcLengthOptions};

let path = arc_length_continuation(&residual_at, &jacobian_at, &tangent, ArcLengthOptions {
    radius, ..Default::default()
})?;
```

## API highlights

| Item | Description |
|------|-------------|
| `newton` | Newton–Raphson solve with `NewtonOptions`. |
| `NewtonOptions` / `NewtonError` | Tolerances, iteration cap, error reporting. |
| `continuation` | Parameter-continuation driver over a residual/jacobian family. |
| `arc_length_continuation` | Arc-length (Ramm) continuation. |
| `ArcLengthOptions` / `ArcLengthError` | Step control and error reporting. |

## Position in the crate stack

```text
tpt-fem-assembly / tpt-fem-sparse ──► tpt-fem-solve (physics-agnostic)
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `newton_scalar` | `cargo run -p tpt-fem-solve --example newton_scalar` | Newton–Raphson root of `u³ - u - 1 = 0` (plastic number ≈ 1.3247). |
| `continuation` | `cargo run -p tpt-fem-solve --example continuation` | Parameter continuation of `u² = λ` from λ = 1 to 4, checking `u = 2`. |
| `arc_length` | `cargo run -p tpt-fem-solve --example arc_length` | Arc-length continuation through the fold of `u³ - 3u = λ`. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
