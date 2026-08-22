# tpt-fem-amr

Adaptive *h*-refinement for the scalar Poisson problem on the unit square,
part of the [tpt-fem](https://github.com/tpt-solutions/tpt-fem) finite-element
core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

The mesh is a **1-irregular quadtree** over `[0,1]²`: every leaf can be split
into four children, and the tree is re-balanced after each refinement round so
adjacent leaves differ by at most one level. Hanging nodes are eliminated at
assembly time with the linear constraint `u_h = (u_a + u_b)/2`, so the Q1
space stays conforming without multi-point-constraint bookkeeping.

The adaptive loop is the classic identify–mark–refine cycle:

1. assemble and solve `−Δu = f` (Q1, 2×2 Gauss),
2. estimate per-element error with a Zienkiewicz–Zhu gradient-recovery
   indicator (`zz_estimates`),
3. mark with Dörfler bulk marking (`theta` fraction of the estimated error),
4. refine, re-balance, repeat until the element budget is reached.

## Usage

```rust
use tpt_fem_amr::{solve_adaptive, AmrOptions};

// Localized source: adaptivity concentrates elements near the peak.
let f = |x: f64, y: f64| (-500.0 * ((x - 0.5).powi(2) + (y - 0.25).powi(2))).exp();
let g = |_: f64, _: f64| 0.0;
let res =
    solve_adaptive(&f, &g, &AmrOptions { max_elements: 800, ..Default::default() }).unwrap();
assert!(res.mesh.len() <= 800);
println!("final error estimate: {}", res.estimated_error);
```

Via the umbrella crate (feature-gated):

```toml
[dependencies]
tpt-fem = { version = "0.1", features = ["amr"] }
```

```rust
use tpt_fem::{solve_adaptive, AmrOptions};
```

Or from the command line:

```sh
cargo run -p tpt-fem-cli -- amr --max-elements 512 --output amr.vtk
```

## Key API

| Item | Description |
|------|-------------|
| `QuadTree` / `CellKey` | The leaf set of a quadtree over `[0,1]²`; `refine`/`balance`. |
| `build_mesh(&QuadTree) -> HangingMesh` | Conforming leaf mesh + hanging-node constraints. |
| `solve_poisson(&HangingMesh, &f, &g)` | One Q1 solve with hanging-node elimination. |
| `zz_estimates(&HangingMesh, &u)` | Zienkiewicz–Zhu per-element error estimates. |
| `solve_adaptive(&f, &g, &AmrOptions)` | Full identify–mark–refine loop → `AdaptiveSolution`. |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
