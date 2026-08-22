# tpt-fem-topopt

SIMP topology optimization (compliance minimisation) built on the
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) elasticity stack.

```rust
use tpt_fem_topopt::{cantilever_problem, simp_optimize, SimpOptions};

let problem = cantilever_problem(32, 16);
let opts = SimpOptions { volfrac: 0.4, ..Default::default() };
let result = simp_optimize(&problem, &opts).unwrap();
// result.densities: one density per element in mesh element order.
```

The optimizer runs the classic SIMP loop: assemble `K(x)` from
`x_e^p`-scaled element stiffness matrices, solve the condensed system,
compute compliance sensitivities, apply a sensitivity density filter over
element centroids, and update the design with an optimality-criteria step
(bisection on the Lagrange multiplier to meet the volume fraction) with
under-relaxation for stability.
