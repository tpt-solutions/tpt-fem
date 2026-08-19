# tpt-fem-contact

Contact mechanics for [tpt-fem](https://github.com/tpt-solutions/tpt-fem) —
the mesh-based finite element core from
[tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate provides unilateral (non-penetration) contact enforcement written
from scratch. Two enforcement strategies operate on a caller-supplied base
stiffness [`Coo`] plus a set of [`ContactConstraint`]s, each expressing a
global DOF that must satisfy `x_dof ≥ lower` (a rigid obstacle located at
`lower` along that DOF's axis — e.g. a node's normal coordinate relative to a
wall):

- `penalty_contact` — adds the diagonal stiffness `ε_N` at the constrained DOF
  and a balancing load `ε_N·lower`, i.e. the energy penalty
  `½ ε_N (min(0, lower − x))²`. It is simple and cheap but, being soft, leaves
  a residual penetration that shrinks only as `ε_N → ∞`.
- `augmented_lagrangian` — iterates the penalty solve while folding the current
  Lagrange multiplier into the load and updating
  `λ ← max(0, λ + ε_N·(lower − x))` each step. This drives the penetration to
  (numerically) zero at finite `ε_N`, recovering the hard-contact reaction.

A minimal from-scratch [`contact_pairs`] helper performs a brute-force
nearest-node pairing between two surfaces (an O(|a|·|b|) scan; a BVH/octree is
a future optimisation and is independent of the constraint enforcement above).

Sign conventions:

- The constraint is `x ≥ lower`. A node penetrating the obstacle has
  `x < lower`, giving a positive penetration (gap)
  `g_N = lower − x > 0`.
- The penalty contact force on the node is `f = ε_N ⟨g_N⟩` (active only while
  `g_N > 0`, i.e. while penetrating).
- The Lagrange multiplier `λ` is the contact reaction and balances the applied
  load at convergence (e.g. for a single DOF pushed by `F` into a wall at
  `lower = 0`, `λ → −F`).

Out of scope / not implemented:

- No surface-to-surface mortar or segment-to-segment contact; constraints are
  expressed on individual global DOFs against a fixed lower bound.
- No friction (the model is frictionless, normal-only contact).
- No curved-obstacle or self-contact search; [`contact_pairs`] is a
  brute-force nearest-node probe, not a full detection pipeline.

## Installation

```toml
[dependencies]
tpt-fem-contact = "0.1"
# The constraint solvers augment a sparse system and return multipliers:
tpt-fem-sparse = "0.1"
```

## Usage

```rust
use tpt_fem_contact::{augmented_lagrangian, ContactConstraint};
use tpt_fem_sparse::Coo;

// One DOF on a spring (k=10) to ground, pushed into a wall at x=0 by F=-4
// (force toward -x). Hard contact should hold x ~ 0 with reaction lambda ~ 4.
let base = Coo { rows: vec![0], cols: vec![0], vals: vec![10.0] };
let load = vec![-4.0];
let con = ContactConstraint { dof: 0, lower: 0.0 };
let (u, lambda) = augmented_lagrangian(&base, &load, &[con], 1e4, 50, 1e-9);
assert!(u[0].abs() < 1e-6);
assert!((lambda[0] - 4.0).abs() < 1e-3);
```

## Examples

- `contact_pairs` — brute-force proximity search between two opposing node
  sets, printing the detected pairs and their gaps.

  ```text
  cargo run -p tpt-fem-contact --example contact_pairs
  ```

- `penalty_contact` — penalty force and stiffness for a prescribed penetration,
  checked against the hand-computed `ε_N·g_N`. The node sits just past the wall.

  ```text
  cargo run -p tpt-fem-contact --example penalty_contact
  ```

- `augmented_lagrangian` — augmented-Lagrangian iteration on the same problem,
  showing the multiplier converge to the wall reaction while the penetration
  collapses far below the pure-penalty result.

  ```text
  cargo run -p tpt-fem-contact --example augmented_lagrangian
  ```

## API highlights

| Item | Description |
|------|-------------|
| `ContactConstraint` | A unilateral `x_dof ≥ lower` constraint on a global DOF. |
| `penalty_contact` | Augments a base `Coo` with `ε_N` stiffness and `ε_N·lower` load; returns `(Coo, Vec<f64>)`. |
| `augmented_lagrangian` | ALM iteration returning `(displacements, multipliers)`. |
| `contact_pairs` | Brute-force nearest-node pairing between two surfaces; returns `(node_a, idx_b, gap)`. |

## Position in the crate stack

```text
tpt-fem-sparse ──► tpt-fem-contact
```

## References

- P. Wriggers, *Computational Contact Mechanics*, 2nd ed., Springer, 2006 —
  penalty and augmented-Lagrangian enforcement of unilateral contact.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
