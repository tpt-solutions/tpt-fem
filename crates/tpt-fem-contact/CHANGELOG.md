# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `Octree` — uniform octree over a fixed point set with pruned
  nearest-neighbour queries (O(log n) per query on evenly distributed
  points).
- `contact_pairs_octree` — drop-in octree-accelerated replacement for the
  brute-force `contact_pairs` surface pairing; results agree up to
  equidistant tie-breaking. This is the BVH/octree spatial search flagged as
  a future optimisation since Phase 12; the constraint enforcement is
  unchanged and remains search-independent.
- Regression tests: octree-vs-brute-force agreement on deterministic
  pseudo-random point clouds, 2-D/duplicate/empty edge cases, and full-pairing
  agreement of `contact_pairs_octree` with `contact_pairs`.

## [0.1.0] - 2026-08-20

### Added

- `ContactConstraint` — unilateral `x_dof ≥ lower` constraint on a global DOF.
- `penalty_contact` — augments a base `Coo` with penalty stiffness `ε_N` and load `ε_N·lower`, returning `(Coo, Vec<f64>)`.
- `augmented_lagrangian` — augmented-Lagrangian iteration (`λ ← max(0, λ + ε_N·(lower − x))`) returning `(Vec<f64> displacements, Vec<f64> multipliers)`.
- `contact_pairs` — brute-force nearest-node pairing between two surfaces, returning `(node_a, idx_b, gap)` tuples.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-contact-0.1.0
