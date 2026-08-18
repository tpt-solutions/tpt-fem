# tpt-fem crates.io publish tracker

Publish order is a topological sort: every crate is published only after all of
its internal `tpt-fem-*` dependencies. Crates are grouped in batches of 5.
Mark each line `[x]` as you publish it.

## Batch 1 (publish first — no internal deps / leaf crates)

- [ ] tpt-fem-quadrature  — no internal deps
- [ ] tpt-fem-element     — depends on: quadrature
- [ ] tpt-fem-mesh        — depends on: element
- [ ] tpt-fem-sparse      — no internal deps
- [ ] tpt-fem-assembly    — depends on: sparse, element, mesh, quadrature

## Batch 2

- [ ] tpt-fem-solve       — depends on: sparse
- [ ] tpt-fem-eigen       — depends on: sparse
- [ ] tpt-fem-thermal     — depends on: assembly, sparse, element, mesh, quadrature
- [ ] tpt-fem-elasticity  — depends on: assembly, sparse, element, mesh, quadrature, eigen
- [ ] tpt-fem-io-vtk      — depends on: mesh

## Batch 3

- [ ] tpt-fem-io-abaqus   — depends on: mesh
- [ ] tpt-fem-io-exodus   — depends on: mesh
- [ ] tpt-fem-mesh-gen    — depends on: mesh
- [ ] tpt-fem            — umbrella crate; depends on: all of the above
- [ ] tpt-fem-cli        — depends on: tpt-fem

## Summary

Total crates: 15
Batches: 3 (5 + 5 + 5)

Dependency notes:
- `tpt-fem-sparse` has NO internal deps (it only depends on external crates), so it
  can be published alongside the early batches even though it is a "parent" of
  assembly/solve/eigen/thermal/elasticity.
- `tpt-fem` must come after every sub-crate it re-exports, since it references them
  by published version. `tpt-fem-cli` comes last as it depends on `tpt-fem`.
- The `io-*` and `mesh-gen` crates only depend on `tpt-fem-mesh`, so they can be
  published any time after `mesh`.
