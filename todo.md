# tpt-fem — Build Todo

> Tracks bootstrap + Phase 1 build-out for the tpt-fem core, per `spec.txt`
> and `tpt-rust-map/registry.toml`. Crates.io publishing is intentionally
> **out of scope** for this pass — crates stop at `status = "git"` in the
> registry, not `"published"`. License for every crate: `MIT OR Apache-2.0`.
> Author: TPT Solutions.

## Phase 0 — Repo Bootstrap

- [x] Copy `template/Cargo.toml` → root `Cargo.toml`, adapt `[workspace.package]`
- [x] Copy `template/rust-toolchain.toml`
- [x] Copy `template/rustfmt.toml`
- [x] Copy `template/deny.toml`
- [x] Copy `template/.github/workflows/ci.yml` (drop the `no_std` job — no
      `no_std = true` crates here yet)
- [x] Copy `template/LICENSE-MIT` and `template/LICENSE-APACHE`
- [x] Create `crates/` directory
- [x] Add a Rust `.gitignore`
- [x] Write root `README.md`
- [x] Adapt `spec.txt` from `template/spec.txt`
- [x] `git init` (local only — no GitHub remote/push)
- [x] Initial commit
- [x] Sanity check: `cargo build` succeeds

## Per-Crate Checklist Template

**Standard crate:**
1. Scaffold `crates/<name>/` (Cargo.toml inheriting workspace fields, `lib.rs`)
2. Wire dependencies (internal `tpt-fem-*` + external wraps)
3. Implement scope
4. Unit tests + doctests
5. Rustdoc (crate-level + public API)
6. `cargo fmt --check` / `cargo clippy --all-targets --all-features -- -D warnings` clean
7. `cargo deny check` clean
8. Update `tpt-rust-map/registry.toml`: `status = "planned"` → `"git"`

**Umbrella crate:** same as tpt-math's umbrella template — Cargo features
gating each constituent re-export, no direct implementation.

---

## Phase 1 — Core

### 1a — tpt-fem-quadrature

*Gauss quadrature rules for reference elements. No internal deps.*

- [x] Scaffold `crates/tpt-fem-quadrature/`
- [x] Implement 1D Gauss-Legendre rules (orders 1-5), tensor-product
      quad/hex rules built from the 1D rule
- [x] Implement fixed low-order triangle rules (1-point degree-1, 3-point
      degree-2, Hammer-Stroud)
- [x] Implement fixed low-order tetrahedron rules (1-point degree-1,
      4-point degree-2, Keast)
- [x] Unit tests: each rule integrates low-degree monomials exactly
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] registry.toml: `tpt-fem-quadrature` → `"git"`

### 1b — tpt-fem-element

*Reference elements + P1 Lagrange shape functions + isoparametric Jacobian
mapping. Depends on: tpt-fem-quadrature.*

- [x] Scaffold `crates/tpt-fem-element/`
- [x] Implement `Line2`, `Tri3`, `Quad4`, `Tet4`, `Hex8` shape functions +
      local (reference-coordinate) derivatives
- [x] Implement isoparametric Jacobian: physical-coordinate derivatives,
      determinant, from node coordinates + local derivatives
- [x] Unit tests: partition of unity (shape functions sum to 1), Jacobian
      of an undistorted reference element is identity/known constant
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] registry.toml: `tpt-fem-element` → `"git"`
- [ ] **Follow-up, not this pass:** P2/quadratic shape functions (Tri6,
      Quad8/9, Tet10, Hex20/27)

### 1c — tpt-fem-mesh

*Mesh data structures, DOF numbering, Gmsh import. No internal deps (see
spec.txt for the planned tpt-geom migration).*

- [x] Scaffold `crates/tpt-fem-mesh/`
- [x] Wire deps: `mshio`
- [x] Implement `Node`/`Element`/`Mesh` types + manual mesh-builder API
- [x] Implement DOF numbering (configurable dofs-per-node)
- [x] Implement Gmsh `.msh` v4.1 import via `mshio`
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] registry.toml: `tpt-fem-mesh` → `"git"`

### 1d — tpt-fem-sparse

*COO/CSR assembly adapter + faer-backed sparse solve. No internal deps.*

- [x] Scaffold `crates/tpt-fem-sparse/`
- [x] Wire deps: `faer`
- [x] Implement COO triplet accumulator (duplicate-summing)
- [x] Implement CSR conversion
- [x] Implement `solve()` dispatching to faer's sparse LU (default backend)
- [x] Unit tests: solve a small known linear system, check against expected
      solution
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] registry.toml: `tpt-fem-sparse` → `"git"`
- [ ] **Follow-up, not this pass:** feature-gated `russell_sparse`
      (SuiteSparse/MUMPS) backend for large-scale problems — verify
      SuiteSparse component licenses individually before enabling

### 1e — tpt-fem (umbrella)

*Re-exports quadrature + element + mesh + sparse behind Cargo features.*

- [x] Scaffold `crates/tpt-fem/`
- [x] Wire optional deps + feature flags for quadrature/element/mesh/sparse
- [x] Re-export each constituent's public API behind its feature
- [x] **End-to-end patch test** (`tests/patch_test.rs`): hand-built 2-3
      element mesh with a textbook-known stiffness matrix and analytical
      solution, driven through element shape functions + quadrature →
      triplet assembly → sparse solve, asserting the result matches the
      analytical solution within tolerance
- [x] Rustdoc documenting the feature matrix
- [x] `cargo fmt` / `clippy` / `deny` clean across feature combinations
- [x] registry.toml: `tpt-fem` → `"git"`

## Final Phase 1 Closeout

- [x] `cargo test --workspace --all-features` passes
- [x] `cargo clippy --workspace --all-targets --all-features -- -D warnings` clean
- [x] `cargo deny check` clean workspace-wide
- [x] Confirm every Phase 1 `tpt-fem-*` entry in `tpt-rust-map/registry.toml`
      reads `status = "git"`

---

## Phase 2 — Assembly + First Physics (not started)

- [ ] `tpt-fem-assembly` — element-to-global scatter, Dirichlet/Neumann/Robin
      boundary conditions. Depends on: tpt-fem-sparse, tpt-fem-element,
      tpt-fem-mesh.
- [ ] `tpt-fem-thermal` — heat conduction / Poisson elements. Depends on:
      tpt-fem-assembly.
- [ ] `tpt-fem-io-vtk` — wrap vtkio for ParaView export. Depends on:
      tpt-fem-mesh.

## Phase 3 — Structural + Nonlinear (not started)

- [ ] `tpt-fem-elasticity` — bar/beam/plane-stress/plane-strain/3D continuum
      linear elasticity. Depends on: tpt-fem-assembly.
- [ ] `tpt-fem-solve` — Newton-Raphson (reuse tpt-math-optimize-general's
      argmin pattern) + arc-length/continuation. Depends on: tpt-fem-assembly.

## Phase 4 — Ecosystem-Gap Crates, Stretch (not started)

- [ ] `tpt-fem-eigen` — sparse shift-invert Lanczos/Arnoldi eigensolver.
      Depends on: tpt-fem-sparse.
- [ ] `tpt-fem-io-exodus` — Exodus II reader/writer. Depends on: tpt-fem-mesh.
- [ ] `tpt-fem-io-abaqus` — Abaqus `.inp` reader/writer. Depends on:
      tpt-fem-mesh.
- [ ] Native 3D quality tet-mesh generation — currently blocked on licensing
      (tritet/TetGen is AGPL). Revisit only if Gmsh-file-import proves
      insufficient for a specific downstream need.
