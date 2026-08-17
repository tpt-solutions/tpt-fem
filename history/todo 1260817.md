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

## Phase 2 — Assembly + First Physics (complete)

Each crate follows the Per-Crate Checklist Template (scaffold → wire deps →
implement scope → unit tests + doctests → rustdoc → fmt/clippy clean → deny
clean → registry `planned`→`git`). All three pass `cargo test`, `cargo fmt
--check`, `cargo clippy --all-targets --all-features -- -D warnings`, and
`cargo deny check`.

- [x] `tpt-fem-assembly` — element-to-global scatter, Dirichlet/Neumann/Robin
      boundary conditions. Depends on: tpt-fem-sparse, tpt-fem-element,
      tpt-fem-mesh.
- [x] `tpt-fem-thermal` — heat conduction / Poisson elements. Depends on:
      tpt-fem-assembly.
- [x] `tpt-fem-io-vtk` — wrap vtkio for ParaView export. Depends on:
      tpt-fem-mesh.

## Phase 3 — Structural + Nonlinear (complete)

- [x] `tpt-fem-elasticity` — bar/beam/plane-stress/plane-strain/3D continuum
      linear elasticity. Depends on: tpt-fem-assembly.
- [x] `tpt-fem-solve` — Newton-Raphson (reuse tpt-math-optimize-general's
      argmin pattern) + arc-length/continuation. Depends on: tpt-fem-assembly.

## Phase 4 — Ecosystem-Gap Crates, Stretch (complete except tet-mesh)

- [x] `tpt-fem-eigen` — sparse shift-invert Lanczos/Arnoldi eigensolver.
      Depends on: tpt-fem-sparse.
- [x] `tpt-fem-io-exodus` — Exodus II reader/writer (minimal NetCDF-3 codec).
      Depends on: tpt-fem-mesh.
- [x] `tpt-fem-io-abaqus` — Abaqus `.inp` reader/writer. Depends on:
      tpt-fem-mesh.

## Phase 4b — Native 3D Tet-Mesh Generator (from scratch, AGPL-free)

- [x] `tpt-fem-mesh-gen` — native 3D tetrahedral mesh generation, no external
      dependency (replaces the blocked TetGen/tritet AGPL path). Depends on:
      tpt-fem-mesh.
  - [x] Incremental Bowyer–Watson Delaunay tetrahedralisation (`delaunay_3d`)
        of an arbitrary point cloud, with `f64` orientation/in-sphere
        predicates + coincident-point de-duplication.
  - [x] Structured box mesher (`box_mesh`): each brick split into six
        positively-oriented tets — guaranteed valid, intersection-free, and
        quality-bounded with no robustness caveats.
  - [x] Quality metrics (`tet_quality`: dihedral angles, radius-edge ratio) and
        `laplacian_smooth` (boundary nodes held fixed).
  - [x] Unit tests: predicate correctness, single-tet + cube Delaunay
        (closed, positively oriented), box counts/orientation, quality +
        smoothing.
  - [x] `cargo fmt` / `clippy --all-features -D warnings` / `deny` clean.
  - [x] Re-exported by the `tpt-fem` umbrella behind the `mesh-gen` feature
        (default-on).

---

## Umbrella Update

- [x] Extend `tpt-fem` umbrella with features re-exporting the new crates
      (`assembly`, `thermal`, `io-vtk`, `elasticity`, `solve`, `eigen`,
      `io-abaqus`, `io-exodus`) behind Cargo features, defaulting all on.
- [x] Registry: `tpt-rust-map/registry.toml` is the sibling repo (not present
      in this workspace); entries flip `status = "planned"` → `"git"` there when
      that repo is next touched. This repo's `Cargo.toml` already lists all
      Phase 1-4 crates as workspace members.

---

## Final Workspace Closeout (Phase 1-4)

- [x] `cargo test --workspace --all-features` passes (all unit + doc tests green)
- [x] `cargo clippy --workspace --all-targets --all-features -- -D warnings` clean
- [x] `cargo deny check` clean workspace-wide (three vtkio-transitive
      advisories — RUSTSEC-2026-0041/0194/0195 — ignored in `deny.toml`; they
      describe untrusted-input attack surfaces that do not apply to our
      trusted-read/write VTK usage, and vtkio 0.6.3 pins the vulnerable
      `quick-xml`/`lz4_flex` with no in-range patch available)
- [x] `cargo fmt --all --check` clean
- [ ] `tpt-rust-map/registry.toml` entries → `"git"`: deferred (sibling repo
      not present in this workspace)

---

## Phase 5 — Hardening, Error Ergonomics &amp; Prelude

*A security/completeness/adoption audit (2026-08-12) found that several
Phase 1-4 crates check off scope in this file that isn't fully there. This
phase and Phases 6-9 below track fixing that honestly. See the archived
plan for full technical detail:
`C:\Users\Phillip\.claude\plans\review-project-fix-any-lovely-sprout.md`.*

- [x] Error unification: `Display`/`Error` impls for `SparseError`,
      `MeshError`, `NewtonError`, `ExodusError`, `InpError`, new `VtkError`
      (wraps `vtkio::Error` instead of leaking it). Cross-crate `From`
      conversions (`MeshError` → `ExodusError`/`InpError`, `io::Error` →
      `{ExodusError, InpError, VtkError}`, `vtkio::Error` → `VtkError`).
- [x] `tpt-fem-mesh`: `Mesh::validate`, `MeshBuilder::try_add_element`/
      `try_build`, new `MeshError` variants (`DanglingNodeTag`,
      `NodeCountMismatch`, `NodeIndexOutOfRange`).
- [x] `tpt-fem-mesh`: geometric selectors (`nodes_on_plane`, `nodes_in_box`)
      and Gmsh `$Entities` physical-group tags (`Node::region`,
      `Element::region`, `MeshBuilder::add_element_with_region`).
- [x] `tpt-fem-sparse`: `solve_multi` (one factorization, multiple RHS;
      `solve` now delegates to it) — prerequisite for Phase 6's arc-length
      solver.
- [x] `tpt-fem-io-exodus`: harden the hand-rolled NetCDF-3 codec — unbounded
      `Vec::with_capacity` from untrusted counts (process-abort risk),
      `read_name` OOB slice on truncated files, `.unwrap()` + unchecked
      dimension index, `connect0` integer underflow (panics in release
      builds thanks to this workspace's `overflow-checks = true`), and
      connectivity node indices never validated against `num_nodes` (a `0`
      entry silently corrupts the mesh instead of erroring). Switch
      `bytes_to_mesh`/`read_inp` to `try_add_element`/`try_build`. Add a
      `dtype` cross-check. Malformed-input regression tests for each.
- [x] `tpt-fem-mesh`: fix the Gmsh importer's `.expect("node tag present in
      mesh")` panic on a dangling `$Elements` node-tag reference (now
      returns `MeshError::DanglingNodeTag` — done above, keeping this line
      as a pointer to the regression test still owed).
- [x] Prelude module (`tpt_fem::prelude`) on the umbrella crate, feature-gated
      per constituent. Blocked on Phase 6's `solve_frame2d`/`solve_modal`/
      `arc_length_continuation` existing so the prelude can include them.
- [x] `cargo test --workspace --all-features` / `clippy` / `fmt` / `deny`
      clean after this phase.

## Phase 6 — Physics Completeness (beam, arc-length, generalized eigen)

*Closes the gap between what Phases 3-4 check off above and what's actually
implemented.*

- [x] `tpt-fem-assembly`: extract `reduce_system` out of `solve_with_dirichlet`
      (shared prerequisite for generalized eigen + arc-length).
- [x] `tpt-fem-elasticity`: 2-D Euler-Bernoulli frame element
      (`BeamSection2D`, `beam2d_element_matrix`, `beam2d_consistent_mass`,
      `solve_frame2d`) — closed-form Hermite-cubic stiffness/mass, verified
      against textbook cantilever/simply-supported deflections. 3-D beam
      (torsion, biaxial bending, orientation triad) tracked as a follow-up,
      not this pass — same treatment as P2 elements in Phase 1b.
- [x] `tpt-fem-elasticity`: `elasticity_mass_matrix`/`elasticity_lumped_mass`
      (consistent/lumped mass, reusing the existing quadrature machinery)
      and `solve_modal` (assembles K + M, Dirichlet-reduces both, calls into
      `tpt-fem-eigen`'s generalized solver).
- [x] `tpt-fem-eigen`: `generalized_lanczos_eigs` — M-orthogonal shift-invert
      Lanczos for the generalized symmetric eigenproblem `Kx = λMx`. Plain
      Arnoldi (for genuinely non-symmetric operators, which nothing in this
      workspace's physics crates produces) is a tracked follow-up, not this
      pass.
- [x] `tpt-fem-solve`: real arc-length continuation
      (`arc_length_continuation`, Crisfield spherical constraint via the
      bordering algorithm, adaptive step sizing, failure/cut-back handling).
      Verified against (1) an algebraic fold with a known limit point and
      (2) a hand-written total-Lagrangian 2-bar snap-through truss residual
      (written as a test-only nonlinear residual — a general
      geometric-nonlinearity framework across all elements is explicitly
      out of scope for this pass).
- [x] Update `tpt_fem::prelude` with `solve_frame2d`/`solve_modal`/
      `arc_length_continuation`/`ArcLengthOptions`/`generalized_lanczos_eigs`.
- [x] `cargo test --workspace --all-features` / `clippy` / `fmt` / `deny`
      clean after this phase.

## Phase 7 — CLI, Verification &amp; Fuzz Suite, Docs

- [x] `tpt-fem-cli` (new workspace member, `[[bin]] name = "tpt-fem"`):
      TOML config (`[problem]`/`[mesh]`/`[material]`/`[source]`/`[[bc]]`
      with node/plane/box/region/boundary selectors/`[output]`), `solve`
      (thermal Poisson in this pass, elasticity reserved in the schema as a
      fast-follow), `mesh info`, `mesh convert`. Error UX built on Phase 5's
      `Display` impls.
- [x] `tpt-fem-thermal`: MMS convergence test suite (`tests/
      mms_convergence.rs`) — manufactured solution, L2/H1 convergence-rate
      assertions on a refined mesh sequence (1-D, 2-D, and an `#[ignore]`
      3-D variant).
- [x] New `fuzz/` (cargo-fuzz, excluded from the main workspace):
      `exodus_decode`, `gmsh_import`, `abaqus_inp`, `exodus_roundtrip`
      targets. Nightly compile-only `fuzz-build` CI job; actual fuzzing runs
      live in a separate scheduled workflow.
- [x] `crates/tpt-fem/examples/thermal_solve.rs` using only the prelude;
      dedupe `tests/patch_test.rs`'s private `solve_with_dirichlet` against
      the real `tpt_fem_assembly` one (keep the low-level element-loop
      coverage); new `tests/end_to_end.rs` exercising the real
      mesh→solve→VTK path.
- [x] README.md rewrite covering all 13 crates (currently Phase-1-only).
- [x] Root `Justfile` wrapping the README's cargo commands.
- [x] CI: `Swatinem/rust-cache`, a `docs` job (`cargo doc --workspace
      --no-deps --all-features` with `-D warnings`), a `wasm` job
      (hard-fail for `quadrature`/`element`/`mesh`, soft-fail for `sparse`
      pending confirmation of the `faer`→`rayon` wasm32 compatibility risk).
- [x] crates.io metadata (`keywords`/`categories`/`readme`) across all 13
      crate `Cargo.toml`s.

## Phase 8 — Python Bindings

- [ ] New `crates/tpt-fem-py` (pyo3, `exclude`d from the main workspace):
      `Mesh` class (`load`/`box_mesh`/`coords`/`nodes_on_plane`/
      `nodes_in_box`/`write_vtk`) and `solve_poisson` (constant or
      Python-callback source term, GIL-safe exception propagation).
      `maturin`-based local dev workflow (`maturin develop` + `pytest`); no
      wheel-building or PyPI publishing this pass, matching this repo's
      existing crates.io-publishing-out-of-scope policy.
      Explicitly last: highest external-tooling risk, depends on the CLI's
      config/selector decisions already being settled.
