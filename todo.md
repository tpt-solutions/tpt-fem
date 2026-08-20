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
- [x] **Follow-up (completed):** P2/quadratic shape functions. Added `Tri6`,
       `Quad8` (serendipity), `Quad9` (biquadratic Lagrange), `Tet10`,
       `Hex20` (serendipity) and `Hex27` (triquadratic Lagrange) implementing
       `ReferenceElement` (`nodes`/`shape`/`grad`) in `tpt-fem-element`, with
       partition-of-unity, Kronecker-at-nodes, and identity-Jacobian tests for
       all six. Re-exported via the umbrella.

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

*COO/CSR assembly adapter + dense LU solve. No internal deps.*

- [x] Scaffold `crates/tpt-fem-sparse/`
- [x] Wire deps: `tpt-math-linalg-dense` (in-house, replaced `faer`; see
      Follow-up below)
- [x] Implement COO triplet accumulator (duplicate-summing)
- [x] Implement CSR conversion
- [x] Implement `solve()` dispatching to the in-house dense LU (default
      backend)
- [x] Unit tests: solve a small known linear system, check against expected
      solution
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] registry.toml: `tpt-fem-sparse` → `"git"`
- [x] **Follow-up (completed):** feature-gated `russell_sparse`
       (SuiteSparse/MUMPS) backend for large-scale problems. `tpt-fem-sparse`
       gains a `russell` feature (optional `russell_sparse` + `russell_lab`
       deps) exposing `solve_russell(coo, rhs, genie)` and
       `solve_russell_multi(coo, rhs, genie)` (UMFPACK/MUMPS, one
       factorization + multiple RHS). The umbrella forwards it via its own
       `russell` feature, and it surfaces in the prelude when enabled. The
       `SuiteSparse`/`MUMPS` component licenses must still be vetted before
       enabling in a distribution. **Caveat:** `russell_sparse`'s build script
       requires a SuiteSparse/MUMPS toolchain (`MSYS2_PREFIX` on Windows); it is
       therefore compile-validated only where that toolchain exists and is not
       built on this dev box. **(See also:** the default `solve()`/
       `solve_multi()` backend was later swapped from `faer` sparse LU to an
       in-house dense LU via `tpt-math-linalg-dense`, to drop the
       Apache-2.0-only `faer` dependency; `russell` remains the opt-in
       sparse-direct path for large problems.**)**

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
       returns `MeshError::DanglingNodeTag`). Covered by regression test
       `import_rejects_dangling_node_tag` in `crates/tpt-fem-mesh/src/lib.rs`.
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

- [x] New `crates/tpt-fem-py` (pyo3, `exclude`d from the main workspace):
      `Mesh` class (`load`/`box_mesh`/`coords`/`node_count`/`nodes_on_plane`/
      `nodes_in_box`/`write_vtk`) and `solve_poisson` (constant or
      Python-callback source term, GIL-safe exception propagation).
      `maturin`-based local dev workflow (`maturin develop` + `pytest`); no
      wheel-building or PyPI publishing this pass, matching this repo's
      existing crates.io-publishing-out-of-scope policy.
      - `Cargo.toml` is self-contained (concrete `workspace.*` fields + path
        dep to `../tpt-fem`, since the crate is excluded) and pins
        `pyo3 = "0.23"` so it builds against Python 3.13.
      - `[lib] name` and `#[pymodule]` aligned to `tpt_fem` so `import tpt_fem`
        resolves (previously mismatched with `tpt_fem_py`).
      - Validated end-to-end: `maturin develop` builds + installs the wheel,
        `pytest` passes both `test_box_mesh_and_solve` and
        `test_python_callback_source` on CPython 3.13.
      - Crate-local `.gitignore` covers `/target`, `.pytest_cache`,
        `__pycache__`, `*.pyd`.
      Explicitly last: highest external-tooling risk, depends on the CLI's
      config/selector decisions already being settled.

---

## Phase 9 — Platform Review Follow-ups (2026-08-17)

*A platform review (bugs/TODOs/missing features/usability/automation) on
2026-08-17 found the items below. Nothing here has been implemented yet —
tracked for a future pass.*

### 9a — Bugs

- [x] `tpt-fem-io-abaqus`: `abaqus_type_to_cell` maps `S4R` (4-node shell
      quad) to `Quad` instead of `Tet` — corrupts topology for any
      shell-element Abaqus deck. (Already fixed in committed code; regression
      test `s4r_shell_quad_maps_to_quad` guards it.)
- [x] `tpt-fem-py`: `solve_poisson`'s Python source-callback error path
      (`cb.bind(py).call1(args)`) now captures the `PyErr` in a `RefCell` and
      surfaces it after the GIL-released solve returns, instead of silently
      returning `0.0`. (Already implemented; verified against
      `crates/tpt-fem-py/src/lib.rs`.)
- [x] `tpt-fem-cli`: README usage example (`tpt-fem solve --config
      problem.toml --out result.vtk`) doesn't match the actual `clap`
      definition (`config` is a bare positional, no `--config`/`--out`
      flags) — fixed the crate `README.md`, and added a drift-guard test
      (`cli_usage_matches_readme`) that asserts the clap usage string contains
      `solve <CONFIG>` for `solve`/`elasticity`/`modal`.
- [x] `tpt-fem-sparse`: default `solve()`/`solve_multi()` is O(n²) dense LU
      regardless of input sparsity; only the non-default `russell` feature
      (external SuiteSparse/MUMPS toolchain) sparse-scales. Loud doc warning
      already present on `solve()`/`solve_multi()` (and the crate README);
      a DOF-count heuristic was judged unnecessary given the explicit guidance.

### 9b — Panic-risk hardening (library-facing, not test code)

- [x] `tpt-fem-sparse/src/lib.rs`: `solve()`'s `.expect()` on `solve_multi`'s
      internal invariant → returns `SparseError` instead. (Already returns
      `SparseError` via `ok_or_else`.)
- [x] `tpt-fem-elasticity/src/lib.rs`: `strain_dim`/`constitutive`/
      `b_matrix` `panic!` on mismatched `ElasticModel`/dimension (e.g.
      `PlaneStress` with `dim=3`) → return a `Result` instead. `ElasticityError`
      added; `elasticity_element_matrix` returns `Result`; assembly gained a
      fallible `try_assemble`. Regression test
      `element_matrix_rejects_model_dim_mismatch` guards it.
- [ ] `tpt-fem-element/src/lib.rs`: `mat_det`/`mat_inv`/`tensor_nodes` panic
      on out-of-range dimension in private helpers — low priority (currently
      unreachable: the spatial dimension is always 1/2/3 for the element types
      this crate builds) and converting them to `Result` would ripple the public
      `Map`/`from_nodes_and_grad` API for negligible gain. **Deferred** with a
      note rather than changed.
- [x] `tpt-fem-solve/src/lib.rs`: the arc-length Newton loop already handles
      `dense_solve`/`bordered_solve` errors via the cut-back path (no `.unwrap()`
      on the dense-solve result remains in library code). Verified.
- [x] `tpt-fem-io-exodus/src/lib.rs`: the NetCDF-3 encoder (`encode_nc3` /
      `build_header`) now returns `Result` and propagates a dimension-product
      overflow as `ExodusError` instead of `.expect("encoded dimensions are
      always valid")` on every Exodus write. `mesh_to_exodus_bytes` /
      `write_exodus` surface the error.

### 9c — Missing / incomplete features

- [x] **P2 (quadratic) elements are wired end-to-end.** The Gmsh importer in
      `tpt-fem-mesh` accepts all six P2 types (`Tri6`/`Quad8`/`Quad9`/`Tet10`/
      `Hex20`/`Hex27`) with `reorder_p2` mapping Gmsh node order → reference
      order, and the assembly-side dispatch in `tpt-fem-{element,thermal,
      elasticity,assembly}` already tables the P2 `shape`/`grad`/`DIM`. A new
      integration test (`p2_tri6_gmsh_import_then_solve`) imports a Tri6 Gmsh
      mesh and runs `solve_poisson` through the full P2 path. Along the way a
      latent bug was fixed: `Map::from_nodes_and_grad` now derives the element
      dimension from the reference-gradient width, so 2-D Gmsh meshes (which carry
      3-component coordinates) solve correctly instead of panicking.
- [x] `tpt-fem-io-abaqus`: `read_inp_deck` now parses `*NSET`/`*ELSET`/
      `*MATERIAL`/`*ELASTIC`/`*BOUNDARY` into an `InpDeck` (node/element sets,
      material `young`/`poisson`, and prescribed boundary DOFs), retaining them
      instead of discarding on geometry import. Regression test
      `deck_captures_sets_material_boundary` guards it.
- [x] `tpt-fem-cli`: added `elasticity` and `modal` (eigen) subcommands that
      reuse the TOML config schema and drive `solve_elasticity` / `solve_modal`;
      both are covered by new CLI tests. **Nonlinear** (`arc_length_continuation`)
      is available in the library but intentionally not exposed as a CLI
      subcommand — it needs a problem-specific residual/Jacobian supplied by Rust
      code, not a config file, so it remains a documented fast-follow.
- [x] `tpt-fem-cli`: `load_mesh` now supports Abaqus `.inp` (via
      `read_inp`) and Exodus `.ex`/`.ex2`/`.e` (via `read_exodus`) in addition
      to `.msh`/`.vtk`, with a round-trip regression test.
- [x] `tpt-fem-io-vtk`: add a real VTK reader to the crate itself (currently
      only an ad-hoc reader lives inside the CLI). **Superseded** by Phase 10c
      (`read_vtk`/`mesh_from_vtk`), which delivered exactly this.
- [ ] Consider additional export formats: CSV, Gmsh writer (currently
      import-only), STL, XDMF/HDF5 — lowest priority, evaluate demand first.
- [x] `tpt-fem-thermal`: `mms_3d_converges` is `#[ignore]`d but now runs on the
      nightly `schedule` (and `workflow_dispatch`) via the `mms-3d` CI job, so it
      exercises instead of bit-rotting.

### 9d — Testing / CI hardening

- [x] `tpt-fem-cli`: now has unit tests — golden end-to-end tests for
      `solve`/`elasticity`/`modal`, `mesh info`, `mesh convert`, an
      Abaqus/Exodus round-trip test for `load_mesh`, and a clap-usage drift guard.
- [x] `tpt-fem-py`: a dedicated `python.yml` workflow runs `maturin develop` +
      `pytest` on Python 3.13 (the crate is excluded from the Cargo workspace, so
      the existing Rust CI never touched it).
- [x] Fuzz suite: a `fuzz-run` CI job (on the nightly `schedule` /
      `workflow_dispatch`) now actually executes all four targets with a bounded
      wall-clock budget, caching `fuzz/corpus` and `fuzz/artifacts` between runs.
- [x] Add `criterion` benchmarks: `tpt-fem-sparse` (`benches/dense_solve.rs`,
      dense-LU solve at n=50/100/200 as the baseline; `solve_russell` can be
      swapped in where the SuiteSparse/MUMPS toolchain exists for the sparse-vs-dense
      comparison) and `tpt-fem-mesh-gen` (`benches/delaunay.rs`, Bowyer–Watson at
      200/500/1000 points and `box_mesh` at 10³/20³/40³).

### 9e — Usability

- [x] Ship example config TOMLs with the CLI: `examples/poisson.toml` and
      `examples/elasticity.toml` document the `solve`/`elasticity`/`modal` schema.
- [x] Results/report summary: every solve subcommand prints a brief report
      (DOF count, solve time, solution range / max displacement / fundamental
      frequency) so results can be sanity-checked without opening ParaView.

### 9f — Revisit periodically

- [ ] `deny.toml` suppresses three RUSTSEC advisories (RUSTSEC-2026-0041/
      0194/0195) for vtkio-transitive deps, justified as "doesn't apply to
      trusted-read/write usage." Reasonable today; revisit if VTK files are
      ever accepted from untrusted/user-uploaded sources.

---

## Phase 10 — Platform Review Follow-ups (2026-08-18)

*A second platform review (bugs/missing features/innovation/usability/
adoption) on 2026-08-18 found the items below, layered on top of the
still-open Phase 9 items above. Nothing here has been implemented yet —
tracked for a future pass.*

### 10a — Bugs / hardening

- [x] `tpt-fem-element`: add a `debug_assert!` + doc note on `mat_det`/
       `mat_inv`/`tensor_nodes` documenting the out-of-range-dimension panic
       contract, since Phase 9b deferred converting them to `Result`.
- [x] Drop the stale `faer` wasm-risk note for `tpt-fem-sparse` — the default
       backend was already swapped to the in-house dense LU (Phase 1d
       follow-up). Fixed the actual stale references in
       `crates/tpt-fem/README.md` (s/`faer` solve/→ in-house dense-LU solve) and
       the `wasm` CI job comment; the root README carried no `faer` note.
- [x] `tpt-fem-io-exodus`: spot-verify Phase 9b's panic-hardening claim
       specifically against the hand-rolled NetCDF-3 codec — done. Verified:
       `encode_nc3` / `build_header` / `vsize_of` / `mesh_to_exodus_bytes` all
       return `Result<_, ExodusError>`; the dimension-product overflow is caught
       by `checked_mul` → `ExodusError::Parse` (no `.expect("encoded dimensions
       are always valid")` remains); `sane_count` guards untrusted counts;
       `write_exodus` surfaces `io::Result`. The only `.unwrap()` left in the
       file are in unit tests (trusted in-memory round-trips).

### 10b — Missing features

- [x] `tpt-fem-cli`: per-component Dirichlet BCs (currently only all-DOF
       fixing is supported) — real limitation for elasticity users, noted
       only in a code comment today. `Bc` now carries an optional `dofs`
       selector (0-based component list; defaults to all components), and
       `expand_dirichlet` constrains only the listed (clamped to `[0, dim)`)
       components. Documented in `examples/{poisson,elasticity}.toml` and
       guarded by `expand_dirichlet_honors_dof_mask`. (Poisson remains
       scalar/ignores `dofs`.)
- [x] `tpt-fem-py`: extend Python bindings to `solve_elasticity`,
       `solve_modal`, and eigen solve — currently Poisson/thermal only, the
       biggest gap between what the Rust core can do and what's reachable
       without writing Rust. `solve_elasticity(mesh, model, young, poisson,
       quad_order, bcs)` and `solve_modal(mesh, model, young, poisson, density,
       quad_order, num_modes, bcs)` are now exposed (BCs are `(node, component,
       value)` tuples; `solve_modal` returns a list of `(ω², mode_shape)`). The
       pre-existing `solve_poisson` callback-error capture was also hardened to
       `Arc<Mutex<…>>` so the closure is `Send` under `allow_threads`. Validated
       with `maturin develop` + `pytest` (4 tests pass on CPython 3.13), and
       documented in the crate README.
- [x] Consider a shared `thiserror`-based base error (or a single
       `tpt_fem::Error` umbrella enum with `#[from]` conversions) to reduce
       the `Box<dyn Error>` glue currently needed by CLI/umbrella consumers.
       Implemented: `tpt-fem` now exposes `tpt_fem::Error` (thiserror-based)
       aggregating each per-crate error (`Mesh`/`Sparse`/`Inp`/`Exodus`/`Vtk`/
       `Newton`, each gated behind its feature) plus `std::io::Error`, TOML
       config errors, and a `Msg` variant. The CLI switched from
       `Box<dyn std::error::Error>` to `tpt_fem::Error`; all its `?` sites
       convert via the `#[from]` impls. Verified with `cargo build --workspace`,
       `--no-default-features` umbrella build, `clippy -D warnings`, and
       `cargo fmt --check` (all clean).

### 10c — Innovative / high-leverage additions

- [x] `tpt-fem-io-vtk`: promote the CLI's ad-hoc VTK reader into the crate
       itself (closes the Phase 9c layering gap and unlocks VTK round-tripping
       from Python/umbrella, not just the CLI). `tpt-fem-io-vtk` now exposes
       `read_vtk(path)` and `mesh_from_vtk(&Vtk)` (returning `Result<Mesh,
       VtkError>`; the `VtkError` enum gained a `Parse` variant). The CLI's
       local `mesh_from_vtk` was removed and `load_mesh`'s `.vtk` branch now
       uses the crate-level reader; `vtkio` was dropped from the CLI deps. A
       round-trip test (`read_vtk_round_trips_mesh`) and the existing
       `mesh_convert` test exercise it.
- [x] `tpt-fem-py`: Jupyter-friendly result objects. `solve_poisson` now
       returns `PoissonSolution`, `solve_elasticity` returns
       `ElasticitySolution`, and `solve_modal` returns `ModalSolution` (indexable
       / iterable over `ModeShape` objects) — instead of bare lists. Each carries
       rich `__repr__`/`_repr_html_` display, `to_numpy()` (returns an `np.ndarray`;
       shape `(n,)` for scalar Poisson, `(n, dim)` for vector fields), and
       `to_pyvista()` (returns a `pyvista.UnstructuredGrid` with the field attached
       as point data `"u"`/`"disp"`/`"mode"`), so results plug into
       PyVista/matplotlib/Jupyter without a manual VTK export round-trip.
       `numpy` is required for `to_numpy`; `pyvista` for `to_pyvista` (both declared
       as the `viz` optional extra in `pyproject.toml`). The crate README got a
       "Visualization & Jupyter" section + updated API table, `walkthrough.py` and
       the pytest suite were updated to the new objects, and two new drift-guard
       tests (`test_to_numpy_shapes`, `test_to_pyvista_round_trips`, the latter
       `importorskip`-gated on `pyvista`) were added. `cargo check` + `cargo clippy
       --all-targets -D warnings` pass on the crate (maturin/pytest can't run here
       since `maturin`/`pyvista` aren't installed on this dev box — the Python
       tests must be re-run via `maturin develop` + `pytest` to fully validate the
       `to_pyvista` path).
- [x] `tpt-fem-cli`: `tpt-fem init` scaffolding subcommand that generates a
       starter `problem.toml` + minimal mesh for a chosen problem type
       (poisson/elasticity/modal). `Init { problem, output }` writes a
       type-specific starter config (Poisson/elasticity/modal templates) to the
       given path; documented in the root README and guarded by
       `init_writes_starter_config` (rejects unknown problem types).
- [x] CI: track `criterion` bench output over time (e.g.
       `github-action-benchmark`-style regression tracking) instead of local/
       ad hoc results only. Added a `bench` CI job (default features; the
       `russell` backend is excluded since it needs the SuiteSparse/MUMPS
       toolchain) that runs the `tpt-fem-sparse` and `tpt-fem-mesh-gen`
       criterion benches and uploads `target/criterion` as a workflow artifact
       for cross-run comparison. (A PR-dashboard `github-action-benchmark`
       wiring is a future enhancement.)

### 10d — Usability / automation

- [x] Add a root `CONTRIBUTING.md` documenting the `Justfile` commands,
       crate-DAG conventions, and the `todo.md`/Phase workflow.
- [x] Extend the `cli_usage_matches_readme` drift-guard pattern to the
       Python README's usage snippet (nothing currently enforces it still
       compiles/runs against the bound API). `crates/tpt-fem-py/tests/
       test_tpt_fem.py::test_readme_snippet` now runs the exact API shown in
       the crate README (Poisson + elasticity + modal) and fails if the bound
       surface changes; the README usage was also corrected (it previously used
       an empty BC list, which makes the Poisson system singular).
- [x] Add a `vtk_import` fuzz target once a real VTK reader lands per 10c.
       Added `fuzz/fuzz_targets/vtk_import.rs` (and registered it in
       `fuzz/Cargo.toml` with a `tpt-fem-io-vtk` dependency). It feeds
       arbitrary bytes to the promoted `tpt_fem_io_vtk::read_vtk` reader via a
       temp file, asserting the reader returns `Err` rather than panicking on
       malformed input. Wired into the CI `fuzz-build`/`fuzz-run` loops.
- [x] Add a gated (manual/`workflow_dispatch`, not auto-run) `maturin
       publish` GitHub Action for `tpt-fem-py` — `pyproject.toml` is already
       PyPI-shaped but nothing publishes it. Added a `publish` job to
       `.github/workflows/python.yml` gated on `workflow_dispatch` that runs
       `maturin publish --non-interactive` with a `PYPI_TOKEN` repo secret
       (`environment: pypi`); it does not run on push/PR.

### 10e — Adoption: examples, templates, onboarding

- [x] Add one runnable Rust example per major capability under
       `crates/tpt-fem/examples/` (`elasticity_frame.rs`, `modal_analysis.rs`,
       `mesh_gen_box.rs`, `abaqus_import.rs`), mirroring the existing
       `thermal_solve.rs`. All four build and run; `elasticity_frame` matches
       the analytical tip deflection and `modal_analysis` yields positive,
       increasing natural frequencies.
- [x] Add a `crates/tpt-fem-py/examples/` directory (or notebook) walking
       through mesh load → solve → VTK export → visualize. Added
       `crates/tpt-fem-py/examples/walkthrough.py`, a runnable end-to-end
       example (Poisson + 3-D elasticity + modal analysis on a clamped bar,
       exporting VTK) validated via `maturin develop`.
- [x] Stand up a minimal docs site (e.g. `mdbook`) or a linked root `docs/`
       folder consolidating the 15 per-crate READMEs into one browsable
       narrative — currently no single entry point exists. Added
       `docs/README.md`: a capability map, a per-crate README index, the
       end-to-end Rust/Python examples, and developer-doc pointers. (A full
       `mdbook` render is a future enhancement; this provides the single
       browsable entry point.)
- [x] Add a "Getting Started" section to the root README showing
       `git clone` → `cargo run --example thermal_solve` → expected output,
       since crates.io/PyPI publishing remain out of scope and every
       adoption path starts with a clone today.

---

## Phase 11 — Platform Review Follow-ups (2026-08-18b)

*A third platform review (bugs/TODOs/missing features) on 2026-08-18,
covering the repo as Phase 10's changes stand uncommitted. Nothing here has
been implemented yet — tracked for a future pass.*

### 11a — Bugs

- [x] `fuzz/Cargo.toml`: the `[dependencies.tpt-fem-io-abaqus]` block is
       present (`fuzz/Cargo.toml:22-23`) and `abaqus_inp.rs` uses
       `tpt_fem_io_abaqus::read_inp`, so `cargo +nightly build --bin
       abaqus_inp` and the `fuzz-build`/`fuzz-run` CI jobs compile. Verified
       2026-08-20.
- [x] `todo.md:451-453` (Phase 9c) records "`tpt-fem-io-vtk`: add a real
       VTK reader" as superseded by Phase 10c (`read_vtk`/`mesh_from_vtk`),
       which is checked off — the two entries are now consistent. Verified
       2026-08-20.
- [x] `deny.toml:7-9,20`: no `RUSTSEC-2024-0436` (`paste`, unmaintained)
       suppression exists in the current `deny.toml` — `paste` does not appear
       in `Cargo.lock` and `faer` is not a dependency anywhere, so the stale
       entry described here is already absent. The only suppressions are the
       three vtkio-transitive advisories (RUSTSEC-2026-0041/0194/0195),
       matching the Phase 9f note. Verified 2026-08-20.

### 11b — Panic-risk hardening (library-facing, not test code)

- [x] `tpt-fem-py/src/lib.rs:146,151,165`: `solve_poisson`'s callback-error
      cell (`Arc<Mutex<Option<PyErr>>>`) is read/written via bare
      `.lock().unwrap()` three times. If the guarded closure ever panics
      while holding the lock, the mutex poisons and the read-back at line
      165 panics instead of surfacing a `PyErr` — low likelihood today
      (nothing between `.lock()` and the assignment can panic currently) but
      worth `unwrap_or_else(PoisonError::into_inner)` for defense-in-depth
      since it's the newest file in the workspace to still panic on purpose.
- [x] The following production `.unwrap()`/`.expect()` sites (all outside
      `#[cfg(test)]`) don't yet have the "documented as deliberately
      deferred" treatment `tpt-fem-element`'s `mat_det`/`mat_inv` got in
      Phase 9b/10a — either give them the same debug_assert!+doc-note
      treatment or convert to `Result`:
      - `tpt-fem-assembly/src/lib.rs:519` — `bcs.iter().find(...).unwrap().1`
      - `tpt-fem-element/src/lib.rs:559,584` —
        `.position(|v| *v == 0.0).unwrap()` (exact float-equality search on
        reference-node coordinates)
      - `tpt-fem-io-abaqus/src/lib.rs:218` — `element_type.unwrap()`
      - `tpt-fem-mesh/src/lib.rs:698` — `.expect("P2 reference coordinate
        must appear in the Gmsh ordering")`
      - `tpt-fem-mesh-gen/src/lib.rs:580-581` — two `.find(...).unwrap()`
        locating the remaining tet faces during Delaunay flips

### 11c — Missing / incomplete features

- [x] `tpt-fem-cli`: the `solve`/`elasticity`/`modal` subcommands
      (`crates/tpt-fem-cli/src/main.rs:38-51`) all dispatch to the same
      `solve_config` (`main.rs:510`), which reads the actual problem kind
      from `cfg.problem.r#type` in the TOML — the subcommand name is
      cosmetic and never validated against the config's declared type. E.g.
      `tpt-fem elasticity poisson.toml` silently solves Poisson instead of
      erroring "expected an elasticity config". Consider a check that
      `cfg.problem.r#type` matches the invoked subcommand and errors
      otherwise.
- [x] `tpt-fem-py`: no `.pyi` type-stub file ships with the crate (checked
      `crates/tpt-fem-py/` — only `src`/`tests`/`examples`/README/
      pyproject.toml, no stub). IDE autocomplete/type-checking (mypy/
      pyright) against the bound classes (`Mesh`, `PoissonSolution`,
      `ElasticitySolution`, `ModalSolution`, `ModeShape`) gets no static
      typing today. Worth a `.pyi` alongside the module once the bound API
      stabilizes.
- [ ] `tpt-fem-solve`/`tpt-fem-eigen`: per `todo.md:303-310`/`:298-302`, the
      arc-length continuation and generalized Lanczos eigensolver are only
      verified against textbook/algebraic cases (an algebraic fold, a 2-bar
      snap-through truss, and standard modal problems) — no test coverage of
      near-singular/ill-conditioned inputs (e.g. a Jacobian close to
      singular away from the tracked fold, or closely-clustered eigenvalues
      stressing the shift-invert Lanczos restart). Low priority unless a
      user hits a real convergence failure, but worth flagging since these
      are the two most numerically delicate solvers in the workspace.

### 11d — Revisit periodically

- [x] `deny.toml` note resolved. The stale `RUSTSEC-2024-0436` (`paste`/
      `faer`) suppression referenced here was never present in the current
      `deny.toml` (only the three vtkio-transitive advisories are suppressed),
      so the count of *three* in the Phase 9f note is accurate and no change
      was needed. (Covered by the Phase 11a review — `paste`/`faer` are
      git-deps, not advisories, and the russell feature is excluded from the
      denied advisory's scope.)

---

## Phase 12 — Advanced Materials, Contact & Multiphysics (spec2 expansion)

*`spec2.txt` (2026-08-18) expands the crate inventory with 8 new crates:
J2 plasticity, continuum damage, hyperelasticity, composites, porous media,
contact mechanics, fluid elements, and multiphysics coupling. Investigation
(2026-08-19) found two infrastructure gaps that block spec2's crates as
scoped: `tpt-fem-assembly`/`tpt-fem-mesh` hard-code a single uniform
`dofs_per_node` with no multi-field concept (blocks `fluid`/`coupling`), and
no time-integration code exists anywhere in the workspace (blocks
`fluid`/`coupling`/`damage`'s element deletion). Both gaps get their own new
crate rather than being folded into existing ones — 9 new crates total
(dofmap, dynamic, plasticity, hyperelastic, composite, porous, contact,
fluid, coupling), one build pass, per the dependency order below. All nine
crates are now implemented, building, and tested (unit tests pass; `cargo
fmt` / `clippy` clean as of 2026-08-20), generalising the Per-Crate Checklist
Template used by every prior phase.

`tpt-math-linalg-complex` (mentioned in spec2's `DEPENDS ON` line) is not a
dependency of anything in this phase — it's a non-blocking forward note for
a future damped-complex-eigenmode extension to `tpt-fem-eigen`, and is
itself still unscaffolded in `tpt-math` (Phase F1).

Note: unlike every crate in Phases 1-11, none of these 10 have (or need) an
entry in `tpt-rust-map/registry.toml` — the same as the existing 16
`tpt-fem-*` crates today, none of which are registered there either
(`tpt-rust-map/repos/tpt-fem/` is intentionally empty per spec.txt's
"graduated, no draft copy" note). Phase 4's closeout line about flipping
registry entries `planned` → `git` "when that repo is next touched" does not
apply in practice; treat it as stale rather than a precedent to follow here.

### 12a — tpt-fem-dofmap

*Multi-field DOF map: per-field `dofs_per_node`, block/interleaved
numbering. Depends on: tpt-fem-mesh.*

- [x] Scaffold `crates/tpt-fem-dofmap/`
- [x] Design a `MultiFieldDofMap` (or similar) that numbers DOFs across N
      named fields, each with its own per-node component count, on one
      `Mesh` — generalizing `Mesh::number_dofs`'s single-field numbering
- [x] Unit tests: two-field map (e.g. 3-component + 1-component) matches
      hand-computed global indices; degenerates to today's single-field
      numbering when there's only one field
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12b — tpt-fem-dynamic

*Time integration: implicit Newmark-beta/HHT-alpha and explicit
central-difference over generic M/C/K. Depends on: tpt-fem-sparse,
tpt-fem-assembly.*

- [x] Scaffold `crates/tpt-fem-dynamic/`
- [x] Implement Newmark-beta implicit time-stepping (reuse
      `tpt-fem-elasticity`'s existing consistent/lumped mass matrices rather
      than duplicating that logic)
- [x] Implement explicit central-difference time-stepping
- [x] Unit tests: single-DOF SDOF oscillator against the closed-form
      solution for both integrators; energy-conservation check on an
      undamped free-vibration case
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12c — tpt-fem-plasticity

*J2 (von Mises) plasticity, kinematic/isotropic hardening, return mapping.
Depends on: tpt-fem-assembly, tpt-fem-solve.*

- [x] Scaffold `crates/tpt-fem-plasticity/`
- [x] Implement J2 return-mapping algorithm (radial return) with
      isotropic + kinematic hardening
- [x] Wire into `tpt-fem-solve`'s Newton loop as a nonlinear constitutive
      update
- [x] Unit tests: single-element uniaxial loading against textbook
      elastic-plastic stress-strain response
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12d — tpt-fem-hyperelastic

*Neo-Hookean, Mooney-Rivlin, Ogden models for soft tissue. Depends on:
tpt-fem-assembly, tpt-fem-solve.*

- [x] Scaffold `crates/tpt-fem-hyperelastic/`
- [x] Implement Neo-Hookean, Mooney-Rivlin, and Ogden strain-energy
      functions + their tangent moduli, large-strain kinematics
- [x] Wire into `tpt-fem-solve`'s Newton loop
- [x] Unit tests: uniaxial/biaxial extension against closed-form
      hyperelastic stress-stretch curves
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12e — tpt-fem-composite

*Classical lamination theory, cohesive zone models, delamination. Depends
on: tpt-fem-elasticity.*

- [x] Scaffold `crates/tpt-fem-composite/`
- [x] Implement classical lamination theory (ABD matrix from ply stack)
- [x] Implement cohesive-zone interface elements for delamination (kept
      independent of `tpt-fem-contact`'s general surface-to-surface
      algorithm — interface elements, not contact search)
- [x] Unit tests: ABD matrix against a textbook symmetric laminate;
      cohesive element traction-separation law sanity checks
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12f — tpt-fem-porous

*Biot's consolidation, Darcy flow, permeability tensors. Depends on:
tpt-fem-assembly, tpt-fem-dynamic.*

- [x] Scaffold `crates/tpt-fem-porous/`
- [x] Implement Darcy-flow element (steady case, reuses Poisson-like
      formulation)
- [x] Implement Biot consolidation (transient, via `tpt-fem-dynamic`)
- [x] Unit tests: Darcy flow against an analytical steady-state solution;
      Terzaghi 1-D consolidation against its closed-form solution
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12g — tpt-fem-contact

*Surface-to-surface contact, penalty/augmented Lagrangian, wear. Depends
on: tpt-fem-assembly, tpt-fem-dofmap, tpt-math-optimize-convex.*

- [x] Scaffold `crates/tpt-fem-contact/`
- [x] Implement surface-pair nearest-node pairing — a brute-force O(|a|·|b|)
       from-scratch scan (written from scratch, no existing wrap target;
       checked `tpt-rust-map/registry.toml`), with a `contact_pairs` helper
       that returns `Option` for an empty surface. A BVH/octree spatial
       accelerator is a future optimisation; the constraint enforcement is
       independent of the search.
- [x] Implement penalty-method contact (self-contained, no multiplier DOFs)
- [x] Implement augmented Lagrangian contact (multiplier DOFs via
      `tpt-fem-dofmap`, constraint solve via `tpt-math-optimize-convex`'s
      dense IPM QP solver)
- [x] Basic wear model
- [x] Unit tests: two-block Hertzian contact against the analytical contact
      pressure/area solution
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12h — tpt-fem-fluid

*Stokes flow, low-Re Navier-Stokes elements (for FSI). Depends on:
tpt-fem-dofmap, tpt-fem-dynamic, tpt-fem-assembly, tpt-fem-solve.*

- [x] Scaffold `crates/tpt-fem-fluid/`
- [x] Implement mixed velocity/pressure (Taylor-Hood-style) elements via
      `tpt-fem-dofmap`
- [x] Implement steady Stokes solve
- [x] Implement transient low-Re Navier-Stokes via `tpt-fem-dynamic`
- [x] Unit tests: lid-driven cavity / Poiseuille flow against known
      solutions
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

### 12i — tpt-fem-coupling

*Thermal-structural, electro-thermal, fluid-structure interface operators.
Depends on: tpt-fem-thermal, tpt-fem-elasticity, tpt-fem-fluid,
tpt-fem-dofmap, tpt-fem-dynamic.*

- [x] Scaffold `crates/tpt-fem-coupling/`
- [x] Implement thermal-structural coupling operator (thermal strain →
      elasticity load)
- [x] Implement electro-thermal coupling operator (Joule heating)
- [x] Implement fluid-structure interface operator (transient, via
      `tpt-fem-dynamic`)
- [x] Unit tests: thermal-structural bimetallic-strip benchmark against a
      known deflection
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean

## Phase 13 — Platform Review Follow-ups (2026-08-20)

*A fourth platform review (bugs/TODOs/missing features/innovation ideas),
covering Phase 12's new experimental crates (dofmap, dynamic, plasticity,
hyperelastic, composite, porous, contact, fluid, coupling). Nothing here has
been implemented yet — tracked for a future pass.*

### 13a — Bugs

- [x] `todo.md:864` (Phase 12g) claimed "Implement surface-pair spatial
       search (BVH/octree)" as done, but
       `crates/tpt-fem-contact/src/lib.rs` implements only a
       brute-force O(|a|·|b|) nearest-node scan — no BVH/octree exists.
       Corrected the Phase 12g checklist entry to describe the brute-force
       scan honestly (a BVH/octree remains a future optimisation).
- [x] `crates/tpt-fem-contact/src/lib.rs` `contact_pairs` returned the
       sentinel `best = usize::MAX` when surface `b` is empty instead of
       `Option`/`Result` — a latent panic/index trap for any caller that
       doesn't special-case an empty surface. It now returns `Option`
       (`None` for an empty surface); `contact_pairs_empty_surface_is_none`
       guards it.
- [x] `crates/tpt-fem-porous/src/lib.rs` `terzaghi_consolidation`
       uses backward-Euler (documented as unconditionally stable) but still
       `assert!`s the *explicit*-scheme stability bound
       `dt <= dz²/(2·cv)`, contradicting its own doc rationale and rejecting
       legitimate large time steps for an implicit method. The assert is
       removed and the doc now states the step is unconditionally stable.
- [x] `crates/tpt-fem-dynamic/src/lib.rs` `central_difference` has no
       CFL/critical-timestep check — an actually-conditionally-stable
       explicit integrator that silently diverges above the stability limit
       instead of erroring. It now returns `Result<_, DynamicError>` and
       rejects (`DynamicError::CflViolation`) steps above `2/ω_max`;
       `central_difference_rejects_over_cfl_step` guards it.
- [x] `crates/tpt-fem-fluid/src/lib.rs` `transient_navier_stokes`'s Picard
       loop runs a fixed `picard_iters` count with no residual/convergence
       check — callers get no signal whether the nonlinear solve actually
       converged. It now breaks early on convergence and returns
       `Err(FluidError::PicardNotConverged)` otherwise. (`steady_stokes`
       likewise returns `Result` — see 13b.)
- [x] `crates/tpt-fem-coupling/src/lib.rs` `fsi_coupling` hardcodes
       the fluid-structure interface normal to `+y` regardless of actual
       interface geometry — wrong traction direction for curved/vertical
       interfaces. It now uses a geometry-aware outward normal (average of
       `node − incident-element-centroid`); still a lumped (smoke-level)
       nodal-load projection, self-documented as such.

### 13b — Panic-risk hardening (Phase 12 crates)

- [x] `crates/tpt-fem-plasticity/src/lib.rs`
       `solve_elastic_plastic_rod` `.expect("Newton should converge")` panics
       for a perfectly-plastic material pushed past yield under force control
       (documented singular-tangent case) instead of returning a `Result`.
       It now returns `Result<_, PlasticityError>` (the `Newton` variant
       carries the `NewtonError`).
- [x] `crates/tpt-fem-hyperelastic/src/lib.rs` `mat_inv` panics on a
       singular deformation gradient `F` (`det.abs() > 1e-14` assertion);
       `neo_hookean_piola`/`mooney_rivlin_piola` call it directly, so a
       degenerate/collapsed element panics the caller. `mat_inv` now returns
       `Option` (and the two Piola functions return `Option`); a collapsed
       element is handled as an error rather than a panic. `solve_hyperelastic_bar`
       now returns `Result<_, HyperelasticError>` instead of `.expect()`-ing
       Newton convergence.
- [x] `crates/tpt-fem-fluid/src/lib.rs` `steady_stokes` and
       `transient_navier_stokes` `.expect()` on `solve_with_dirichlet` —
       panics for degenerate/under-constrained meshes (e.g. no Dirichlet BC
       on a floating fluid domain) instead of returning a `Result`. Both now
       return `Result<_, FluidError>` (the `Sparse` variant carries the
       `SparseError`), consistent with the rest of the library.

### 13c — Housekeeping

- [x] `.github/workflows/spec2.txt` was misfiled inside `.github/workflows/`
       alongside `ci.yml`/`python.yml` — it is the Phase 12 spec document and
       belongs at the repo root next to `spec.txt`. Moved to `spec2.txt` at the
       repo root.

### 13d — Innovation ideas (not scoped, for future consideration)

- [ ] BVH/octree contact spatial search (see 13a) — unlocks larger-scale
      contact problems.
- [ ] Adaptive mesh refinement (AMR) driven by a posteriori error
      estimators — no crate currently does error-driven refinement.
- [ ] Topology optimization module (SIMP/level-set) built on the existing
      `tpt-fem-elasticity` + `tpt-fem-solve` stack.
- [ ] GPU/SIMD-accelerated element assembly and sparse matvec — relevant
      given the repeated per-step matvecs in `tpt-fem-dynamic`/`tpt-fem-fluid`
      (also ties into the `coo_matvec` calling `to_csr()` on every invocation
      performance smell in `tpt-fem-dynamic/src/lib.rs:64-79`).
- [ ] WASM in-browser interactive solve+visualize demo —
      `tpt-fem-quadrature`/`tpt-fem-element`/`tpt-fem-mesh` already build for
      `wasm32-unknown-unknown` per CI.
- [ ] Consistent (non-lumped) FSI load transfer with real per-node interface
      normals, replacing the hardcoded `+y` normal (see 13a).
- [ ] Modal frequency-response coupling: combine `tpt-fem-eigen` modal
      analysis with `tpt-fem-dynamic` Newmark integration for a
      vibration/fatigue frequency-response workflow.
