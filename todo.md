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
- [ ] `tpt-fem-io-vtk`: add a real VTK reader to the crate itself (currently
      only an ad-hoc reader lives inside the CLI). **Deferred** — low priority;
      the CLI's reader already covers the `.msh`→`.vtk` round-trip need.
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
