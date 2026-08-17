# Contributing to tpt-fem

Thanks for your interest in improving the tpt-fem workspace. This document
covers the local workflow, the crate dependency DAG, and the phased `todo.md`
process we use to track work.

## Local workflow

We use a [`Justfile`](Justfile) to wrap the cargo commands. Install
[`just`](https://github.com/casey/just) and run:

```bash
just                 # list recipes
just build           # cargo build  --workspace --all-features
just test            # cargo test   --workspace --all-features
just fmt             # cargo fmt    --all --check
just fmt-fix         # cargo fmt    --all
just clippy          # cargo clippy --workspace --all-targets --all-features -- -D warnings
just deny            # cargo deny   check
just verify          # fmt + clippy + deny + test
```

All four gates (`fmt`, `clippy`, `deny`, `test`) must be green before a change
is considered complete. CI runs the same `just verify` plus a `wasm32` build of
the pure-math crates and a nightly `cargo doc --all-features -D warnings`.

### Other useful recipes

```bash
just example-thermal # run the prelude-only Poisson example
just mms-3d          # run the (heavy, #[ignore]d) 3-D MMS convergence test
just cli             # build the tpt-fem-cli driver
```

### Python bindings

`tpt-fem-py` is excluded from the Cargo workspace (it builds against a fixed
`pyo3 = "0.23"` and Python 3.13). Develop it out-of-tree:

```bash
cd crates/tpt-fem-py
maturin develop
pytest
```

## Crate dependency DAG

The workspace is a layered DAG. Lower layers never depend on higher layers; a
crate may only depend on crates listed above it:

```text
tpt-fem-quadrature   (no internal deps)
tpt-fem-mesh         (mshio; no internal deps)
tpt-fem-element      (tpt-fem-quadrature)
tpt-fem-sparse       (tpt-math-linalg-dense; no internal deps)
tpt-fem-assembly     (tpt-fem-sparse, tpt-fem-element, tpt-fem-mesh)
tpt-fem-thermal      (tpt-fem-assembly)
tpt-fem-io-vtk       (tpt-fem-mesh, vtkio)
tpt-fem-elasticity   (tpt-fem-assembly)
tpt-fem-solve        (tpt-fem-assembly)
tpt-fem-eigen        (tpt-fem-sparse)
tpt-fem-io-abaqus    (tpt-fem-mesh)
tpt-fem-io-exodus    (tpt-fem-mesh)
tpt-fem-mesh-gen     (tpt-fem-mesh)
tpt-fem-cli          (tpt-fem umbrella)
tpt-fem             (umbrella: every constituent behind a feature)
tpt-fem-py          (excluded; depends on tpt-fem)
```

When adding a dependency, keep this layering intact. Cross-layer cycles are not
allowed.

## Umbrella feature flags

The `tpt-fem` umbrella re-exports each constituent behind a Cargo feature
(`quadrature`, `element`, `mesh`, `sparse`, `assembly`, `thermal`, `io-vtk`,
`elasticity`, `solve`, `eigen`, `io-abaqus`, `io-exodus`, `mesh-gen`). All are
enabled by default. New public API that should be reachable from
`tpt_fem::prelude::*` must be added to the umbrella's re-export list behind the
appropriate feature.

## The `todo.md` / Phase workflow

`todo.md` is the single source of truth for what is done and what is planned. It
is organized into Phases (bootstrap, core, assembly/physics, hardening, and the
per-review follow-up phases). Conventions:

- Every crate change follows the **Per-Crate Checklist Template** in `todo.md`
  (scaffold → wire deps → implement → unit/doctests → rustdoc → fmt/clippy →
  deny → flip registry status).
- Tasks are tracked as `- [ ]` (open) and `- [x]` (done). Prefer editing the
  checklist item in place over rewriting it.
- When a task is intentionally deferred (e.g. a known low-priority gap), mark it
  `- [ ]` **Deferred** with a one-line note rather than silently dropping it, so
  the rationale survives for the next review.
- When a follow-up review (Phase 9/10/...) adds new work, append a new phase
  section rather than mutating completed ones.

### Adding a new crate

1. Scaffold `crates/<name>/` with a `Cargo.toml` inheriting `[workspace.package]`
   and `lib.rs` (or `main.rs`).
2. Add the crate to `members` in the root `Cargo.toml` (exclude it instead if it
   needs a divergent toolchain, like `tpt-fem-py` or `fuzz`).
3. Add it to the umbrella crate's features if users should reach it via
   `tpt_fem::prelude`.
4. Write the crate README, run the `just verify` gates, and add the
   `planned` → `git` registry flip (in the sibling `tpt-rust-map` repo).

## License

Contributions are accepted under the same dual `MIT OR Apache-2.0` license as
the project.
