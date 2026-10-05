# crates.io publish tracker

Tracks `cargo publish` passes across the workspace, in dependency order.
Out of scope for crates.io: `tpt-fem-py` (PyPI via maturin), `tpt-fem-wasm` and
`tpt-fem-capi` (both `publish = false`, dev-only, excluded from the workspace).

- [Release 0.2.0 (next round)](#release-020--next-round) — **plan, not started**
- [Release 0.1.0 (done)](#release-010--done)

---

## Release 0.2.0 — next round

Baseline: the 0.1.0 publish (commit `9a5da2f`, 2026-08-23). Everything since
(`git log 9a5da2f..HEAD`) is unreleased.

### 1. Scope decisions

| Item | Decision | Why |
|------|----------|-----|
| `tpt-fem-wasm` (new) | stay unpublished | `publish = false`; browser demo, built with wasm-pack |
| `tpt-fem-capi` (new) | stay unpublished (for now) | `publish = false`, "dev-only this pass". To publish later: remove the flag, add `keywords`/`categories`, depend on the released `tpt-fem` 0.2, and make sure the committed cbindgen header is in the package (`cargo package --list`) |
| `tpt-fem-py` | PyPI only | bump its `tpt-fem` pin to 0.2.0 so it builds against the release |

**Versioning rule (0.x semver):** breaking API change → minor bump (0.1 → 0.2);
additive/fix → patch (0.1.1). A dependency bumped to 0.2.0 can't be satisfied by
`^0.1.0`, so every crate that depends on one must be bumped and republished too —
otherwise the graph ends up with two copies of e.g. `tpt-fem-mesh` and the types
don't match. Simplest policy: **every republished dependent goes to 0.2.0.**

### 2. What changed and the bump it needs

Root causes of the cascade (from CHANGELOGs / API diff against `9a5da2f`):

| Crate | Change | Bump |
|-------|--------|------|
| tpt-fem-quadrature | Keast4 tet table bug fix; new `try_gauss_legendre*` + `QuadratureError` | **0.1.1** (additive + fix) |
| tpt-fem-eigen | Lanczos twice-is-enough reorthogonalization (behavioural fix only) | **0.1.1** |
| tpt-fem-element | `line_rule`/`quad_rule`/`hex_rule` now return `Result` | **0.2.0** breaking |
| tpt-fem-assembly | `apply_neumann_order`/`apply_robin_order` return `Result` | **0.2.0** breaking |
| tpt-fem-thermal | `poisson_element_matrix`/`poisson_source_vector` return `Result` | **0.2.0** breaking |
| tpt-fem-elasticity | `elasticity_body_vector`/`_mass_matrix`/`_lumped_mass` return `Result` | **0.2.0** breaking |
| tpt-fem-contact | new `ContactError`; augmented-Lagrangian returns `Result`; octree search | **0.2.0** breaking |
| tpt-fem-coupling | `fsi_interface_loads` path now returns `Result<_, CouplingError>` | **0.2.0** breaking — ⚠ no CHANGELOG entry yet |
| tpt-fem-io-exodus | `mesh_to_exodus_bytes` returns `Result` | **0.2.0** breaking |
| tpt-fem-cli | new `topopt` subcommand (enables `tpt-fem/topopt`) | **0.2.0** (follows `tpt-fem`) |

Not changed and **not republished**: tpt-fem-sparse, tpt-fem-solve, tpt-fem-amr,
tpt-fem-composite (no internal deps on a bumped crate, no source changes).

### 3. Publish order & status

Topological order (checked with `cargo metadata`). Wait for each crate to appear
in the index (`cargo search <name>` / crates.io page) before publishing its
dependents; dry-runs of dependents only pass once their deps are live.

| # | Crate | 0.1.0 → | Status | Notes |
|---|-------|---------|--------|-------|
| 1 | tpt-fem-quadrature | 0.1.1 | ⬜ | |
| 2 | tpt-fem-eigen | 0.1.1 | ⬜ | |
| 3 | tpt-fem-element | 0.2.0 | ⬜ | dep `tpt-fem-quadrature = "0.1.1"` (needs `try_gauss_legendre`) |
| 4 | tpt-fem-mesh | 0.2.0 | ⬜ | follows element |
| 5 | tpt-fem-dofmap | 0.2.0 | ⬜ | follows mesh |
| 6 | tpt-fem-mesh-gen | 0.2.0 | ⬜ | follows mesh |
| 7 | tpt-fem-io-vtk | 0.2.0 | ⬜ | follows mesh |
| 8 | tpt-fem-io-abaqus | 0.2.0 | ⬜ | follows mesh |
| 9 | tpt-fem-io-exodus | 0.2.0 | ⬜ | breaking + follows mesh |
| 10 | tpt-fem-assembly | 0.2.0 | ⬜ | breaking |
| 11 | tpt-fem-hyperelastic | 0.2.0 | ⬜ | follows mesh |
| 12 | tpt-fem-plasticity | 0.2.0 | ⬜ | follows mesh |
| 13 | tpt-fem-elasticity | 0.2.0 | ⬜ | breaking |
| 14 | tpt-fem-thermal | 0.2.0 | ⬜ | breaking |
| 15 | tpt-fem-contact | 0.2.0 | ⬜ | breaking |
| 16 | tpt-fem-dynamic | 0.2.0 | ⬜ | follows assembly/elasticity |
| 17 | tpt-fem-fluid | 0.2.0 | ⬜ | follows assembly/dynamic |
| 18 | tpt-fem-modal | 0.2.0 | ⬜ | follows dynamic |
| 19 | tpt-fem-porous | 0.2.0 | ⬜ | follows assembly/dynamic |
| 20 | tpt-fem-topopt | 0.2.0 | ⬜ | follows elasticity |
| 21 | tpt-fem-coupling | 0.2.0 | ⬜ | breaking |
| 22 | tpt-fem | 0.2.0 | ⬜ | umbrella; bump pins for every crate above |
| 23 | tpt-fem-cli | 0.2.0 | ⬜ | last |
| — | tpt-fem-sparse, -solve, -amr, -composite | 0.1.0 | ➖ unchanged | skip |
| — | tpt-fem-py | 0.2.0 | 🚫 PyPI | build/upload via maturin after crates.io is done |
| — | tpt-fem-wasm, tpt-fem-capi | 0.1.0 | 🚫 `publish = false` | not published |

> Before bumping, double-check the "unchanged" four with
> `git diff 9a5da2f HEAD -- crates/tpt-fem-{sparse,solve,amr,composite}` (the
> `git diff --stat` currently shows no source changes for them).

### 4. Pre-flight checklist

- [ ] **Resolve stale CHANGELOG sections.** `contact`, `eigen` and `io-exodus`
      already had `[Unreleased]` sections in the 0.1.0 commit. Work out which
      entries actually shipped in 0.1.0 (compare against the crates.io source,
      e.g. docs.rs "source" view) and which are new; file them under `[0.1.0]` or
      the new version accordingly.
- [ ] Add the missing `tpt-fem-coupling` CHANGELOG entry (Result-returning FSI path).
- [ ] Rename `[Unreleased]` → `[<version>] - <date>` in each republished crate; add
      a short "dependency bump" line for follow-only crates (mesh, dofmap, io-*,
      etc.), and cli's `topopt` entry.
- [ ] Bump `version` in each republished crate's `Cargo.toml`.
- [ ] Bump the matching entries in root `Cargo.toml` `[workspace.dependencies]`
      (all currently `"0.1.0"`), including `tpt-fem-quadrature = "0.1.1"` /
      `tpt-fem-eigen = "0.1.1"`.
- [ ] Bump `tpt-fem` pins in `tpt-fem-py`, `tpt-fem-capi`, `tpt-fem-wasm`
      (`tpt-fem-topopt` pin in wasm) — these live outside the workspace, so
      `cargo` won't flag them.
- [ ] Update version snippets in `README.md` and per-crate READMEs (`grep -rn '0\.1\.0'`).
- [ ] `cargo update -w` to refresh `Cargo.lock`.
- [ ] Remove or merge the stale `crates-publish-order.md` (all boxes unchecked,
      duplicates this file).

### 5. Verification gates (all green before the first `cargo publish`)

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo test -p tpt-fem --test manifest_drift   # version pins consistent
cargo deny check
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo semver-checks check-release              # optional; confirms breaking vs patch
# non-workspace crates still build against the bumped pins:
cargo build --manifest-path crates/tpt-fem-capi/Cargo.toml
cargo build --manifest-path crates/tpt-fem-wasm/Cargo.toml --target wasm32-unknown-unknown
cargo check --manifest-path crates/tpt-fem-py/Cargo.toml
```

Per crate, before publishing: `cargo package -p <crate> --list` (nothing stray —
no logs/`.vtk` files) and `cargo publish -p <crate> --dry-run` (dependents can only
dry-run after their deps are live).

### 6. Publish

```sh
cargo publish -p tpt-fem-quadrature     # then follow table order in §3
```

- Updates to **existing** crates aren't subject to the new-crate rate limit, so no
  cooldown waits are expected (the 0.1.0 pass hit them only for brand-new crates).
- If a crate fails (e.g. metadata validation like the keyword-length rejection in
  0.1.0), fix, commit, bump nothing already published, and resume from that row.
- Mark each row ✅ as it lands. A published version can't be overwritten — fix
  forward with a patch bump (`cargo yank` only if truly broken).

### 7. Post-release

- [ ] `git tag v0.2.0` (+ per-crate tags for the 0.1.1 patches if wanted) and push;
      create a GitHub release from the CHANGELOGs.
- [ ] Add fresh empty `[Unreleased]` sections to each CHANGELOG.
- [ ] Build and upload `tpt-fem-py` to PyPI (maturin) at 0.2.0.
- [ ] Confirm docs.rs builds succeeded for all 23 crates.
- [ ] Add a short note to the 0.1.x users: the `Result`-returning API migration
      (quadrature order outside `1..=5` no longer panics).

### Open questions

1. Should `tpt-fem-capi` go on crates.io (it's useful mainly as a C library, so
   GitHub releases with prebuilt `.dll`/`.so`/header may fit better)?
2. Is the `tpt-fem-coupling` change intended as breaking? (Assumed yes.)
3. Lockstep 0.2.0 vs. leaving follow-only crates (dofmap, io-vtk, mesh-gen…) to
   the minimum — lockstep is assumed here; it's simpler and avoids duplicate
   `tpt-fem-mesh` versions in the graph.

---

## Release 0.1.0 — done

First-ever `cargo publish` pass (2026-08), in dependency order.

| # | Crate | Status | Notes |
|---|-------|--------|-------|
| 1 | tpt-fem-quadrature | ✅ published | |
| 2 | tpt-fem-sparse | ✅ published | |
| 3 | tpt-fem-composite | ✅ published | |
| 4 | tpt-fem-element | ✅ published | |
| 5 | tpt-fem-eigen | ✅ published | |
| 6 | tpt-fem-solve | ✅ published | |
| 7 | tpt-fem-amr | ✅ published | hit crates.io new-crate rate limit once, retried after cooldown |
| 8 | tpt-fem-mesh | ✅ published | hit rate limit once, auto-retried |
| 9 | tpt-fem-dofmap | ✅ published | |
| 10 | tpt-fem-mesh-gen | ✅ published | |
| 11 | tpt-fem-io-vtk | ✅ published | |
| 12 | tpt-fem-io-abaqus | ✅ published | |
| 13 | tpt-fem-io-exodus | ✅ published | |
| 14 | tpt-fem-assembly | ✅ published | |
| 15 | tpt-fem-hyperelastic | ✅ published | |
| 16 | tpt-fem-plasticity | ✅ published | |
| 17 | tpt-fem-elasticity | ✅ published | |
| 18 | tpt-fem-contact | ✅ published | |
| 19 | tpt-fem-thermal | ✅ published | |
| 20 | tpt-fem-dynamic | ✅ published | |
| 21 | tpt-fem-fluid | ✅ published | |
| 22 | tpt-fem-modal | ✅ published | |
| 23 | tpt-fem-porous | ✅ published | |
| 24 | tpt-fem-topopt | ✅ published | initial attempt rejected by crates.io for the 21-char keyword `topology-optimization`; fixed to `topology-opt` in df2c3fa |
| 25 | tpt-fem-coupling | ✅ published | |
| 26 | tpt-fem | ✅ published | umbrella crate |
| 27 | tpt-fem-cli | ✅ published | |
| — | tpt-fem-py | 🚫 not published | PyPI/maturin, not crates.io |

All crates were v0.1.0. crates.io enforces a new-crate creation rate limit
(~1 new crate per ~10 min once the initial burst allowance is used), which the
publish script waited out automatically between crates.
