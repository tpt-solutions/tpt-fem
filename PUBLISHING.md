# crates.io publish tracker

Tracks `cargo publish` passes across the workspace, in dependency order.
Out of scope for crates.io: `tpt-fem-py` (PyPI via maturin), `tpt-fem-wasm` and
`tpt-fem-capi` (both `publish = false`, dev-only, excluded from the workspace).

- [Release 0.2.0 (next round)](#release-020--next-round) — **plan, not started**
- [Release 0.1.0 (done)](#release-010--done)

---

## Release 0.2.0 — next round

Baseline: the 0.1.0 publish. **Verified against crates.io (2026-10-05):** the
published 0.1.0 tarballs of all 27 crates are source-identical (ignoring line
endings) to commit `9a5da2f`, so `git diff 9a5da2f..HEAD` is exactly what is
unreleased. Pre-release code work is tracked in `todo.md` Phase 16.

### 1. Scope decisions

| Item | Decision | Why |
|------|----------|-----|
| `tpt-fem-wasm` (new) | stay unpublished | `publish = false`; browser demo, built with wasm-pack; CI job `wasm-demo` |
| `tpt-fem-capi` (new) | stay unpublished | `publish = false`; shipped as prebuilt libs + header via `.github/workflows/capi-release.yml` (tag `capi-v*`), which suits a C library better than crates.io |
| `tpt-fem-py` | PyPI only | bump its `tpt-fem` pin to the released version; `maturin publish` workflow already exists |
| MSRV | **1.85** | the old `rust-version = "1.75"` could not build the dependency tree; now enforced by the CI `msrv` job |

**Versioning rule (0.x semver).** Breaking public-API change → minor bump
(0.1 → 0.2). Additive change or fix → patch (0.1 → 0.1.1). A crate whose own
source did not change is **not** republished unless a type from a
breaking-bumped dependency appears in its public API — a dependency pinned
`^0.1.0` simply resolves to the published release. A crate whose only change is
bumping an internal pin (no public API change) gets a patch release so the
registry matches the repo.

> `cargo semver-checks` is run in CI as an advisory job but **cannot see
> return-type changes** (the main break in this release), so the bump table
> below was derived by diffing public signatures by hand
> (`git diff 9a5da2f..HEAD`). Re-derive it if more API changes land.

### 2. Bump table (derived from the code, not guessed)

| Crate | Why | New version |
|-------|-----|-------------|
| tpt-fem-quadrature | Keast4 table fix; new `try_gauss_legendre*`, `QuadratureError` (additive) | **0.1.1** |
| tpt-fem-sparse | new `solve_cg` (additive) | **0.1.1** |
| tpt-fem-eigen | `lanczos_eigs` → `Result` (breaking) + Lanczos fixes | **0.2.0** |
| tpt-fem-composite | internal `unwrap` removal only | **0.1.1** |
| tpt-fem-element | `line/quad/hex_rule` → `Result` (breaking) | **0.2.0** |
| tpt-fem-mesh | pin `tpt-fem-element` 0.2 (private use only); new `to_msh_string`, `nodal_csv`, `ExportError` (additive) | **0.1.1** |
| tpt-fem-mesh-gen | internal `unwrap` removal only | **0.1.1** |
| tpt-fem-io-abaqus | untrusted-input `unwrap` fix | **0.1.1** |
| tpt-fem-io-exodus | `connect*` suffix now rejected, not defaulted | **0.1.1** |
| tpt-fem-assembly | `apply_neumann_order`/`apply_robin_order` → `Result` (breaking) | **0.2.0** |
| tpt-fem-elasticity | mass/body-vector fns → `Result` (breaking) | **0.2.0** |
| tpt-fem-thermal | `poisson_*` → `Result` (breaking); new transient heat solver | **0.2.0** |
| tpt-fem-contact | `augmented_lagrangian` → `Result`, new `ContactError` (breaking) | **0.2.0** |
| tpt-fem-dynamic | `newmark` → `Result` (breaking) | **0.2.0** |
| tpt-fem-fluid | `stokes_dofmap`, `transient_stokes` → `Result`; new error variants (breaking) | **0.2.0** |
| tpt-fem-modal | `modal_superposition` → `Result` (breaking) | **0.2.0** |
| tpt-fem-porous | `PorousError`; `solve_darcy`/`terzaghi_consolidation` signatures (breaking) | **0.2.0** |
| tpt-fem-topopt | pin bumps only (no public API change) | **0.1.1** |
| tpt-fem-coupling | `fsi_interface_loads` → `Result` (breaking) | **0.2.0** |
| tpt-fem | umbrella re-exports all of the above; `Error` gains variants | **0.2.0** |
| tpt-fem-cli | new `topopt` subcommand, `mesh convert` formats; follows umbrella | **0.2.0** |
| tpt-fem-amr, -dofmap, -hyperelastic, -io-vtk, -plasticity, -solve | no source change, no breaking dependency in their public API | unchanged (skip) |

21 crates are republished (the first plan, which bumped every dependent, was
23). Migration notes for users: [`docs/MIGRATING-0.2.md`](docs/MIGRATING-0.2.md).

### 3. Publish order & status

Topological order. Wait for each crate to appear in the index
(`cargo search <name>` / crates.io page) before publishing its dependents.

| # | Crate | → | Status |
|---|-------|---|--------|
| 1 | tpt-fem-quadrature | 0.1.1 | ⬜ |
| 2 | tpt-fem-sparse | 0.1.1 | ⬜ |
| 3 | tpt-fem-eigen | 0.2.0 | ⬜ |
| 4 | tpt-fem-composite | 0.1.1 | ⬜ |
| 5 | tpt-fem-element | 0.2.0 | ⬜ |
| 6 | tpt-fem-mesh | 0.1.1 | ⬜ |
| 7 | tpt-fem-mesh-gen | 0.1.1 | ⬜ |
| 8 | tpt-fem-io-abaqus | 0.1.1 | ⬜ |
| 9 | tpt-fem-io-exodus | 0.1.1 | ⬜ |
| 10 | tpt-fem-assembly | 0.2.0 | ⬜ |
| 11 | tpt-fem-elasticity | 0.2.0 | ⬜ |
| 12 | tpt-fem-thermal | 0.2.0 | ⬜ |
| 13 | tpt-fem-contact | 0.2.0 | ⬜ |
| 14 | tpt-fem-dynamic | 0.2.0 | ⬜ |
| 15 | tpt-fem-fluid | 0.2.0 | ⬜ |
| 16 | tpt-fem-modal | 0.2.0 | ⬜ |
| 17 | tpt-fem-porous | 0.2.0 | ⬜ |
| 18 | tpt-fem-topopt | 0.1.1 | ⬜ |
| 19 | tpt-fem-coupling | 0.2.0 | ⬜ |
| 20 | tpt-fem | 0.2.0 | ⬜ |
| 21 | tpt-fem-cli | 0.2.0 | ⬜ |
| — | tpt-fem-py | 0.2.0 | 🚫 PyPI (`maturin publish`, after crates.io) |
| — | tpt-fem-wasm, tpt-fem-capi | 0.1.0 | 🚫 `publish = false` |

### 4. Pre-flight checklist

Code work (Phase 16 of `todo.md`) is done; what remains is mechanical:

- [ ] Rename each republished crate's `[Unreleased]` → `[<version>] - <date>`
      in its `CHANGELOG.md` and add the `[x.y.z]:` link at the bottom.
- [ ] Bump `version` in each republished crate's `Cargo.toml` per §2.
- [ ] Bump the matching entries in root `Cargo.toml` `[workspace.dependencies]`
      (quadrature/sparse/mesh... patch pins may stay `"0.1.0"`, which still
      resolves to the new patch; breaking crates **must** move to `"0.2.0"`,
      and `tpt-fem-element`'s dependency on quadrature should be `"0.1.1"`
      because it needs `try_gauss_legendre`).
- [ ] Update the `tpt-fem*` pins in `tpt-fem-py`, `tpt-fem-capi`, `tpt-fem-wasm`
      (checked by `just pins` / the CI `pins` job).
- [ ] Update version snippets in `README.md` and per-crate READMEs
      (`grep -rn '0\.1\.0'`).
- [ ] `cargo update -w` and commit the refreshed `Cargo.lock` files.
- [ ] Delete or fold in the stale `crates-publish-order.md`.
- [ ] Run `just release-check` (below) and the manual `fuzz` workflow once.

### 5. Verification gates

`just release-check` runs everything that can run locally: `fmt`, `clippy
-D warnings`, `cargo deny`, `cargo test`, the pin check, an MSRV (1.85) build,
`cargo package --list` for every publishable crate, rustdoc with warnings
denied, the manifest-drift test, and builds of capi / py / wasm. CI runs the same
set plus `cargo semver-checks` (advisory), `cargo machete`, and (on PRs) the
CHANGELOG guard.

Per crate, just before publishing: `cargo publish -p <crate> --dry-run`
(dependents only dry-run once their dependencies are live on crates.io).

### 6. Publish

```sh
cargo publish -p tpt-fem-quadrature     # then follow the table in §3
```

- Updates to **existing** crates are not subject to the new-crate rate limit, so
  no cooldown waits are expected.
- If a crate fails (e.g. metadata validation like the 0.1.0 keyword-length
  rejection), fix, commit, and resume from that row; mark each row ✅ as it lands.
- A published version can't be overwritten — fix forward with a patch bump
  (`cargo yank` only if truly broken).

### 7. Post-release

- [ ] `git tag v0.2.0` (+ per-crate tags if wanted), push, and create a GitHub
      release from the CHANGELOGs (link `docs/MIGRATING-0.2.md`).
- [ ] Add fresh empty `[Unreleased]` sections to each CHANGELOG.
- [ ] `maturin publish` `tpt-fem-py` to PyPI (workflow_dispatch `publish` job).
- [ ] Tag `capi-v0.2.0` to produce the prebuilt C libraries.
- [ ] Confirm docs.rs builds for all 21 republished crates.

### Open items deliberately left for a later release

- CLI / Python exposure of `dynamic`, `plasticity`, `hyperelastic`, `fluid`,
  `porous`, `contact`, `coupling`, `composite` (library-only today).
- GPU/SIMD assembly; additional export formats beyond CSV/Gmsh (STL, XDMF).
- A native sparse *direct* solver (the new `solve_cg` covers SPD systems only).

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
