# crates.io publish tracker — 0.1.0 first release

Tracks the first-ever `cargo publish` pass across the workspace, in dependency order.
`tpt-fem-py` is excluded (PyO3/maturin bindings, not part of the cargo workspace,
shipped to PyPI separately, not crates.io).

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
| 24 | tpt-fem-topopt | ✅ published | initial attempt rejected by `crates.io` for the 21-char keyword `topology-optimization`; fixed to `topology-opt` in df2c3fa |
| 25 | tpt-fem-coupling | ✅ published | |
| 26 | tpt-fem | ✅ published | umbrella crate |
| 27 | tpt-fem-cli | ✅ published | |
| — | tpt-fem-py | 🚫 not published | excluded from this pass (PyPI/maturin, not crates.io) |

All crates are v0.1.0 for this first release. crates.io enforces a new-crate
creation rate limit (~1 new crate per ~10 min once the initial burst allowance
is used), which the publish script waits out automatically between crates.
