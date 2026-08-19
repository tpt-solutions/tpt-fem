# tpt-fem-plasticity

J2 (von Mises) plasticity for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

Provides the radial-return (closest-point projection) algorithm for associative
J2 plasticity with combined isotropic + kinematic hardening, operating on a
single integration point in Voigt `6`-vector form, together with a force-
controlled 1-D rod solver wired into [`tpt-fem-solve`]'s Newton loop.

## Overview

The crate implements the small-strain, rate-independent J2 (von Mises)
elastoplasticity model with associative flow. It is an *incremental* flow-theory
model: at each integration point the trial stress is the elastic predictor built
from the elastic strain, and if it falls outside the von Mises cylinder the
deviatoric stress is projected back onto the yield surface along its normal (the
radial-return map):

```text
elastic predictor :  σ_tr = C : (ε − εᵖ)
relative stress   :  η = dev(σ_tr) − α
yield function    :  f = √(3/2)·‖η‖ − σ_y(ε̄ᵖ) ,   σ_y(ε̄ᵖ) = σ_y⁰ + H_iso·ε̄ᵖ
consistency       :  Δγ = f / [ 2G·√(3/2) + √(2/3)·(H_iso + H_kin) ]
radial return     :  dev(σ) ← dev(σ_tr) − 2G·Δγ·n̂ ,   n̂ = η/‖η‖
flow rule         :  εᵖ ← εᵖ + Δγ·n̂ ,   α ← α + (2/3)·H_kin·Δγ·n̂
                     ε̄ᵖ ← ε̄ᵖ + √(2/3)·Δγ
```

Because the plastic strain `εᵖ` is carried in [`PlasticState`], unloading from a
plastic state is elastic (slope `E`, or `G` in shear) and reversed cycles behave
correctly — kinematic hardening reproduces the Bauschinger effect. Note that
kinematic hardening enters the consistency denominator as well as translating the
surface: shifting the yield surface consumes part of the trial overshoot even
though it leaves the radius unchanged.

Combined hardening is supported: `H_iso` grows the yield radius and `H_kin`
translates it via the back stress `α`. The 3-D update is available through
[`j2_return_mapping`]; a scalar uniaxial-*stress* analogue [`plastic_1d`] (and its
consistent tangent [`uniaxial_tangent_1d`], `E` elastic and `E·H/(E + H)` plastic)
drives the 1-D [`solve_elastic_plastic_rod`] rod solver, which carries a
per-element equivalent-plastic-strain state through [`tpt-fem-solve`]'s Newton
iteration.

> **Uniaxial stress is not `ε = (ε_x, −ν ε_x, −ν ε_x, 0, 0, 0)` past yield.** That
> strain state is the uniaxial-stress state only while the response is elastic:
> plastic flow is incompressible, so holding the lateral strains at the elastic
> contraction develops spurious lateral stress and an artificial `E/3` hardening
> slope. Use [`plastic_1d`] for a genuine uniaxial-stress path, or solve for the
> lateral strain that makes `σ_y = σ_z = 0` (the `uniaxial_sweep` example does
> both, side by side).

Voigt / tensor conventions:

- Strain and stress are stored as `6`-vectors in the order
  `[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]` (and the stress analogue), using the
  **engineering-shear** convention (no `1/2` factor on `γ`; the shear terms
  enter the Frobenius norm with the `2(γ²)` Voigt doubling). The flow rule
  therefore adds `2·Δγ·n̂` to the shear components of `εᵖ`, so that
  `σ = C:(ε − εᵖ)` reproduces the returned stress exactly.
- The plastic strain `εᵖ`, back stress `α`, and equivalent plastic strain `ε̄ᵖ`
  live in [`PlasticState`].
- Units: stress in Pa, strain dimensionless, Young's modulus and the hardening
  moduli in Pa, Poisson's ratio dimensionless.

Out of scope / not implemented: large-deformation (the return mapping is
small-strain), rate-dependence / viscoplasticity, non-associative flow,
anisotropic or pressure-dependent yield, and thermal coupling. The scalar
[`plastic_1d`] law threads only the (non-negative) `ε̄ᵖ` between calls, so it
reconstructs the signed plastic strain as `εᵖ = ε̄ᵖ·sign(ε)` — exact for monotonic
(proportional) loading, which is the regime the rod solver drives. For reversed
cycles use [`j2_return_mapping`] and carry the full [`PlasticState`].

## Installation

```toml
[dependencies]
tpt-fem-plasticity = "0.1"
# siblings a caller realistically needs alongside it:
tpt-fem-mesh = "0.1"      # build the rod mesh
tpt-fem-solve = "0.1"     # (re-exported use) non-linear Newton loop
tpt-fem-sparse = "0.1"    # sparse assembly fed to the solver
```

## Usage

```rust
use tpt_fem_plasticity::{j2_return_mapping, PlasticState, PlasticityParams};

// Uniaxial-stress strain (lateral strains free): ε_y = ε_z = −ν ε_x.
let p = PlasticityParams::steel();
let strain = [1e-3, -p.poisson * 1e-3, -p.poisson * 1e-3, 0.0, 0.0, 0.0];
let (stress, state) = j2_return_mapping(&p, &strain, &PlasticState::new());
// Elastic: σ_x = E ε_x.
assert!((stress[0] - p.young * 1e-3).abs() < 1.0);

// Past yield (strain-driven) the trial stress is projected back onto the
// surface and `state.plastic_strain` records the irrecoverable part:
let ep = 5e-3; // > σ_y/E ≈ 1.25e-3 for steel()
let big = [ep, -p.poisson * ep, -p.poisson * ep, 0.0, 0.0, 0.0];
let (s2, state2) = j2_return_mapping(&p, &big, &PlasticState::new());
assert!(state2.eps_eq > 0.0, "this strain is plastic");
assert!(
    (s2[0] / p.young - ep).abs() > (stress[0] / p.young - 1e-3).abs(),
    "the elastic strain is smaller than the applied strain"
);
```

A force-controlled rod:

```rust
use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_plasticity::{solve_elastic_plastic_rod, PlasticityParams};

let mut b = MeshBuilder::new();
let n0 = b.add_node(vec![0.0]);
let n1 = b.add_node(vec![1.0]);
b.add_element(CellType::Line, vec![n0, n1]);
let mesh = b.build();
let p = PlasticityParams::steel();
let history = solve_elastic_plastic_rod(&mesh, 1.0, &p, &[1e8, 2e8]);
let u_free = history.last().unwrap()[1]; // displacement of the loaded node
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `uniaxial_sweep` | `cargo run -p tpt-fem-plasticity --example uniaxial_sweep` | Load / unload / reload sweep through `j2_return_mapping`, checked against the closed-form `σ(ε)`. |
| `isotropic_vs_kinematic` | `cargo run -p tpt-fem-plasticity --example isotropic_vs_kinematic` | Isotropic vs kinematic hardening with the Bauschinger effect on reversed loading. |
| `elastic_stiffness_demo` | `cargo run -p tpt-fem-plasticity --example elastic_stiffness_demo` | Voigt `6×6` sanity check: `σ_x = E ε_x`, `τ_xy = μ γ_xy`. |
| `rod` | `cargo run -p tpt-fem-plasticity --example rod` | Force-controlled 1-D rod driven through yield. |

## API highlights

| Item | Description |
|------|-------------|
| `PlasticityParams` | Material parameters (`young`, `poisson`, `yield_stress`, `iso_hardening`, `kin_hardening`); `steel()` preset. |
| `PlasticState` | Per-integration-point state: `plastic_strain` (Voigt `6`), `back_stress` (Voigt `6`) and `eps_eq` (`ε̄ᵖ`). |
| `elastic_stiffness` | Isotropic elastic `C` as a Voigt `6×6` from `E` and `ν`. |
| `j2_return_mapping` | 3-D incremental J2 radial-return: `(total_strain [6], state) → (stress [6], state)` via the elastic predictor `σ_tr = C:(ε − εᵖ)`. |
| `plastic_1d` | Scalar uniaxial-*stress* law: `(total_strain, prev_eps_eq) → (σ, eps_eq, plastic)`. |
| `solve_elastic_plastic_rod` | Force-controlled 1-D rod via Newton, carrying `ε̄ᵖ` per element. |
| `uniaxial_tangent_1d` | Consistent tangent `dσ/dε` for the 1-D law (`E` elastic, `E·H/(E+H)` plastic). |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-sparse ──► tpt-fem-solve ──► tpt-fem-plasticity ──► tpt-fem (optional)
```

## Limitations

- Small-strain only; the return mapping is not valid for finite rotations or
  large stretches.
- No rate dependence, no viscoplasticity, no non-associative or anisotropic
  yield, no thermal coupling.
- The 1-D rod law threads only `ε̄ᵖ` between elements, so within a single
  [`plastic_1d`] call the kinematic back stress starts from zero; for the full
  3-D back-stress evolution use [`j2_return_mapping`].

## References

- Simo, J. C. & Hughes, T. J. R., *Computational Inelasticity* (Springer,
  1998) — the J2 radial-return (closest-point projection) algorithm.
- Belytschko, T., Liu, W. K. & Moran, B., *Nonlinear Finite Elements for
  Continua and Structures* (Wiley, 2000) — von Mises plasticity and return
  mapping.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
