# tpt-fem-hyperelastic

Hyperelastic constitutive models for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

Provides the first Piola–Kirchhoff stress `P = ∂Ψ/∂F` for the incompressible
[Neo-Hookean](https://en.wikipedia.org/wiki/Neo-Hookean_solid),
[Mooney–Rivlin](https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid), and
[Ogden](https://en.wikipedia.org/wiki/Ogden_(material_model)) strain-energy
functions, plus a 1-D bar solver wired into [`tpt-fem-solve`]'s Newton loop.

## Overview

The crate evaluates the first Piola–Kirchhoff stress (the stress measure
conjugate to the deformation gradient `F`) for three classic **incompressible**
hyperelastic models. Each model is a strain-energy density `Ψ` of the deviatoric
stretches; incompressibility (`J = det F = 1`) is enforced by a user-supplied
hydrostatic pressure `p`:

```text
Neo-Hookean   :  Ψ = μ/2 (I₁ − 3)                       P = μ F − p F^{-T}
Mooney-Rivlin :  Ψ = c₁(I₁−3) + c₂(I₂−3)                P = 2c₁F + 2c₂(I₁F − C F) − p F^{-T}
Ogden         :  Ψ = Σ μᵢ/αᵢ (λ₁^{αᵢ}+λ₂^{αᵢ}+λ₃^{αᵢ}−3)  P = Σ μᵢαᵢ λᵢ^{αᵢ−1} nᵢnᵢᵀ − p F^{-T}
```

where `C = FᵀF`, `I₁ = tr(C)`, `λᵢ` are the principal stretches and `nᵢ` the
principal directions (the columns of the rotation `N`). The 3-D routines operate
on a `3×3` deformation gradient; a scalar 1-D incompressible neo-Hookean law
`P = μ(F − F⁻²)` and a [`solve_hyperelastic_bar`] Newton solver are also
provided.

Conventions:

- `F` is a row-major `3×3` matrix (`Mat3`); `P` and `C` share the same layout.
- The models are **incompressible**: the caller supplies the pressure `p`
  (e.g. for traction-free lateral faces in uniaxial stretch,
  `p = μ λ⁻¹` for Neo-Hookean). The crate does not solve for `p` itself.
- Units: `μ`, `c₁`, `c₂`, `μᵢ` in Pa; stretches dimensionless; `P` in Pa.

Out of scope / not implemented: compressible (decoupled-volumetric) forms,
automatic pressure / Lagrange-multiplier enforcement, anisotropy (e.g. fibre
reinforcement), and any model beyond the three listed.

## Installation

```toml
[dependencies]
tpt-fem-hyperelastic = "0.1"
# siblings a caller realistically needs alongside it:
tpt-fem-mesh = "0.1"      # build the bar mesh
tpt-fem-solve = "0.1"     # non-linear Newton loop
tpt-fem-sparse = "0.1"    # sparse assembly fed to the solver
```

## Usage

```rust
use tpt_fem_hyperelastic::neo_hookean_piola;

// Uniaxial stretch λ with traction-free sides (p = μ λ⁻¹ for incompressible
// neo-Hookean): the axial PK1 stress is μ(λ − λ⁻²).
let lam = 1.5;
let mut f = [[0.0; 3]; 3];
f[0][0] = lam;
f[1][1] = lam.powf(-0.5);
f[2][2] = lam.powf(-0.5);
let p = 100.0 * lam.powf(-1.0);
let pk = neo_hookean_piola(&f, 100.0, p);
let want = 100.0 * (lam - lam.powf(-2.0));
assert!((pk[0][0] - want).abs() / want < 1e-12);
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `neo_hookean_uniaxial` | `cargo run -p tpt-fem-hyperelastic --example neo_hookean_uniaxial` | Neo-Hookean uniaxial stretch sweep vs the closed-form nominal stress. |
| `model_comparison` | `cargo run -p tpt-fem-hyperelastic --example model_comparison` | Neo-Hookean vs Mooney-Rivlin vs Ogden stress-stretch curves on one deformation. |
| `bar` | `cargo run -p tpt-fem-hyperelastic --example bar` | `solve_hyperelastic_bar` large-stretch Newton solve. |
| `kinematics` | `cargo run -p tpt-fem-hyperelastic --example kinematics` | `Mat3` demo: deformation gradient → `J`, `C`, invariants. |

## API highlights

| Item | Description |
|------|-------------|
| `Mat3` | Row-major `3×3` matrix type. |
| `mat_mul` / `mat_det` / `mat_inv` / `mat_transpose` | Small `3×3` linear-algebra helpers. |
| `neo_hookean_piola` | First PK stress `P = μ F − p F^{-T}`. |
| `mooney_rivlin_piola` | First PK stress for Mooney–Rivlin. |
| `ogden_piola` / `OgdenTerm` | First PK stress for the (multi-term) Ogden model from principal stretches/directions. |
| `neo_hookean_1d` | Scalar incompressible neo-Hookean nominal stress `μ(F − F⁻²)`. |
| `solve_hyperelastic_bar` | 1-D bar Newton solve to a target end stretch. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-sparse ──► tpt-fem-solve ──► tpt-fem-hyperelastic ──► tpt-fem (optional)
```

## References

- Ogden, R. W., "Large Deformation Isotropic Elasticity — On the Correlation
  of Theory and Experiment for Incompressible Rubberlike Solids", *Proceedings
  of the Royal Society A* **326**, 565–584 (1972) — the Ogden model.
- Holzapfel, G. A., *Nonlinear Solid Mechanics* (Wiley, 2000) — Neo-Hookean and
  Mooney–Rivlin models and the first Piola–Kirchhoff stress.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
