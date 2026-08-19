# tpt-fem-composite

Classical lamination theory and cohesive-zone delamination models for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

This crate implements the standard first-order (classical) laminated-plate
theory used to compute the effective in-plane stiffness of a stacked fibre
composite. [`laminate_abd`] assembles the `6×6` extensional–bending–coupling
matrix from a bottom-to-top stack of orthotropic [`Ply`] layers:

```text
{N}   [A  B] {ε⁰}
{M} = [B  D] {κ}
```

where `{N}`, `{M}` are the in-plane force and moment resultants, `{ε⁰}` the
mid-surface strains, `{κ}` the curvatures, and `A`/`B`/`D` the
extensional / coupling / bending stiffness blocks (units `N/m`, `N`, `N·m`
respectively for a plate of unit width). A stack that is mirrored about the
mid-plane (e.g. `[0/90]s`) yields `B = 0`: extension and bending are decoupled.

[`CohesiveLaw`] provides a bilinear traction–separation model for
interface (cohesive-zone) delamination elements. The normal traction rises
linearly to a peak `σ̄` at `δc`, then drops linearly to zero at `δf` so the
enclosed area equals the fracture toughness `Gc`. It is intentionally kept
independent of `tpt-fem-contact`'s surface-to-surface algorithm.

Sign / unit conventions:

- Elastic moduli `E₁`, `E₂`, shear `G₁₂` in pascals (Pa); ply `thickness` in
  metres (m).
- `angle_deg` is the fibre orientation relative to the laminate x-axis, in
  degrees, using the standard engineering (counter-clockwise) convention.
- The laminate normal `z` is measured from the mid-plane; `z = −h/2` at the
  bottom of the stack and `+h/2` at the top.
- All ABD entries follow the row-major `6×6` ordering
  `[A B; B D]` with the in-plane indices `0..3` and curvature indices `3..6`.

Out of scope / not implemented:

- No full three-dimensional shell or plate finite-element discretisation (this
  crate produces the constitutive matrix only).
- No geometrically nonlinear or higher-order shear-deformation theory (only
  classical lamination theory).
- No mixed-mode (I/II/III) cohesive laws; [`CohesiveLaw`] is a
  mode-independent scalar (single effective opening) bilinear model.
- No temperature / moisture or viscoelastic effects.

## Installation

```toml
[dependencies]
tpt-fem-composite = "0.1"
# (optional) feed the ABD matrix into an elasticity solver, e.g.:
# tpt-fem-elasticity = "0.1"
```

## Usage

```rust
use tpt_fem_composite::{laminate_abd, Ply};

// A single 1 mm isotropic ply at 0°: A11 should equal Q11 * t.
let ply = Ply {
    e1: 70e9,
    e2: 70e9,
    nu12: 0.3,
    g12: 70e9 / 2.6,
    thickness: 1e-3,
    angle_deg: 0.0,
};
let abd = laminate_abd(&[ply]);
let nu = 0.3;
let q11 = 70e9 / (1.0 - nu * nu);
assert!((abd[0][0] - q11 * 1e-3).abs() / (q11 * 1e-3) < 1e-9);
```

A bilinear cohesive law is built from the toughness `Gc`, peak traction and a
prescribed critical opening:

```rust
use tpt_fem_composite::CohesiveLaw;

let law = CohesiveLaw::from_toughness(1.0, 10.0, 0.1);
assert!((law.traction(law.critical_opening) - 10.0).abs() < 1e-9);
assert!((law.toughness() - 1.0).abs() < 1e-9);
```

## Examples

- `laminate_abd` — builds the ABD matrix of a symmetric cross-ply `[0/90]s`
  stack, prints the `A`/`B`/`D` blocks, and asserts `B ≈ 0`.

  ```text
  cargo run -p tpt-fem-composite --example laminate_abd
  ```

- `unsymmetric_laminate` — builds an unsymmetric `[0/90]` stack and confirms
  the extension–bending coupling block `B` is non-zero.

  ```text
  cargo run -p tpt-fem-composite --example unsymmetric_laminate
  ```

- `cohesive_law` — sweeps the bilinear traction–separation curve and asserts
  the numerically integrated area recovers the fracture toughness `Gc`.

  ```text
  cargo run -p tpt-fem-composite --example cohesive_law
  ```

## API highlights

| Item | Description |
|------|-------------|
| `Ply` | A single orthotropic ply (moduli, Poisson ratio, thickness, orientation). |
| `laminate_abd` | Assembles the `6×6` `A`/`B`/`D` stiffness of a ply stack (row-major). |
| `CohesiveLaw` | Bilinear traction–separation law for a delamination interface. |
| `CohesiveLaw::from_toughness` | Builds the law from `Gc`, peak traction and `δc`; `δf` is derived. |
| `CohesiveLaw::traction` | Traction for an effective opening (monotonic loading branch). |
| `CohesiveLaw::toughness` | Fracture toughness (area under the traction–separation curve). |

## Position in the crate stack

```text
tpt-fem-elasticity ──► tpt-fem-composite
```

## References

- R. M. Jones, *Mechanics of Composite Materials*, 2nd ed., Taylor & Francis,
  1999 — classical lamination theory and the `A`/`B`/`D` construction.

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
