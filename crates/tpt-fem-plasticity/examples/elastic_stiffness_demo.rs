//! Voigt-convention sanity check for `elastic_stiffness`.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-plasticity --example elastic_stiffness_demo
//! ```
//!
//! Builds the isotropic 6×6 Voigt stiffness `C` and verifies two things: a
//! uniaxial-stress strain produces `σ_x = E ε_x` with traction-free lateral
//! faces (`σ_y = σ_z = 0`), and a pure engineering-shear strain `γ` produces
//! `τ = μ γ` (no factor-of-two surprise), confirming the Voigt convention used
//! throughout the crate.

use tpt_fem_plasticity::{elastic_stiffness, PlasticityParams};

fn main() {
    let p = PlasticityParams::steel();
    let e = p.young;
    let nu = p.poisson;
    let mu = e / (2.0 * (1.0 + nu));
    let c = elastic_stiffness(e, nu);

    // Matrix-vector product C·v for a Voigt 6-vector.
    let mv = |c: &[[f64; 6]; 6], v: &[f64; 6]| -> [f64; 6] {
        let mut o = [0.0; 6];
        for i in 0..6 {
            for j in 0..6 {
                o[i] += c[i][j] * v[j];
            }
        }
        o
    };

    // 1) Uniaxial-stress strain: ε = [ε_x, −ν ε_x, −ν ε_x, 0, 0, 0].
    let eps = 1e-3;
    let strain = [eps, -nu * eps, -nu * eps, 0.0, 0.0, 0.0];
    let s = mv(&c, &strain);
    println!(
        "Uniaxial strain ε_x = {:.3e}  ->  σ_x = {:.6e} Pa  (E ε_x = {:.6e})",
        eps,
        s[0],
        e * eps
    );
    println!(
        "  σ_y = {:.3e}, σ_z = {:.3e}  (should be ~0, traction-free lateral faces)",
        s[1], s[2]
    );
    assert!((s[0] - e * eps).abs() < 1.0, "axial stress wrong");
    assert!(
        s[1].abs() < 1.0 && s[2].abs() < 1.0,
        "lateral faces not traction-free"
    );

    // 2) Pure engineering shear γ on the xy component.
    let g = 2e-3;
    let shear = [0.0, 0.0, 0.0, g, 0.0, 0.0];
    let t = mv(&c, &shear);
    println!(
        "\nPure shear γ_xy = {:.3e}  ->  τ_xy = {:.6e} Pa  (μ γ = {:.6e})",
        g,
        t[3],
        mu * g
    );
    assert!((t[3] - mu * g).abs() < 1e-6, "shear stress wrong");
    assert!(
        t[0].abs() < 1e-9
            && t[1].abs() < 1e-9
            && t[2].abs() < 1e-9
            && t[4].abs() < 1e-9
            && t[5].abs() < 1e-9,
        "spurious tractions from pure shear"
    );

    println!(
        "\nOK: Voigt 6×6 stiffness — σ_x = E ε_x, τ_xy = μ γ_xy, lateral faces traction-free."
    );
}
