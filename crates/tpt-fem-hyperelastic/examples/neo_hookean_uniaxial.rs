//! Uniaxial extension of an incompressible Neo-Hookean solid.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-hyperelastic --example neo_hookean_uniaxial
//! ```
//!
//! Sweeps the axial stretch `λ` and compares the first Piola–Kirchhoff (nominal)
//! axial stress computed by `neo_hookean_piola` against the closed form
//! `P = μ(λ − λ⁻²)` for traction-free lateral faces, where the pressure is
//! `p = μ λ⁻¹`.

use tpt_fem_hyperelastic::neo_hookean_piola;

fn main() {
    let mu = 100.0;
    println!("Incompressible Neo-Hookean, μ = {mu} Pa");
    println!(
        "{:<8} {:<16} {:<16} {:<8}",
        "λ", "P_num (Pa)", "P_closed (Pa)", "err"
    );
    println!("{}", "-".repeat(46));
    for &lam in &[1.1_f64, 1.3, 1.5, 1.8, 2.0, 0.8] {
        let mut f = [[0.0; 3]; 3];
        f[0][0] = lam;
        f[1][1] = lam.powf(-0.5);
        f[2][2] = lam.powf(-0.5);
        let p = mu * lam.powf(-1.0); // traction-free lateral faces
        let pk = neo_hookean_piola(&f, mu, p);
        let closed = mu * (lam - lam.powf(-2.0));
        let err = (pk[0][0] - closed).abs() / closed.abs();
        println!(
            "{:<8.3} {:<16.6e} {:<16.6e} {:<8.2e}",
            lam, pk[0][0], closed, err
        );
        assert!(err < 1e-12, "Neo-Hookean mismatch at λ={lam}");
        // Lateral tractions vanish.
        assert!(pk[1][1].abs() < 1e-9 && pk[2][2].abs() < 1e-9);
    }
    println!("\nOK: P matches μ(λ − λ⁻²); lateral faces traction-free.");
}
