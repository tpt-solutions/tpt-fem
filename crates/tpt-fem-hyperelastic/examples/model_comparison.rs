//! Stress–stretch comparison of Neo-Hookean, Mooney–Rivlin and Ogden.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-hyperelastic --example model_comparison
//! ```
//!
//! On the same incompressible uniaxial stretch `λ`, the three models' principal
//! nominal stress
//! ```text
//!   Neo-Hookean : S = μ(λ − λ⁻²)
//!   Mooney-Rivlin: S = 2(λ − λ⁻²)(c₁ + c₂ λ⁻¹)
//!   Ogden       : S = Σ μᵢ(λ^{αᵢ−1} − λ^{−αᵢ/2 − 1})
//! ```
//! are compared, each verified against its closed form.

use tpt_fem_hyperelastic::{mooney_rivlin_piola, neo_hookean_piola, ogden_piola, OgdenTerm};

fn identity() -> [[f64; 3]; 3] {
    let mut f = [[0.0; 3]; 3];
    for i in 0..3 {
        f[i][i] = 1.0;
    }
    f
}

fn main() {
    let mu = 100.0;
    let c1 = 80.0;
    let c2 = 20.0;
    let terms = [
        OgdenTerm {
            mu: 60.0,
            alpha: 1.3,
        },
        OgdenTerm {
            mu: 40.0,
            alpha: -2.0,
        },
    ];
    println!("μ = {mu}, c1 = {c1}, c2 = {c2}, Ogden[ (60,1.3), (40,-2.0) ]");
    println!(
        "{:<7} {:<14} {:<14} {:<14}",
        "λ", "Neo-Hookean", "Mooney-Rivlin", "Ogden"
    );
    println!("{}", "-".repeat(48));
    for &lam in &[1.1_f64, 1.3, 1.5, 1.8, 2.0] {
        // Deformation gradient and traction-free pressures.
        let mut f = identity();
        f[0][0] = lam;
        f[1][1] = lam.powf(-0.5);
        f[2][2] = lam.powf(-0.5);

        let p_nh = mu * lam.powf(-1.0);
        let s_nh = neo_hookean_piola(&f, mu, p_nh).expect("F is invertible")[0][0];
        let closed_nh = mu * (lam - lam.powf(-2.0));

        let p_mr = 2.0 * c1 * lam.powf(-1.0) + 2.0 * c2 * (lam + lam.powf(-2.0));
        let s_mr = mooney_rivlin_piola(&f, c1, c2, p_mr).expect("F is invertible")[0][0];
        let closed_mr = 2.0 * (lam - lam.powf(-2.0)) * (c1 + c2 * lam.powf(-1.0));

        let lm = lam.powf(-0.5);
        let lambdas = [lam, lm, lm];
        let p_og: f64 = terms
            .iter()
            .map(|t| t.mu * t.alpha * lambdas[1].powf(t.alpha))
            .sum();
        let s_og = ogden_piola(&[lam, lm, lm], &identity(), &terms, p_og)[0][0];
        let closed_og: f64 = terms
            .iter()
            .map(|t| t.mu * t.alpha * (lam.powf(t.alpha - 1.0) - lam.powf(-(t.alpha / 2.0 + 1.0))))
            .sum();

        println!(
            "{:<7.3} {:<14.6e} {:<14.6e} {:<14.6e}",
            lam, s_nh, s_mr, s_og
        );
        assert!((s_nh - closed_nh).abs() / closed_nh.abs() < 1e-12);
        assert!((s_mr - closed_mr).abs() / closed_mr.abs() < 1e-12);
        assert!((s_og - closed_og).abs() / closed_og.abs() < 1e-12);
    }
    println!("\nOK: all three models match their closed-form principal nominal stress.");
}
