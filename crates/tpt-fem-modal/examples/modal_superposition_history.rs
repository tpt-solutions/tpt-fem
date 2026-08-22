//! Free-vibration time history of a 2-DOF spring–mass chain computed by modal
//! superposition: each independent modal equation is integrated with Newmark's
//! scheme, then the modes are recombined into physical displacements. The
//! result is compared against the analytic superposition of the two normal
//! modes.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-modal --example modal_superposition_history
//! ```

use tpt_fem_dynamic::NewmarkOptions;
use tpt_fem_modal::modal_analysis;
use tpt_fem_sparse::Coo;

fn main() {
    // Same chain as the frequency-response example:
    //   M = I,  K = [[2, -1], [-1, 1]],  ω₁ ≈ 0.618, ω₂ ≈ 1.618.
    let m = Coo {
        rows: vec![0, 1],
        cols: vec![0, 1],
        vals: vec![1.0, 1.0],
    };
    let k = Coo {
        rows: vec![0, 0, 1, 1],
        cols: vec![0, 1, 0, 1],
        vals: vec![2.0, -1.0, -1.0, 1.0],
    };

    let data = modal_analysis(&k, &m, 0.0, 2, 8, 0.0).unwrap();

    // Initial conditions: pull mass 1 by u₀ = (1, 0), release from rest,
    // no external load — free vibration.
    let u0 = [1.0_f64, 0.0];
    let opts = NewmarkOptions {
        dt: 0.002,
        beta: 0.25,
        gamma: 0.5,
    };
    let nsteps = 1000; // simulate t ∈ [0, 2]
    let hist = data.modal_superposition(&u0, &[0.0, 0.0], &|_| vec![0.0; 2], &opts, nsteps);

    // Analytic free vibration: u(t) = Σ_i φ_i · (φᵀ M u₀ / m_i) · cos(ω_i t).
    // The modes are unnormalized, so project explicitly with the stored modal
    // masses. Here M = I for simplicity.
    let proj: Vec<f64> = data
        .modes
        .iter()
        .zip(&data.modal_mass)
        .map(|(phi, mi)| {
            let mu0: f64 = phi.iter().zip(&u0).map(|(p, v)| p * v).sum();
            mu0 / mi
        })
        .collect();


    println!("Free vibration of the 2-DOF chain (ζ = 0):");
    println!("  {:>8}  {:>14}  {:>14}", "t", "u1 numeric", "u1 analytic");
    for step in [0usize, 250, 500, 750, 1000] {
        let (t, u) = &hist[step];
        let analytic: f64 = (0..2)
            .map(|i| data.modes[i][0] * proj[i] * (data.omega[i] * t).cos())
            .sum();
        let u1 = u[0];
        println!("  {t:>8.3}  {u1:>14.8}  {analytic:>14.8}");
        assert!(
            (u1 - analytic).abs() < 5e-3,
            "mismatch at t={t}: {u1} vs {analytic}"
        );
    }

    // Energy sanity: undamped free vibration conserves energy, so the response
    // amplitude cannot grow.
    let max_abs_u1 = hist
        .iter()
        .map(|(_, u)| u[0].abs())
        .fold(0.0_f64, f64::max);
    assert!(max_abs_u1 < 1.0 + 1e-6 + 2.0); // generous bound; no blow-up
    println!("\nmax |u₁(t)| = {:.6} over t ∈ [0, 2]", max_abs_u1);
    println!("Modal superposition matches the closed form within Newmark accuracy.");
}
