//! Isotropic vs kinematic hardening and the Bauschinger effect.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem-plasticity --example isotropic_vs_kinematic
//! ```
//!
//! A material point is taken through a *reversed* load cycle in *pure shear*:
//! sheared past yield, then sheared back the other way. Pure shear is chosen
//! because it is purely deviatoric — the volumetric response is identically zero,
//! so the strain path is exact at every stage and no lateral-strain solve is
//! needed (unlike a uniaxial-stress path, see the `uniaxial_sweep` example).
//!
//! With **isotropic** hardening the von Mises surface grows about the origin, so
//! the reverse yield stress has the *same magnitude* as the forward one. With
//! **kinematic** hardening the surface keeps its radius but translates with the
//! back stress `α`, so reverse yield happens much earlier — the Bauschinger
//! effect.
//!
//! Recall that in pure shear `q = √3·τ`, so yield begins at `τ_y = σ_y⁰/√3`.

use tpt_fem_plasticity::{j2_return_mapping, PlasticState, PlasticityParams};

/// Voigt shear-doubled norm of a deviatoric 6-vector.
fn dev_norm(s: &[f64; 6]) -> f64 {
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
        .sqrt()
}

/// von Mises equivalent stress of the *relative* stress `η = dev(σ) − α`.
fn relative_q(sig: &[f64; 6], back: &[f64]) -> f64 {
    let mean = (sig[0] + sig[1] + sig[2]) / 3.0;
    let mut eta = *sig;
    for i in 0..3 {
        eta[i] -= mean;
    }
    for i in 0..6 {
        eta[i] -= back[i];
    }
    (1.5_f64).sqrt() * dev_norm(&eta)
}

/// Pure-shear strain state with engineering shear `γ_xy`.
fn shear(gamma: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, gamma, 0.0, 0.0]
}

struct Cycle {
    /// Shear stress at the end of forward loading.
    forward_tau: f64,
    /// Accumulated equivalent plastic strain at load reversal.
    eps_eq_at_reversal: f64,
    /// Shear stress at which plastic flow resumes in reverse.
    reverse_tau: f64,
    /// Radius of the yield surface (in `q`) when reverse flow starts.
    reverse_radius: f64,
    /// `q` of the relative stress when reverse flow starts.
    reverse_q: f64,
}

/// Shear forward to `gamma_max`, then reverse, and report the reverse-yield point.
fn run_cycle(p: &PlasticityParams, gamma_max: f64) -> Cycle {
    let steps = 400;
    let mut state = PlasticState::new();
    let mut forward_tau = 0.0;

    // Forward loading, in increments so the plastic history accumulates.
    for k in 1..=steps {
        let g = gamma_max * k as f64 / steps as f64;
        let (sig, next) = j2_return_mapping(p, &shear(g), &state);
        forward_tau = sig[3];
        state = next;
    }
    let eps_eq_at_reversal = state.eps_eq;

    // Reverse loading: sweep gamma back down through zero into negative shear.
    let g_end = -2.5 * gamma_max;
    for k in 1..=steps * 4 {
        let g = gamma_max + (g_end - gamma_max) * k as f64 / (steps * 4) as f64;
        let before = state.eps_eq;
        let (sig, next) = j2_return_mapping(p, &shear(g), &state);
        if next.eps_eq > before + 1e-18 {
            // First step where flow resumes: this is the reverse yield point. The
            // consistency condition holds for the *updated* state, so measure the
            // relative stress and the surface radius against `next`.
            return Cycle {
                forward_tau,
                eps_eq_at_reversal,
                reverse_tau: sig[3],
                reverse_radius: p.yield_stress + p.iso_hardening * next.eps_eq,
                reverse_q: relative_q(&sig, &next.back_stress),
            };
        }
        state = next;
    }
    panic!("reverse yield was never reached");
}

fn main() {
    let h = 5e9; // hardening modulus, E/40
    let base = PlasticityParams::steel();
    let tau_y = base.yield_stress / 3.0_f64.sqrt();
    let mu = base.young / (2.0 * (1.0 + base.poisson));
    let gamma_y = tau_y / mu;
    let gamma_max = 6.0 * gamma_y;

    println!("Reversed pure-shear cycle: isotropic vs kinematic hardening");
    println!(
        "  E = {:.3e} Pa, nu = {:.2}, G = {:.3e} Pa, sigma_y = {:.3e} Pa",
        base.young, base.poisson, mu, base.yield_stress
    );
    println!("  tau_y = sigma_y/sqrt(3) = {tau_y:.6e} Pa, gamma_y = {gamma_y:.4e}");
    println!("  hardening modulus H = {h:.3e} Pa, gamma_max = {gamma_max:.4e}\n");

    let iso = PlasticityParams {
        iso_hardening: h,
        kin_hardening: 0.0,
        ..base
    };
    let kin = PlasticityParams {
        iso_hardening: 0.0,
        kin_hardening: h,
        ..base
    };

    let c_iso = run_cycle(&iso, gamma_max);
    let c_kin = run_cycle(&kin, gamma_max);

    println!(
        "  {:<12} {:<18} {:<18} {:<16}",
        "model", "forward tau (Pa)", "reverse tau (Pa)", "eps_p_bar"
    );
    println!("  {}", "-".repeat(66));
    println!(
        "  {:<12} {:<18.6e} {:<18.6e} {:<16.4e}",
        "isotropic", c_iso.forward_tau, c_iso.reverse_tau, c_iso.eps_eq_at_reversal
    );
    println!(
        "  {:<12} {:<18.6e} {:<18.6e} {:<16.4e}",
        "kinematic", c_kin.forward_tau, c_kin.reverse_tau, c_kin.eps_eq_at_reversal
    );

    // Isotropic hardening: the surface grows about the origin, so reverse yield
    // mirrors the forward stress.
    let iso_sym =
        (c_iso.reverse_tau.abs() - c_iso.forward_tau.abs()).abs() / c_iso.forward_tau.abs();
    println!(
        "\n  isotropic |reverse tau| / |forward tau| = {:.6}",
        c_iso.reverse_tau.abs() / c_iso.forward_tau.abs()
    );
    assert!(
        iso_sym < 5e-3,
        "isotropic hardening must yield symmetrically, off by {iso_sym:.3e}"
    );

    // Kinematic hardening: the radius is unchanged (still sigma_y0), the centre
    // has moved, so reverse yield occurs 2*tau_y earlier in stress.
    println!(
        "  kinematic reverse q = {:.6e} Pa vs unchanged radius sigma_y = {:.6e} Pa",
        c_kin.reverse_q, c_kin.reverse_radius
    );
    assert!(
        (c_kin.reverse_q - c_kin.reverse_radius).abs() / c_kin.reverse_radius < 1e-6,
        "kinematic yield surface must keep its radius"
    );
    let want_reverse = c_kin.forward_tau - 2.0 * tau_y;
    println!(
        "  kinematic reverse tau = {:.6e} Pa vs forward tau - 2*tau_y = {:.6e} Pa",
        c_kin.reverse_tau, want_reverse
    );
    assert!(
        (c_kin.reverse_tau - want_reverse).abs() / tau_y < 5e-2,
        "kinematic reverse yield should sit 2*tau_y below the forward stress"
    );

    // The Bauschinger effect itself.
    assert!(
        c_kin.reverse_tau.abs() < c_iso.reverse_tau.abs(),
        "Bauschinger effect missing: |tau_kin| {:.3e} >= |tau_iso| {:.3e}",
        c_kin.reverse_tau.abs(),
        c_iso.reverse_tau.abs()
    );
    println!(
        "\nOK: Bauschinger effect - kinematic reverse yield |tau| = {:.3e} Pa is well",
        c_kin.reverse_tau.abs()
    );
    println!(
        "    below the isotropic {:.3e} Pa, and each model keeps its own yield surface.",
        c_iso.reverse_tau.abs()
    );
}
