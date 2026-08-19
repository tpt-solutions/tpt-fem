//! Uniaxial load / unload / reload through the 1-D return-mapping law.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem-plasticity --example uniaxial_sweep
//! ```
//!
//! Drives a single material point along a uniaxial-*stress* path with
//! [`plastic_1d`] and checks it against the closed-form elastic-plastic response
//! for two materials:
//!
//! * perfect plasticity — the stress plateaus at `σ_y⁰`;
//! * linear isotropic hardening — `σ = σ_y⁰ + H ε̄ᵖ` with `ε = σ/E + ε̄ᵖ`.
//!
//! The final section shows *why* the scalar law is not simply the 3-D law
//! evaluated at `ε = (ε_x, −ν ε_x, −ν ε_x, 0, 0, 0)`: that strain state is the
//! uniaxial-stress state only while the response is elastic.

use tpt_fem_plasticity::{j2_return_mapping, plastic_1d, PlasticState, PlasticityParams};

/// Closed-form uniaxial-stress response for linear isotropic hardening.
///
/// Combining `σ = σ_y⁰ + H ε̄ᵖ` (the grown yield stress) with the strain split
/// `ε = σ/E + ε̄ᵖ` gives `ε̄ᵖ = (E ε − σ_y⁰)/(E + H)`.
fn closed_form(eps: f64, e: f64, sig_y: f64, h: f64) -> f64 {
    if e * eps <= sig_y {
        e * eps
    } else {
        let epsp = (e * eps - sig_y) / (e + h);
        sig_y + h * epsp
    }
}

/// Monotonic tension sweep, printed and checked against the closed form.
fn sweep(label: &str, p: &PlasticityParams) {
    let e = p.young;
    let h = p.iso_hardening + p.kin_hardening;
    let ey = p.yield_stress / e;

    println!("\n{label}");
    println!(
        "  E = {:.3e} Pa, sigma_y = {:.3e} Pa, H = {:.3e} Pa, eps_y = {:.4e}",
        e, p.yield_stress, h, ey
    );
    println!(
        "  {:<12} {:<16} {:<16} {:<14} {:<8}",
        "eps", "sigma (Pa)", "sigma_closed", "eps_p_bar", "regime"
    );
    println!("  {}", "-".repeat(70));

    for k in 1..=8 {
        let eps = 0.5 * ey * k as f64;
        // Each point is evaluated from a virgin state, so this traces the
        // monotonic loading curve.
        let (sig, epsp, plastic) = plastic_1d(p, eps, 0.0);
        let want = closed_form(eps, e, p.yield_stress, h);
        println!(
            "  {:<12.4e} {:<16.6e} {:<16.6e} {:<14.4e} {:<8}",
            eps,
            sig,
            want,
            epsp,
            if plastic { "plastic" } else { "elastic" }
        );
        assert!(
            (sig - want).abs() / want.abs() < 1e-12,
            "stress mismatch at eps={eps:.3e}: {sig:.6e} vs {want:.6e}"
        );
    }
}

/// Load past yield, unload elastically, then reload — threading `eps_p_bar`
/// between steps so the history is remembered.
fn load_unload_reload(p: &PlasticityParams) {
    let e = p.young;
    let ey = p.yield_stress / e;
    println!("\nLoad / unload / reload (history threaded through eps_p_bar)");
    println!(
        "  {:<12} {:<16} {:<14} {:<8}",
        "eps", "sigma (Pa)", "eps_p_bar", "regime"
    );
    println!("  {}", "-".repeat(56));

    // Up to 3 eps_y, back down to 1 eps_y, then up to 4 eps_y.
    let path: Vec<f64> = (1..=6)
        .map(|k| 0.5 * ey * k as f64)
        .chain((2..=5).rev().map(|k| 0.5 * ey * k as f64))
        .chain((3..=8).map(|k| 0.5 * ey * k as f64))
        .collect();

    let mut epsp = 0.0;
    let mut peak_sig = 0.0_f64;
    let mut unload_slope_checked = false;
    let mut prev: Option<(f64, f64)> = None;

    for eps in path {
        let (sig, next, plastic) = plastic_1d(p, eps, epsp);
        println!(
            "  {:<12.4e} {:<16.6e} {:<14.4e} {:<8}",
            eps,
            sig,
            next,
            if plastic { "plastic" } else { "elastic" }
        );

        // Unloading from a plastic state must be elastic with slope E.
        if let Some((peps, psig)) = prev {
            if eps < peps && !plastic {
                let slope = (sig - psig) / (eps - peps);
                assert!(
                    (slope - e).abs() / e < 1e-9,
                    "elastic unloading slope should be E, got {slope:.6e}"
                );
                unload_slope_checked = true;
            }
        }
        // The accumulated plastic strain never decreases.
        assert!(next >= epsp - 1e-18, "eps_p_bar must be monotone");
        peak_sig = peak_sig.max(sig);
        epsp = next;
        prev = Some((eps, sig));
    }
    assert!(unload_slope_checked, "the path never exercised unloading");
    println!("  peak stress = {peak_sig:.6e} Pa, final eps_p_bar = {epsp:.4e}");
}

/// Show that prescribing the *elastic* lateral contraction in the 3-D law is not
/// uniaxial stress once the point yields.
fn uniaxial_stress_vs_prescribed_lateral(p: &PlasticityParams) {
    let eps = 4.0 * p.yield_stress / p.young; // well past yield
    let nu = p.poisson;

    // (a) Prescribe eps_lat = -nu*eps (valid only in the elastic range).
    let (s_pre, _) = j2_return_mapping(
        p,
        &[eps, -nu * eps, -nu * eps, 0.0, 0.0, 0.0],
        &PlasticState::new(),
    );

    // (b) Solve for the lateral strain that actually makes sigma_y = sigma_z = 0.
    let lateral_stress = |lat: f64| {
        let (s, _) = j2_return_mapping(p, &[eps, lat, lat, 0.0, 0.0, 0.0], &PlasticState::new());
        s[1]
    };
    let (mut a, mut b) = (-nu * eps, -0.5 * eps);
    let (mut fa, mut fb) = (lateral_stress(a), lateral_stress(b));
    for _ in 0..200 {
        if fb.abs() < 1e-6 || (fb - fa).abs() < 1e-30 {
            break;
        }
        let c = b - fb * (b - a) / (fb - fa);
        (a, fa) = (b, fb);
        (b, fb) = (c, lateral_stress(c));
    }
    let (s_true, _) = j2_return_mapping(p, &[eps, b, b, 0.0, 0.0, 0.0], &PlasticState::new());
    let (s_1d, _, _) = plastic_1d(p, eps, 0.0);

    println!("\nWhy the scalar law is not the 3-D law at eps_lat = -nu*eps");
    println!("  axial strain eps_x            = {eps:.4e}");
    println!(
        "  (a) eps_lat = -nu*eps_x       -> sigma_x = {:.6e} Pa, sigma_y = {:.6e} Pa",
        s_pre[0], s_pre[1]
    );
    println!(
        "  (b) eps_lat solved for sigma_y = 0 -> eps_lat = {:.4e}, sigma_x = {:.6e} Pa",
        b, s_true[0]
    );
    println!("  1-D law plastic_1d            -> sigma_x = {s_1d:.6e} Pa");
    println!(
        "  incompressible plastic flow drives eps_lat from {:.4e} towards {:.4e}",
        -nu * eps,
        -0.5 * eps
    );

    // (a) develops spurious lateral stress; (b) agrees with the scalar law.
    assert!(
        s_pre[1].abs() > 0.01 * p.yield_stress,
        "prescribing the elastic contraction should leave lateral stress"
    );
    assert!(
        (s_true[0] - s_1d).abs() / s_1d.abs() < 1e-6,
        "genuine uniaxial stress must match the 1-D law: {:.6e} vs {:.6e}",
        s_true[0],
        s_1d
    );
}

fn main() {
    println!("J2 uniaxial-stress return mapping");

    let perfect = PlasticityParams::steel();
    sweep(
        "Perfect plasticity (H = 0): stress plateaus at sigma_y",
        &perfect,
    );

    let hardening = PlasticityParams {
        iso_hardening: 4e9, // E/50
        ..PlasticityParams::steel()
    };
    sweep("Linear isotropic hardening (H = E/50)", &hardening);

    load_unload_reload(&hardening);
    uniaxial_stress_vs_prescribed_lateral(&perfect);

    println!("\nOK: the 1-D law matches the closed form, unloading is elastic with");
    println!("    slope E, and genuine uniaxial stress in 3-D reproduces the 1-D law.");
}
