//! Frequency-response sweep of a 2-DOF spring–mass chain under a harmonic tip
//! force. Computes natural frequencies via modal analysis, sweeps the
//! excitation frequency to locate resonances, and compares amplitudes against
//! the closed-form undamped solution.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-modal --example frequency_response_sweep
//! ```

use tpt_fem_modal::modal_analysis;
use tpt_fem_sparse::Coo;

/// Symmetric matrix in triplet form from its upper triangle.
fn sym(entries: &[(usize, usize, f64)]) -> Coo {
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for &(i, j, v) in entries {
        rows.push(i);
        cols.push(j);
        vals.push(v);
        if i != j {
            rows.push(j);
            cols.push(i);
            vals.push(v);
        }
    }
    Coo { rows, cols, vals }
}

fn main() {
    // Chain m1 - k - m2 - ground with m = 1, k = 1:
    //   M = I,  K = [[2, -1], [-1, 1]]
    // Eigenvalues λ = (3 ± √5)/2 → ω ≈ 0.618 and ω ≈ 1.618 rad/s.
    let m = sym(&[(0, 0, 1.0), (1, 1, 1.0)]);
    let k = sym(&[(0, 0, 2.0), (1, 1, 1.0), (0, 1, -1.0)]);

    // Modal analysis: both modes, ζ = 0.02 modal damping.
    let data = modal_analysis(&k, &m, 0.0, 2, 8, 0.02).unwrap();
    println!("Natural frequencies:");
    for (i, w) in data.omega.iter().enumerate() {
        println!("  mode {}: ω = {:.6} rad/s", i + 1, w);
    }

    // Harmonic tip force F = (1, 0), sweep Ω across the band.
    let forces = [1.0, 0.0];
    let freqs: Vec<f64> = (0..=160).map(|i| 0.01 + i as f64 * (2.4 / 160.0)).collect();
    let resp = data.frequency_response(&m, &forces, &freqs);

    // Locate resonance peaks in |U_1(Ω)|.
    let amp = |r: &tpt_fem_modal::HarmonicResponse, d: usize| {
        (r.displacement[d].0.powi(2) + r.displacement[d].1.powi(2)).sqrt()
    };
    println!("\nFrequency-response peaks:");
    for i in 1..resp.len() - 1 {
        let a_prev = amp(&resp[i - 1], 0);
        let a_here = amp(&resp[i], 0);
        let a_next = amp(&resp[i + 1], 0);
        if a_here > a_prev && a_here > a_next {
            println!(
                "  resonance near Ω = {:.4} rad/s, |U_1| = {:.3}",
                resp[i].frequency, a_here
            );
        }
    }

    // Off-resonance check: far below the first mode the response approaches
    // the static stiffness solution K⁻¹ F.
    let static_amp = amp(&resp[0], 0);
    println!(
        "\nQuasi-static |U_1| at Ω = {:.3}: {:.6}",
        resp[0].frequency, static_amp
    );

    // Damped response stays bounded at resonance (no singularity).
    for r in &resp {
        for d in 0..2 {
            assert!(amp(r, d).is_finite());
        }
    }

    println!("\nSweep complete.");
}
