//! Modal analysis and frequency-response for `tpt-fem`.
//!
//! Combines `tpt-fem-eigen`'s generalized eigensolver with
//! `tpt-fem-dynamic`'s Newmark integrator to support two standard vibration
//! workflows:
//!
//! * [`modal_analysis`] + [`ModalData::frequency_response`] — the steady-state
//!   harmonic response `u(Ω)` of a structure under a sinusoidal load
//!   `f(t) = Re(F·e^{iΩt})`, resolved through a truncated set of eigenmodes
//!   (the classic modal-superposition / frequency-response workflow used for
//!   vibration qualification and fatigue screening).
//! * [`ModalData::modal_superposition`] — the *time-domain* complement: the same
//!   eigenbasis reduces the second-order system to independent single-DOF
//!   modal equations, each integrated with `tpt-fem-dynamic`'s Newmark
//!   scheme, then recombined into the physical response `u(t)`.
//!
//! ```
//! use tpt_fem_modal::{modal_analysis, ModalData};
//! use tpt_fem_sparse::Coo;
//!
//! // SDOF: M=1, K=4 (ω=2). Modal analysis recovers the single mode.
//! let m = Coo { rows: vec![0], cols: vec![0], vals: vec![1.0] };
//! let k = Coo { rows: vec![0], cols: vec![0], vals: vec![4.0] };
//! let data = modal_analysis(&k, &m, 0.0, 1, 8, 0.0).unwrap();
//! assert!((data.omega[0] - 2.0).abs() < 1e-6);
//! // Harmonic load F=1 at Ω=1: amplitude |U| = 1/|4 - 1²| = 1/3 (undamped).
//! let resp = data.frequency_response(&m, &[1.0], &[1.0]);
//! let amp = (resp[0].displacement[0].0.powi(2) + resp[0].displacement[0].1.powi(2)).sqrt();
//! assert!((amp - 1.0 / 3.0).abs() < 1e-4);
//! ```

use tpt_fem_dynamic::{newmark, NewmarkOptions};
use tpt_fem_eigen::{generalized_lanczos_eigs, matvec};
use tpt_fem_sparse::{Coo, SparseError};

/// A truncated modal basis of a second-order system `M·ü + C·u̇ + K·u = f(t)`.
///
/// Produced by [`modal_analysis`]. The modes are the eigenvectors of the
/// generalized eigenproblem `K φ = ω² M φ`; they are *not* normalized, so the
/// per-mode modal mass `φᵀ M φ` is stored alongside for use in the reduced
/// equations.
pub struct ModalData {
    /// Natural circular frequencies `ω_i = sqrt(λ_i)` (rad/s), ascending.
    pub omega: Vec<f64>,
    /// Mode shapes `φ_i` (one `Vec<f64>` of length = number of DOFs, per mode).
    pub modes: Vec<Vec<f64>>,
    /// Modal mass `m_i = φ_iᵀ M φ_i` for each mode.
    pub modal_mass: Vec<f64>,
    /// Modal damping ratio `ζ_i` assumed for each mode (Rayleigh/proportional
    /// damping is implied: `c_i = 2 ζ_i ω_i m_i`).
    pub zeta: Vec<f64>,
}

/// Steady-state harmonic response at a single excitation frequency.
pub struct HarmonicResponse {
    /// Excitation frequency `Ω` (rad/s).
    pub frequency: f64,
    /// Complex displacement amplitude `U` such that the physical response is
    /// `Re(U·e^{iΩt})`. Each entry is `(real, imag)` for one global DOF.
    pub displacement: Vec<(f64, f64)>,
}

/// Solve the generalized eigenproblem `K φ = λ M φ` for the `num` smallest
/// eigenpairs, returning the modal basis (frequencies, shapes, modal masses)
/// with a uniform modal damping ratio `zeta`.
///
/// `sigma` is the Lanczos shift (use `0.0` for the smallest modes; a positive
/// shift biases toward modes near `σ`). `lanczos_dim` is the Krylov dimension
/// (must exceed `num`; 2·`num` is a safe default).
pub fn modal_analysis(
    k: &Coo,
    m: &Coo,
    sigma: f64,
    num: usize,
    lanczos_dim: usize,
    zeta: f64,
) -> Result<ModalData, SparseError> {
    let pairs = generalized_lanczos_eigs(k, m, sigma, num, lanczos_dim)?;
    let mut omega = Vec::with_capacity(pairs.len());
    let mut modes = Vec::with_capacity(pairs.len());
    let mut modal_mass = Vec::with_capacity(pairs.len());
    for (lam, phi) in pairs {
        // λ should be ≥ 0 for a stable (positive-definite K) system.
        let w = if lam > 0.0 { lam.sqrt() } else { 0.0 };
        omega.push(w);
        let mp = matvec(m, &phi);
        let mi = phi.iter().zip(mp.iter()).map(|(a, b)| a * b).sum::<f64>();
        modes.push(phi);
        modal_mass.push(mi);
    }
    let zeta = vec![zeta; omega.len()];
    Ok(ModalData {
        omega,
        modes,
        modal_mass,
        zeta,
    })
}

impl ModalData {
    /// Steady-state harmonic response under a sinusoidal load with complex
    /// amplitude `F` (real part supplied; imaginary part taken as zero).
    ///
    /// For each excitation frequency `Ω ∈ freqs`, returns the complex
    /// displacement amplitude `U(Ω)` via modal superposition:
    /// `U = Σ_i φ_i · F_i / (k_i − Ω² m_i + i·Ω·c_i)`, where `F_i = φ_iᵀ F` is the
    /// modal load, `k_i = ω_i² m_i`, and `c_i = 2 ζ_i ω_i m_i`.
    pub fn frequency_response(
        &self,
        _m: &Coo,
        force: &[f64],
        freqs: &[f64],
    ) -> Vec<HarmonicResponse> {
        let n = force.len();
        let nm = self.omega.len();
        let mut out = Vec::with_capacity(freqs.len());
        for &omega_f in freqs {
            let mut re = vec![0.0; n];
            let mut im = vec![0.0; n];
            for i in 0..nm {
                let wi = self.omega[i];
                let mi = self.modal_mass[i];
                let fi = self.modes[i]
                    .iter()
                    .zip(force)
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                let ci = 2.0 * self.zeta[i] * wi * mi;
                let ki = wi * wi * mi;
                let dr = ki - omega_f * omega_f * mi;
                let di = omega_f * ci;
                let denom = dr * dr + di * di;
                let qr = fi * dr / denom;
                let qi = -fi * di / denom;
                for (d, &phi) in self.modes[i].iter().enumerate() {
                    re[d] += phi * qr;
                    im[d] += phi * qi;
                }
            }
            out.push(HarmonicResponse {
                frequency: omega_f,
                displacement: re.into_iter().zip(im).collect(),
            });
        }
        out
    }

    /// Time-domain response via modal superposition.
    ///
    /// Projects the initial state `(u0, v0)` and the load `f(t)` onto the modal
    /// basis, integrates each independent single-DOF modal equation with
    /// [`tpt_fem_dynamic`]'s Newmark scheme, and recombines into the physical
    /// displacement history `u(t)`.
    pub fn modal_superposition<F>(
        &self,
        u0: &[f64],
        v0: &[f64],
        f: &F,
        opts: &NewmarkOptions,
        nsteps: usize,
    ) -> Vec<(f64, Vec<f64>)>
    where
        F: Fn(f64) -> Vec<f64>,
    {
        let n = u0.len();
        let nm = self.omega.len();
        // Initial modal coordinates (M-orthogonal projection).
        let mut q0 = vec![0.0; nm];
        let mut qd0 = vec![0.0; nm];
        for i in 0..nm {
            let mi = self.modal_mass[i];
            let mu0 = self.modes[i]
                .iter()
                .zip(u0)
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let mv0 = self.modes[i]
                .iter()
                .zip(v0)
                .map(|(a, b)| a * b)
                .sum::<f64>();
            q0[i] = mu0 / mi;
            qd0[i] = mv0 / mi;
        }
        // Integrate each modal coordinate independently.
        let mut q_hist = Vec::with_capacity(nm);
        for i in 0..nm {
            let mi = self.modal_mass[i];
            let ci = 2.0 * self.zeta[i] * self.omega[i] * mi;
            let ki = self.omega[i] * self.omega[i] * mi;
            let mm = Coo {
                rows: vec![0],
                cols: vec![0],
                vals: vec![mi],
            };
            let cc = Coo {
                rows: vec![0],
                cols: vec![0],
                vals: vec![ci],
            };
            let kk = Coo {
                rows: vec![0],
                cols: vec![0],
                vals: vec![ki],
            };
            let fi = |t: f64| {
                let fv = f(t);
                vec![self.modes[i]
                    .iter()
                    .zip(&fv)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()]
            };
            let hist = newmark(&mm, &cc, &kk, &[q0[i]], &[qd0[i]], fi, opts, nsteps);
            q_hist.push(hist);
        }
        // Recombine into the physical response.
        let nsteps_out = q_hist[0].len();
        let mut out = Vec::with_capacity(nsteps_out);
        for step in 0..nsteps_out {
            let t = q_hist[0][step].0;
            let mut u = vec![0.0; n];
            for i in 0..nm {
                let qi = q_hist[i][step].1[0];
                for d in 0..n {
                    u[d] += self.modes[i][d] * qi;
                }
            }
            out.push((t, u));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_dynamic::NewmarkOptions;
    use tpt_fem_sparse::Coo;

    fn sdof(mass: f64, stiff: f64) -> (Coo, Coo) {
        let m = Coo {
            rows: vec![0],
            cols: vec![0],
            vals: vec![mass],
        };
        let k = Coo {
            rows: vec![0],
            cols: vec![0],
            vals: vec![stiff],
        };
        (m, k)
    }

    #[test]
    fn sdof_frequency_response_closed_form() {
        // M=1, K=4 (ω=2). Force amplitude 1 at Ω=1 -> |U| = 1/|4 - 1| = 1/3.
        let (m, k) = sdof(1.0, 4.0);
        let data = modal_analysis(&k, &m, 0.0, 1, 8, 0.0).unwrap();
        let resp = data.frequency_response(&m, &[1.0], &[1.0]);
        let (re, im) = resp[0].displacement[0];
        let amp = (re * re + im * im).sqrt();
        assert!((amp - 1.0 / 3.0).abs() < 1e-4, "amp = {amp}");
    }

    #[test]
    fn sdof_modal_superposition_matches_newmark() {
        // Free vibration of the single mode: u(t) = cos(2t) from u0=1, v0=0.
        // One mode spans the full (1-DOF) space, so modal superposition must
        // reproduce the closed form to Newmark accuracy.
        let (m, k) = sdof(1.0, 4.0);
        let data = modal_analysis(&k, &m, 0.0, 1, 8, 0.0).unwrap();
        let opts = NewmarkOptions {
            dt: 0.005,
            beta: 0.25,
            gamma: 0.5,
        };
        let nsteps = 400; // t = 2.0
        let hist = data.modal_superposition(&[1.0], &[0.0], &|_| vec![0.0], &opts, nsteps);
        let (t, u) = &hist[nsteps];
        let want = (2.0 * t).cos();
        assert!((u[0] - want).abs() < 5e-3, "got {} want {}", u[0], want);
    }
}
