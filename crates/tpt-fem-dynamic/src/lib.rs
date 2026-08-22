//! Time integration for `tpt-fem`.
//!
//! Provides two standard schemes for the generic second-order system
//!
//! ```text
//! M·ü + C·v + K·u = f(t)
//! ```
//!
//! assembled from any physics crate:
//!
//! * [`newmark`] — implicit [Newmark-beta](https://en.wikipedia.org/wiki/Newmark-beta_method)
//!   (defaults to the unconditionally-stable average-acceleration rule), and
//! * [`central_difference`] — explicit central-difference, using the lumped
//!   (diagonal) mass for a cheap `M⁻¹`.
//!
//! Both reuse the consistent/lumped mass matrices produced by
//! [`tpt-fem-elasticity`] rather than re-deriving them.
//!
//! ```
//! use tpt_fem_dynamic::{central_difference, newmark, NewmarkOptions};
//! use tpt_fem_sparse::Coo;
//!
//! // SDOF oscillator: M=1, K=4 (ω=2). Free vibration from u0=1, v0=0.
//! let m = Coo { rows: vec![0], cols: vec![0], vals: vec![1.0] };
//! let k = Coo { rows: vec![0], cols: vec![0], vals: vec![4.0] };
//! let c = Coo::new();
//! let nsteps = 200;
//! let opts = NewmarkOptions { dt: 0.01, beta: 0.25, gamma: 0.5 };
//! let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);
//! let (t, u) = hist[nsteps].clone();
//! // Closed form u(t) = cos(ω t), ω = 2.
//! let want = (2.0 * t).cos();
//! assert!((u[0] - want).abs() < 1e-3, "got {} want {}", u[0], want);
//! ```

use tpt_fem_sparse::{solve, Coo};

/// Errors returned by the time-integration routines.
#[derive(Debug)]
pub enum DynamicError {
    /// The explicit central-difference step exceeds the critical (CFL) timestep
    /// `2/ω_max`, where `ω_max` is the largest natural frequency of the
    /// (lumped-mass) `M⁻¹K` system. Integration above this limit is unstable and
    /// diverges, so the step is rejected rather than silently producing
    /// meaningless displacements.
    CflViolation { dt: f64, critical: f64 },
    /// The underlying sparse solve (e.g. the modal eigen solve in
    /// [`modal_frequency_response`]) failed.
    Sparse(tpt_fem_sparse::SparseError),
    /// A caller-supplied parameter is invalid (e.g. a negative damping ratio,
    /// or a stiffness that is not positive-definite on the free DOFs).
    InvalidInput(String),
}

impl std::fmt::Display for DynamicError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DynamicError::CflViolation { dt, critical } => write!(
                f,
                "central-difference timestep {dt:e} exceeds the critical (CFL) limit {critical:e}"
            ),
            DynamicError::Sparse(e) => write!(f, "sparse solve failed: {e}"),
            DynamicError::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
        }
    }
}

impl std::error::Error for DynamicError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DynamicError::Sparse(e) => Some(e),
            _ => None,
        }
    }
}

impl From<tpt_fem_sparse::SparseError> for DynamicError {
    fn from(e: tpt_fem_sparse::SparseError) -> Self {
        DynamicError::Sparse(e)
    }
}

/// Return a new [`Coo`] equal to `coo` scaled by `s`.
pub fn coo_scale(coo: &Coo, s: f64) -> Coo {
    let mut out = Coo::new();
    for i in 0..coo.rows.len() {
        out.push(coo.rows[i], coo.cols[i], coo.vals[i] * s);
    }
    out
}

/// Return a new [`Coo`] equal to `a + b` (duplicate entries summed on collapse).
pub fn coo_add(a: &Coo, b: &Coo) -> Coo {
    let mut out = Coo::new();
    for i in 0..a.rows.len() {
        out.push(a.rows[i], a.cols[i], a.vals[i]);
    }
    for i in 0..b.rows.len() {
        out.push(b.rows[i], b.cols[i], b.vals[i]);
    }
    out
}

/// Matrix–vector product `y = A x` for a [`Coo`] matrix.
///
/// The output has the same length as `x` (a square `n×n` matrix times an
/// `n`-vector); an empty matrix contributes the zero vector, which makes this
/// safe to call with a zero/placeholder damping or mass matrix of any dimension.
/// Matrix–vector product `y = A x` for a [`Coo`] matrix.
///
/// The output has the same length as `x` (a square `n×n` matrix times an
/// `n`-vector); an empty matrix contributes the zero vector, which makes this
/// safe to call with a zero/placeholder damping or mass matrix of any dimension.
///
/// # Cost note
///
/// This converts `coo` to CSR on **every** call. Inside time-stepping or
/// iterative loops, convert once with [`Coo::to_csr`] and use
/// [`Csr::matvec`](tpt_fem_sparse::Csr::matvec) instead — that is what the
/// integrators in this crate do internally.
pub fn coo_matvec(coo: &Coo, x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if coo.rows.is_empty() {
        return vec![0.0; n];
    }
    coo.to_csr().matvec(x)
}

/// Options for the [`newmark`] integrator.
#[derive(Clone, Copy, Debug)]
pub struct NewmarkOptions {
    /// Time step.
    pub dt: f64,
    /// Newmark `β` (default `0.25` = average acceleration).
    pub beta: f64,
    /// Newmark `γ` (default `0.5`).
    pub gamma: f64,
}

impl Default for NewmarkOptions {
    fn default() -> Self {
        NewmarkOptions {
            dt: 0.01,
            beta: 0.25,
            gamma: 0.5,
        }
    }
}

/// Implicit Newmark-beta time integration of `M·ü + C·v + K·u = f(t)`.
///
/// `f(t)` returns the global load vector at time `t`. Returns the displacement
/// history as `(t, u)` pairs for steps `0..=nsteps` (step 0 is the initial
/// state).
pub fn newmark(
    mass: &Coo,
    damping: &Coo,
    stiffness: &Coo,
    u0: &[f64],
    v0: &[f64],
    f: impl Fn(f64) -> Vec<f64>,
    opts: &NewmarkOptions,
    nsteps: usize,
) -> Vec<(f64, Vec<f64>)> {
    let n = u0.len();
    let dt = opts.dt;
    let b = opts.beta;
    let g = opts.gamma;

    // Convert the constant operators to CSR once; every step below does two
    // matvecs against damping and stiffness, and re-converting per call would
    // dominate the runtime on anything but toy problems.
    let csr_d = damping.to_csr();
    let csr_k = stiffness.to_csr();

    // Initial acceleration a0 = M⁻¹ (f0 - C v0 - K u0).
    let f0 = f(0.0);
    let r0: Vec<f64> = (0..n)
        .map(|i| f0[i] - csr_d.matvec(v0)[i] - csr_k.matvec(u0)[i])
        .collect();
    let a0 = solve(mass, &r0).expect("mass matrix must be invertible");

    // Effective stiffness is constant; build it once.
    let k_hat = coo_add(
        stiffness,
        &coo_add(
            &coo_scale(mass, 1.0 / (b * dt * dt)),
            &coo_scale(damping, g / (b * dt)),
        ),
    );
    let csr_m = mass.to_csr();

    let mut history = Vec::with_capacity(nsteps + 1);
    let mut u = u0.to_vec();
    let mut v = v0.to_vec();
    let mut a = a0;
    history.push((0.0, u.clone()));

    for step in 1..=nsteps {
        let t = step as f64 * dt;
        let ft = f(t);

        let p_m: Vec<f64> = (0..n)
            .map(|i| u[i] / (b * dt * dt) + v[i] / (b * dt) + (0.5 / b - 1.0) * a[i])
            .collect();
        let p_c: Vec<f64> = (0..n)
            .map(|i| g / (b * dt) * u[i] + (g / b - 1.0) * v[i] + dt * (g / (2.0 * b) - 1.0) * a[i])
            .collect();
        let rhs: Vec<f64> = (0..n)
            .map(|i| ft[i] + csr_m.matvec(&p_m)[i] + csr_d.matvec(&p_c)[i])
            .collect();

        let u_new = solve(&k_hat, &rhs).expect("effective stiffness must be invertible");

        let a_new: Vec<f64> = (0..n)
            .map(|i| (u_new[i] - u[i]) / (b * dt * dt) - v[i] / (b * dt) - (0.5 / b - 1.0) * a[i])
            .collect();
        let v_new: Vec<f64> = (0..n)
            .map(|i| v[i] + dt * ((1.0 - g) * a[i] + g * a_new[i]))
            .collect();

        u = u_new;
        v = v_new;
        a = a_new;
        history.push((t, u.clone()));
    }
    history
}

/// Options for the [`central_difference`] integrator.
#[derive(Clone, Copy, Debug)]
pub struct CentralOptions {
    /// Time step.
    pub dt: f64,
}

impl Default for CentralOptions {
    fn default() -> Self {
        CentralOptions { dt: 0.01 }
    }
}

/// Explicit central-difference time integration of `M·ü + C·v + K·u = f(t)`.
///
/// The (possibly consistent) mass matrix is internally lumped to a diagonal so
/// that `M⁻¹` is a cheap component-wise division — the standard choice for an
/// explicit scheme. Returns the displacement history `(t, u)` for
/// `0..=nsteps`, or [`DynamicError::CflViolation`] if the step exceeds the
/// conditionally-stable critical timestep.
pub fn central_difference(
    mass: &Coo,
    damping: &Coo,
    stiffness: &Coo,
    u0: &[f64],
    v0: &[f64],
    f: impl Fn(f64) -> Vec<f64>,
    opts: &CentralOptions,
    nsteps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, DynamicError> {
    let n = u0.len();
    let dt = opts.dt;
    let csr_m = mass.to_csr();
    let mut mdiag = vec![0.0; n];
    for r in 0..n {
        let mut s = 0.0;
        for c in csr_m.row_ptrs[r]..csr_m.row_ptrs[r + 1] {
            s += csr_m.values[c];
        }
        mdiag[r] = s;
    }
    let inv = |x: &[f64]| (0..n).map(|i| x[i] / mdiag[i]).collect::<Vec<_>>();

    // Explicit central-difference is only conditionally stable: the step must
    // not exceed the critical (CFL) timestep `2/ω_max`, with `ω_max` the
    // largest natural frequency of `M⁻¹K`. A conservative guard uses the
    // infinity norm `‖M⁻¹K‖_∞ = max_i ( Σ_j |K_ij| / m_i ) ≥ ρ(M⁻¹K)`, giving a
    // *lower* bound on the true critical timestep; any step above it is
    // therefore definitely unstable and is rejected rather than silently
    // diverging.
    let mut rho_max = 0.0_f64;
    let mut has_mass = false;
    let csr_k = stiffness.to_csr();
    for i in 0..n {
        if mdiag[i] > 0.0 {
            let mut row_sum = 0.0_f64;
            for idx in csr_k.row_ptrs[i]..csr_k.row_ptrs[i + 1] {
                row_sum += csr_k.values[idx].abs();
            }
            rho_max = rho_max.max(row_sum / mdiag[i]);
            has_mass = true;
        }
    }
    let dt_safe = if has_mass && rho_max > 0.0 {
        2.0 / rho_max.sqrt()
    } else {
        f64::INFINITY
    };
    if dt > dt_safe {
        return Err(DynamicError::CflViolation {
            dt,
            critical: dt_safe,
        });
    }

    let f0 = f(0.0);
    let csr_d = damping.to_csr();
    let a0: Vec<f64> = inv(&(0..n)
        .map(|i| f0[i] - csr_d.matvec(v0)[i] - csr_k.matvec(u0)[i])
        .collect::<Vec<_>>());

    let mut history = Vec::with_capacity(nsteps + 1);
    let mut u = u0.to_vec();
    let mut v_half: Vec<f64> = (0..n).map(|i| v0[i] - 0.5 * dt * a0[i]).collect();
    history.push((0.0, u.clone()));

    for step in 1..=nsteps {
        let t = step as f64 * dt;
        let ft = f(t);
        let a: Vec<f64> = inv(&(0..n)
            .map(|i| ft[i] - csr_d.matvec(&v_half)[i] - csr_k.matvec(&u)[i])
            .collect::<Vec<_>>());
        let mut v_next = v_half.clone();
        for i in 0..n {
            v_next[i] += dt * a[i];
        }
        let mut u_next = u.clone();
        for i in 0..n {
            u_next[i] += dt * v_next[i];
        }
        u = u_next;
        v_half = v_next;
        history.push((t, u.clone()));
    }
    Ok(history)
}

/// Modal frequency-response (harmonic) analysis result.
///
/// `displacement_real[k][d]` / `displacement_imag[k][d]` hold the real and
/// imaginary parts of the steady-state displacement amplitude at DOF `d` when
/// the structure is excited by the harmonic force `f·cos(Ω t)` with
/// `Ω = 2π·frequencies[k]` and a uniform modal damping ratio `zeta`.
#[derive(Clone, Debug)]
pub struct ModalFrequencyResponse {
    /// Excitation frequencies in Hz (same order as the caller's slice).
    pub frequencies: Vec<f64>,
    /// Natural frequencies of the retained modes, in Hz.
    pub mode_frequencies_hz: Vec<f64>,
    /// Real part of the displacement amplitude per excitation frequency.
    pub displacement_real: Vec<Vec<f64>>,
    /// Imaginary part of the displacement amplitude per excitation frequency.
    pub displacement_imag: Vec<Vec<f64>>,
}

impl ModalFrequencyResponse {
    /// Magnitude `|u_d(Ω)|` of the displacement amplitude at DOF `d` for
    /// excitation index `k`.
    pub fn magnitude(&self, k: usize, d: usize) -> f64 {
        let re = self.displacement_real[k][d];
        let im = self.displacement_imag[k][d];
        (re * re + im * im).sqrt()
    }
}

/// Frequency-response workflow combining [`tpt-fem-eigen`] modal analysis with
/// modal superposition — the frequency-domain counterpart of stepping
/// [`newmark`] to steady state for every excitation frequency at once.
///
/// The generalized eigenproblem `K φ = ω² M φ` is solved for the `n_modes`
/// lowest modes (M-orthonormal shift-invert Lanczos). With the classical
/// modal-damping assumption (mode `i` damped at `2 ζ ω_i`), the steady-state
/// response to an in-phase harmonic force `f cos(Ω t)` is the modal sum
///
/// ```text
/// u(Ω) = Σᵢ φᵢ (φᵢᵀ f) / (ωᵢ² − Ω² + 2i ζ ωᵢ Ω)
/// ```
///
/// evaluated here in real/imaginary parts without needing a complex type.
/// Frequencies are given in Hz; internally they are converted to rad/s.
///
/// Modes beyond `n_modes` are truncated, so results are accurate only when
/// the response is dominated by the lowest modes (add static-residual or
/// missing-mass correction for more demanding cases — not implemented here).
pub fn modal_frequency_response(
    k: &Coo,
    m: &Coo,
    force: &[f64],
    zeta: f64,
    n_modes: usize,
    frequencies_hz: &[f64],
    lanczos_dim: usize,
) -> Result<ModalFrequencyResponse, DynamicError> {
    if zeta.is_nan() || zeta < 0.0 {
        return Err(DynamicError::InvalidInput(
            "zeta must be non-negative".into(),
        ));
    }
    let pairs = tpt_fem_eigen::generalized_lanczos_eigs(k, m, 0.0, n_modes, lanczos_dim)
        .map_err(DynamicError::Sparse)?;
    let mut mode_freqs = Vec::with_capacity(pairs.len());
    let mut modes: Vec<&Vec<f64>> = Vec::with_capacity(pairs.len());
    let mut modal_force = Vec::with_capacity(pairs.len());
    for (lam, phi) in &pairs {
        if *lam <= 0.0 {
            return Err(DynamicError::InvalidInput(format!(
                "eigen solver returned a non-positive eigenvalue {lam}; \
                 check that K is positive-definite on the free DOFs"
            )));
        }
        mode_freqs.push(lam.sqrt() / (2.0 * std::f64::consts::PI));
        modes.push(phi);
        let q: f64 = phi.iter().zip(force).map(|(a, b)| a * b).sum();
        modal_force.push(q);
    }

    let mut out = ModalFrequencyResponse {
        frequencies: frequencies_hz.to_vec(),
        mode_frequencies_hz: mode_freqs,
        displacement_real: Vec::with_capacity(frequencies_hz.len()),
        displacement_imag: Vec::with_capacity(frequencies_hz.len()),
    };
    for &f_hz in frequencies_hz {
        let omega = std::f64::consts::TAU * f_hz;
        let mut ur = vec![0.0; force.len()];
        let mut ui = vec![0.0; force.len()];
        for (i, phi) in modes.iter().enumerate() {
            let w2 = pairs[i].0;
            let denom_r = w2 - omega * omega;
            let denom_i = 2.0 * zeta * w2.sqrt() * omega;
            let denom = denom_r * denom_r + denom_i * denom_i;
            // u += phi · q_i · conj(denom)/|denom|²
            let scale_r = modal_force[i] * denom_r / denom;
            let scale_i = -modal_force[i] * denom_i / denom;
            for (d, val) in phi.iter().enumerate() {
                ur[d] += val * scale_r;
                ui[d] += val * scale_i;
            }
        }
        out.displacement_real.push(ur);
        out.displacement_imag.push(ui);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sdof(mass: f64, stiff: f64) -> (Coo, Coo, Coo) {
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
        let c = Coo::new();
        (m, c, k)
    }

    #[test]
    fn newmark_matches_closed_form_sdof() {
        // ω = 2, u(t) = cos(2t) from u0=1, v0=0.
        let (m, c, k) = sdof(1.0, 4.0);
        let opts = NewmarkOptions {
            dt: 0.005,
            beta: 0.25,
            gamma: 0.5,
        };
        let nsteps = 400; // t = 2.0
        let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, nsteps);
        let (t, u) = hist[nsteps].clone();
        let want = (2.0 * t).cos();
        assert!((u[0] - want).abs() < 5e-3, "got {} want {}", u[0], want);
    }

    #[test]
    fn central_difference_rejects_over_cfl_step() {
        // ω = 2, critical dt = 2/ω = 1.0. A step of 2.0 is well above it and
        // must be rejected rather than silently diverging.
        let (m, c, k) = sdof(1.0, 4.0);
        let opts = CentralOptions { dt: 2.0 };
        let res = central_difference(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, 10);
        assert!(res.is_err(), "over-CFL step must be rejected");
    }

    #[test]
    fn newmark_conserves_energy_undamped() {
        // Undamped free vibration: total energy E = ½ v² + ½ k u² is constant.
        let (m, c, k) = sdof(2.0, 8.0); // ω = 2
        let opts = NewmarkOptions {
            dt: 0.002,
            beta: 0.25,
            gamma: 0.5,
        };
        let hist = newmark(&m, &c, &k, &[1.0], &[0.0], |_| vec![0.0], &opts, 500);
        let dt = opts.dt;
        let e0 = 0.5 * 8.0 * 1.0_f64.powi(2); // u0=1, v0=0 -> E0 = ½ k
        for w in 1..hist.len() - 1 {
            let v_mid = (hist[w + 1].1[0] - hist[w - 1].1[0]) / (2.0 * dt);
            let e = 0.5 * 2.0 * v_mid * v_mid + 0.5 * 8.0 * hist[w].1[0] * hist[w].1[0];
            assert!((e - e0).abs() < 1e-2, "energy drift: {} vs {}", e, e0);
        }
    }

    fn diag2() -> (Coo, Coo) {
        let mut k = Coo::new();
        let mut m = Coo::new();
        k.push(0, 0, 4.0);
        m.push(0, 0, 1.0);
        (k, m)
    }

    #[test]
    fn modal_frf_sdof_static_and_resonance() {
        // SDOF: K=4, M=1 (ω=2 rad/s), ζ=0.05, unit force.
        // At Ω=0 the response is the static deflection 1/K = 0.25 (real).
        // At resonance Ω=ω the magnitude is 1/(2ζK) = 2.5 and purely
        // imaginary (90° phase lag), negative imaginary part.
        let (k, m) = diag2();
        let f_res = 2.0 / std::f64::consts::TAU;
        let resp = modal_frequency_response(&k, &m, &[1.0], 0.05, 1, &[0.0, f_res], 1).unwrap();
        assert!((resp.displacement_real[0][0] - 0.25).abs() < 1e-12);
        assert!(resp.displacement_imag[0][0].abs() < 1e-12);
        let mag = resp.magnitude(1, 0);
        assert!((mag - 2.5).abs() < 1e-6, "resonance magnitude {mag}");
        assert!((resp.displacement_real[1][0]).abs() < 1e-8);
        assert!(
            (resp.displacement_imag[1][0] + 2.5).abs() < 1e-4,
            "imag {}",
            resp.displacement_imag[1][0]
        );
    }

    #[test]
    fn modal_frf_two_dof_matches_closed_form() {
        // K = [[5,1],[1,5]], M = I has exact modes at ω² = 4, 6 with
        // φ₁ = [1,1]/√2 and φ₂ = [1,-1]/√2. With force f = [1,0] and
        // ζ = 0.02, evaluate the two-mode sum in closed form at Ω = 1 rad/s
        // and compare against the solver.
        let mut k = Coo::new();
        let mut m = Coo::new();
        k.push(0, 0, 5.0);
        k.push(1, 1, 5.0);
        k.push(0, 1, 1.0);
        k.push(1, 0, 1.0);
        m.push(0, 0, 1.0);
        m.push(1, 1, 1.0);
        let s2 = 1.0 / 2.0_f64.sqrt();
        let zeta = 0.02;
        let omega = 1.0_f64;
        // Closed form. Eigenvalue 6 pairs with [1,1]/√2 ("A"),
        // eigenvalue 4 pairs with [1,-1]/√2 ("B").
        let q = s2; // phiᵀ [1,0] for both modes
        let d = |w2: f64| {
            let r = w2 - omega * omega;
            let i = 2.0 * zeta * w2.sqrt() * omega;
            (r, i)
        };
        let (ra, ia) = d(6.0);
        let (rb, ib) = d(4.0);
        let ma = ra * ra + ia * ia;
        let mb = rb * rb + ib * ib;
        let expect = [
            [
                s2 * (q * ra / ma + q * rb / mb),
                -s2 * (q * ia / ma + q * ib / mb),
            ],
            [
                s2 * (q * ra / ma - q * rb / mb),
                -s2 * (q * ia / ma - q * ib / mb),
            ],
        ];
        let freqs = [omega / std::f64::consts::TAU];
        let resp = modal_frequency_response(&k, &m, &[1.0, 0.0], zeta, 2, &freqs, 2).unwrap();
        for kdof in 0..2 {
            assert!(
                (resp.displacement_real[0][kdof] - expect[kdof][0]).abs() < 1e-10,
                "real dof{kdof}: {} vs {}",
                resp.displacement_real[0][kdof],
                expect[kdof][0]
            );
            assert!(
                (resp.displacement_imag[0][kdof] - expect[kdof][1]).abs() < 1e-10,
                "imag dof{kdof}: {} vs {}",
                resp.displacement_imag[0][kdof],
                expect[kdof][1]
            );
        }
        // Mode frequencies reported back in Hz.
        assert!(
            (resp.mode_frequencies_hz[0] - 4.0_f64.sqrt() / std::f64::consts::TAU).abs() < 1e-8
        );
    }

    #[test]
    fn modal_frf_rejects_bad_input() {
        let (k, m) = diag2();
        assert!(modal_frequency_response(&k, &m, &[1.0], -0.1, 1, &[0.0], 1).is_err());
    }
}
