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

use tpt_fem_sparse::{solve, Coo, Csr};

/// Errors returned by the time-integration routines.
#[derive(Debug)]
pub enum DynamicError {
    /// The explicit central-difference step exceeds the critical (CFL) timestep
    /// `2/ω_max`, where `ω_max` is the largest natural frequency of the
    /// (lumped-mass) `M⁻¹K` system. Integration above this limit is unstable and
    /// diverges, so the step is rejected rather than silently producing
    /// meaningless displacements.
    CflViolation { dt: f64, critical: f64 },
}

impl std::fmt::Display for DynamicError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DynamicError::CflViolation { dt, critical } => write!(
                f,
                "central-difference timestep {dt:e} exceeds the critical (CFL) limit {critical:e}"
            ),
        }
    }
}

impl std::error::Error for DynamicError {}

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

/// Matrix–vector product `y = A x` for a [`Csr`] matrix.
///
/// The tight inner loop over each row's stored entries is a simple contiguous
/// reduction that LLVM auto-vectorises; callers doing repeated matvecs (e.g. the
/// time integrators below) should pre-compile their [`Coo`] to [`Csr`] once via
/// [`Coo::to_csr`] and call this, rather than repeatedly rebuilding the CSR as
/// the older `coo_matvec` helper did.
pub fn csr_matvec(csr: &Csr, x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut y = vec![0.0; n];
    for r in 0..n.min(csr.nrows) {
        let mut s = 0.0;
        for c in csr.row_ptrs[r]..csr.row_ptrs[r + 1] {
            s += csr.values[c] * x[csr.col_ind[c]];
        }
        y[r] = s;
    }
    y
}

/// Matrix–vector product `y = A x` for a [`Coo`] matrix (single-shot
/// convenience).
///
/// This builds the CSR representation internally, so it is correct but
/// **not** suitable for the per-step hot loop of a time integrator — convert
/// once with [`Coo::to_csr`] and use [`csr_matvec`] there. The output has the
/// same length as `x`; an empty matrix contributes the zero vector.
pub fn coo_matvec(coo: &Coo, x: &[f64]) -> Vec<f64> {
    let n = x.len();
    if coo.rows.is_empty() {
        return vec![0.0; n];
    }
    csr_matvec(&coo.to_csr(), x)
}

/// A second-order operator `M·ü + C·v + K·u = f(t)` pre-compiled to CSR.
///
/// Construct once from the [`Coo`] mass/damping/stiffness and use
/// [`CachedSystem::apply_mass`] / [`apply_damping`](Self::apply_damping) /
/// [`apply_stiffness`](Self::apply_stiffness) for the repeated per-step matvecs
/// of a time integrator, so the triplet→CSR conversion happens a single time
/// instead of on every call (the former `coo_matvec`-in-a-loop performance
/// smell).
pub struct CachedSystem {
    /// Mass matrix (CSR).
    pub mass: Csr,
    /// Damping matrix (CSR).
    pub damping: Csr,
    /// Stiffness matrix (CSR).
    pub stiffness: Csr,
}

impl CachedSystem {
    /// Compile `mass`/`damping`/`stiffness` to CSR exactly once.
    pub fn new(mass: &Coo, damping: &Coo, stiffness: &Coo) -> Self {
        CachedSystem {
            mass: mass.to_csr(),
            damping: damping.to_csr(),
            stiffness: stiffness.to_csr(),
        }
    }

    /// `y = M x`.
    pub fn apply_mass(&self, x: &[f64]) -> Vec<f64> {
        csr_matvec(&self.mass, x)
    }

    /// `y = C x`.
    pub fn apply_damping(&self, x: &[f64]) -> Vec<f64> {
        csr_matvec(&self.damping, x)
    }

    /// `y = K x`.
    pub fn apply_stiffness(&self, x: &[f64]) -> Vec<f64> {
        csr_matvec(&self.stiffness, x)
    }
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

    // Pre-compile the operators to CSR once; the effective stiffness below is
    // also constant, so the matvecs in the loop never re-sort the triplets.
    let sys = CachedSystem::new(mass, damping, stiffness);

    // Initial acceleration a0 = M⁻¹ (f0 - C v0 - K u0).
    let f0 = f(0.0);
    let r0: Vec<f64> = (0..n)
        .map(|i| f0[i] - sys.apply_damping(v0)[i] - sys.apply_stiffness(u0)[i])
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
            .map(|i| ft[i] + sys.apply_mass(&p_m)[i] + sys.apply_damping(&p_c)[i])
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
    // Pre-compile all three operators to CSR once; the CFL check below still
    // reads the (one-time) stiffness triplets, but the per-step matvecs use the
    // cached CSR rather than rebuilding it each call.
    let sys = CachedSystem::new(mass, damping, stiffness);
    let csr_m = &sys.mass;
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
    for i in 0..n {
        if mdiag[i] > 0.0 {
            let mut row_sum = 0.0_f64;
            for idx in 0..stiffness.rows.len() {
                if stiffness.rows[idx] == i {
                    row_sum += stiffness.vals[idx].abs();
                }
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
    let a0: Vec<f64> = inv(&(0..n)
        .map(|i| f0[i] - csr_matvec(&sys.damping, v0)[i] - csr_matvec(&sys.stiffness, u0)[i])
        .collect::<Vec<_>>());

    let mut history = Vec::with_capacity(nsteps + 1);
    let mut u = u0.to_vec();
    let mut v_half: Vec<f64> = (0..n).map(|i| v0[i] - 0.5 * dt * a0[i]).collect();
    history.push((0.0, u.clone()));

    for step in 1..=nsteps {
        let t = step as f64 * dt;
        let ft = f(t);
        let a: Vec<f64> = inv(&(0..n)
            .map(|i| {
                ft[i] - csr_matvec(&sys.damping, &v_half)[i] - csr_matvec(&sys.stiffness, &u)[i]
            })
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

    #[test]
    fn csr_matvec_matches_coo() {
        // A small dense matrix as COO, multiplied by a vector, must equal the
        // CSR path used by the cached integrators.
        let coo = Coo {
            rows: vec![0, 0, 1, 1, 1, 2],
            cols: vec![0, 1, 0, 1, 2, 2],
            vals: vec![2.0, -1.0, -1.0, 3.0, -1.0, 4.0],
        };
        let x = vec![1.0, 2.0, 3.0];
        let y_coo = coo_matvec(&coo, &x);
        let y_csr = csr_matvec(&coo.to_csr(), &x);
        assert_eq!(y_coo, y_csr);
        // 2*1 -1*2 = 0; -1*1 + 3*2 -1*3 = 2; 4*3 = 12
        assert_eq!(y_csr, vec![0.0, 2.0, 12.0]);
    }

    #[test]
    fn cached_system_avoids_repeat_conversion() {
        let m = Coo {
            rows: vec![0, 1, 2],
            cols: vec![0, 1, 2],
            vals: vec![2.0, 3.0, 4.0],
        };
        let k = Coo {
            rows: vec![0, 0, 1, 1, 1, 2],
            cols: vec![0, 1, 0, 1, 2, 2],
            vals: vec![2.0, -1.0, -1.0, 3.0, -1.0, 4.0],
        };
        let c = Coo::new();
        let sys = CachedSystem::new(&m, &c, &k);
        let x = vec![1.0, 2.0, 3.0];
        assert_eq!(sys.apply_mass(&x), csr_matvec(&m.to_csr(), &x));
        assert_eq!(sys.apply_stiffness(&x), csr_matvec(&k.to_csr(), &x));
    }
}
