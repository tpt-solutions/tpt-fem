//! J2 (von Mises) plasticity for `tpt-fem`.
//!
//! Provides the radial-return (closest-point projection) algorithm for
//! associative J2 plasticity with combined isotropic + kinematic hardening,
//! operating on a single integration point in Voigt `6`-vector form. The
//! constitutive update is wired into [`tpt-fem-solve`]'s Newton loop by
//! [`solve_elastic_plastic_rod`], a force-controlled 1-D rod solver that
//! carries per-element plastic state.
//!
//! ```
//! use tpt_fem_plasticity::{j2_return_mapping, PlasticityParams, PlasticState};
//!
//! // Elastic step: stress is the isotropic elastic response.
//! let p = PlasticityParams::steel();
//! let s0 = PlasticState::new();
//! let (stress, _) = j2_return_mapping(&p, &[1e-4, -0.3e-4, -0.3e-4, 0.0, 0.0, 0.0], &s0);
//! let e = 200e9;
//! assert!((stress[0] - e * 1e-4).abs() < 1.0);
//! ```

/// Material parameters for J2 plasticity.
#[derive(Clone, Copy, Debug)]
pub struct PlasticityParams {
    /// Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Initial yield stress `σ_y⁰`.
    pub yield_stress: f64,
    /// Isotropic hardening modulus `H_iso`.
    pub iso_hardening: f64,
    /// Kinematic hardening modulus `H_kin`.
    pub kin_hardening: f64,
}

impl PlasticityParams {
    /// A representative structural steel: `E = 200 GPa`, `ν = 0.3`,
    /// `σ_y = 250 MPa`, perfect plasticity.
    pub fn steel() -> Self {
        PlasticityParams {
            young: 200e9,
            poisson: 0.3,
            yield_stress: 250e6,
            iso_hardening: 0.0,
            kin_hardening: 0.0,
        }
    }
}

/// Per-integration-point plastic state (Voigt `6`-vectors).
#[derive(Clone, Debug)]
pub struct PlasticState {
    /// Plastic strain `εᵖ` (Voigt `6`, *engineering* shear `γᵖ = 2εᵖ_ij`, to match
    /// the convention of [`elastic_stiffness`]).
    ///
    /// This is what makes the update an incremental flow-theory model: the trial
    /// stress is the elastic predictor `C:(ε − εᵖ)`, so unloading from a plastic
    /// state is elastic and reversed cycles behave correctly.
    pub plastic_strain: Vec<f64>,
    /// Back stress `α` (deviatoric, Voigt `6`).
    pub back_stress: Vec<f64>,
    /// Equivalent plastic strain `ε̄ᵖ`.
    pub eps_eq: f64,
}

impl PlasticState {
    /// A fresh, virgin state (no plastic strain, no hardening, no back stress).
    pub fn new() -> Self {
        PlasticState {
            plastic_strain: vec![0.0; 6],
            back_stress: vec![0.0; 6],
            eps_eq: 0.0,
        }
    }
}

impl Default for PlasticState {
    fn default() -> Self {
        Self::new()
    }
}

/// Isotropic elastic stiffness (Voigt `6×6`) for `E`, `ν`.
pub fn elastic_stiffness(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let mu = e / (2.0 * (1.0 + nu));
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mut c = [[0.0; 6]; 6];
    let diag = [lam + 2.0 * mu, lam + 2.0 * mu, lam + 2.0 * mu];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = if i == j { diag[i] } else { lam };
        }
    }
    for i in 3..6 {
        c[i][i] = mu;
    }
    c
}

/// Deviatoric part of a Voigt `6` symmetric tensor (shear components unchanged).
fn deviatoric(s: &[f64; 6]) -> [f64; 6] {
    let mean = (s[0] + s[1] + s[2]) / 3.0;
    let mut d = *s;
    d[0] -= mean;
    d[1] -= mean;
    d[2] -= mean;
    d
}

/// Frobenius norm of the deviatoric part using the Voigt shear doubling:
/// `||dev s||² = s₀²+s₁²+s₂² + 2(s₃²+s₄²+s₅²)`.
fn dev_norm(s: &[f64; 6]) -> f64 {
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
        .sqrt()
}

fn matvec6(c: &[[f64; 6]; 6], v: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0; 6];
    for i in 0..6 {
        let mut s = 0.0;
        for j in 0..6 {
            s += c[i][j] * v[j];
        }
        out[i] = s;
    }
    out
}

/// J2 radial-return update for a single integration point.
///
/// Given the *total* strain (Voigt `6`) and the previous state, returns the
/// updated stress (Voigt `6`) and the new plastic state. The trial stress is the
/// elastic predictor `σ_tr = C:(ε − εᵖ)`; when it lies inside the yield surface
/// the step is elastic, otherwise the deviatoric stress is projected back onto
/// the surface along its normal (closest-point / radial return) and the plastic
/// strain, back stress, and equivalent plastic strain are all advanced.
///
/// ```text
/// σ_tr = C:(ε − εᵖ)
/// η    = dev(σ_tr) − α
/// f    = √(3/2)·‖η‖ − (σ_y⁰ + H_iso ε̄ᵖ)
/// Δγ   = f / [ 2G·√(3/2) + √(2/3)·(H_iso + H_kin) ]      (when f > 0)
/// ```
///
/// Because `εᵖ` is carried in [`PlasticState`], unloading from a plastic state is
/// elastic and reversed cycles reproduce the Bauschinger effect under kinematic
/// hardening.
pub fn j2_return_mapping(
    p: &PlasticityParams,
    total_strain: &[f64; 6],
    state: &PlasticState,
) -> ([f64; 6], PlasticState) {
    let c = elastic_stiffness(p.young, p.poisson);
    // Elastic predictor from the *elastic* strain ε − εᵖ.
    let mut eps_el = [0.0; 6];
    for i in 0..6 {
        eps_el[i] = total_strain[i] - state.plastic_strain[i];
    }
    let sig_tr = matvec6(&c, &eps_el);
    let s_tr = deviatoric(&sig_tr);
    let mut eta = s_tr;
    for i in 0..6 {
        eta[i] -= state.back_stress[i];
    }
    let eta_norm = dev_norm(&eta);
    let q = (1.5_f64).sqrt() * eta_norm;
    let sig_y = p.yield_stress + p.iso_hardening * state.eps_eq;
    let f = q - sig_y;

    if f <= 0.0 || eta_norm < 1e-14 {
        return (sig_tr, state.clone());
    }

    let mu = p.young / (2.0 * (1.0 + p.poisson));
    // Radial-return increment from `f = √(3/2)·‖η‖ − σ_y(ε̄ᵖ)`.
    //
    // The deviatoric update is `s ← s − 2G·Δγ·n̂` and the back stress moves as
    // `α ← α + (2/3)·H_kin·Δγ·n̂`, so the *relative* stress contracts by
    // `(2G + (2/3)H_kin)·Δγ` while the radius grows by `H_iso·√(2/3)·Δγ`
    // (because `ε̄ᵖ ← ε̄ᵖ + √(2/3)·Δγ`). Enforcing `f = 0` after the update gives
    //
    //   Δγ = f / [ 2G·√(3/2) + √(2/3)·(H_iso + H_kin) ]
    //
    // Kinematic hardening therefore *does* enter the denominator: translating the
    // yield surface consumes part of the trial overshoot even though it leaves
    // the radius unchanged.
    let denom =
        2.0 * mu * (1.5_f64).sqrt() + (2.0 / 3.0_f64).sqrt() * (p.iso_hardening + p.kin_hardening);
    let dgamma = f / denom;

    let mut n = eta;
    for i in 0..6 {
        n[i] /= eta_norm;
    }
    let mut s = s_tr;
    for i in 0..6 {
        s[i] -= 2.0 * mu * dgamma * n[i];
    }
    let mean = (sig_tr[0] + sig_tr[1] + sig_tr[2]) / 3.0;
    let mut sig = s;
    sig[0] += mean;
    sig[1] += mean;
    sig[2] += mean;

    let mut alpha = state.back_stress.clone();
    for i in 0..6 {
        alpha[i] += (2.0 / 3.0) * p.kin_hardening * dgamma * n[i];
    }
    // Associative flow rule `Δεᵖ = Δγ·n̂`. The Voigt strain vector carries
    // *engineering* shear (`γ = 2ε_ij`), so the shear components pick up a
    // factor of two; this makes `σ = C:(ε − εᵖ)` reproduce the returned stress
    // exactly (see `stress_equals_elastic_predictor_of_plastic_strain`).
    let mut eps_p = state.plastic_strain.clone();
    for i in 0..6 {
        let voigt = if i < 3 { 1.0 } else { 2.0 };
        eps_p[i] += voigt * dgamma * n[i];
    }
    let eps_eq = state.eps_eq + (2.0 / 3.0_f64).sqrt() * dgamma;
    let new_state = PlasticState {
        plastic_strain: eps_p,
        back_stress: alpha,
        eps_eq,
    };
    (sig, new_state)
}

/// 1-D (uniaxial-*stress*) return-mapping law with combined isotropic +
/// kinematic hardening — the scalar analogue of [`j2_return_mapping`].
///
/// Returns `(σ, ε̄ᵖ, plastic)`, where `plastic` flags whether this step produced
/// plastic flow. Used by the rod solver's element constitutive update.
///
/// The update is the textbook 1-D radial return:
///
/// ```text
/// σ_tr = E (ε − εᵖ)              trial (elastic predictor)
/// η    = σ_tr − α,   α = H_kin εᵖ    relative stress (back stress α)
/// f    = |η| − (σ_y⁰ + H_iso ε̄ᵖ)     yield function
/// Δγ   = f / (E + H_iso + H_kin)     consistency parameter (f > 0)
/// σ    = σ_tr − E Δγ sign(η)
/// ```
///
/// so a perfectly-plastic material (`H_iso = H_kin = 0`) gives the exact stress
/// plateau `σ = σ_y⁰`, and linear isotropic hardening gives
/// `σ = σ_y⁰ + H_iso ε̄ᵖ` with `ε = σ/E + ε̄ᵖ`.
///
/// # Uniaxial stress vs uniaxial strain
///
/// This is *not* [`j2_return_mapping`] evaluated at
/// `ε = (ε_x, −ν ε_x, −ν ε_x, 0, 0, 0)`. That strain state is the uniaxial-stress
/// state only while the response is elastic: once the point yields, plastic flow
/// is incompressible, so holding the lateral strains at the elastic contraction
/// `−ν ε_x` develops spurious lateral stress and an artificial `E/3` hardening
/// slope. Enforcing genuine uniaxial stress in the 3-D law requires solving for
/// the lateral strain that makes `σ_y = σ_z = 0` (see the
/// `uniaxial_stress_consistency` test and the `uniaxial_sweep` example).
///
/// # Limitations
///
/// Only the *equivalent* plastic strain `ε̄ᵖ ≥ 0` is threaded between calls, not
/// its sign or an independent back stress. The signed plastic strain is therefore
/// reconstructed as `εᵖ = ε̄ᵖ·sign(ε)`, which is exact for monotonic
/// (proportional) loading — the regime the rod solver drives. For a full
/// reversed-cycle history with kinematic hardening, drive
/// [`j2_return_mapping`] and carry the complete [`PlasticState`] instead.
pub fn plastic_1d(
    params: &PlasticityParams,
    total_strain: f64,
    prev_eps_eq: f64,
) -> (f64, f64, bool) {
    let e = params.young;
    // Signed plastic strain, reconstructed under the monotonic-loading
    // assumption documented above.
    let sign_load = if total_strain >= 0.0 { 1.0 } else { -1.0 };
    let eps_p = prev_eps_eq * sign_load;

    let sigma_tr = e * (total_strain - eps_p);
    let alpha = params.kin_hardening * eps_p;
    let sig_y = params.yield_stress + params.iso_hardening * prev_eps_eq;
    let eta = sigma_tr - alpha;
    let f = eta.abs() - sig_y;

    if f <= 0.0 {
        return (sigma_tr, prev_eps_eq, false);
    }

    let h = params.iso_hardening + params.kin_hardening;
    let dgamma = f / (e + h);
    let sigma = sigma_tr - e * dgamma * eta.signum();
    (sigma, prev_eps_eq + dgamma, true)
}

/// Consistent tangent `dσ/dε` for [`plastic_1d`]: the elastic modulus `E` before
/// yield, and the combined-hardening slope `E·H/(E + H)` with
/// `H = H_iso + H_kin` once plastic.
///
/// Note that a perfectly-plastic material (`H = 0`) has a **zero** plastic
/// tangent, so a *force*-controlled solve cannot pass the yield load — the
/// problem is genuinely singular there. Use a non-zero hardening modulus (or
/// displacement control) past yield.
pub fn uniaxial_tangent_1d(params: &PlasticityParams, total_strain: f64, prev_eps_eq: f64) -> f64 {
    let (_, _, plastic) = plastic_1d(params, total_strain, prev_eps_eq);
    if !plastic {
        params.young
    } else {
        let h = params.iso_hardening + params.kin_hardening;
        params.young * h / (params.young + h)
    }
}

/// Errors returned by the plasticity solvers.
#[derive(Debug)]
pub enum PlasticityError {
    /// The coupled Newton iteration (per load step) did not converge. For a
    /// perfectly-plastic material under force control this is a genuinely
    /// singular tangent — the problem has no bounded solution past yield, so a
    /// hard panic is replaced by a propagated error.
    Newton(tpt_fem_solve::NewtonError),
}

impl std::fmt::Display for PlasticityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlasticityError::Newton(e) => {
                write!(f, "elastic-plastic rod Newton iteration failed: {e}")
            }
        }
    }
}

impl std::error::Error for PlasticityError {}

impl From<tpt_fem_solve::NewtonError> for PlasticityError {
    fn from(e: tpt_fem_solve::NewtonError) -> Self {
        PlasticityError::Newton(e)
    }
}

/// Solve a 1-D rod (built from `Line2` elements, one axial DOF per node) under a
/// sequence of *force* loads `P` applied at the last node, with the first node
/// fixed. Each load is solved with [`tpt-fem-solve`]'s Newton iteration using a
/// per-element 1-D J2 (radial-return) plastic state and the consistent tangent.
/// Returns the displacement vector at each load step, or a [`PlasticityError`]
/// if the Newton iteration fails to converge (e.g. a perfectly-plastic rod
/// pushed past yield under force control).
///
/// Because the loading is force-controlled, the material must harden
/// (`iso_hardening + kin_hardening > 0`) for loads beyond the yield force: a
/// perfectly-plastic rod has a zero tangent stiffness there and the problem has
/// no bounded solution.
pub fn solve_elastic_plastic_rod(
    mesh: &tpt_fem_mesh::Mesh,
    area: f64,
    params: &PlasticityParams,
    loads: &[f64],
) -> Result<Vec<Vec<f64>>, PlasticityError> {
    use std::cell::RefCell;
    use std::rc::Rc;

    let nnodes = mesh.node_count();
    let ndof = nnodes;
    let mut elem_len = Vec::with_capacity(mesh.elements.len());
    for elem in &mesh.elements {
        let a = mesh.node_coords(elem.nodes[0])[0];
        let b = mesh.node_coords(elem.nodes[1])[0];
        elem_len.push((b - a).abs());
    }
    // Previous-converged equivalent plastic strain per element, read immutably by
    // the Newton closures and committed after each converged load step.
    let base: Rc<RefCell<Vec<f64>>> = Rc::new(RefCell::new(vec![0.0; mesh.elements.len()]));

    let assemble = |u: &[f64]| -> (Vec<f64>, tpt_fem_sparse::Coo, Vec<f64>) {
        let b = base.borrow();
        let mut fint = vec![0.0; ndof];
        let mut coo = tpt_fem_sparse::Coo::new();
        let mut next = vec![0.0; mesh.elements.len()];
        for (e, elem) in mesh.elements.iter().enumerate() {
            let n0 = elem.nodes[0];
            let n1 = elem.nodes[1];
            let l = elem_len[e];
            let e_axial = (u[n1] - u[n0]) / l;
            let (sigma, eps_eq, _plastic) = plastic_1d(params, e_axial, b[e]);
            next[e] = eps_eq;
            let axial_force = sigma * area;
            fint[n0] -= axial_force;
            fint[n1] += axial_force;
            let et = uniaxial_tangent_1d(params, e_axial, b[e]);
            let ke = et * area / l;
            coo.push(n0, n0, ke);
            coo.push(n0, n1, -ke);
            coo.push(n1, n0, -ke);
            coo.push(n1, n1, ke);
        }
        (fint, coo, next)
    };

    let mut u = vec![0.0; ndof];
    // `NewtonOptions::tol` is an *absolute* residual-force norm, so it has to be
    // scaled to the problem: the internal force is evaluated as `σ·A`, which
    // carries a rounding floor of order `|σ·A|·ε_machine`. Tie the tolerance to
    // the largest applied load so the same code converges for a rod measured in
    // newtons and one measured in meganewtons.
    let load_scale = loads
        .iter()
        .fold(0.0_f64, |acc, p| acc.max(p.abs()))
        .max(1.0);
    let opts = tpt_fem_solve::NewtonOptions {
        tol: 1e-9 * load_scale,
        ..Default::default()
    };
    let mut history = Vec::with_capacity(loads.len());
    for &p_load in loads {
        let dirichlet = [(0usize, 0.0)];
        let u_sol = tpt_fem_solve::newton(
            &u,
            |u| {
                let (mut fint, _, _) = assemble(u);
                // Prescribed DOF (node 0) carries the reaction; zero its residual
                // so the Newton norm measures only free DOFs.
                fint[0] = 0.0;
                fint[nnodes - 1] -= p_load;
                fint
            },
            |u| {
                let (_, coo, _) = assemble(u);
                coo
            },
            &dirichlet,
            &opts,
        )?;
        let (_, _, next) = assemble(&u_sol);
        *base.borrow_mut() = next;
        u = u_sol;
        history.push(u.clone());
    }
    Ok(history)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    /// Uniaxial-*stress* strain vector: lateral strains free (`ε_y = ε_z = -ν ε_x`).
    fn uniaxial_stress_strain(epsx: f64, nu: f64) -> [f64; 6] {
        [epsx, -nu * epsx, -nu * epsx, 0.0, 0.0, 0.0]
    }

    fn rod(n_elements: usize, length: f64) -> tpt_fem_mesh::Mesh {
        let mut b = MeshBuilder::new();
        let mut prev = b.add_node(vec![0.0]);
        for i in 1..=n_elements {
            let node = b.add_node(vec![length * i as f64 / n_elements as f64]);
            b.add_element(CellType::Line, vec![prev, node]);
            prev = node;
        }
        b.build()
    }

    #[test]
    fn elastic_stress_is_linear() {
        let p = PlasticityParams::steel();
        let e = 1e-4;
        let strain = uniaxial_stress_strain(e, p.poisson);
        let (s, st) = j2_return_mapping(&p, &strain, &PlasticState::new());
        // Uniaxial stress: σ_x = E ε_x, σ_y = σ_z = 0.
        assert!((s[0] - p.young * e).abs() < 1.0);
        assert!(s[1].abs() < 1.0);
        assert_eq!(st.eps_eq, 0.0);
    }

    #[test]
    fn plastic_step_lands_on_yield_surface() {
        let p = PlasticityParams::steel();
        let strain = uniaxial_stress_strain(0.01, p.poisson); // far past yield
        let (s, st) = j2_return_mapping(&p, &strain, &PlasticState::new());
        let q = (1.5_f64).sqrt() * dev_norm(&deviatoric(&s));
        assert!((q - p.yield_stress).abs() / p.yield_stress < 1e-9, "q={q}");
        assert!(st.eps_eq > 0.0);
    }

    #[test]
    fn combined_hardening_lands_on_yield_surface() {
        // With kinematic hardening the yield surface *translates*, so the
        // consistency condition must be checked on the relative stress
        // `η = dev(σ) − α`, whose radius is `σ_y⁰ + H_iso ε̄ᵖ`.
        for (h_iso, h_kin) in [(0.0, 5e9), (3e9, 5e9), (5e9, 0.0)] {
            let p = PlasticityParams {
                iso_hardening: h_iso,
                kin_hardening: h_kin,
                ..PlasticityParams::steel()
            };
            let strain = uniaxial_stress_strain(0.01, p.poisson);
            let (s, st) = j2_return_mapping(&p, &strain, &PlasticState::new());
            let mut eta = deviatoric(&s);
            for i in 0..6 {
                eta[i] -= st.back_stress[i];
            }
            let q = (1.5_f64).sqrt() * dev_norm(&eta);
            let radius = p.yield_stress + h_iso * st.eps_eq;
            assert!(
                (q - radius).abs() / radius < 1e-9,
                "H_iso={h_iso:.1e} H_kin={h_kin:.1e}: q={q:.6e} vs radius={radius:.6e}"
            );
        }
    }

    /// Solve for the lateral strain that makes `σ_y = σ_z = 0` — genuine uniaxial
    /// *stress* — in the 3-D return-mapping law, and return the axial stress.
    ///
    /// Needed because `ε_lat = −ν ε_x` is the uniaxial-stress state only while the
    /// response is elastic; plastic flow is incompressible, so past yield the
    /// lateral strain relaxes towards `−ε_x/2`.
    fn uniaxial_stress_axial(p: &PlasticityParams, eps_x: f64) -> f64 {
        let lateral_stress = |lat: f64| {
            let (s, _) =
                j2_return_mapping(p, &[eps_x, lat, lat, 0.0, 0.0, 0.0], &PlasticState::new());
            s[1]
        };
        // Secant iteration between the elastic and the fully-plastic contraction.
        let mut a = -p.poisson * eps_x;
        let mut b = -0.5 * eps_x;
        let (mut fa, mut fb) = (lateral_stress(a), lateral_stress(b));
        for _ in 0..200 {
            if fb.abs() < 1e-6 || (fb - fa).abs() < 1e-30 {
                break;
            }
            let c = b - fb * (b - a) / (fb - fa);
            (a, fa) = (b, fb);
            (b, fb) = (c, lateral_stress(c));
        }
        let (s, _) = j2_return_mapping(p, &[eps_x, b, b, 0.0, 0.0, 0.0], &PlasticState::new());
        s[0]
    }

    #[test]
    fn stress_equals_elastic_predictor_of_plastic_strain() {
        // The returned stress must equal C:(ε − εᵖ) for the *updated* εᵖ — this
        // pins down the Voigt engineering-shear factor in the flow rule.
        let p = PlasticityParams {
            iso_hardening: 3e9,
            kin_hardening: 2e9,
            ..PlasticityParams::steel()
        };
        let c = elastic_stiffness(p.young, p.poisson);
        // A general (non-proportional) strain state with shear.
        let strain = [4e-3, -1e-3, -5e-4, 2e-3, 1e-3, -8e-4];
        let (sig, st) = j2_return_mapping(&p, &strain, &PlasticState::new());
        assert!(st.eps_eq > 0.0, "the step must be plastic");
        let mut eps_el = [0.0; 6];
        for i in 0..6 {
            eps_el[i] = strain[i] - st.plastic_strain[i];
        }
        let recomputed = matvec6(&c, &eps_el);
        for i in 0..6 {
            assert!(
                (recomputed[i] - sig[i]).abs() <= 1e-6 * p.yield_stress,
                "component {i}: C:(ε−εᵖ) = {:.6e} but stress = {:.6e}",
                recomputed[i],
                sig[i]
            );
        }
        // Plastic flow is incompressible: tr(εᵖ) = 0.
        let tr = st.plastic_strain[0] + st.plastic_strain[1] + st.plastic_strain[2];
        assert!(
            tr.abs() < 1e-18,
            "plastic strain must be deviatoric, tr={tr:.3e}"
        );
    }

    #[test]
    fn unloading_from_a_plastic_state_is_elastic() {
        // Load past yield in pure shear, then unload: the stress must retrace the
        // elastic slope G and accumulate no further plastic strain.
        let p = PlasticityParams {
            iso_hardening: 5e9,
            ..PlasticityParams::steel()
        };
        let mu = p.young / (2.0 * (1.0 + p.poisson));
        let gamma_y = (p.yield_stress / 3.0_f64.sqrt()) / mu;

        let mut state = PlasticState::new();
        let g_max = 4.0 * gamma_y;
        for k in 1..=100 {
            let g = g_max * k as f64 / 100.0;
            let (_, next) = j2_return_mapping(&p, &[0.0, 0.0, 0.0, g, 0.0, 0.0], &state);
            state = next;
        }
        let eps_eq_peak = state.eps_eq;
        assert!(eps_eq_peak > 0.0);

        // Unload by a small decrement: elastic, slope G.
        let (sig_a, _) = j2_return_mapping(&p, &[0.0, 0.0, 0.0, g_max, 0.0, 0.0], &state);
        let dg = -0.1 * gamma_y;
        let (sig_b, st_b) = j2_return_mapping(&p, &[0.0, 0.0, 0.0, g_max + dg, 0.0, 0.0], &state);
        assert!(
            (st_b.eps_eq - eps_eq_peak).abs() < 1e-18,
            "unloading must not accumulate plastic strain"
        );
        let slope = (sig_b[3] - sig_a[3]) / dg;
        assert!(
            (slope - mu).abs() / mu < 1e-9,
            "unloading shear slope should be G = {mu:.6e}, got {slope:.6e}"
        );
    }

    #[test]
    fn plastic_1d_matches_closed_form() {
        // Perfect plasticity: the axial stress plateaus at σ_y⁰.
        let p = PlasticityParams::steel();
        let ey = p.yield_stress / p.young;
        for &mult in &[1.5_f64, 2.0, 5.0] {
            let (s, eps_eq, plastic) = plastic_1d(&p, ey * mult, 0.0);
            assert!(plastic);
            assert!(
                (s - p.yield_stress).abs() / p.yield_stress < 1e-12,
                "perfect plasticity must plateau at σ_y, got {s:.6e}"
            );
            // ε = σ/E + ε̄ᵖ.
            assert!((ey * mult - (s / p.young + eps_eq)).abs() < 1e-15);
        }

        // Linear isotropic hardening: σ = σ_y⁰ + H ε̄ᵖ with ε = σ/E + ε̄ᵖ, so
        // ε̄ᵖ = (E ε − σ_y⁰)/(E + H).
        let h = 4e9;
        let ph = PlasticityParams {
            iso_hardening: h,
            ..PlasticityParams::steel()
        };
        for &mult in &[1.5_f64, 3.0, 10.0] {
            let eps = ey * mult;
            let (s, eps_eq, _) = plastic_1d(&ph, eps, 0.0);
            let want_epsp = (ph.young * eps - ph.yield_stress) / (ph.young + h);
            let want_sig = ph.yield_stress + h * want_epsp;
            assert!(
                (eps_eq - want_epsp).abs() / want_epsp < 1e-12,
                "ε̄ᵖ mismatch"
            );
            assert!((s - want_sig).abs() / want_sig < 1e-12, "σ mismatch");
        }
    }

    #[test]
    fn uniaxial_tangent_matches_hardening_slope() {
        let h = 4e9;
        let p = PlasticityParams {
            iso_hardening: h,
            ..PlasticityParams::steel()
        };
        let ey = p.yield_stress / p.young;
        // Elastic branch.
        assert_eq!(uniaxial_tangent_1d(&p, 0.5 * ey, 0.0), p.young);
        // Plastic branch: the combined-hardening slope E·H/(E + H).
        let want = p.young * h / (p.young + h);
        let got = uniaxial_tangent_1d(&p, 2.0 * ey, 0.0);
        assert!((got - want).abs() / want < 1e-12, "{got} vs {want}");
        // Cross-check against a finite difference of the stress law itself.
        let d = 1e-7;
        let sp = plastic_1d(&p, 2.0 * ey + d, 0.0).0;
        let sm = plastic_1d(&p, 2.0 * ey - d, 0.0).0;
        let fd = (sp - sm) / (2.0 * d);
        assert!((fd - want).abs() / want < 1e-6, "fd {fd} vs {want}");
    }

    #[test]
    fn uniaxial_stress_consistency_3d_vs_1d() {
        // The 3-D law under *genuine* uniaxial stress must reproduce the scalar
        // law, for both perfect plasticity and linear isotropic hardening.
        for p in [
            PlasticityParams::steel(),
            PlasticityParams {
                iso_hardening: 4e9,
                ..PlasticityParams::steel()
            },
        ] {
            let ey = p.yield_stress / p.young;
            for &mult in &[0.5_f64, 1.5, 2.0, 4.0] {
                let eps = ey * mult;
                let s3d = uniaxial_stress_axial(&p, eps);
                let s1d = plastic_1d(&p, eps, 0.0).0;
                assert!(
                    (s3d - s1d).abs() / s1d.abs() < 1e-6,
                    "3-D uniaxial-stress {s3d:.6e} vs 1-D {s1d:.6e} at ε={eps:.3e}"
                );
            }
        }
    }

    #[test]
    fn rod_reaction_matches_single_point() {
        // A 1-element rod: applied force P should give axial stress σ = P/A and,
        // for P below yield, strain ε = P/(EA).
        let mesh = rod(1, 1.0);
        let a = 1.0;
        let p = PlasticityParams::steel();
        let p_load = 200e6 * a * 0.5; // half-yield force -> elastic
        let history =
            solve_elastic_plastic_rod(&mesh, a, &p, &[p_load]).expect("rod solve must converge");
        let u = &history[0];
        let eps = (u[1] - u[0]) / 1.0;
        let want_eps = p_load / (p.young * a);
        assert!(
            (eps - want_eps).abs() / want_eps < 1e-6,
            "eps {} vs {}",
            eps,
            want_eps
        );
    }

    #[test]
    fn rod_past_yield_matches_hardening_closed_form() {
        // Force control past yield needs a hardening material (a perfectly
        // plastic rod has zero tangent stiffness there).
        let h = 4e9;
        let p = PlasticityParams {
            iso_hardening: h,
            ..PlasticityParams::steel()
        };
        let area = 2e-4;
        let mesh = rod(4, 1.0);
        let p_y = p.yield_stress * area;
        let load = 1.3 * p_y;
        let history =
            solve_elastic_plastic_rod(&mesh, area, &p, &[load]).expect("rod solve must converge");
        let u_end = *history[0].last().unwrap();

        // Uniform bar: σ = P/A, ε̄ᵖ = (σ − σ_y⁰)/H, ε = σ/E + ε̄ᵖ, u = ε·L.
        let sigma = load / area;
        let epsp = (sigma - p.yield_stress) / h;
        let want = sigma / p.young + epsp;
        assert!(
            (u_end - want).abs() / want < 1e-8,
            "u_end {u_end:.6e} vs closed form {want:.6e}"
        );
    }
}
