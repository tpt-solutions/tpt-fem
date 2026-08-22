//! Topology optimization for 2-D linear elasticity.
//!
//! This crate provides a self-contained minimum-compliance (stiffness
//! maximization) optimizer using the **Solid Isotropic Material with
//! Penalization (SIMP)** model, a density (sensitivity) filter, and the
//! **optimality-criteria (OC)** update. It sits on top of the existing
//! `tpt-fem` stack:
//!
//! * [`tpt_fem_element::Quad4`] — bilinear quadrilateral shape functions and
//!   reference-coordinate derivatives.
//! * [`tpt_fem_quadrature`] — Gauss-Legendre quadrature.
//! * [`tpt_fem_sparse::Coo`] — triplet accumulation.
//! * [`tpt_fem_assembly::solve_with_dirichlet`] — essential-boundary-condition
//!   solve for each design iteration.
//!
//! The optimizer minimizes the structural compliance `c = uᵀ K(ρ) u` subject to
//! a volume fraction `vol_frac`, where the Young's modulus of each element is
//! interpolated as `E(ρ) = E_min + ρᵖ (E₀ − E_min)` (`ρ ∈ [0, 1]`, `p ≥ 1` the
//! SIMP penalty). This is the classic Bendsøe/Sigmund formulation; see
//! `crates/tpt-fem-topopt/README.md` for the mathematical setup.
//!
//! ```
//! use tpt_fem_topopt::{cantilever_load, topopt_simp, TopOptParams, Grid};
//!
//! // Small cantilever: 20×10 grid of unit squares, 50% volume fraction.
//! let grid = Grid::new(20, 10, 1.0);
//! let (f, bcs) = cantilever_load(&grid, 1.0);
//! let params = TopOptParams {
//!     grid: grid.clone(),
//!     e0: 1.0,
//!     nu: 0.3,
//!     vol_frac: 0.5,
//!     penal: 3.0,
//!     filter_radius: 2.0,
//!     max_iter: 40,
//!     move_limit: 0.2,
//! };
//! let res = topopt_simp(&params, &f, &bcs).unwrap();
//! // The optimized design is lighter (compliance is lower) than the uniform
//! // starting point, while still consuming exactly `vol_frac` of the domain.
//! assert!(res.compliance.last().unwrap() < &res.compliance[0]);
//! let used: f64 = res.densities.iter().sum();
//! assert!((used - 0.5 * (grid.n_elem() as f64)).abs() < 1e-6);
//! ```

use tpt_fem_assembly::solve_with_dirichlet;
use tpt_fem_element::{Quad4, ReferenceElement};
use tpt_fem_quadrature::gauss_legendre;
use tpt_fem_sparse::{Coo, SparseError};

/// Minimum (void) Young's modulus, expressed as a fraction of the solid
/// modulus `E₀`. A small non-zero void stiffness keeps the global system
/// non-singular when an element is driven to `ρ ≈ 0`.
const E_MIN_FACTOR: f64 = 1e-3;
/// Lower bound on a design density. The SIMP model is only physically
/// meaningful for `ρ > 0`; `RHO_MIN` prevents a fully-void (and singular)
/// element while still allowing near-void material.
const RHO_MIN: f64 = 1e-3;
/// Upper bound on a design density (a solid element).
const RHO_MAX: f64 = 1.0;
/// Damping exponent in the optimality-criteria move (`η = 0.5`).
const OC_DAMP: f64 = 0.5;

/// A structured grid of `nx·ny` nodes and `(nx−1)·(ny−1)` axis-aligned
/// `Quad4` elements of side `h`.
///
/// Node `n = j·nx + i` sits at `(i·h, j·h)`. Element `(ei, ej)` (with
/// `ei ∈ [0, nx−1)`, `ej ∈ [0, ny−1)`) connects four nodes in the `Quad4`
/// reference order `[-1,-1] → [1,-1] → [1,1] → [-1,1]` (bottom-left,
/// bottom-right, top-right, top-left).
#[derive(Clone)]
pub struct Grid {
    /// Number of nodes along the `x` axis.
    pub nx: usize,
    /// Number of nodes along the `y` axis.
    pub ny: usize,
    /// Element side length.
    pub h: f64,
    /// Node coordinates, `coords[j·nx + i] = [i·h, j·h]`.
    pub coords: Vec<[f64; 2]>,
    /// Element → four node indices (Quad4 reference order).
    pub elems: Vec<[usize; 4]>,
    /// Element centroids in physical coordinates.
    pub centers: Vec<[f64; 2]>,
}

impl Grid {
    /// Build an `nx·ny` node grid of `h`-sized `Quad4` elements.
    pub fn new(nx: usize, ny: usize, h: f64) -> Self {
        let mut coords = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                coords.push([i as f64 * h, j as f64 * h]);
            }
        }
        let mut elems = Vec::with_capacity((nx - 1) * (ny - 1));
        let mut centers = Vec::with_capacity((nx - 1) * (ny - 1));
        for ej in 0..ny - 1 {
            for ei in 0..nx - 1 {
                let n0 = ej * nx + ei;
                let n1 = ej * nx + (ei + 1);
                let n2 = (ej + 1) * nx + (ei + 1);
                let n3 = (ej + 1) * nx + ei;
                elems.push([n0, n1, n2, n3]);
                centers.push([(ei as f64 + 0.5) * h, (ej as f64 + 0.5) * h]);
            }
        }
        Grid {
            nx,
            ny,
            h,
            coords,
            elems,
            centers,
        }
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.coords.len()
    }

    /// Number of elements.
    pub fn n_elem(&self) -> usize {
        self.elems.len()
    }
}

/// Optimizer configuration for [`topopt_simp`].
pub struct TopOptParams {
    /// Structured mesh the design lives on.
    pub grid: Grid,
    /// Solid (void-free) Young's modulus `E₀`.
    pub e0: f64,
    /// Poisson ratio (plane-stress constitutive model).
    pub nu: f64,
    /// Target volume fraction `vol_frac ∈ (0, 1]` of solid material.
    pub vol_frac: f64,
    /// SIMP penalty exponent `p` (≥ 1; typically 3).
    pub penal: f64,
    /// Sensitivity-filter radius in physical units.
    pub filter_radius: f64,
    /// Maximum number of OC iterations.
    pub max_iter: usize,
    /// Per-iteration density move limit (e.g. `0.2`).
    pub move_limit: f64,
}

/// Outcome of a [`topopt_simp`] run.
pub struct TopOptResult {
    /// Final design density of every element (`Vec` length = `grid.n_elem()`).
    pub densities: Vec<f64>,
    /// Structural compliance `c = uᵀ K u` at each iteration (index 0 is the
    /// uniform starting point). Strictly non-increasing under OC.
    pub compliance: Vec<f64>,
    /// Number of OC iterations actually performed.
    pub iterations: usize,
}

/// Assemble the `8×8` `Quad4` plane-stress element stiffness matrix for unit
/// Young's modulus (`E = 1`) and Poisson ratio `nu` on a square of side `h`.
///
/// Because the element is an axis-aligned square, the Jacobian to physical
/// coordinates is constant (`J = (h/2)·I`), so a single `2×2` Gauss rule
/// integrates the stiffness exactly.
fn element_stiffness_unit(nu: f64, h: f64) -> Vec<Vec<f64>> {
    // Plane-stress constitutive matrix for E = 1.
    let c = 1.0 / (1.0 - nu * nu);
    let d = [
        [c, c * nu, 0.0],
        [c * nu, c, 0.0],
        [0.0, 0.0, c * (1.0 - nu) / 2.0],
    ];
    let g = gauss_legendre(2);
    let jdet = h * h / 4.0;
    let scale = 2.0 / h; // d(·)/dx = (2/h)·d(·)/dξ for an h-sided square
    let mut ke = vec![vec![0.0; 8]; 8];

    for a in 0..g.points.len() {
        for b in 0..g.points.len() {
            let xi = [g.points[a], g.points[b]];
            let w = g.weights[a] * g.weights[b];
            let dndxi = Quad4::grad(&xi); // [node][dξ, dη]
                                          // Physical-gradient B matrix (3×8).
            let mut bm = vec![vec![0.0; 8]; 3];
            for n in 0..4 {
                let nx_ = scale * dndxi[n][0];
                let ny_ = scale * dndxi[n][1];
                bm[0][2 * n] = nx_;
                bm[1][2 * n + 1] = ny_;
                bm[2][2 * n] = ny_;
                bm[2][2 * n + 1] = nx_;
            }
            // ke += w·|J|·Bᵀ D B
            let mut db = vec![vec![0.0; 8]; 3];
            for i in 0..3 {
                for j in 0..8 {
                    let mut s = 0.0;
                    for k in 0..3 {
                        s += d[i][k] * bm[k][j];
                    }
                    db[i][j] = s;
                }
            }
            for i in 0..8 {
                for j in 0..8 {
                    let mut s = 0.0;
                    for k in 0..3 {
                        s += bm[k][i] * db[k][j];
                    }
                    ke[i][j] += w * jdet * s;
                }
            }
        }
    }
    ke
}

/// Assemble the global stiffness `Coo` for the current design densities.
fn assemble(grid: &Grid, densities: &[f64], ke0: &[Vec<f64>], e0: f64, penal: f64) -> Coo {
    let emin = E_MIN_FACTOR * e0;
    let mut coo = Coo::new();
    for (e, &nodes) in grid.elems.iter().enumerate() {
        let e_rho = emin + densities[e].powf(penal) * (e0 - emin);
        // Global DOF ordering: [2·n0, 2·n0+1, 2·n1, 2·n1+1, …].
        let mut gdof = [0usize; 8];
        for (k, &n) in nodes.iter().enumerate() {
            gdof[2 * k] = 2 * n;
            gdof[2 * k + 1] = 2 * n + 1;
        }
        for i in 0..8 {
            for j in 0..8 {
                coo.push(gdof[i], gdof[j], e_rho * ke0[i][j]);
            }
        }
    }
    coo
}

/// `uₑᵀ K₀ uₑ` for element `e` (used in the compliance sensitivity).
fn element_quad_form(ke0: &[Vec<f64>], u: &[f64], nodes: &[usize; 4]) -> f64 {
    let mut ue = [0.0f64; 8];
    for (k, &n) in nodes.iter().enumerate() {
        ue[2 * k] = u[2 * n];
        ue[2 * k + 1] = u[2 * n + 1];
    }
    let mut s = 0.0;
    for i in 0..8 {
        for j in 0..8 {
            s += ue[i] * ke0[i][j] * ue[j];
        }
    }
    s
}

/// Filter a raw per-element sensitivity through the density-weighted kernel
/// `w_ef = max(0, R − |c_e − c_f|)`, normalising by the weight sum. This is the
/// standard (Sigmund) sensitivity filter that removes checkerboard artefacts.
fn filter_sensitivity(grid: &Grid, raw: &[f64], radius: f64) -> Vec<f64> {
    let r2 = radius * radius;
    let mut out = vec![0.0; raw.len()];
    for e in 0..raw.len() {
        let ce = grid.centers[e];
        let mut num = 0.0;
        let mut den = 0.0;
        for f in 0..raw.len() {
            let dx = ce[0] - grid.centers[f][0];
            let dy = ce[1] - grid.centers[f][1];
            let d2 = dx * dx + dy * dy;
            if d2 <= r2 {
                let w = radius - d2.sqrt();
                num += w * raw[f];
                den += w;
            }
        }
        out[e] = if den > 0.0 { num / den } else { raw[e] };
    }
    out
}

/// One optimality-criteria update: choose `λ` by bisection so the resulting
/// densities satisfy the volume constraint, then move each density toward the
/// OC point `ρ·(B_e / λ)^{η}` within the per-step move limit.
fn oc_update(densities: &[f64], sens: &[f64], vol_frac: f64, move_limit: f64) -> Vec<f64> {
    let n = densities.len();
    let target = vol_frac * (n as f64);
    // B_e = −∂c/∂ρ_e (positive for a compliance-minimizing move).
    let b: Vec<f64> = sens.iter().map(|&s| -s).collect();

    // Bisection on λ. Larger λ → smaller ρ_new (monotonic in λ).
    let mut lo = 1e-9;
    let mut hi = 1e9;
    let mut lambda = 1.0;
    for _ in 0..60 {
        lambda = 0.5 * (lo + hi);
        let mut sum = 0.0;
        for e in 0..n {
            let ratio = if b[e] > 0.0 {
                (b[e] / lambda).powf(OC_DAMP)
            } else {
                0.0
            };
            let cand = densities[e] * ratio;
            let clamped = cand.clamp(
                (densities[e] - move_limit).max(RHO_MIN),
                (densities[e] + move_limit).min(RHO_MAX),
            );
            sum += clamped;
        }
        if sum > target {
            lo = lambda; // need to shrink densities → larger λ
        } else {
            hi = lambda;
        }
    }

    let mut out = vec![0.0; n];
    for e in 0..n {
        let ratio = if b[e] > 0.0 {
            (b[e] / lambda).powf(OC_DAMP)
        } else {
            0.0
        };
        let cand = densities[e] * ratio;
        out[e] = cand.clamp(
            (densities[e] - move_limit).max(RHO_MIN),
            (densities[e] + move_limit).min(RHO_MAX),
        );
    }
    out
}

/// Run the SIMP minimum-compliance optimization.
///
/// `f` is the global load vector (one entry per nodal DOF, `2·n_nodes` long)
/// and `bcs` the Dirichlet conditions as `(global_dof, value)` pairs. Returns
/// the optimized densities, the compliance history, and the iteration count.
///
/// The volume constraint `Σ ρ_e = vol_frac·n_elem` is enforced to bisection
/// tolerance on every iteration, so `densities.iter().sum()` equals the target
/// (minus a sub-`1e-6` residual).
pub fn topopt_simp(
    params: &TopOptParams,
    f: &[f64],
    bcs: &[(usize, f64)],
) -> Result<TopOptResult, SparseError> {
    let grid = &params.grid;
    let ne = grid.n_elem();
    let ke0 = element_stiffness_unit(params.nu, grid.h);
    let emin = E_MIN_FACTOR * params.e0;

    let mut densities = vec![params.vol_frac; ne];
    let mut compliance = Vec::with_capacity(params.max_iter + 1);

    // Initial (uniform) compliance.
    {
        let coo = assemble(grid, &densities, &ke0, params.e0, params.penal);
        let u = solve_with_dirichlet(&coo, f, bcs)?;
        compliance.push(dot(f, &u));
    }

    let mut iter = 0;
    for _ in 0..params.max_iter {
        iter += 1;
        let coo = assemble(grid, &densities, &ke0, params.e0, params.penal);
        let u = solve_with_dirichlet(&coo, f, bcs)?;
        compliance.push(dot(f, &u));

        // Raw sensitivity: ∂c/∂ρ_e = −p·ρ^{p−1}·(E₀−E_min)·uₑᵀ K₀ uₑ.
        let mut raw = vec![0.0; ne];
        for e in 0..ne {
            let ue_k0 = element_quad_form(&ke0, &u, &grid.elems[e]);
            raw[e] =
                -params.penal * densities[e].powf(params.penal - 1.0) * (params.e0 - emin) * ue_k0;
        }
        let sens = filter_sensitivity(grid, &raw, params.filter_radius);
        densities = oc_update(&densities, &sens, params.vol_frac, params.move_limit);
    }

    Ok(TopOptResult {
        densities,
        compliance,
        iterations: iter,
    })
}

/// `Σ_i a_i·b_i`.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Build the load vector and Dirichlet conditions for a classic cantilever:
/// the entire left edge (`x = 0`) is clamped, and a downward point load is
/// applied at the bottom-right node.
///
/// Returns `(f, bcs)` where `f` has length `2·grid.n_nodes()` and `bcs` fixes
/// every DOF on the left edge to zero.
pub fn cantilever_load(grid: &Grid, load: f64) -> (Vec<f64>, Vec<(usize, f64)>) {
    let ndof = 2 * grid.n_nodes();
    let mut f = vec![0.0; ndof];
    let mut bcs = Vec::new();
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            if i == 0 {
                let n = j * grid.nx + i;
                bcs.push((2 * n, 0.0));
                bcs.push((2 * n + 1, 0.0));
            }
        }
    }
    // Bottom-right node (i = nx−1, j = 0) gets a downward (−y) load.
    let lr = grid.nx - 1;
    f[2 * lr + 1] = -load;
    (f, bcs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cantilever_lowers_compliance_and_keeps_volume() {
        let grid = Grid::new(20, 10, 1.0);
        let (f, bcs) = cantilever_load(&grid, 1.0);
        let params = TopOptParams {
            grid: grid.clone(),
            e0: 1.0,
            nu: 0.3,
            vol_frac: 0.5,
            penal: 3.0,
            filter_radius: 2.0,
            max_iter: 40,
            move_limit: 0.2,
        };
        let res = topopt_simp(&params, &f, &bcs).unwrap();
        assert_eq!(res.densities.len(), grid.n_elem());
        // Optimization strictly reduces compliance from the uniform start.
        assert!(
            res.compliance.last().unwrap() < res.compliance.first().unwrap(),
            "compliance did not decrease: {} -> {}",
            res.compliance[0],
            res.compliance.last().unwrap()
        );
        // Volume constraint honoured to bisection tolerance.
        let used: f64 = res.densities.iter().sum();
        assert!((used - 0.5 * (grid.n_elem() as f64)).abs() < 1e-6);
        // Compliance is monotonically non-increasing under OC.
        for w in res.compliance.windows(2) {
            assert!(w[1] <= w[0] + 1e-9, "compliance increased at an iteration");
        }
    }

    #[test]
    fn full_volume_stays_solid() {
        // With vol_frac = 1 the optimizer has no material to remove; densities
        // remain solid and the solve must succeed.
        let grid = Grid::new(16, 8, 1.0);
        let (f, bcs) = cantilever_load(&grid, 1.0);
        let params = TopOptParams {
            grid: grid.clone(),
            e0: 1.0,
            nu: 0.3,
            vol_frac: 1.0,
            penal: 3.0,
            filter_radius: 2.0,
            max_iter: 10,
            move_limit: 0.2,
        };
        let res = topopt_simp(&params, &f, &bcs).unwrap();
        let used: f64 = res.densities.iter().sum();
        assert!((used - grid.n_elem() as f64).abs() < 1e-6);
        for &rho in &res.densities {
            assert!((rho - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn denser_initial_compliance_exceeds_optimized() {
        // Sanity: a 40%-volume optimized design is stiffer (lower compliance)
        // than a 40%-volume uniform plate.
        let grid = Grid::new(24, 12, 1.0);
        let (f, bcs) = cantilever_load(&grid, 1.0);
        let params = TopOptParams {
            grid: grid.clone(),
            e0: 1.0,
            nu: 0.3,
            vol_frac: 0.4,
            penal: 3.0,
            filter_radius: 2.0,
            max_iter: 30,
            move_limit: 0.2,
        };
        let res = topopt_simp(&params, &f, &bcs).unwrap();
        assert!(res.compliance.last().unwrap() < &res.compliance[0]);
    }
}
