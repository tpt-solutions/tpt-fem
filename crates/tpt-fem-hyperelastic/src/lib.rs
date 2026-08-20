//! Hyperelastic constitutive models for `tpt-fem`.
//!
//! Provides the first Piola–Kirchhoff stress `P = ∂Ψ/∂F` for the
//! incompressible [Neo-Hookean](https://en.wikipedia.org/wiki/Neo-Hookean_solid),
//! [Mooney–Rivlin](https://en.wikipedia.org/wiki/Mooney%E2%80%93Rivlin_solid),
//! and [Ogden](https://en.wikipedia.org/wiki/Ogden_(material_model)) strain-energy
//! functions, plus a 1-D bar solver wired into [`tpt-fem-solve`]'s Newton loop.
//!
//! ```
//! use tpt_fem_hyperelastic::neo_hookean_piola;
//!
//! // Uniaxial stretch λ with traction-free sides (p = μ λ⁻¹ for incompressible
//! // neo-Hookean): the axial PK1 stress is μ(λ - λ⁻²).
//! let lam = 1.5;
//! let mut f: [[f64; 3]; 3] = [[0.0; 3]; 3];
//! f[0][0] = lam;
//! f[1][1] = lam.powf(-0.5);
//! f[2][2] = lam.powf(-0.5);
//! let p = 100.0 * lam.powf(-1.0);
//! let pk = neo_hookean_piola(&f, 100.0, p).expect("F is invertible");
//! let want = 100.0 * (lam - lam.powf(-2.0));
//! assert!((pk[0][0] - want).abs() / want < 1e-12);
//! ```

/// A `3×3` matrix (row-major).
pub type Mat3 = [[f64; 3]; 3];

/// Matrix product `a·b`.
pub fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}

/// Determinant of a `3×3` matrix.
pub fn mat_det(m: &Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Inverse of a `3×3` matrix. Returns `None` if the matrix is (numerically)
/// singular rather than panicking, so callers that may supply a degenerate
/// deformation gradient (collapsed element) can handle it as an error.
pub fn mat_inv(m: &Mat3) -> Option<Mat3> {
    let det = mat_det(m);
    if det.abs() <= 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    let mut inv = [[0.0; 3]; 3];
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;
    Some(inv)
}

/// Transpose of a `3×3` matrix.
pub fn mat_transpose(m: &Mat3) -> Mat3 {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

/// First Piola–Kirchhoff stress for the incompressible **Neo-Hookean** model,
/// `P = μ F − p F^{-T}`. Returns `None` if `F` is singular.
pub fn neo_hookean_piola(f: &Mat3, mu: f64, pressure: f64) -> Option<Mat3> {
    let finv = mat_inv(f)?;
    let finvt = mat_transpose(&finv);
    let mut p = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] = mu * f[i][j] - pressure * finvt[i][j];
        }
    }
    Some(p)
}

/// First Piola–Kirchhoff stress for the incompressible **Mooney–Rivlin** model,
/// `P = 2 c₁ F + 2 c₂ (I₁ F − C F) − p F^{-T}` where `C = FᵀF` and
/// `I₁ = tr(C)`. Returns `None` if `F` is singular.
pub fn mooney_rivlin_piola(f: &Mat3, c1: f64, c2: f64, pressure: f64) -> Option<Mat3> {
    let ft = mat_transpose(f);
    let c = mat_mul(&ft, f);
    let i1 = c[0][0] + c[1][1] + c[2][2];
    let cf = mat_mul(&c, f);
    let finv = mat_inv(f)?;
    let finvt = mat_transpose(&finv);
    let mut p = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] =
                2.0 * c1 * f[i][j] + 2.0 * c2 * (i1 * f[i][j] - cf[i][j]) - pressure * finvt[i][j];
        }
    }
    Some(p)
}

/// A single Ogden term `(μ_i, α_i)`.
#[derive(Clone, Copy, Debug)]
pub struct OgdenTerm {
    /// Shear-like modulus `μ_i`.
    pub mu: f64,
    /// Exponent `α_i`.
    pub alpha: f64,
}

/// First Piola–Kirchhoff stress for the incompressible **Ogden** model,
/// `P = Σ_i μ_i α_i λ_i^{α_i−1} n_i n_iᵀ − p F^{-T}`, evaluated from the
/// principal stretches `λ` (length 3) and the corresponding principal directions
/// (columns of the rotation `n`).
pub fn ogden_piola(
    lambdas: &[f64; 3],
    directions: &Mat3,
    terms: &[OgdenTerm],
    pressure: f64,
) -> Mat3 {
    let finv = {
        // F = N diag(λ) Nᵀ ; F^{-T} = N diag(λ⁻¹) Nᵀ.
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += directions[i][k] * (1.0 / lambdas[k]) * directions[j][k];
                }
                m[i][j] = s;
            }
        }
        m
    };
    let mut p = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                let pk = terms
                    .iter()
                    .map(|t| t.mu * t.alpha * lambdas[k].powf(t.alpha - 1.0))
                    .sum::<f64>();
                s += pk * directions[i][k] * directions[j][k];
            }
            p[i][j] = s - pressure * finv[i][j];
        }
    }
    p
}

/// 1-D incompressible neo-Hookean nominal (PK1) stress for a stretch `F = λ`:
/// `P = μ (F − F⁻²)`.
pub fn neo_hookean_1d(f: f64, mu: f64) -> f64 {
    mu * (f - f.powi(-2))
}

/// Errors returned by the hyperelastic solvers.
#[derive(Debug)]
pub enum HyperelasticError {
    /// The Newton iteration on the (large-strain) residual did not converge.
    Newton(tpt_fem_solve::NewtonError),
}

impl std::fmt::Display for HyperelasticError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HyperelasticError::Newton(e) => {
                write!(f, "hyperelastic bar Newton iteration failed: {e}")
            }
        }
    }
}

impl std::error::Error for HyperelasticError {}

impl From<tpt_fem_solve::NewtonError> for HyperelasticError {
    fn from(e: tpt_fem_solve::NewtonError) -> Self {
        HyperelasticError::Newton(e)
    }
}

/// Solve a 1-D bar (built from `Line2` elements, one axial DOF per node) under a
/// target end stretch using [`tpt-fem-solve`]'s Newton loop with the
/// incompressible neo-Hookean response. Returns the displacement field whose
/// rightmost node has been displaced to achieve the prescribed stretch, or a
/// [`HyperelasticError`] if the Newton iteration fails to converge.
pub fn solve_hyperelastic_bar(
    mesh: &tpt_fem_mesh::Mesh,
    area: f64,
    mu: f64,
    target_stretch: f64,
) -> Result<Vec<f64>, HyperelasticError> {
    let nnodes = mesh.node_count();
    let mut elem_len = Vec::with_capacity(mesh.elements.len());
    for elem in &mesh.elements {
        let a = mesh.node_coords(elem.nodes[0])[0];
        let b = mesh.node_coords(elem.nodes[1])[0];
        elem_len.push((b - a).abs());
    }
    let total_len: f64 = elem_len.iter().sum();
    let target_disp = total_len * (target_stretch - 1.0);

    let dirichlet = [(0usize, 0.0), (nnodes - 1, target_disp)];
    tpt_fem_solve::newton(
        &vec![0.0; nnodes],
        |u| {
            let mut fint = vec![0.0; nnodes];
            for (e, elem) in mesh.elements.iter().enumerate() {
                let n0 = elem.nodes[0];
                let n1 = elem.nodes[1];
                let l0 = elem_len[e];
                let l = l0 + (u[n1] - u[n0]);
                let f = l / l0; // axial stretch
                let p = neo_hookean_1d(f, mu) * area;
                fint[n0] -= p;
                fint[n1] += p;
            }
            // Prescribed (fixed) DOFs carry a nonzero reaction; zero their
            // residual so the Newton convergence norm measures only free DOFs.
            fint[0] = 0.0;
            fint[nnodes - 1] = 0.0;
            fint
        },
        |u| {
            let mut coo = tpt_fem_sparse::Coo::new();
            for (e, elem) in mesh.elements.iter().enumerate() {
                let n0 = elem.nodes[0];
                let n1 = elem.nodes[1];
                let l0 = elem_len[e];
                let l = l0 + (u[n1] - u[n0]);
                let f = l / l0;
                // dP/du = dP/df * df/du ; dP/df = μ(1 + 2 f⁻³); df/du_{n1}=1/l0
                let dpdf = mu * (1.0 + 2.0 * f.powi(-3));
                let ke = dpdf * area / l0;
                coo.push(n0, n0, ke);
                coo.push(n0, n1, -ke);
                coo.push(n1, n0, -ke);
                coo.push(n1, n1, ke);
            }
            coo
        },
        &dirichlet,
        &tpt_fem_solve::NewtonOptions::default(),
    )
    .map_err(HyperelasticError::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    fn identity() -> Mat3 {
        let mut f = [[0.0; 3]; 3];
        for i in 0..3 {
            f[i][i] = 1.0;
        }
        f
    }

    #[test]
    fn neo_hookean_uniaxial() {
        for &lam in &[1.1, 1.5, 2.0, 0.8] {
            let mut f = identity();
            f[0][0] = lam;
            f[1][1] = lam.powf(-0.5);
            f[2][2] = lam.powf(-0.5);
            let mu = 100.0;
            let p = mu * lam.powf(-1.0); // traction-free sides
            let pk = neo_hookean_piola(&f, mu, p).expect("F is invertible");
            let want = mu * (lam - lam.powf(-2.0));
            assert!(
                (pk[0][0] - want).abs() / want < 1e-12,
                "lam={lam} got {} want {}",
                pk[0][0],
                want
            );
            // Side tractions vanish.
            assert!(pk[1][1].abs() < 1e-9);
            assert!(pk[2][2].abs() < 1e-9);
        }
    }

    #[test]
    fn mooney_rivlin_uniaxial() {
        // Principal nominal stress: S1 = 2(λ - λ⁻²)(c1 + c2 λ⁻¹).
        for &lam in &[1.2, 1.6, 2.0] {
            let mut f = identity();
            f[0][0] = lam;
            f[1][1] = lam.powf(-0.5);
            f[2][2] = lam.powf(-0.5);
            let c1 = 80.0;
            let c2 = 20.0;
            let _mu = 2.0 * (c1 + c2);
            // Traction-free sides give p = 2 c1 λ⁻¹ + 2 c2 (λ + λ⁻²).
            let p = 2.0 * c1 * lam.powf(-1.0) + 2.0 * c2 * (lam + lam.powf(-2.0));
            let pk = mooney_rivlin_piola(&f, c1, c2, p).expect("F is invertible");
            let want = 2.0 * (lam - lam.powf(-2.0)) * (c1 + c2 * lam.powf(-1.0));
            assert!(
                (pk[0][0] - want).abs() / want < 1e-12,
                "lam={lam} got {} want {}",
                pk[0][0],
                want
            );
        }
    }

    #[test]
    fn ogden_uniaxial() {
        // S1 = Σ μ_i (λ^{α_i-1} - λ^{(α_i/2)-1}) for traction-free incompressible.
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
        for lam in [1.2_f64, 1.7_f64, 2.2_f64] {
            let n = identity(); // principal directions = axes
            let lm = lam.powf(-0.5);
            let lambdas = [lam, lm, lm];
            // Traction-free sides: p = Σ μ_i α_i λ₂^{α_i} with λ₂ = λ^{-1/2}.
            let p: f64 = terms
                .iter()
                .map(|t| t.mu * t.alpha * lambdas[1].powf(t.alpha))
                .sum();
            let pk = ogden_piola(&lambdas, &n, &terms, p);
            let want: f64 = terms
                .iter()
                .map(|t| {
                    t.mu * t.alpha * (lam.powf(t.alpha - 1.0) - lam.powf(-(t.alpha / 2.0 + 1.0)))
                })
                .sum();
            assert!(
                (pk[0][0] - want).abs() / want < 1e-12,
                "lam={lam} got {} want {}",
                pk[0][0],
                want
            );
        }
    }

    #[test]
    fn bar_solve_matches_stretch() {
        let mut b = MeshBuilder::new();
        let mut prev = b.add_node(vec![0.0]);
        for i in 1..=4 {
            let node = b.add_node(vec![1.0 * i as f64 / 4.0]);
            b.add_element(CellType::Line, vec![prev, node]);
            prev = node;
        }
        let mesh = b.build();
        let u = solve_hyperelastic_bar(&mesh, 1.0, 100.0, 1.5)
            .expect("hyperelastic bar solve must converge");
        // Rightmost node displacement = total length * 0.5 = 0.5.
        assert!((u[4] - 0.5).abs() < 1e-6);
        // Linear stretch -> uniform displacement field.
        assert!((u[2] - 0.25).abs() < 1e-6);
    }
}
