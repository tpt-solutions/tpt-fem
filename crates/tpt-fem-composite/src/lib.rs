//! Classical lamination theory and cohesive-zone delamination models for
//! `tpt-fem`.
//!
//! [`laminate_abd`] builds the `6×6` extensional–bending–coupling (`A`/`B`/`D`)
//! stiffness matrix of a stacked plate from individual [`Ply`] properties and
//! orientations. A symmetric stacking sequence yields `B = 0` (no
//! extension–bending coupling), which the unit tests verify.
//!
//! [`CohesiveLaw`] provides a bilinear traction–separation model for interface
//! (cohesive-zone) elements, kept independent of `tpt-fem-contact`'s
//! surface-to-surface algorithm.
//!
//! ```
//! use tpt_fem_composite::{laminate_abd, Ply};
//!
//! // A single 1 mm isotropic ply at 0°: A11 should equal Q11 * t.
//! let ply = Ply { e1: 70e9, e2: 70e9, nu12: 0.3, g12: 70e9 / 2.6, thickness: 1e-3, angle_deg: 0.0 };
//! let abd = laminate_abd(&[ply]);
//! let nu = 0.3;
//! let q11 = 70e9 / (1.0 - nu * nu);
//! assert!((abd[0][0] - q11 * 1e-3).abs() / (q11 * 1e-3) < 1e-9);
//! ```

/// A single orthotropic ply in the laminate stack.
#[derive(Clone, Copy, Debug)]
pub struct Ply {
    /// Longitudinal Young's modulus `E₁`.
    pub e1: f64,
    /// Transverse Young's modulus `E₂`.
    pub e2: f64,
    /// Major Poisson's ratio `ν₁₂`.
    pub nu12: f64,
    /// In-plane shear modulus `G₁₂`.
    pub g12: f64,
    /// Ply thickness (along the laminate normal).
    pub thickness: f64,
    /// Fibre orientation relative to the laminate x-axis, in degrees.
    pub angle_deg: f64,
}

/// Reduced plane-stress stiffness `[Q11, Q12, Q22, Q66]` in material axes.
fn q_material(p: &Ply) -> [f64; 4] {
    let nu21 = p.nu12 * p.e2 / p.e1;
    let denom = 1.0 - p.nu12 * nu21;
    let q11 = p.e1 / denom;
    let q22 = p.e2 / denom;
    let q12 = p.nu12 * p.e2 / denom;
    let q66 = p.g12;
    [q11, q12, q22, q66]
}

/// Transform a material-axis `[Q11,Q12,Q22,Q66]` into laminate-axis `Q̄` as a
/// `6`-vector `[Q̄11, Q̄12, Q̄22, Q̄16, Q̄26, Q̄66]`.
fn q_bar(q: [f64; 4], theta_deg: f64) -> [f64; 6] {
    let (q11, q12, q22, q66) = (q[0], q[1], q[2], q[3]);
    let t = theta_deg.to_radians();
    let c = t.cos();
    let s = t.sin();
    let c2 = c * c;
    let s2 = s * s;
    let c4 = c2 * c2;
    let s4 = s2 * s2;
    let let_ = c2 * s2;
    let qbar11 = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * let_ + q22 * s4;
    let qbar22 = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * let_ + q22 * c4;
    let qbar12 = (q11 + q22 - 4.0 * q66) * let_ + q12 * (s4 + c4);
    let qbar16 = (q11 - q12 - 2.0 * q66) * s * c2 + (q22 - q12 - 2.0 * q66) * (-s2 * c);
    let qbar26 = (q11 - q12 - 2.0 * q66) * s2 * c + (q22 - q12 - 2.0 * q66) * (-s * c2);
    let qbar66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * let_ + q66 * (s4 + c4);
    [qbar11, qbar12, qbar22, qbar16, qbar26, qbar66]
}

/// Build the `6×6` ABD matrix of a laminate from its (bottom-to-top) ply stack.
///
/// Returns a row-major `6×6` matrix `[A, B; B, D]` where the in-plane force /
/// moment resultants relate to mid-surface strain / curvature by
///
/// ```text
/// {N}   [A  B] {ε⁰}
/// {M} = [B  D] {κ}
/// ```
///
/// A symmetric stack (plies mirrored about the mid-plane) produces `B = 0`.
pub fn laminate_abd(plies: &[Ply]) -> [[f64; 6]; 6] {
    let total: f64 = plies.iter().map(|p| p.thickness).sum();
    // Ply interface z-levels, from bottom (-h/2) to top (+h/2).
    let mut z = vec![-total / 2.0];
    for p in plies {
        z.push(z.last().unwrap() + p.thickness);
    }

    // Per-ply Q̄ mapped into the 3×3 laminate constitutive form
    // [[Q̄11, Q̄12, Q̄16], [Q̄12, Q̄22, Q̄26], [Q̄16, Q̄26, Q̄66]].
    let mut a = [[0.0; 3]; 3];
    let mut b = [[0.0; 3]; 3];
    let mut d = [[0.0; 3]; 3];
    for (i, p) in plies.iter().enumerate() {
        let zb = z[i];
        let zt = z[i + 1];
        let q = q_bar(q_material(p), p.angle_deg);
        // q = [Q̄11, Q̄12, Q̄22, Q̄16, Q̄26, Q̄66]
        let m = [[q[0], q[1], q[3]], [q[1], q[2], q[4]], [q[3], q[4], q[5]]];
        let dz = zt - zb;
        let dz2 = (zt * zt - zb * zb) * 0.5;
        let dz3 = (zt * zt * zt - zb * zb * zb) / 3.0;
        for r in 0..3 {
            for c in 0..3 {
                a[r][c] += m[r][c] * dz;
                b[r][c] += m[r][c] * dz2;
                d[r][c] += m[r][c] * dz3;
            }
        }
    }

    let mut out = [[0.0; 6]; 6];
    for r in 0..3 {
        for c in 0..3 {
            out[r][c] = a[r][c];
            out[r][c + 3] = b[r][c];
            out[r + 3][c] = b[r][c];
            out[r + 3][c + 3] = d[r][c];
        }
    }
    out
}

/// A bilinear cohesive (traction–separation) law for delamination interfaces.
///
/// The normal traction `T = σ̄ · δ / δc` rises linearly from `0` to the peak
/// `σ̄` over `[0, δc)` then drops linearly to `0` at `δf`, so that the enclosed
/// area equals the fracture toughness `Gc`.
#[derive(Clone, Copy, Debug)]
pub struct CohesiveLaw {
    /// Peak traction `σ̄` (mode-independent scalar for this bilinear model).
    pub peak_traction: f64,
    /// Critical separation `δc` where the peak is reached.
    pub critical_opening: f64,
    /// Final separation `δf` where the traction returns to zero.
    pub final_opening: f64,
}

impl CohesiveLaw {
    /// Construct from the fracture toughness `gc`, peak traction `peak`, and a
    /// prescribed `critical_opening` `dc`. `final_opening` is derived so the
    /// triangle area equals `gc`: `δf = δc + 2·gc/σ̄ - δc = 2·gc/σ̄`.
    pub fn from_toughness(gc: f64, peak: f64, dc: f64) -> Self {
        let df = 2.0 * gc / peak;
        CohesiveLaw {
            peak_traction: peak,
            critical_opening: dc,
            final_opening: df,
        }
    }

    /// Traction for an effective opening `delta` (monotonic loading branch).
    pub fn traction(&self, delta: f64) -> f64 {
        let d = delta.abs();
        if d <= 0.0 {
            0.0
        } else if d <= self.critical_opening {
            self.peak_traction * d / self.critical_opening
        } else if d < self.final_opening {
            self.peak_traction * (self.final_opening - d)
                / (self.final_opening - self.critical_opening)
        } else {
            0.0
        }
    }

    /// Fracture toughness (area under the traction–separation curve).
    pub fn toughness(&self) -> f64 {
        0.5 * self.peak_traction * self.final_opening
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_laminate_has_zero_coupling() {
        // [0/90]s symmetric stack should give B = 0.
        let ply = Ply {
            e1: 140e9,
            e2: 10e9,
            nu12: 0.3,
            g12: 5e9,
            thickness: 0.125e-3,
            angle_deg: 0.0,
        };
        let p90 = Ply {
            angle_deg: 90.0,
            ..ply
        };
        let stack = [ply, p90, p90, ply]; // symmetric about mid-plane
        let abd = laminate_abd(&stack);
        for i in 0..3 {
            for j in 3..6 {
                assert!(abd[i][j].abs() < 1e-6, "B[{i}][{j}] = {}", abd[i][j]);
                assert!(abd[j][i].abs() < 1e-6, "B[{j}][{i}] = {}", abd[j][i]);
            }
        }
    }

    #[test]
    fn single_isotropic_ply_a11() {
        let nu = 0.3;
        let e = 70e9;
        let ply = Ply {
            e1: e,
            e2: e,
            nu12: nu,
            g12: e / 2.6,
            thickness: 1e-3,
            angle_deg: 0.0,
        };
        let abd = laminate_abd(&[ply]);
        let q11 = e / (1.0 - nu * nu);
        assert!((abd[0][0] - q11 * 1e-3).abs() / (q11 * 1e-3) < 1e-9);
    }

    #[test]
    fn cohesive_law_sanity() {
        // Peak at δc, zero after δf, area equals Gc.
        let law = CohesiveLaw::from_toughness(1.0, 10.0, 0.1);
        assert!((law.traction(0.0)).abs() < 1e-12);
        assert!((law.traction(0.1) - 10.0).abs() < 1e-9);
        // Just past peak it must be lower than the peak.
        assert!(law.traction(0.15) < 10.0);
        assert!((law.traction(law.final_opening)).abs() < 1e-9);
        assert!(law.traction(law.final_opening * 2.0) == 0.0);
        assert!((law.toughness() - 1.0).abs() < 1e-9);
    }
}
