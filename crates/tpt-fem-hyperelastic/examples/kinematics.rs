//! `Mat3` matrix-kernel demo: deformation gradient → J, C, invariants.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-hyperelastic --example kinematics
//! ```
//!
//! Exercises the small linear-algebra helpers (`mat_mul`, `mat_det`,
//! `mat_inv`, `mat_transpose`) on a deformation gradient `F`, checking the
//! determinant `J = det F`, the right Cauchy–Green tensor `C = FᵀF`, and the
//! principal invariants against closed forms.

use tpt_fem_hyperelastic::{mat_det, mat_inv, mat_mul, mat_transpose, Mat3};

fn main() {
    // Uniform stretch λ along x, plus a simple shear in the xy plane.
    let lam = 2.0;
    let mut f: Mat3 = [[0.0; 3]; 3];
    f[0][0] = lam;
    f[1][1] = 1.0;
    f[2][2] = 1.0;
    f[0][1] = 0.3; // shear

    let j = mat_det(&f);
    let ft = mat_transpose(&f);
    let c = mat_mul(&ft, &f);
    let i1 = c[0][0] + c[1][1] + c[2][2];
    let i3 = j * j;

    // Closed forms: J = λ (a shear preserves volume); I1 = λ² + 2 + shear².
    let j_closed = lam;
    let i1_closed = lam * lam + 2.0 + 0.3 * 0.3;

    println!("F =");
    for r in 0..3 {
        println!("  [{:.3}, {:.3}, {:.3}]", f[r][0], f[r][1], f[r][2]);
    }
    println!("J = det F        = {:.6}  (closed {:.6})", j, j_closed);
    println!("I1 = tr(C)       = {:.6}  (closed {:.6})", i1, i1_closed);
    println!(
        "I3 = J²          = {:.6}  (closed {:.6})",
        i3,
        j_closed * j_closed
    );

    let finv = mat_inv(&f).expect("F is invertible");
    let prod = mat_mul(&f, &finv);
    println!("F·F⁻¹ =");
    for r in 0..3 {
        println!(
            "  [{:.3}, {:.3}, {:.3}]",
            prod[r][0], prod[r][1], prod[r][2]
        );
    }

    assert!((j - j_closed).abs() < 1e-9);
    assert!((i1 - i1_closed).abs() < 1e-9);
    assert!((i3 - j_closed * j_closed).abs() < 1e-9);
    // F·F⁻¹ = I.
    for r in 0..3 {
        for cc in 0..3 {
            let want = if r == cc { 1.0 } else { 0.0 };
            assert!((prod[r][cc] - want).abs() < 1e-9, "F·F⁻¹ off at ({r},{cc})");
        }
    }
    println!("\nOK: determinant, invariants, and inverse verified.");
}
