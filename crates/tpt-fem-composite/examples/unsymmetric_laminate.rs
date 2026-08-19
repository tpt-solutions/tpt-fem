//! Demonstrates that an unsymmetric laminate exhibits non-zero
//! extension-bending coupling B in its ABD matrix (a `[0/90]` stack is not
//! mirrored about the mid-plane, so stretching couples to bending).
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-composite --example unsymmetric_laminate
//! ```

use tpt_fem_composite::{laminate_abd, Ply};

fn main() {
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

    // [0/90]: bottom-to-top 0, 90 -> NOT symmetric about the mid-plane.
    let stack = [ply, p90];
    let abd = laminate_abd(&stack);

    let b = block(&abd, 0, 3);
    let mut max_b = 0.0_f64;
    for r in 0..3 {
        for c in 0..3 {
            max_b = max_b.max(b[r][c].abs());
        }
    }

    println!("Unsymmetric [0/90] laminate (2 plies)");
    println!("B (coupling, N):");
    for r in 0..3 {
        println!("  {:14.6e}  {:14.6e}  {:14.6e}", b[r][0], b[r][1], b[r][2]);
    }
    println!("max |B| = {:.3e} N", max_b);
    assert!(
        max_b > 1.0,
        "expected a significant coupling block B for an unsymmetric laminate"
    );
    println!("OK: non-zero coupling block B confirms extension-bending coupling.");
}

fn block(m: &[[f64; 6]; 6], r0: usize, c0: usize) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            out[r][c] = m[r0 + r][c0 + c];
        }
    }
    out
}
