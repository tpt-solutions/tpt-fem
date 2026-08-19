//! Demonstrates classical lamination theory for a symmetric cross-ply
//! `[0/90]s` laminate: builds the 6x6 ABD matrix and verifies that the
//! extension-bending coupling block B is (numerically) zero.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-composite --example laminate_abd
//! ```

use tpt_fem_composite::{laminate_abd, Ply};

fn main() {
    // Material: a typical unidirectional graphite/epoxy (Jones, 1999).
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

    // [0/90]s: bottom-to-top 0, 90, 90, 0 (symmetric about mid-plane).
    let stack = [ply, p90, p90, ply];
    let abd = laminate_abd(&stack);

    let a = block(&abd, 0, 0);
    let b = block(&abd, 0, 3);
    let d = block(&abd, 3, 3);

    println!(
        "Symmetric cross-ply [0/90]s  (4 plies, t = {:.3} mm each)",
        ply.thickness * 1e3
    );
    println!();
    print_block("A  (extensional, N/m)", &a);
    print_block("B  (coupling, N)", &b);
    print_block("D  (bending, N.m)", &d);

    // Symmetric stacking => no extension-bending coupling.
    let mut max_b = 0.0_f64;
    for r in 0..3 {
        for c in 0..3 {
            max_b = max_b.max(b[r][c].abs());
        }
    }
    println!("max |B| = {:.3e} N", max_b);
    assert!(max_b < 1e-3, "expected B ~ 0 for a symmetric laminate");
    println!("OK: coupling block B is zero (symmetric laminate).");
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

fn print_block(title: &str, b: &[[f64; 3]; 3]) {
    println!("{title}");
    for r in 0..3 {
        println!("  {:14.6e}  {:14.6e}  {:14.6e}", b[r][0], b[r][1], b[r][2]);
    }
    println!();
}
