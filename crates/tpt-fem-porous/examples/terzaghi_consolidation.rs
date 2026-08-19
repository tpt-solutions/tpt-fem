//! 1-D Terzaghi consolidation of a saturated column resolved with the FEM and
//! checked against the closed-form solution for the average degree of
//! consolidation U(Tv) under one-way drainage:
//!
//!   U(Tv) = 1 - Σ_{m=1,3,5,...} 8/((2m-1)² π²) exp(-(2m-1)² π² Tv / 4)
//!
//! where Tv = cv t / H² is the dimensionless time factor.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-porous --example terzaghi_consolidation
//! ```

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_fem_porous::terzaghi_consolidation;

fn line_mesh(n: usize, length: f64) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut prev = b.add_node(vec![0.0]);
    for i in 1..=n {
        let node = b.add_node(vec![length * i as f64 / n as f64]);
        b.add_element(CellType::Line, vec![prev, node]);
        prev = node;
    }
    b.build()
}

/// Closed-form average degree of consolidation (one-way drainage).
fn degree_of_consolidation(tv: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let mut u = 0.0_f64;
    for m in 1..=200 {
        let mm = 2 * m - 1; // odd integers
        let mmf = mm as f64;
        u += 8.0 / (mmf * mmf * pi * pi) * (-(mmf * mmf * pi * pi * tv / 4.0)).exp();
    }
    1.0 - u
}

fn main() {
    let n = 20;
    let h = 1.0;
    let mesh = line_mesh(n, h);
    let q0 = 100.0;
    let ev = 1e6;
    let cv = 0.1;
    let dt = 0.005;
    // Largest target time factor -> total time.
    let tv_max = 1.0;
    let total = tv_max * h * h / cv; // = 10.0
    let hist = terzaghi_consolidation(&mesh, q0, cv, ev, total, dt);

    let s_inf = q0 * h / ev; // closed-form drained settlement
    let targets = [0.2_f64, 0.5, 1.0];

    println!("Terzaghi 1-D consolidation (H={h}, cv={cv}, Ev={:e})", ev);
    println!(
        "  {:>6} {:>10} {:>12} {:>12} {:>8}",
        "Tv", "t", "U(num)", "U(analyt)", "err"
    );
    let mut worst = 0.0_f64;
    for &tv in &targets {
        let t_target = tv * h * h / cv;
        let rec = hist
            .iter()
            .min_by(|a, b| {
                (a.0 - t_target)
                    .abs()
                    .partial_cmp(&(b.0 - t_target).abs())
                    .unwrap()
            })
            .unwrap();
        let u_num = rec.1 / s_inf; // degree of consolidation
        let u_an = degree_of_consolidation(tv);
        let err = (u_num - u_an).abs();
        worst = worst.max(err);
        println!(
            "  {:6.2} {:10.4} {:12.6} {:12.6} {:8.2e}",
            tv, rec.0, u_num, u_an, err
        );
    }
    println!();
    println!("worst |U(num)-U(analyt)| = {:.3e}", worst);
    assert!(
        worst < 0.05,
        "FE degree of consolidation must match the series"
    );
    println!("OK: FE settlement history matches the Terzaghi closed form.");
}
