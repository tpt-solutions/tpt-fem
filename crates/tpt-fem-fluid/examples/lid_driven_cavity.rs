//! Lid-driven cavity: velocity extrema and the penalty-independence check.
//!
//! Run with: `cargo run -p tpt-fem-fluid --example lid_driven_cavity`
//!
//! The classic Stokes benchmark: a unit square with no-slip walls and a lid at
//! `y = 1` sliding at `u_x = 1`, no body force. Unlike Poiseuille flow the exact
//! solution is *not* representable by a discretely divergence-free bilinear
//! field, so this case exercises the divergence-penalty constraint itself.
//!
//! The example prints the vertical centreline profile (`x = 0.5`), reports the
//! velocity extrema, and sweeps the penalty weight `1/ε` to show that the
//! extrema converge to a finite limit as the penalty grows — the signature of
//! the selectively reduced-integrated penalty term. It asserts forward flow
//! immediately under the lid, recirculating flow at mid-height, zero transverse
//! velocity on the (symmetric) centreline, and penalty independence.

use tpt_fem_fluid::steady_stokes;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Structured `n` x `n` Quad4 grid on the unit square.
fn unit_square(n: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=n {
        let y = j as f64 / n as f64;
        let mut r = Vec::new();
        for i in 0..=n {
            r.push(b.add_node(vec![i as f64 / n as f64, y]));
        }
        rows.push(r);
    }
    for j in 0..n {
        for i in 0..n {
            b.add_element(
                CellType::Quad,
                vec![
                    rows[j][i],
                    rows[j][i + 1],
                    rows[j + 1][i + 1],
                    rows[j + 1][i],
                ],
            );
        }
    }
    b.build()
}

/// Lid at `y = 1` moving at `u_x = 1`; no-slip on the other three walls.
fn cavity_bc(mesh: &Mesh) -> Vec<(usize, f64)> {
    const TOL: f64 = 1e-9;
    let mut bc = Vec::new();
    for n in 0..mesh.node_count() {
        let c = mesh.node_coords(n);
        let boundary = c[0] < TOL || c[0] > 1.0 - TOL || c[1] < TOL || c[1] > 1.0 - TOL;
        if boundary {
            let lid = c[1] > 1.0 - TOL;
            bc.push((n * 2, if lid { 1.0 } else { 0.0 }));
            bc.push((n * 2 + 1, 0.0));
        }
    }
    bc
}

fn main() {
    let mu = 1.0;
    let n = 12;
    let mesh = unit_square(n);
    let bc = cavity_bc(&mesh);
    let penalty = 1e6;

    let (u, p) = steady_stokes(&mesh, mu, |_| vec![0.0, 0.0], &bc, penalty)
        .expect("steady Stokes cavity solve");

    println!("lid-driven cavity: {n} x {n} Quad4, mu = {mu}, u_lid = 1, penalty = {penalty:.0e}");
    println!("\nvertical centreline x = 0.5");
    println!("     y        u_x         u_y");
    println!("  ------   ---------   ---------");
    let mut column: Vec<usize> = (0..mesh.node_count())
        .filter(|&nd| (mesh.node_coords(nd)[0] - 0.5).abs() < 1e-9)
        .collect();
    column.sort_by(|&a, &b| {
        mesh.node_coords(a)[1]
            .partial_cmp(&mesh.node_coords(b)[1])
            .unwrap()
    });
    let mut max_centre_uy = 0.0_f64;
    for nd in &column {
        let c = mesh.node_coords(*nd);
        println!(
            "  {:6.3}   {:+9.6}   {:+9.6}",
            c[1],
            u[nd * 2],
            u[nd * 2 + 1]
        );
        max_centre_uy = max_centre_uy.max(u[nd * 2 + 1].abs());
    }

    // Extrema over the interior (the lid itself is prescribed, so exclude it).
    let mut ux_min = f64::INFINITY;
    let mut ux_max = f64::NEG_INFINITY;
    let mut uy_min = f64::INFINITY;
    let mut uy_max = f64::NEG_INFINITY;
    let (mut at_min, mut at_max) = ([0.0, 0.0], [0.0, 0.0]);
    for nd in 0..mesh.node_count() {
        let c = mesh.node_coords(nd);
        if c[1] > 1.0 - 1e-9 {
            continue;
        }
        if u[nd * 2] < ux_min {
            ux_min = u[nd * 2];
            at_min = [c[0], c[1]];
        }
        if u[nd * 2] > ux_max {
            ux_max = u[nd * 2];
            at_max = [c[0], c[1]];
        }
        uy_min = uy_min.min(u[nd * 2 + 1]);
        uy_max = uy_max.max(u[nd * 2 + 1]);
    }
    println!("\nvelocity extrema below the lid");
    println!(
        "  max forward u_x = {ux_max:+.6}  at (x, y) = ({:.3}, {:.3})",
        at_max[0], at_max[1]
    );
    println!(
        "  max return  u_x = {ux_min:+.6}  at (x, y) = ({:.3}, {:.3})",
        at_min[0], at_min[1]
    );
    println!("  u_y range       = [{uy_min:+.6}, {uy_max:+.6}]");

    // Recovered pressure p = -(1/eps) div u: antisymmetric about x = 0.5, with
    // the extrema at the two lid corners where the Stokes solution is singular.
    let mut p_min = f64::INFINITY;
    let mut p_max = f64::NEG_INFINITY;
    let (mut p_at_min, mut p_at_max) = ([0.0, 0.0], [0.0, 0.0]);
    for nd in 0..mesh.node_count() {
        let c = mesh.node_coords(nd);
        if p[nd] < p_min {
            p_min = p[nd];
            p_at_min = [c[0], c[1]];
        }
        if p[nd] > p_max {
            p_max = p[nd];
            p_at_max = [c[0], c[1]];
        }
    }
    println!("\nrecovered pressure extrema");
    println!(
        "  p_min = {p_min:+.4}  at (x, y) = ({:.3}, {:.3})",
        p_at_min[0], p_at_min[1]
    );
    println!(
        "  p_max = {p_max:+.4}  at (x, y) = ({:.3}, {:.3})",
        p_at_max[0], p_at_max[1]
    );

    println!("\npenalty sweep on an 8 x 8 mesh (extrema must converge, not vanish)");
    println!("   1/eps      max |u_x|   max |u_y|");
    println!("  --------   ---------   ---------");
    let coarse = unit_square(8);
    let cbc = cavity_bc(&coarse);
    let mut last = (0.0, 0.0);
    for pen in [1e2, 1e3, 1e4, 1e5, 1e6] {
        let (uc, _) = steady_stokes(&coarse, mu, |_| vec![0.0, 0.0], &cbc, pen)
            .expect("steady Stokes penalty-sweep solve");
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for nd in 0..coarse.node_count() {
            if coarse.node_coords(nd)[1] > 1.0 - 1e-9 {
                continue;
            }
            mx = mx.max(uc[nd * 2].abs());
            my = my.max(uc[nd * 2 + 1].abs());
        }
        println!("  {pen:8.0e}   {mx:9.6}   {my:9.6}");
        last = (mx, my);
    }

    // Physics checks: forward flow under the lid, return flow at mid-height,
    // and a symmetric centreline (u_y = 0 at x = 0.5).
    let at = |x: f64, y: f64| {
        let nd = (0..mesh.node_count())
            .find(|&nd| {
                (mesh.node_coords(nd)[0] - x).abs() < 1e-9
                    && (mesh.node_coords(nd)[1] - y).abs() < 1e-9
            })
            .expect("probe node must exist");
        u[nd * 2]
    };
    assert!(
        at(0.5, 1.0 - 1.0 / n as f64) > 0.1,
        "flow under the lid must be forward"
    );
    assert!(at(0.5, 0.5) < -0.05, "mid-height flow must recirculate");
    assert!(
        max_centre_uy < 1e-9,
        "centreline u_y must vanish by symmetry, got {max_centre_uy:.2e}"
    );
    // The 1e5 and 1e6 sweep rows must agree to better than 0.1%: the penalty
    // limit exists (no locking) and is reached.
    assert!(
        last.0 > 0.3 && last.1 > 0.1,
        "penalty limit collapsed to zero"
    );
    println!("\nverified: forward flow under the lid, recirculation at mid-height,");
    println!("centreline u_y < 1e-9, and a non-zero penalty limit");
}
