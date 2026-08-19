//! Bimetal-style thermal bending from a through-thickness temperature gradient.
//!
//! Run with: `cargo run -p tpt-fem-coupling --example thermal_bimetal_strip`
//!
//! A plane-stress cantilever strip (length `L`, thickness `h`) is heated with a
//! linear temperature profile through the thickness, `dT(y)` running from
//! `-dT0` at the bottom face to `+dT0` at the top face. `thermal_structural`
//! turns that into an initial strain `eps_th = alpha * dT(y) * I`, whose
//! antisymmetric part is a pure bending load. Beam theory gives the thermal
//! curvature and cantilever tip deflection
//!
//! ```text
//! kappa = -alpha * (dT_top - dT_bot) / h,     v(x) = kappa * x^2 / 2
//! ```
//!
//! so the hotter (longer) top fibre makes the strip curl downwards, `v < 0`.
//!
//! The example prints the deflected centreline against `kappa x^2 / 2`, sweeps
//! `dT0` to confirm the response is linear in the gradient, and asserts the sign
//! and magnitude of the tip deflection. With 20 x 4 bilinear (`Quad4`) elements
//! the tip deflection lands about 1% above the beam-theory value (the clamped
//! end restrains the transverse thermal strain), so 5% is asserted.

use tpt_fem_coupling::thermal_structural;
use tpt_fem_elasticity::ElasticModel;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Structured `nx` x `ny` Quad4 strip of length `l` and thickness `h`.
fn strip(l: f64, h: f64, nx: usize, ny: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=ny {
        let y = h * j as f64 / ny as f64;
        let mut r = Vec::new();
        for i in 0..=nx {
            r.push(b.add_node(vec![l * i as f64 / nx as f64, y]));
        }
        rows.push(r);
    }
    for j in 0..ny {
        for i in 0..nx {
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

/// Solve one gradient and return `(tip deflection, analytic tip deflection)`.
fn bend(mesh: &Mesh, l: f64, h: f64, alpha: f64, dt0: f64, verbose: bool) -> (f64, f64) {
    // Linear profile: dT(y) = dt0 * (2y/h - 1), i.e. -dt0 at y=0, +dt0 at y=h.
    let temp: Vec<f64> = (0..mesh.node_count())
        .map(|n| dt0 * (2.0 * mesh.node_coords(n)[1] / h - 1.0))
        .collect();
    // Clamp the left edge in both components.
    let mut fixed = Vec::new();
    for n in 0..mesh.node_count() {
        if mesh.node_coords(n)[0] < 1e-9 {
            fixed.push((n * 2, 0.0));
            fixed.push((n * 2 + 1, 0.0));
        }
    }
    let u = thermal_structural(
        mesh,
        ElasticModel::PlaneStress,
        70.0e9, // aluminium, Pa
        0.33,
        alpha,
        &temp,
        &fixed,
    )
    .expect("thermal-structural solve");

    let kappa = -alpha * (2.0 * dt0) / h;
    // Mid-surface row (y = h/2) exists whenever ny is even.
    let mut mid: Vec<usize> = (0..mesh.node_count())
        .filter(|&n| (mesh.node_coords(n)[1] - h / 2.0).abs() < 1e-12)
        .collect();
    mid.sort_by(|&a, &b| {
        mesh.node_coords(a)[0]
            .partial_cmp(&mesh.node_coords(b)[0])
            .unwrap()
    });
    if verbose {
        println!("      x        v (FEM)      kappa x^2/2     ratio");
        println!("   -------   -----------   -------------   -------");
        for &n in &mid {
            let x = mesh.node_coords(n)[0];
            let beam = kappa * x * x / 2.0;
            let ratio = if beam.abs() > 0.0 {
                u[n * 2 + 1] / beam
            } else {
                1.0
            };
            println!(
                "   {x:7.3}   {:11.3e}   {beam:13.3e}   {ratio:7.4}",
                u[n * 2 + 1]
            );
        }
    }
    let tip = *mid.last().expect("mid-surface row");
    (u[tip * 2 + 1], kappa * l * l / 2.0)
}

fn main() {
    let l = 1.0; // length, m
    let h = 0.1; // thickness, m
    let (nx, ny) = (20, 4);
    let alpha = 2.3e-5; // 1/K, aluminium
    let dt0 = 20.0; // K at each face (gradient = 40 K across h)

    let mesh = strip(l, h, nx, ny);
    println!("plane-stress strip: L = {l}, h = {h}, {nx} x {ny} Quad4");
    println!("alpha = {alpha:.3e} 1/K, dT = -{dt0} K (bottom) .. +{dt0} K (top)");
    println!(
        "thermal curvature kappa = -alpha*(dT_top - dT_bot)/h = {:.4e} 1/m\n",
        -alpha * 2.0 * dt0 / h
    );
    let (tip, beam_tip) = bend(&mesh, l, h, alpha, dt0, true);
    println!("\n   tip v (FEM)       = {tip:.6e}");
    println!("   tip kappa L^2 / 2 = {beam_tip:.6e}");
    println!("   FEM / beam theory = {:.4}", tip / beam_tip);

    println!("\ngradient sweep (response must be linear in dT0)");
    println!("    dT0 [K]     tip v (FEM)     beam theory     ratio");
    println!("   ---------   -------------   -------------   -------");
    for d in [5.0, 10.0, 20.0, 40.0] {
        let (t, b) = bend(&mesh, l, h, alpha, d, false);
        println!("   {d:9.1}   {t:13.4e}   {b:13.4e}   {:7.4}", t / b);
        assert!(t < 0.0, "hotter top face must curl the strip downwards");
        assert!(
            (t / b - tip / beam_tip).abs() < 1e-6,
            "response must be linear in the gradient"
        );
    }

    assert!(tip < 0.0, "tip must deflect downwards, got {tip:e}");
    assert!(
        (tip / beam_tip - 1.0).abs() < 0.05,
        "tip deflection within 5% of beam theory, ratio = {:.4}",
        tip / beam_tip
    );
    println!("\nverified: downward curl, linear in the gradient, and within 5% of");
    println!("the beam-theory tip deflection kappa L^2 / 2");
}
