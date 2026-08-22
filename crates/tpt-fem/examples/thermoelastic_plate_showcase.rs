//! Coupled thermo-elastic pipeline: conduction -> thermal expansion -> VTK.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example thermoelastic_plate_showcase --features coupling
//! ```
//!
//! Stage 1 solves steady heat conduction on a plane-stress strip (Quad4) with
//! `solve_poisson`: hot at the left edge, cold at the right edge, no source.
//! The exact solution is the linear field `T(x) = T_hot * (1 - x/L)`; bilinear
//! quadrilaterals reproduce linear fields exactly, so the FEM temperature must
//! match it to machine precision.
//!
//! Stage 2 feeds that per-node temperature into `thermal_structural` as a
//! `delta_T` and clamps the left edge. Two thermal loads are run:
//!
//! * a **uniform** rise `dt0`, whose free-expansion tip elongation is
//!   `alpha * dt0 * L` (the clamp restrains lateral Poisson contraction near
//!   the fixed edge, so FEM lands slightly below);
//! * the **computed field** from stage 1, whose free-expansion elongation is
//!   `alpha * int_0^L T(x) dx = alpha * T_hot * L / 2`.
//!
//! A sweep over load scales confirms both responses are exactly linear, and
//! everything is exported to ParaView via `tpt-fem-io-vtk`.

use std::path::Path;

use tpt_fem_coupling::thermal_structural;
use tpt_fem_elasticity::ElasticModel;
use tpt_fem_io_vtk::{write_vtk_with_data, PointData};
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_fem_thermal::solve_poisson;

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

/// Mid-surface nodes sorted by x; returns `(ux, uy)` at the tip.
fn midplane_tip(mesh: &Mesh, u: &[f64], h: f64) -> (f64, f64) {
    let mut m: Vec<usize> = (0..mesh.node_count())
        .filter(|&n| (mesh.node_coords(n)[1] - h / 2.0).abs() < 1e-12)
        .collect();
    m.sort_by(|&a, &b| {
        mesh.node_coords(a)[0]
            .partial_cmp(&mesh.node_coords(b)[0])
            .unwrap()
    });
    let tip = *m.last().unwrap();
    (u[tip * 2], u[tip * 2 + 1])
}

fn main() {
    let l = 1.0;
    let h = 0.1;
    let (nx, ny) = (40usize, 4usize);
    let (e, nu) = (70.0e9, 0.33); // aluminium
    let alpha = 2.3e-5; // 1/K
    let k_cond = 200.0; // W/(m K)
    let t_hot = 100.0; // K above reference at x = 0

    let mesh = strip(l, h, nx, ny);

    // ---- Stage 1: steady conduction ------------------------------------
    let mut bcs = Vec::new();
    for n in 0..mesh.node_count() {
        let x = mesh.node_coords(n)[0];
        if x < 1e-9 {
            bcs.push((n, t_hot));
        } else if (x - l).abs() < 1e-9 {
            bcs.push((n, 0.0));
        }
    }
    let temp = solve_poisson(&mesh, k_cond, 2, |_| 0.0, &bcs, None, None).expect("heat solve");

    // Linear fields are in the Quad4 space: verify against T_hot*(1 - x/L).
    let mut max_err = 0.0_f64;
    for n in 0..mesh.node_count() {
        let x = mesh.node_coords(n)[0];
        max_err = max_err.max((temp[n] - t_hot * (1.0 - x / l)).abs());
    }
    println!("stage 1: steady conduction on {nx} x {ny} Quad4 strip");
    println!(
        "  T in [{:.4}, {:.4}] K, max |T - linear| = {max_err:.3e} K",
        0.0,
        t_hot
    );
    assert!(max_err < 1e-8, "linear temperature field must be exact");

    // Clamp the left edge in both components.
    let fixed: Vec<(usize, f64)> = (0..mesh.node_count())
        .filter(|&n| mesh.node_coords(n)[0] < 1e-9)
        .flat_map(|n| [(n * 2, 0.0), (n * 2 + 1, 0.0)])
        .collect();

    let solve_case = |dt: &Vec<f64>| -> Vec<f64> {
        thermal_structural(&mesh, ElasticModel::PlaneStress, e, nu, alpha, dt, &fixed)
            .expect("thermal-structural solve")
    };

    // ---- Stage 2: thermal expansion ------------------------------------
    let dt_field = temp.clone();
    let dt0 = 50.0;
    let dt_uniform = vec![dt0; mesh.node_count()];

    let u_uni = solve_case(&dt_uniform);
    let (ux_uni, uy_uni) = midplane_tip(&mesh, &u_uni, h);
    let uni_exact = alpha * dt0 * l;
    println!("\nstage 2a: uniform dT = {dt0} K");
    println!(
        "  tip ux = {ux_uni:.6e} vs alpha*dT*L = {uni_exact:.6e} (ratio {:.4})",
        ux_uni / uni_exact
    );
    assert!(
        (0.85..=1.02).contains(&(ux_uni / uni_exact)),
        "uniform expansion must approach the free-expansion value"
    );
    assert!(uy_uni.abs() < 1e-9, "mid-plane stays straight");

    let u_field = solve_case(&dt_field);
    let (ux_fld, _) = midplane_tip(&mesh, &u_field, h);
    let fld_exact = alpha * t_hot * l / 2.0;
    println!("\nstage 2b: computed conduction field");
    println!(
        "  tip ux = {ux_fld:.6e} vs alpha*int(T)dx = {fld_exact:.6e} (ratio {:.4})",
        ux_fld / fld_exact
    );
    assert!(
        (0.7..=1.05).contains(&(ux_fld / fld_exact)),
        "field case should be near the integral free-expansion estimate"
    );

    // ---- Linearity sweep -------------------------------------------------
    println!("\nload-scale sweep (response must be exactly linear)");
    for s in [0.25, 0.5, 2.0] {
        let a = solve_case(&dt_uniform.iter().map(|t| s * t).collect::<Vec<f64>>());
        let b = solve_case(&dt_field.iter().map(|t| s * t).collect::<Vec<f64>>());
        let ra = midplane_tip(&mesh, &a, h).0 / ux_uni;
        let rb = midplane_tip(&mesh, &b, h).0 / ux_fld;
        println!("  scale {s:5}: uniform ratio = {ra:.10}, field ratio = {rb:.10}");
        assert!((ra - s).abs() < 1e-9 && (rb - s).abs() < 1e-9);
    }

    // ---- Export ----------------------------------------------------------
    let ux: Vec<f64> = (0..mesh.node_count()).map(|n| u_field[n * 2]).collect();
    let uy: Vec<f64> = (0..mesh.node_count()).map(|n| u_field[n * 2 + 1]).collect();
    let path = Path::new("thermoelastic_plate.vtk");
    write_vtk_with_data(
        &mesh,
        &[
            PointData::new("temperature", temp),
            PointData::new("ux", ux),
            PointData::new("uy", uy),
        ],
        path,
    )
    .expect("write vtk");
    println!("\nwrote {}", path.display());
    println!("OK: conduction -> expansion pipeline verified and exported.");
}

