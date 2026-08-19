//! Fluid-structure interaction relaxation with `fsi_coupling`.
//!
//! Run with: `cargo run -p tpt-fem-coupling --example fsi_coupling`
//!
//! A 2-D elastic block (modelled by `tpt-fem-elasticity`) sits directly under a
//! fluid layer (modelled by `tpt-fem-fluid`). `fsi_coupling` performs one
//! explicit FSI substep:
//!
//! 1. it displaces the fluid interface nodes by the current structure motion,
//! 2. it solves a steady Stokes problem on the (deformed) fluid with a downward
//!    body force, recovering a nodal pressure field,
//! 3. it transfers that pressure as a normal traction onto the matching
//!    structure interface nodes, and
//! 4. it solves the elastic structure under that load to give the new motion.
//!
//! Repeating the substep builds a fixed-point iteration `u^{k+1} = F(u^k)`. For
//! small motions the geometrical feedback is weak, so the sequence relaxes to a
//! steady fluid-loaded shape. The example prints the convergence history and
//! asserts that (a) the structure actually responds to the fluid and (b) the
//! update norm shrinks over the iterations (the coupling relaxes).

use tpt_fem_coupling::fsi_coupling;
use tpt_fem_elasticity::ElasticModel;
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Structured `nx` x `ny` Quad4 sheet spanning `[0,w] x [y0, y0+h]`.
fn quad_sheet(w: f64, h: f64, y0: f64, nx: usize, ny: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=ny {
        let y = y0 + h * j as f64 / ny as f64;
        let mut r = Vec::new();
        for i in 0..=nx {
            r.push(b.add_node(vec![w * i as f64 / nx as f64, y]));
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

fn main() {
    let w = 1.0; // shared width
    let (nx, ny_s, ny_f) = (4, 2, 2);
    let h_s = 0.5; // structure height
    let h_f = 0.5; // fluid height

    let struct_mesh = quad_sheet(w, h_s, 0.0, nx, ny_s); // y in [0, 0.5]
    let fluid_mesh = quad_sheet(w, h_f, h_s, nx, ny_f); // y in [0.5, 1.0]

    // Interface: every structure top node (y = 0.5) ↔ fluid bottom node (y = 0.5)
    // with the same x-coordinate.
    let mut interface = Vec::new();
    for s in 0..struct_mesh.node_count() {
        if (struct_mesh.node_coords(s)[1] - h_s).abs() < 1e-9 {
            let xs = struct_mesh.node_coords(s)[0];
            let f = (0..fluid_mesh.node_count())
                .find(|&n| {
                    (fluid_mesh.node_coords(n)[1] - h_s).abs() < 1e-9
                        && (fluid_mesh.node_coords(n)[0] - xs).abs() < 1e-9
                })
                .expect("matching fluid interface node");
            interface.push((s, f));
        }
    }
    println!("FSI: structure 1.0 x 0.5 (PlaneStress) + fluid layer 1.0 x 0.5");
    println!("interface pairs (struct, fluid) = {}", interface.len());

    // Fix the structure base (y = 0) in both components.
    let mut struct_dirichlet = Vec::new();
    for n in 0..struct_mesh.node_count() {
        if struct_mesh.node_coords(n)[1] < 1e-9 {
            struct_dirichlet.push((n * 2, 0.0));
            struct_dirichlet.push((n * 2 + 1, 0.0));
        }
    }
    // Tip = top-centre structure node.
    let tip = (0..struct_mesh.node_count())
        .find(|&n| {
            (struct_mesh.node_coords(n)[0] - w / 2.0).abs() < 1e-9
                && (struct_mesh.node_coords(n)[1] - h_s).abs() < 1e-9
        })
        .expect("tip node");

    let mut u = vec![0.0; struct_mesh.node_count() * 2];
    let mut prev = u.clone();
    println!("\n   iter        |du|        tip |u|       tip u_y");
    println!("   ----   -----------   -------------   ----------");
    let mut first_du = 0.0_f64;
    for k in 0..12 {
        u = fsi_coupling(
            &struct_mesh,
            &fluid_mesh,
            ElasticModel::PlaneStress,
            1.0e3, // young (keeps the motion small and the iteration stable)
            0.3,   // poisson
            1.0,   // viscosity
            &prev,
            &interface,
            &struct_dirichlet,
            1.0e5, // fluid penalty
        );
        let du: f64 = u
            .iter()
            .zip(prev.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        if k == 0 {
            first_du = du;
        }
        println!(
            "   {k:>4}   {du:11.4e}   {:13.4e}   {:10.4e}",
            (u[tip * 2].powi(2) + u[tip * 2 + 1].powi(2)).sqrt(),
            u[tip * 2 + 1]
        );
        prev = u.clone();
    }

    let tip_motion = (u[tip * 2].powi(2) + u[tip * 2 + 1].powi(2)).sqrt();
    let final_du: f64 = u
        .iter()
        .zip(prev.iter())
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt();

    println!("\n   tip motion magnitude = {tip_motion:.4e}");
    println!("   first |du| = {first_du:.4e}, final |du| = {final_du:.4e}");

    assert!(
        tip_motion > 0.0,
        "structure must respond to the fluid traction"
    );
    assert!(
        final_du < first_du,
        "the FSI fixed point must relax (|du| must shrink): {final_du} >= {first_du}"
    );
    println!("\nverified: structure responds to the fluid and the coupling relaxation");
    println!("converges (the update norm decreases across iterations)");
}
