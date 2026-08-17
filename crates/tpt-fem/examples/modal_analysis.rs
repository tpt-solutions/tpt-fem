//! Modal analysis: natural frequencies of a 3-D cantilever block.
//!
//! Builds a structured tetrahedral mesh of a 3-D block with [`box_mesh`], clamps
//! one face, and solves the generalized eigenproblem `K Φ = ω² M Φ` with
//! [`solve_modal`] (from the umbrella prelude). Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example modal_analysis
//! ```

use tpt_fem::prelude::*;

fn main() {
    // A slender 3-D block: long in x, square (and short) in y/z.
    let mesh = box_mesh([0.0, 0.0, 0.0], [2.0, 0.2, 0.2], [16, 2, 2]);

    // Material: steel-ish (E, ν, ρ).
    let young = 200.0e9;
    let poisson = 0.3;
    let density = 7800.0;

    // Clamp the x = 0 face: fix all three displacement DOFs of every node there.
    let mut dirichlet = Vec::new();
    for n in mesh.nodes_on_plane(0, 0.0, 1e-9) {
        for c in 0..3 {
            dirichlet.push((n * 3 + c, 0.0));
        }
    }

    let modes = solve_modal(
        &mesh,
        ElasticModel::Continuum3D,
        young,
        poisson,
        density,
        2,
        4,
        &dirichlet,
    )
    .expect("modal solve");

    println!("Cantilever block — fundamental natural frequencies:");
    for (i, (lam, phi)) in modes.iter().enumerate() {
        // `lam` is ω²; report the frequency in Hz (ω / 2π).
        let omega = lam.sqrt();
        println!(
            "  mode {}: ω² = {:.3e}, f = {:.3e} Hz",
            i + 1,
            lam,
            omega / (2.0 * std::f64::consts::PI)
        );

        // Export the first mode shape for ParaView.
        if i == 0 {
            let mut mag = vec![0.0; mesh.node_count()];
            for n in 0..mesh.node_count() {
                let mut s = 0.0;
                for c in 0..3 {
                    s += phi[n * 3 + c].powi(2);
                }
                mag[n] = s.sqrt();
            }
            write_vtk_with_data(&mesh, &[PointData::new("mode1", mag)], "modal_analysis.vtk")
                .expect("write vtk");
            println!("Wrote modal_analysis.vtk (mode-1 shape)");
        }
    }
}
