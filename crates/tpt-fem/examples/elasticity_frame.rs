//! 2-D Euler–Bernoulli frame analysis: a cantilever beam under a tip load.
//!
//! Builds a 2-D `Line` frame mesh by hand, fixes the left end, applies a
//! transverse tip load, and solves with [`solve_frame2d`] (from the umbrella
//! prelude). Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example elasticity_frame
//! ```

use tpt_fem::prelude::*;

fn main() {
    // A horizontal cantilever of length L with N line segments.
    let n_seg = 10usize;
    let length = 2.0;
    let dx = length / n_seg as f64;

    let mut b = MeshBuilder::new();
    let mut nodes = Vec::with_capacity(n_seg + 1);
    for i in 0..=n_seg {
        nodes.push(b.add_node(vec![i as f64 * dx, 0.0]));
    }
    for i in 0..n_seg {
        b.add_element(CellType::Line, vec![nodes[i], nodes[i + 1]]);
    }
    let mesh = b.build();

    // Steel-ish beam section: EA = 200, EI = 10, mass per length = 1.
    let section = BeamSection2D::from_ea_ei(200.0, 10.0, 1.0);

    let tip = nodes[n_seg];
    let tip_load = 1.0;

    // Concentrated transverse load at the free tip; nothing else loaded.
    let loads = |node: usize, _coords: &[f64]| -> [f64; 3] {
        if node == tip {
            [0.0, -tip_load, 0.0]
        } else {
            [0.0; 3]
        }
    };

    // Clamp the root: fix all three DOFs (axial, transverse, rotation).
    let dirichlet = vec![(0, 0.0), (1, 0.0), (2, 0.0)];

    let u = solve_frame2d(&mesh, &section, loads, &dirichlet).expect("frame solve");

    // Tip transverse deflection (DOF index `tip * 3 + 1`).
    let tip_defl = u[tip * 3 + 1];
    println!(
        "Cantilever: {} segments, tip transverse deflection = {:.6} (analytical ≈ -P L³/3EI = {:.6})",
        n_seg,
        tip_defl,
        -tip_load * length.powi(3) / (3.0 * section.ei)
    );

    // Per-node displacement magnitude for ParaView.
    let mut mag = vec![0.0; mesh.node_count()];
    for n in 0..mesh.node_count() {
        let mut s = 0.0;
        for c in 0..3 {
            s += u[n * 3 + c].powi(2);
        }
        mag[n] = s.sqrt();
    }
    write_vtk_with_data(
        &mesh,
        &[PointData::new("disp_mag", mag)],
        "elasticity_frame.vtk",
    )
    .expect("write vtk");
    println!("Wrote elasticity_frame.vtk");
}
