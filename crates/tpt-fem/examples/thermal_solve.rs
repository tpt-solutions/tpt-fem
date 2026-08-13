//! Solve a 3D Poisson problem end-to-end and export the result to VTK.
//!
//! This example uses *only* the umbrella prelude: every type and function it
//! needs (`box_mesh`, `solve_poisson`, `PointData`, `write_vtk_with_data`)
//! arrives through `tpt_fem::prelude::*`. Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example thermal_solve
//! ```
//!
//! It writes `thermal_solve.vtk` (open in ParaView) containing the temperature
//! field.

use std::path::Path;

use tpt_fem::prelude::*;

fn main() {
    // Structured 3-D tetrahedral mesh of the unit cube.
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [6, 6, 6]);

    // Dirichlet u = 0 on every boundary face (six axis-aligned planes).
    let mut bcs = Vec::new();
    for axis in 0..3 {
        for &coord in &[0.0, 1.0] {
            for nid in mesh.nodes_on_plane(axis, coord, 1e-9) {
                bcs.push((nid, 0.0));
            }
        }
    }

    // -∇·(k ∇u) = 1, with constant conductivity k = 1.
    let u = solve_poisson(&mesh, 1.0, 2, |_| 1.0, &bcs, None, None).expect("solve");

    let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!(
        "Poisson solution on {} tets: u in [{:.6}, {:.6}]",
        mesh.element_count(),
        umin,
        umax
    );

    let path = Path::new("thermal_solve.vtk");
    write_vtk_with_data(&mesh, &[PointData::new("temperature", u)], path).expect("write vtk");
    println!("Wrote {}", path.display());
}
