//! Native 3-D mesh generation: structured box meshing and Delaunay tetrahedralisation.
//!
//! Demonstrates [`box_mesh`] (guaranteed-valid structured tets) and
//! [`delaunay_3d`] (incremental Bowyer–Watson on an arbitrary point cloud), both
//! from `tpt-fem-mesh-gen` via the umbrella prelude. Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example mesh_gen_box
//! ```

use tpt_fem::prelude::*;

fn main() {
    // 1. Structured box mesh: the unit cube split into six tets per brick.
    let box1 = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 4, 4]);
    println!(
        "box_mesh: {} nodes, {} tets",
        box1.node_count(),
        box1.element_count()
    );
    write_vtk(&box1, "mesh_gen_box_structured.vtk").expect("write structured");

    // 2. Delaunay tetrahedralisation of a random-ish point cloud.
    let points: Vec<Point3> = {
        // Deterministic lattice-with-jitter so the example is reproducible.
        let mut v = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    let (i, j, k) = (i as f64, j as f64, k as f64);
                    v.push([
                        i + 0.15 * (j * 0.7 + k * 0.3).sin(),
                        j + 0.15 * (k * 0.7 + i * 0.3).cos(),
                        k + 0.15 * (i * 0.7 + j * 0.3).sin(),
                    ]);
                }
            }
        }
        v
    };
    let del = delaunay_3d(&points);
    println!(
        "delaunay_3d: {} input points -> {} tets",
        points.len(),
        del.element_count()
    );
    write_vtk(&del, "mesh_gen_box_delaunay.vtk").expect("write delaunay");

    println!("Wrote mesh_gen_box_structured.vtk and mesh_gen_box_delaunay.vtk");
}
