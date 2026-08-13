//! Real end-to-end `mesh → solve → VTK` path through the public umbrella API.
//!
//! Unlike `patch_test.rs` (which keeps a low-level hand-rolled element loop and
//! Dirichlet condensation for coverage), this test drives the *public* physics
//! entry point `solve_poisson` and exports the result through `tpt_fem_io_vtk`,
//! asserting both the numerical answer and that the file was written.

use tpt_fem::prelude::*;

#[test]
fn mesh_solve_vtk_round_trip() {
    // 2-D unit square, four triangles meeting at a centre node.
    let mut b = MeshBuilder::new();
    let c00 = b.add_node(vec![0.0, 0.0]);
    let c10 = b.add_node(vec![1.0, 0.0]);
    let c01 = b.add_node(vec![0.0, 1.0]);
    let c11 = b.add_node(vec![1.0, 1.0]);
    let mid = b.add_node(vec![0.5, 0.5]);
    b.add_element(CellType::Tri, vec![c00, c10, mid]);
    b.add_element(CellType::Tri, vec![c10, c11, mid]);
    b.add_element(CellType::Tri, vec![c11, c01, mid]);
    b.add_element(CellType::Tri, vec![c01, c00, mid]);
    let mesh = b.build();

    // -∇²u = 0 with a harmonic linear Dirichlet field u = x + y. The P1
    // interpolant of a linear function is exact, so the centre node recovers
    // u = 1.0 to machine precision.
    let u = solve_poisson(
        &mesh,
        1.0,
        2,
        |_| 0.0,
        &[(c00, 0.0), (c10, 1.0), (c01, 1.0), (c11, 2.0)],
        None,
        None,
    )
    .expect("solve");

    assert!((u[mid] - 1.0).abs() < 1e-10, "got {}", u[mid]);
    assert!((u[c00] - 0.0).abs() < 1e-12);
    assert!((u[c10] - 1.0).abs() < 1e-12);
    assert!((u[c01] - 1.0).abs() < 1e-12);
    assert!((u[c11] - 2.0).abs() < 1e-12);

    // Export through the real VTK writer and confirm the file was produced.
    let path = std::env::temp_dir().join("tpt_fem_end_to_end.vtk");
    write_vtk_with_data(&mesh, &[PointData::new("u", u)], &path).expect("write vtk");
    let meta = std::fs::metadata(&path).expect("vtk file exists");
    assert!(meta.len() > 0);
    let _ = std::fs::remove_file(&path);
}
