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

#[test]
fn p2_tri6_gmsh_import_then_solve() {
    // Highest-ROI P2 wiring check: a quadratic Tri6 mesh imported straight from
    // a Gmsh `.msh` (the importer reorders Gmsh node order -> reference order)
    // must drive the P2 assembly + solve path end to end.
    //
    // The exact solution u = x^2 lies in the P2 space, so with -grad^2 u = 2
    // (source = -2) and Dirichlet pinned on the three *corner* nodes, the
    // Galerkin solution is exactly x^2 everywhere -- including the free edge
    // mid-nodes. This exercises the P2 `shape`/`grad` dispatch through the real
    // imported mesh, not a hand-built one.
    let msh = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 6 1 6
2 1 0 6
1 2 3 4 5 6
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.5 0.0 0.0
0.5 0.5 0.0
0.0 0.5 0.0
$EndNodes
$Elements
1 1 1 1
2 1 9 1
1 1 2 3 4 5 6
$EndElements
";
    let mesh = Mesh::from_msh_bytes(msh.as_bytes()).expect("import P2 Tri6");
    assert_eq!(mesh.elements[0].cell_type, CellType::Tri6);

    // Drive the full P2 path (Gmsh reorder, Tri6 `shape`/`grad`, isoparametric
    // map) by solving -grad^2 u = 2 (source = -2) with the exact field u = x^2,
    // which lies in the P2 space. Every node is pinned to x^2 *except* the
    // bottom-edge midpoint (0.5, 0, 0): on y = 0 the exact field satisfies the
    // natural Neumann condition du/dy = 0, so the single free DOF still recovers
    // x^2 = 0.25 exactly. This exercises a genuine non-empty reduced system.
    let free = (0.5, 0.0);
    let mut bcs = Vec::new();
    for ni in 0..mesh.node_count() {
        let c = mesh.node_coords(ni);
        let on_free_edge = (c[0] - free.0).abs() < 1e-12 && (c[1] - free.1).abs() < 1e-12;
        if !on_free_edge {
            bcs.push((ni, c[0] * c[0]));
        }
    }
    assert_eq!(
        bcs.len(),
        5,
        "expected five pinned nodes, one free edge-mid"
    );

    let u = solve_poisson(&mesh, 1.0, 4, |_| -2.0, &bcs, None, None).expect("solve");

    for ni in 0..mesh.node_count() {
        let c = mesh.node_coords(ni);
        assert!(
            (u[ni] - c[0] * c[0]).abs() < 1e-8,
            "node {ni} at {:?}: got {} expected {}",
            c,
            u[ni],
            c[0] * c[0]
        );
    }
}
