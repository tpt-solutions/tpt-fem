use tpt_fem_mesh_gen::{all_positively_oriented, box_mesh, laplacian_smooth, tet_quality};

fn main() {
    let mut mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3, 3, 3]);

    let q0 = tet_quality(&mesh);
    println!(
        "min dihedral angle (deg) = {:.2}",
        q0.min_dihedral.to_degrees()
    );
    println!("max radius-edge ratio    = {:.2}", q0.max_radius_edge);
    assert!(q0.max_radius_edge.is_finite());
    assert!(q0.min_dihedral > 0.0);
    assert!(all_positively_oriented(&mesh));

    // A few Laplacian-smoothing passes move only interior nodes.
    let worst_disp = laplacian_smooth(&mut mesh, 3);
    println!(
        "worst node displacement after smoothing = {:.3e}",
        worst_disp
    );

    let q1 = tet_quality(&mesh);
    assert!(q1.max_radius_edge.is_finite());
    assert!(all_positively_oriented(&mesh));
    // The corner node stays pinned at the origin.
    assert_eq!(mesh.node_coords(0)[0], 0.0);
}
