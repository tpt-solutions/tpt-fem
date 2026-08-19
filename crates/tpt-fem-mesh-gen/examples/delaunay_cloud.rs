use tpt_fem_mesh_gen::{all_positively_oriented, delaunay_3d};

fn main() {
    // Eight corners of the unit cube: Delaunay tetrahedralises its convex hull.
    let pts = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    let mesh = delaunay_3d(&pts);

    println!("node_count   = {} (expected 8)", mesh.node_count());
    println!("element_count = {} (>= 5 for a cube)", mesh.element_count());

    assert_eq!(mesh.node_count(), 8);
    assert!(mesh.element_count() >= 5);
    assert!(all_positively_oriented(&mesh));
}
