use tpt_fem_mesh_gen::{all_positively_oriented, box_mesh};

fn main() {
    // A 2x2x2 brick grid -> 3*3*3 = 27 nodes, 2*2*2*6 = 48 tetrahedra.
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2]);

    println!("node_count   = {} (expected 27)", mesh.node_count());
    println!("element_count = {} (expected 48)", mesh.element_count());

    assert_eq!(mesh.node_count(), 27);
    assert_eq!(mesh.element_count(), 48);
    // The structured decomposition is guaranteed valid and right-handed.
    assert!(all_positively_oriented(&mesh));
}
