use tpt_fem_mesh::MeshBuilder;

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0, 0.0]);
    let n3 = b.add_node(vec![1.0, 1.0, 1.0]);
    let mesh = b.build();

    let mut on_x0 = mesh.nodes_on_plane(0, 0.0, 1e-9);
    on_x0.sort_unstable();
    println!("nodes on plane x = 0: {:?}", on_x0);
    assert_eq!(on_x0, vec![n0, n2]);

    let mut in_box = mesh.nodes_in_box([0.0, 0.0, 0.0], [1.0, 1.0, 0.5]);
    in_box.sort_unstable();
    println!("nodes in box z <= 0.5: {:?}", in_box);
    assert_eq!(in_box, vec![n0, n1, n2]);
    assert!(!in_box.contains(&n3));
}
