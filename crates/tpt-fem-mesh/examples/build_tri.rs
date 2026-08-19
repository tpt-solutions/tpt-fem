use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    println!("node_count   = {}", mesh.node_count());
    println!("element_count = {}", mesh.element_count());
    println!("node 1 coords = {:?}", mesh.node_coords(n1));

    assert_eq!(mesh.node_count(), 3);
    assert_eq!(mesh.element_count(), 1);
    assert_eq!(mesh.elements[0].cell_type, CellType::Tri);
    assert_eq!(mesh.node_coords(n1), &[1.0, 0.0]);
}
