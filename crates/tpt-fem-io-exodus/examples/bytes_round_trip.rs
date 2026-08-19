use tpt_fem_io_exodus::{bytes_to_mesh, mesh_to_exodus_bytes};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    // In-memory encode/decode round-trip (no file needed).
    let bytes = mesh_to_exodus_bytes(&mesh).unwrap();
    assert_eq!(&bytes[0..4], b"CDF\x01");

    let parsed = bytes_to_mesh(&bytes).unwrap();
    println!("decoded node_count   = {}", parsed.node_count());
    println!("decoded element_count = {}", parsed.element_count());
    assert_eq!(parsed.node_count(), 3);
    assert_eq!(parsed.element_count(), 1);
    assert_eq!(parsed.elements[0].cell_type, CellType::Tri);
    assert!((parsed.node_coords(1)[0] - 1.0).abs() < 1e-6);
}
