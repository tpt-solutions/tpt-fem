use tpt_fem_io_exodus::{read_exodus, write_exodus};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    let path = std::env::temp_dir().join("tpt_fem_io_exodus_demo.ex2");
    write_exodus(&mesh, &path).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    println!("wrote {} bytes", bytes.len());
    // NetCDF-3 classic / CDF-1 magic number.
    assert_eq!(&bytes[0..4], b"CDF\x01");

    let imported = read_exodus(&path).unwrap();
    println!("node_count   = {}", imported.node_count());
    println!("element_count = {}", imported.element_count());
    assert_eq!(imported.node_count(), 3);
    assert_eq!(imported.element_count(), 1);
    assert_eq!(imported.elements[0].cell_type, CellType::Tri);
    assert_eq!(imported.elements[0].nodes, vec![0, 1, 2]);
    let _ = std::fs::remove_file(&path);
}
