use tpt_fem_io_vtk::PointData;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    let path = std::env::temp_dir().join("tpt_fem_io_vtk_roundtrip.vtk");
    tpt_fem_io_vtk::write_vtk_with_data(&mesh, &[PointData::new("u", vec![0.0, 1.0, 1.0])], &path)
        .unwrap();

    // Read the mesh back and confirm the topology survives the round-trip.
    let imported = tpt_fem_io_vtk::read_vtk(&path).unwrap();
    println!("imported node_count   = {}", imported.node_count());
    println!("imported element_count = {}", imported.element_count());
    assert_eq!(imported.node_count(), 3);
    assert_eq!(imported.element_count(), 1);
    assert_eq!(imported.elements[0].cell_type, CellType::Tri);
    assert_eq!(imported.elements[0].nodes, vec![0, 1, 2]);
    let _ = std::fs::remove_file(&path);
}
