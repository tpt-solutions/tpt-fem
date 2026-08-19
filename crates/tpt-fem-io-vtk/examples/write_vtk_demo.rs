use tpt_fem_io_vtk::PointData;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    // Build an in-memory VTK model with a per-node field and inspect it.
    let vtk = tpt_fem_io_vtk::mesh_to_vtk(&mesh, &[PointData::new("u", vec![0.0, 0.5, 1.0])]);
    let vtkio::model::DataSet::UnstructuredGrid { pieces, .. } = &vtk.data else {
        panic!("expected an unstructured grid");
    };
    let vtkio::model::Piece::Inline(p) = &pieces[0] else {
        panic!("expected an inline piece");
    };
    println!("vtk points = {}", p.num_points());
    println!("vtk cells  = {}", p.cells.num_cells());
    assert_eq!(p.num_points(), 3);
    assert_eq!(p.cells.num_cells(), 1);

    // Write an ASCII .vtk file and confirm it carries the expected markers.
    let path = std::env::temp_dir().join("tpt_fem_io_vtk_demo.vtk");
    tpt_fem_io_vtk::write_vtk_ascii(&mesh, &path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    println!("wrote {} bytes to {}", text.len(), path.display());
    assert!(text.contains("POINTS"));
    assert!(text.contains("CELLS"));
    assert!(text.contains("CELL_TYPES"));
    let _ = std::fs::remove_file(&path);
}
