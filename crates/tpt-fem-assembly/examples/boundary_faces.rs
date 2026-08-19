use tpt_fem_assembly::boundary_faces;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    // Unit square split by the diagonal into two triangles.
    let mut b = MeshBuilder::new();
    let n00 = b.add_node(vec![0.0, 0.0]);
    let n10 = b.add_node(vec![1.0, 0.0]);
    let n01 = b.add_node(vec![0.0, 1.0]);
    let n11 = b.add_node(vec![1.0, 1.0]);
    b.add_element(CellType::Tri, vec![n00, n10, n11]);
    b.add_element(CellType::Tri, vec![n00, n11, n01]);
    let mesh = b.build();

    // The square has 4 edges, each belonging to exactly one triangle.
    let bf = boundary_faces(&mesh);
    assert_eq!(bf.len(), 4);
    println!("Boundary faces found: {}", bf.len());
}
