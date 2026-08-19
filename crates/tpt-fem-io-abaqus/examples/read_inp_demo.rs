use tpt_fem_io_abaqus::read_inp;
use tpt_fem_mesh::CellType;

fn main() {
    let deck = "\
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 0.0, 1.0, 0.0
4, 1.0, 1.0, 0.0
*ELEMENT, TYPE=CPS4
1, 1, 2, 4, 3
";
    let mesh = read_inp(deck).unwrap();
    println!("node_count   = {}", mesh.node_count());
    println!("element_count = {}", mesh.element_count());
    println!("first node coords = {:?}", mesh.node_coords(0));

    assert_eq!(mesh.node_count(), 4);
    assert_eq!(mesh.element_count(), 1);
    assert_eq!(mesh.elements[0].cell_type, CellType::Quad);
    // Connectivity maps back to the same relative ordering.
    assert_eq!(mesh.elements[0].nodes, vec![0, 1, 3, 2]);
    assert_eq!(mesh.node_coords(0), &[0.0, 0.0, 0.0]);
}
