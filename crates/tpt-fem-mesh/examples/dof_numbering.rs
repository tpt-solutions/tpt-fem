use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![n0, n1, n2]);
    let mesh = b.build();

    // One DOF per node: 3 nodes -> 3 DOFs, contiguously numbered.
    let dof1 = mesh.number_dofs(1);
    println!("ndof (1 dof/node) = {}", dof1.ndof);
    assert_eq!(dof1.ndof, 3);
    assert_eq!(dof1.dof(n2, 0), 2);

    // Three DOFs per node: 3 nodes -> 9 DOFs; node 1 owns 3..5.
    let dof3 = mesh.number_dofs(3);
    println!("ndof (3 dof/node) = {}", dof3.ndof);
    assert_eq!(dof3.ndof, 9);
    assert_eq!(dof3.dof(n1, 0), 3);
    assert_eq!(dof3.dof(n1, 2), 5);
}
