use tpt_fem_assembly::{assemble, solve_with_dirichlet};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    // Two-element line on [0,1] with unit conductivity.
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0]);
    let n1 = b.add_node(vec![0.5]);
    let n2 = b.add_node(vec![1.0]);
    b.add_element(CellType::Line, vec![n0, n1]);
    b.add_element(CellType::Line, vec![n1, n2]);
    let mesh = b.build();

    // Element stiffness for unit conductivity, length 0.5: 2*[[1,-1],[-1,1]].
    let coo = assemble(&mesh, 1, |_, _| vec![vec![2.0, -2.0], vec![-2.0, 2.0]]);

    // Fix u(0) = 0, u(1) = 1; the P1 solution is linear, u(0.5) = 0.5.
    let u = solve_with_dirichlet(&coo, &[0.0, 0.0, 0.0], &[(n0, 0.0), (n2, 1.0)]).unwrap();
    assert!((u[n0] - 0.0).abs() < 1e-12);
    assert!((u[n1] - 0.5).abs() < 1e-12);
    assert!((u[n2] - 1.0).abs() < 1e-12);
    println!("Nodal solution u = {:?}", u);
}
