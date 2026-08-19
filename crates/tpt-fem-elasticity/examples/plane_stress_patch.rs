//! Plane-stress patch test (linear-displacement reproduction).
//!
//! A unit square split into four triangles around a centre node is loaded with
//! the linear displacement field `u = (0.1 x, 0.2 y)` on the four corners. The
//! free centre node must reproduce `u = (0.05, 0.1)` exactly, exercising the
//! full `solve_elasticity` assemble-and-solve path.

use tpt_fem_elasticity::{solve_elasticity, ElasticModel};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let c00 = b.add_node(vec![0.0, 0.0]);
    let c10 = b.add_node(vec![1.0, 0.0]);
    let c01 = b.add_node(vec![0.0, 1.0]);
    let c11 = b.add_node(vec![1.0, 1.0]);
    let mid = b.add_node(vec![0.5, 0.5]);
    b.add_element(CellType::Tri, vec![c00, c10, mid]);
    b.add_element(CellType::Tri, vec![c10, c11, mid]);
    b.add_element(CellType::Tri, vec![c11, c01, mid]);
    b.add_element(CellType::Tri, vec![c01, c00, mid]);
    let mesh = b.build();

    let mut bcs = Vec::new();
    for (node, x, y) in [
        (c00, 0.0, 0.0),
        (c10, 0.1, 0.0),
        (c01, 0.0, 0.2),
        (c11, 0.1, 0.2),
    ] {
        bcs.push((node * 2, x));
        bcs.push((node * 2 + 1, y));
    }

    let u = solve_elasticity(
        &mesh,
        ElasticModel::PlaneStress,
        1.0,
        0.3,
        2,
        |_| vec![0.0, 0.0],
        &bcs,
    )
    .expect("plane-stress solve");

    println!(
        "Centre node displacement: u = ({}, {})",
        u[mid * 2],
        u[mid * 2 + 1]
    );
    assert!((u[mid * 2] - 0.05).abs() < 1e-9, "got {}", u[mid * 2]);
    assert!(
        (u[mid * 2 + 1] - 0.1).abs() < 1e-9,
        "got {}",
        u[mid * 2 + 1]
    );
    println!("OK: centre node reproduces the linear field exactly");
}
