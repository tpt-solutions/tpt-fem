//! Single axial-bar element stiffness check.
//!
//! For a 1-D bar `[0, 1]` with axial stiffness `EA = 1` the element stiffness is
//! known in closed form: `K = (EA/L) [[1, -1], [-1, 1]] = [[1, -1], [-1, 1]]`.

use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0]);
    let n1 = b.add_node(vec![1.0]);
    b.add_element(CellType::Line, vec![n0, n1]);
    let mesh = b.build();

    let k = elasticity_element_matrix(&mesh, 0, ElasticModel::BarAxial, 1.0, 0.0, 1)
        .expect("bar element matrix");

    println!("Bar (EA=1, L=1) element stiffness:");
    for row in &k {
        println!("  {:8.4} {:8.4}", row[0], row[1]);
    }

    assert!((k[0][0] - 1.0).abs() < 1e-12, "k00 = {}", k[0][0]);
    assert!((k[1][1] - 1.0).abs() < 1e-12, "k11 = {}", k[1][1]);
    assert!((k[0][1] + 1.0).abs() < 1e-12, "k01 = {}", k[0][1]);
    assert!((k[1][0] + 1.0).abs() < 1e-12, "k10 = {}", k[1][0]);
    println!("OK: matches analytic EA/L [[1,-1],[-1,1]]");
}
