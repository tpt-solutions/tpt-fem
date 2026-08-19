//! Axial-bar element stiffness via the umbrella re-exports.
//!
//! Builds a single `Line2` bar through `tpt_fem::prelude` (re-exporting
//! `tpt-fem-mesh` and `tpt-fem-elasticity`) and checks the closed-form element
//! stiffness `K = (EA/L) [[1, -1], [-1, 1]]` for `EA = 1`, `L = 1`.

use tpt_fem::prelude::*;

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
