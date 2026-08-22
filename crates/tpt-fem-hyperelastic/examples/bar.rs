//! Large-stretch Newton solve of a 1-D incompressible neo-Hookean bar.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-hyperelastic --example bar
//! ```
//!
//! `solve_hyperelastic_bar` drives a 4-element bar to a target end stretch with
//! `tpt-fem-solve`'s Newton loop. The free-end displacement and per-element
//! stretch are printed and checked against the closed form.

use tpt_fem_hyperelastic::solve_hyperelastic_bar;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let mut prev = b.add_node(vec![0.0]);
    for i in 1..=4 {
        let node = b.add_node(vec![1.0 * i as f64 / 4.0]);
        b.add_element(CellType::Line, vec![prev, node]);
        prev = node;
    }
    let mesh = b.build();

    let mu = 100.0;
    let area = 1.0;
    let target = 1.5;
    let u = solve_hyperelastic_bar(&mesh, area, mu, target)
        .expect("hyperelastic bar solve must converge");

    let total_len = 1.0; // four elements of length 0.25
    let target_disp = total_len * (target - 1.0);
    println!("4-element bar, μ = {mu}, target stretch λ = {target}");
    println!("target end displacement = {:.4} m", target_disp);
    println!(
        "node displacements: {:?}",
        u.iter().map(|x| format!("{:.4}", x)).collect::<Vec<_>>()
    );

    // Rightmost node (index 4) should reach the target displacement.
    assert!(
        (u[4] - target_disp).abs() < 1e-6,
        "end displacement {}",
        u[4]
    );
    // Linear stretch -> uniform displacement, so node 2 sits at half the end disp.
    assert!((u[2] - target_disp / 2.0).abs() < 1e-6);
    // Per-element stretch should equal the target (uniform bar).
    for e in 0..4 {
        let l0 = 0.25;
        let l = l0 + (u[e + 1] - u[e]);
        let stretch = l / l0;
        assert!(
            (stretch - target).abs() < 1e-6,
            "element {e} stretch {stretch}"
        );
    }
    println!(
        "\nOK: Newton converges; end displacement and per-element stretch match λ = {target}."
    );
}
