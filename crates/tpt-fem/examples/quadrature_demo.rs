//! Quadrature + reference-element sanity check via the umbrella prelude.
//!
//! Uses `tpt_fem::prelude` (which re-exports `tpt-fem-quadrature` and
//! `tpt-fem-element`) to verify two hand-derived facts about the reference
//! triangle: the quadrature weights sum to its area (`½`), and the Lagrange
//! shape functions satisfy partition of unity (`Σ Nᵢ = 1`). It also integrates
//! the coordinate `x` over the reference triangle, whose exact value is `1/6`.

use tpt_fem::prelude::*;

fn main() {
    let nodes = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let rule = triangle(TriangleRule::Degree2);

    // Weights sum to the area of the reference triangle (1/2).
    let wsum: f64 = rule.weights.iter().sum();
    println!("Σ weights = {} (reference-triangle area = 0.5)", wsum);

    // Partition of unity + integrate x over the reference triangle.
    let mut int_x = 0.0;
    for (p, &w) in rule.points.iter().zip(&rule.weights) {
        let n = Tri3::shape(p);
        let s: f64 = n.iter().sum();
        assert!((s - 1.0).abs() < 1e-12, "partition of unity = {s}");
        let x = n[0] * nodes[0][0] + n[1] * nodes[1][0] + n[2] * nodes[2][0];
        int_x += w * x;
    }

    println!("∫ x dA over reference triangle = {} (exact = 1/6)", int_x);
    assert!((wsum - 0.5).abs() < 1e-12, "weight sum = {wsum}");
    assert!((int_x - 1.0 / 6.0).abs() < 1e-12, "int_x = {int_x}");
    println!("OK: quadrature weights = area, shape functions partition of unity");
}
