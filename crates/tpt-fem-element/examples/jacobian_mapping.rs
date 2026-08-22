//! Demonstrates the isoparametric [`Map`]: building the Jacobian from physical
//! node coordinates, checking its determinant (the integration measure), and
//! mapping reference gradients to physical space.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-element --example jacobian_mapping
//! ```

use tpt_fem_element::{quad_rule, Map, Quad4, ReferenceElement};

fn main() {
    // Map the reference square [-1,1]² onto a sheared parallelogram spanned by
    // the edge vectors a = (4, 0) and b = (1, 2) anchored at (1, 1):
    //
    //   x(xi) = (1,1) + (xi+1)/2 · a + (eta+1)/2 · b
    //
    // The map is affine, so J is constant and det J = |a × b| = 8 everywhere;
    // a bilinear Quad4 reproduces it exactly.
    let a = [4.0_f64, 0.0];
    let b = [1.0_f64, 2.0];
    let origin = [1.0_f64, 1.0];
    let physical: Vec<Vec<f64>> = Quad4::nodes()
        .iter()
        .map(|n| {
            [
                origin[0] + (n[0] + 1.0) / 2.0 * a[0] + (n[1] + 1.0) / 2.0 * b[0],
                origin[1] + (n[0] + 1.0) / 2.0 * a[1] + (n[1] + 1.0) / 2.0 * b[1],
            ]
            .to_vec()
        })
        .collect();

    println!("Physical Quad4 nodes:");
    for (i, x) in physical.iter().enumerate() {
        println!("  X_{i} = ({:.3}, {:.3})", x[0], x[1]);
    }

    // Build the Jacobian at the reference centre.
    let centre = [0.0_f64, 0.0];
    let grad = Quad4::grad(&centre);
    let m = Map::from_nodes_and_grad(&physical, &grad);

    println!("\nMap at the reference centre:");
    println!(
        "  J     = [[{:.4}, {:.4}], [{:.4}, {:.4}]]",
        m.jacobian[0], m.jacobian[1], m.jacobian[2], m.jacobian[3]
    );
    println!("  det J = {:.6}", m.determinant);
    // dx/dxi = a/2 = (2, 0) and dx/deta = b/2 = (1/2, 1), so
    // det J = 2·1 − 0·(1/2) = 2 everywhere on this affine map.
    assert!(
        (m.determinant - 2.0).abs() < 1e-12,
        "affine map must have constant det J = 2"
    );

    // Map a reference gradient dN/dxi = (1, 0) into physical coordinates:
    // dN/dx = (dN/dxi) · J⁻¹.
    let local = [1.0_f64, 0.0];
    let physical_grad = m.physical_grad(&local);
    println!(
        "\ndN/dxi = (1,0)  ->  dN/dx = ({:.6}, {:.6})",
        physical_grad[0], physical_grad[1]
    );

    // Integration measure in action: ∫_Ω 1 dΩ over the physical element equals
    // Σ_q w_q · |det J(q)|, i.e. the parallelogram's area |a × b| = 8.
    let rule = quad_rule(2);
    let mut area = 0.0;
    for (p, w) in rule.points.iter().zip(&rule.weights) {
        let g = Quad4::grad(p);
        area += w * Map::from_nodes_and_grad(&physical, &g).determinant;
    }
    println!("\nArea by 2×2 Gauss integration = {area:.12} (expected 8.0)");
    assert!((area - 8.0).abs() < 1e-12);

    println!("\nAll Jacobian checks passed.");
}
