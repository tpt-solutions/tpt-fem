use tpt_fem_quadrature::{tetrahedron, TetrahedronRule};

fn main() {
    // Reference tetrahedron has volume 1/6.
    let t1 = tetrahedron(TetrahedronRule::Degree1);
    assert!((t1.weight_sum() - 1.0 / 6.0).abs() < 1e-12);
    // ∫_T x dV = 1!0!0!/(4!) = 1/24 (linear moment, exact for degree-1 rule).
    let xint: f64 = t1
        .weights
        .iter()
        .zip(&t1.points)
        .map(|(w, p)| w * p[0])
        .sum();
    assert!((xint - 1.0 / 24.0).abs() < 1e-12);
    println!(
        "∫ x over reference tetrahedron ≈ {xint} (exact {})",
        1.0 / 24.0
    );

    // Degree2 rule is exact for quadratics: ∫_T x^2 dV = 2!/(5!) = 1/60.
    let t2 = tetrahedron(TetrahedronRule::Degree2);
    assert!((t2.weight_sum() - 1.0 / 6.0).abs() < 1e-12);
    let x2: f64 = t2
        .weights
        .iter()
        .zip(&t2.points)
        .map(|(w, p)| w * p[0].powi(2))
        .sum();
    assert!((x2 - 1.0 / 60.0).abs() < 1e-12);
    println!(
        "∫ x^2 over reference tetrahedron ≈ {x2} (exact {})",
        1.0 / 60.0
    );
}
