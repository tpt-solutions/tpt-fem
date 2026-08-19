use tpt_fem_quadrature::{gauss_legendre, tensor_square, triangle, TriangleRule};

fn main() {
    // Tensor-product square rule from an order-2 1-D rule; exact up to degree 3.
    let q = tensor_square(&gauss_legendre(2));
    assert!((q.weight_sum() - 4.0).abs() < 1e-12);
    // ∫_{[-1,1]^2} x^2 y^2 dx dy = (2/3)*(2/3) = 4/9.
    let val: f64 = q
        .weights
        .iter()
        .zip(&q.points)
        .map(|(w, p)| w * p[0].powi(2) * p[1].powi(2))
        .sum();
    assert!((val - 4.0 / 9.0).abs() < 1e-12);
    println!("∫ x^2 y^2 over square ≈ {val} (exact {})", 4.0 / 9.0);

    // Reference-triangle rule Degree2 is exact for quadratics.
    let t = triangle(TriangleRule::Degree2);
    assert!((t.weight_sum() - 0.5).abs() < 1e-12);
    // ∫_T x y dA = 1!1!/(4!) = 1/24 (closed form for the reference triangle).
    let tval: f64 = t
        .weights
        .iter()
        .zip(&t.points)
        .map(|(w, p)| w * p[0] * p[1])
        .sum();
    assert!((tval - 1.0 / 24.0).abs() < 1e-12);
    println!(
        "∫ x y over reference triangle ≈ {tval} (exact {})",
        1.0 / 24.0
    );
}
