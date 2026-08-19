use tpt_fem_quadrature::{gauss_legendre, gauss_legendre_unit};

fn main() {
    // Order-3 Gauss-Legendre on [-1,1] is exact for degree 5.
    let r = gauss_legendre(3);
    assert!((r.weight_sum() - 2.0).abs() < 1e-12);

    // ∫_{-1}^{1} x^5 dx = 0 (odd integrand).
    let odd: f64 = r
        .weights
        .iter()
        .zip(&r.points)
        .map(|(w, x)| w * x.powi(5))
        .sum();
    assert!(odd.abs() < 1e-12);
    println!("∫ x^5 dx over [-1,1] ≈ {odd} (exact 0)");

    // ∫_{-1}^{1} x^4 dx = 2/5 (even, degree 4 < 5, exact).
    let even: f64 = r
        .weights
        .iter()
        .zip(&r.points)
        .map(|(w, x)| w * x.powi(4))
        .sum();
    assert!((even - 2.0 / 5.0).abs() < 1e-12);
    println!("∫ x^4 dx over [-1,1] ≈ {even} (exact {})", 2.0 / 5.0);

    // Unit interval [0,1], order 3 exact for degree 5: ∫_0^1 x^5 dx = 1/6.
    let ru = gauss_legendre_unit(3);
    assert!((ru.weight_sum() - 1.0).abs() < 1e-12);
    let pu: f64 = ru
        .weights
        .iter()
        .zip(&ru.points)
        .map(|(w, x)| w * x.powi(5))
        .sum();
    assert!((pu - 1.0 / 6.0).abs() < 1e-12);
    println!("∫ x^5 dx over [0,1] ≈ {pu} (exact {})", 1.0 / 6.0);
}
