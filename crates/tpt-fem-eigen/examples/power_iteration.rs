use tpt_fem_eigen::{matvec, power_iteration};
use tpt_fem_sparse::Coo;

fn main() {
    // Symmetric matrix [[2, 1], [1, 2]] has eigenvalues 1 and 3.
    let mut c = Coo::new();
    c.push(0, 0, 2.0);
    c.push(0, 1, 1.0);
    c.push(1, 0, 1.0);
    c.push(1, 1, 2.0);

    let (lam, v) = power_iteration(&c, 200, 1e-12);
    assert!((lam - 3.0).abs() < 1e-8, "got {lam}");

    // Eigenvector residual |A v - λ v| should be ~ 0.
    let av = matvec(&c, &v);
    let res = av
        .iter()
        .zip(&v)
        .map(|(a, b)| a - lam * b)
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    assert!(res < 1e-6, "eigenvector residual {res}");

    println!("Dominant eigenvalue λ = {lam}");
}
