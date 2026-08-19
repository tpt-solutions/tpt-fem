use tpt_fem_eigen::{inverse_iteration, matvec};
use tpt_fem_sparse::Coo;

fn main() {
    // Symmetric matrix [[2, 1], [1, 2]] has eigenvalues 1 and 3.
    let mut c = Coo::new();
    c.push(0, 0, 2.0);
    c.push(0, 1, 1.0);
    c.push(1, 0, 1.0);
    c.push(1, 1, 2.0);

    // Shift 0 -> eigenpair nearest 0, i.e. the smallest eigenvalue 1.
    let (lam, v) = inverse_iteration(&c, 0.0, 200, 1e-12).unwrap();
    assert!((lam - 1.0).abs() < 1e-8, "got {lam}");

    let av = matvec(&c, &v);
    let res = av
        .iter()
        .zip(&v)
        .map(|(a, b)| a - lam * b)
        .map(|x| x * x)
        .sum::<f64>()
        .sqrt();
    assert!(res < 1e-6, "eigenvector residual {res}");

    println!("Eigenvalue nearest 0: λ = {lam}");
}
