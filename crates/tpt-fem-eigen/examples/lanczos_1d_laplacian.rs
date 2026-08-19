use tpt_fem_eigen::{lanczos_eigs, EigWhich};
use tpt_fem_sparse::Coo;

fn main() {
    // Discrete 1-D Laplacian (n = 10) tridiagonal matrix.
    let n = 10;
    let mut c = Coo::new();
    for i in 0..n {
        c.push(i, i, 2.0);
        if i + 1 < n {
            c.push(i, i + 1, -1.0);
            c.push(i + 1, i, -1.0);
        }
    }

    // Closed-form smallest eigenvalue: 2 - 2 cos(π/(n+1)).
    let smallest = lanczos_eigs(&c, 1, EigWhich::Smallest, n);
    let expected = 2.0 - 2.0 * (std::f64::consts::PI / (n as f64 + 1.0)).cos();
    assert!(
        (smallest[0].0 - expected).abs() < 1e-6,
        "got {}",
        smallest[0].0
    );

    println!(
        "Smallest eigenvalue ≈ {} (exact {})",
        smallest[0].0, expected
    );
}
