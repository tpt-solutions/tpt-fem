use tpt_fem_solve::{newton, NewtonOptions};
use tpt_fem_sparse::Coo;

fn main() {
    // Find the real root of u^3 - u - 1 = 0, starting from u = 1.
    let u = newton(
        &[1.0],
        |u| vec![u[0] * u[0] * u[0] - u[0] - 1.0],
        |u| {
            let mut c = Coo::new();
            c.push(0, 0, 3.0 * u[0] * u[0] - 1.0);
            c
        },
        &[],
        &NewtonOptions::default(),
    )
    .unwrap();

    // The unique real root is the plastic number ≈ 1.324717957.
    assert!((u[0] - 1.324717957).abs() < 1e-8);
    println!("Root of u^3 - u - 1 = 0: u = {}", u[0]);
}
