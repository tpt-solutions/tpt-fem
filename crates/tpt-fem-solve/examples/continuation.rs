use tpt_fem_solve::{continuation, NewtonOptions};
use tpt_fem_sparse::Coo;

fn main() {
    // Trace u^2 = λ from λ = 1 to λ = 4, warm-starting each Newton step.
    let res = continuation(
        &[1.0],
        1.0,
        4.0,
        8,
        &[],
        |u, lam| vec![u[0] * u[0] - lam],
        |u, _lam| {
            let mut c = Coo::new();
            c.push(0, 0, 2.0 * u[0]);
            c
        },
        &NewtonOptions::default(),
    )
    .unwrap();

    let (lam, u) = res.last().unwrap();
    assert!(((*lam) - 4.0).abs() < 1e-12);
    assert!((u[0] - 2.0).abs() < 1e-8, "got {}", u[0]);
    println!("At λ = {lam} the continuation solved u = {}", u[0]);
}
