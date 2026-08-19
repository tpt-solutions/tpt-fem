use tpt_fem_solve::{arc_length_continuation, ArcLengthOptions};
use tpt_fem_sparse::Coo;

fn main() {
    // Cubic fold: R(u, λ) = u^3 - 3u - λ. The equilibrium λ = u^3 - 3u has a
    // limit point at u = 1 (λ = -2). Start on the right branch and trace
    // downward; the arc-length method must pass the fold.
    let res = |u: &[f64], lam: f64| vec![u[0] * u[0] * u[0] - 3.0 * u[0] - lam];
    let jac = |u: &[f64], _lam: f64| {
        let mut c = Coo::new();
        c.push(0, 0, 3.0 * u[0] * u[0] - 3.0);
        c
    };
    let opts = ArcLengthOptions {
        initial_arc_length: 0.15,
        target_arc_length: 0.15,
        max_arc_length: 0.5,
        max_steps: 800,
        initial_direction: -1.0,
        ..Default::default()
    };

    let trace = arc_length_continuation(&[2.1038], 3.0, &[], &[1.0], res, jac, &opts).unwrap();
    for (lam, u) in trace.iter().skip(1) {
        let r = u[0] * u[0] * u[0] - 3.0 * u[0] - lam;
        assert!(r.abs() < 1e-5, "off-curve: {r} at λ={lam}");
    }
    let min_lam = trace.iter().map(|(l, _)| *l).fold(f64::INFINITY, f64::min);
    assert!(
        min_lam <= -1.7,
        "did not reach the fold region, min λ = {min_lam}"
    );
    println!("Arc-length trace reached minimum λ = {min_lam} (fold at -2)");
}
