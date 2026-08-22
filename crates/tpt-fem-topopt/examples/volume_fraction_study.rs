//! Volume-fraction study: how much stiffness survives as material is removed?
//! Runs the SIMP optimizer on the same cantilever at several volume fractions
//! and reports the compliance penalty of each budget.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-topopt --example volume_fraction_study
//! ```

use tpt_fem_topopt::{cantilever_load, topopt_simp, Grid, TopOptParams};

fn main() {
    let grid = Grid::new(21, 11, 1.0);
    let (f, bcs) = cantilever_load(&grid, 1.0);

    println!("Cantilever compliance vs. allowed volume ({} elements):\n", grid.n_elem());
    println!("  {:>8}  {:>14}  {:>10}", "vol_frac", "compliance", "vs. solid");

    let mut baseline = None;
    for &vf in &[0.9_f64, 0.7, 0.5, 0.35, 0.25] {
        let params = TopOptParams {
            grid: grid.clone(),
            e0: 1.0,
            nu: 0.3,
            vol_frac: vf,
            penal: 3.0,
            filter_radius: 2.0,
            max_iter: 20,
            move_limit: 0.2,
        };
        let res = topopt_simp(&params, &f, &bcs).unwrap();

        // Volume constraint honoured exactly.
        let used: f64 = res.densities.iter().sum();
        assert!((used - vf * grid.n_elem() as f64).abs() < 1e-6);

        let c = res.compliance.last().unwrap();
        if baseline.is_none() {
            baseline = Some(*c);
        }
        let rel = c / baseline.unwrap();
        println!("  {vf:>8.2}  {c:>14.6e}  {rel:>9.2}x");

        // Removing material can never lower the optimal compliance.
        assert!(*c >= baseline.unwrap() * (1.0 - 1e-9));
    }

    // Monotonicity sanity check across the whole run: compliance histories are
    // non-increasing under the optimality-criteria update.
    let params = TopOptParams {
        grid: grid.clone(),
        e0: 1.0,
        nu: 0.3,
        vol_frac: 0.5,
        penal: 3.0,
        filter_radius: 2.0,
        max_iter: 20,
        move_limit: 0.2,
    };
    let res = topopt_simp(&params, &f, &bcs).unwrap();
    for w in res.compliance.windows(2) {
        assert!(w[1] <= w[0] + 1e-9, "compliance increased at an iteration");
    }

    println!("\nStudy complete.");
}
