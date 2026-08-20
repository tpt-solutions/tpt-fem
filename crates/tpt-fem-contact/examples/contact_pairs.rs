//! Nearest-node proximity search between two opposing surfaces using the
//! brute-force `contact_pairs` helper.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-contact --example contact_pairs
//! ```

use tpt_fem_contact::contact_pairs;

fn main() {
    // Surface A: nodes along x = 0, y = 0,1,2,...
    // Surface B: nodes along x = 0.25, y = 0,1,2,... (parallel, offset gap).
    let a: Vec<(usize, Vec<f64>)> = (0..5).map(|i| (i, vec![0.0, i as f64, 0.0])).collect();
    let b: Vec<(usize, Vec<f64>)> = (0..5)
        .map(|i| (100 + i, vec![0.25, i as f64, 0.0]))
        .collect();

    let pairs = contact_pairs(&a, &b);

    println!("Proximity search: surface A (x=0) vs surface B (x=0.25)");
    println!("  {:>6} {:>6} {:>10}", "node A", "node B", "gap");
    for (na, best) in &pairs {
        let (ib, gap) = best.expect("surface B is non-empty");
        let nb = b[ib].0;
        println!("  {:6} {:6} {:10.4}", na, nb, gap);
        assert!((gap - 0.25).abs() < 1e-12, "gap must equal the x-offset");
    }
    println!();
    println!("OK: every node of A is paired with its parallel neighbour on B.");
}
