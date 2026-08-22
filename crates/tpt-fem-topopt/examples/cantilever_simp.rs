//! Classic SIMP topology optimization of a tip-loaded cantilever plate.
//! Optimizes a small grid, prints the compliance history, and renders the
//! final density design as ASCII art.
//!
//! Run with:
//!
//! ```sh
//! cargo run -p tpt-fem-topopt --example cantilever_simp
//! ```

use tpt_fem_topopt::{cantilever_load, topopt_simp, Grid, TopOptParams};

fn main() {
    // 26 × 13 grid of unit squares, clamped left edge, downward tip load.
    let grid = Grid::new(26, 13, 1.0);
    let (f, bcs) = cantilever_load(&grid, 1.0);

    let params = TopOptParams {
        grid: grid.clone(),
        e0: 1.0,           // solid Young's modulus
        nu: 0.3,           // Poisson's ratio
        vol_frac: 0.5,     // keep 50% of the material
        penal: 3.0,        // SIMP penalty
        filter_radius: 2.4,// sensitivity filter radius (in elements)
        max_iter: 30,
        move_limit: 0.2,
    };

    let res = topopt_simp(&params, &f, &bcs).unwrap();

    println!("SIMP cantilever, {} elements, target volume fraction {}",
        grid.n_elem(), params.vol_frac);
    println!("\nCompliance history:");
    for (it, c) in res.compliance.iter().enumerate() {
        println!("  iter {:>2}: c = {:.6e}", it + 1, c);
    }
    assert!(res.compliance.last().unwrap() < &res.compliance[0]);
    let used: f64 = res.densities.iter().sum();
    assert!((used - params.vol_frac * grid.n_elem() as f64).abs() < 1e-6);

    // ASCII rendering of the optimized densities (top row printed first).
    println!("\nOptimized design ('#' = solid, '.' = void):");
    let nx = grid.nx - 1;
    let ny = grid.ny - 1;
    for ej in (0..ny).rev() {
        let row: String = (0..nx)
            .map(|ei| {
                if res.densities[ej * nx + ei] > 0.5 { '#' } else { '.' }
            })
            .collect();
        println!("  {row}");
    }

    println!("\nFinal compliance: {:.6e}", res.compliance.last().unwrap());
}
