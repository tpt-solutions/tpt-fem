//! Joule heating: `joule_source` plus a steady `electro_thermal` conduction solve.
//!
//! Run with: `cargo run -p tpt-fem-coupling --example joule_heating`
//!
//! `joule_source(sigma, |E|)` returns the volumetric Ohmic dissipation
//! `q = sigma |E|^2` (W/m^3), and `electro_thermal` feeds that constant source
//! into the steady heat-conduction (Poisson) operator of `tpt-fem-thermal`:
//!
//! ```text
//! -k grad^2 T = q = sigma |E|^2
//! ```
//!
//! For a bar of length `L` held at `T = 0` at both ends the closed form is
//!
//! ```text
//! T(x) = q x (L - x) / (2k),      T_max = q L^2 / (8k)   at x = L/2
//! ```
//!
//! The example prints the nodal temperature rise against that parabola, sweeps
//! the field strength to show the quadratic `T ~ |E|^2` scaling, and asserts both
//! the profile and the scaling.

use tpt_fem_coupling::{electro_thermal, joule_source};
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Uniform 1-D bar of length `length` with `nel` `Line` cells.
fn bar(length: f64, nel: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut prev = b.add_node(vec![0.0]);
    for i in 1..=nel {
        let n = b.add_node(vec![length * i as f64 / nel as f64]);
        b.add_element(CellType::Line, vec![prev, n]);
        prev = n;
    }
    b.build()
}

fn main() {
    let length = 0.02; // 20 mm copper bar
    let nel = 8;
    let k = 400.0; // W/(m K), copper thermal conductivity
    let sigma = 5.8e7; // S/m, copper electrical conductivity
    let e_mag = 0.05; // V/m field magnitude

    let q = joule_source(sigma, e_mag);
    println!("joule_source(sigma = {sigma:.3e} S/m, |E| = {e_mag} V/m) = {q:.4e} W/m^3");
    println!("bar: L = {length} m, {nel} Line cells, k = {k} W/(m K), T = 0 at both ends");
    println!("closed form T(x) = q x (L - x) / (2k), T_max = q L^2 / (8k)\n");

    let mesh = bar(length, nel);
    let last = mesh.node_count() - 1;
    let bc = [(0usize, 0.0), (last, 0.0)];
    let t = electro_thermal(&mesh, k, sigma, e_mag, &bc).expect("electro-thermal solve");

    println!("   node       x [mm]      T (FEM)      T (exact)      error");
    println!("   ----   ----------   -----------   -----------   ---------");
    let mut worst = 0.0_f64;
    for n in 0..mesh.node_count() {
        let x = mesh.node_coords(n)[0];
        let exact = q * x * (length - x) / (2.0 * k);
        println!(
            "   {n:>4}   {:10.3}   {:11.4}   {exact:11.4}   {:9.2e}",
            x * 1000.0,
            t[n],
            (t[n] - exact).abs()
        );
        worst = worst.max((t[n] - exact).abs());
    }
    let t_max_exact = q * length * length / (8.0 * k);
    let t_max = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("\n   peak temperature rise = {t_max:.4} K  (exact {t_max_exact:.4} K)");
    println!("   max nodal error       = {worst:.2e} K");

    println!("\nfield sweep: q = sigma |E|^2 and T_max = q L^2 / (8k), so T_max ~ |E|^2");
    println!("    |E| [V/m]        q [W/m^3]     T_max (FEM)   T_max (exact)   T/T_ref");
    println!("   ----------   --------------   -------------   -------------   -------");
    let mut t_ref = 0.0;
    for (i, e) in [0.025, 0.05, 0.1, 0.2].iter().enumerate() {
        let qe = joule_source(sigma, *e);
        let te = electro_thermal(&mesh, k, sigma, *e, &bc).expect("electro-thermal solve");
        let peak = te.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exact = qe * length * length / (8.0 * k);
        if i == 0 {
            t_ref = peak;
        }
        println!(
            "   {e:10.3}   {qe:14.4e}   {peak:13.4}   {exact:13.4}   {:7.2}",
            peak / t_ref
        );
        assert!(
            (peak - exact).abs() / exact < 1e-9,
            "|E| = {e}: {peak} vs {exact}"
        );
    }
    // Doubling |E| must quadruple the temperature rise.
    let t1 = electro_thermal(&mesh, k, sigma, 0.05, &bc).unwrap();
    let t2 = electro_thermal(&mesh, k, sigma, 0.10, &bc).unwrap();
    let ratio = t2.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        / t1.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("\n   T_max(2|E|) / T_max(|E|) = {ratio:.6}  (exact 4)");

    assert!(
        (t_max - t_max_exact).abs() / t_max_exact < 1e-9,
        "peak temperature {t_max} != {t_max_exact}"
    );
    assert!(
        (ratio - 4.0).abs() < 1e-9,
        "quadratic scaling broken: {ratio}"
    );
    println!("\nverified: nodal temperatures match q x (L-x) / (2k) to < 1e-9 relative");
    println!("and the rise scales as |E|^2");
}
