//! Free thermal expansion of a heated bar versus `dL = alpha * dT * L`.
//!
//! Run with: `cargo run -p tpt-fem-coupling --example thermal_bar_expansion`
//!
//! `thermal_structural` applies the thermal strain `eps_th = alpha * dT * I` as an
//! initial strain: the element load is `f_e = K_e u0` with
//! `u0(x) = alpha * dT * (x - x_centroid)`, which is exact for a uniform
//! temperature rise. An axial bar (`ElasticModel::BarAxial`, `Line` cells) that is
//! held only at one end is therefore free to expand, and the closed form is
//!
//! ```text
//! u(x) = alpha * dT * (x - x_fixed),      dL = alpha * dT * L
//! ```
//!
//! independent of the stiffness `EA`. The example prints the nodal displacement
//! against that closed form, sweeps `dT` and the bar length, and asserts the tip
//! extension to a relative tolerance of 1e-9.

use tpt_fem_coupling::thermal_structural;
use tpt_fem_elasticity::ElasticModel;
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
    let length = 2.0;
    let nel = 8;
    let ea = 2.1e5; // BarAxial reads `young` as E*A
    let alpha = 1.2e-5; // 1/K, mild steel
    let dt = 80.0; // K temperature rise

    let mesh = bar(length, nel);
    let temp = vec![dt; mesh.node_count()];
    let fixed = [(0usize, 0.0)]; // node 0, x-component (dim = 1)

    let u = thermal_structural(
        &mesh,
        ElasticModel::BarAxial,
        ea,
        0.3, // Poisson's ratio is unused by BarAxial (D = [EA])
        alpha,
        &temp,
        &fixed,
    )
    .expect("thermal-structural solve");

    println!(
        "axial bar: L = {length}, {nel} Line cells, EA = {ea:.3e}, alpha = {alpha:.3e} 1/K, dT = {dt} K"
    );
    println!("free expansion is stress free, so u(x) = alpha * dT * x\n");
    println!("   node       x       u (FEM)      u (exact)      error");
    println!("   ----   -------   -----------   -----------   ---------");
    let mut worst = 0.0_f64;
    for n in 0..mesh.node_count() {
        let x = mesh.node_coords(n)[0];
        let exact = alpha * dt * x;
        println!(
            "   {n:>4}   {x:7.3}   {:11.3e}   {exact:11.3e}   {:9.2e}",
            u[n],
            (u[n] - exact).abs()
        );
        worst = worst.max((u[n] - exact).abs());
    }
    let tip = u[mesh.node_count() - 1];
    let exact_tip = alpha * dt * length;
    println!("\n   tip extension dL = {tip:.6e}  (exact {exact_tip:.6e})");
    println!(
        "   axial strain     = {:.6e}  (alpha*dT = {:.6e})",
        tip / length,
        alpha * dt
    );
    println!("   max nodal error  = {worst:.2e}");

    println!("\ntemperature sweep (L = {length}, dL = alpha * dT * L)");
    println!("    dT [K]      dL (FEM)      dL (exact)    rel error");
    println!("   --------   -----------   ------------   ----------");
    for dt_s in [10.0, 50.0, 100.0, 250.0] {
        let temp = vec![dt_s; mesh.node_count()];
        let us = thermal_structural(&mesh, ElasticModel::BarAxial, ea, 0.3, alpha, &temp, &fixed)
            .expect("thermal-structural solve");
        let got = us[mesh.node_count() - 1];
        let want = alpha * dt_s * length;
        println!(
            "   {dt_s:8.1}   {got:11.4e}   {want:12.4e}   {:10.2e}",
            (got - want).abs() / want
        );
        assert!(
            (got - want).abs() / want < 1e-9,
            "dT = {dt_s}: {got} vs {want}"
        );
    }

    println!("\nlength sweep (dT = {dt} K, stiffness independent)");
    println!("      L      cells      dL (FEM)      dL (exact)    rel error");
    println!("   -----   -------   -----------   ------------   ----------");
    for (l, ne) in [(0.5, 2usize), (1.0, 4), (5.0, 10)] {
        let m = bar(l, ne);
        let temp = vec![dt; m.node_count()];
        let us = thermal_structural(&m, ElasticModel::BarAxial, ea, 0.3, alpha, &temp, &fixed)
            .expect("thermal-structural solve");
        let got = us[m.node_count() - 1];
        let want = alpha * dt * l;
        println!(
            "   {l:5.1}   {ne:>7}   {got:11.4e}   {want:12.4e}   {:10.2e}",
            (got - want).abs() / want
        );
        assert!((got - want).abs() / want < 1e-9, "L = {l}: {got} vs {want}");
    }

    assert!(
        (tip - exact_tip).abs() / exact_tip < 1e-9,
        "tip extension {tip} != {exact_tip}"
    );
    println!("\nverified: every free-expansion case matches alpha*dT*L to < 1e-9 relative");
}
