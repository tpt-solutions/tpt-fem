//! Force-controlled elastic-plastic rod solved by `solve_elastic_plastic_rod`.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem-plasticity --example rod
//! ```
//!
//! A steel rod (length 1 m, area 2 cm²) discretised into `Line2` elements is
//! loaded through a sequence of increasing axial forces, with the plastic state
//! carried between load steps. The end displacement is checked against the
//! closed-form uniform-bar response at every step:
//!
//! ```text
//! sigma = P/A,   eps_p_bar = max(0, (sigma - sigma_y)/H),   u = (sigma/E + eps_p_bar)*L
//! ```
//!
//! Note the material *must* harden for force control to pass the yield load: a
//! perfectly plastic rod has zero tangent stiffness there, so the yield load is a
//! genuine limit point with no bounded solution beyond it.

use tpt_fem_mesh::{CellType, MeshBuilder};
use tpt_fem_plasticity::{solve_elastic_plastic_rod, PlasticityParams};

/// A rod of `n` `Line2` elements spanning `[0, length]`, one axial DOF per node.
fn rod(n: usize, length: f64) -> tpt_fem_mesh::Mesh {
    let mut b = MeshBuilder::new();
    let mut prev = b.add_node(vec![0.0]);
    for i in 1..=n {
        let node = b.add_node(vec![length * i as f64 / n as f64]);
        b.add_element(CellType::Line, vec![prev, node]);
        prev = node;
    }
    b.build()
}

fn main() {
    let length = 1.0;
    let area = 2e-4; // 2 cm^2
    let h = 4e9; // isotropic hardening modulus, E/50
    let p = PlasticityParams {
        iso_hardening: h,
        ..PlasticityParams::steel()
    };

    let mesh = rod(4, length);
    let p_y = p.yield_stress * area;
    let u_y = (p.yield_stress / p.young) * length;

    println!("Force-controlled elastic-plastic rod");
    println!(
        "  L = {length} m, A = {area:.1e} m^2, {} elements",
        mesh.element_count()
    );
    println!(
        "  E = {:.3e} Pa, sigma_y = {:.3e} Pa, H = {:.3e} Pa",
        p.young, p.yield_stress, h
    );
    println!("  yield load P_y = {p_y:.6e} N, yield displacement u_y = {u_y:.6e} m\n");

    // Ramp from well below yield to 1.75x yield, in one continuous history.
    let loads: Vec<f64> = (1..=7).map(|k| 0.25 * k as f64 * p_y).collect();
    let history =
        solve_elastic_plastic_rod(&mesh, area, &p, &loads).expect("rod solve must converge");

    println!(
        "  {:<10} {:<14} {:<16} {:<16} {:<8}",
        "P/P_y", "sigma (Pa)", "u_end (m)", "u_closed (m)", "regime"
    );
    println!("  {}", "-".repeat(68));

    let mut prev_u = 0.0;
    for (load, u) in loads.iter().zip(&history) {
        let u_end = *u.last().unwrap();

        // Closed form for a uniform bar carrying sigma = P/A.
        let sigma = load / area;
        let epsp = ((sigma - p.yield_stress) / h).max(0.0);
        let u_closed = (sigma / p.young + epsp) * length;
        let regime = if epsp > 0.0 { "plastic" } else { "elastic" };

        println!(
            "  {:<10.2} {:<14.6e} {:<16.6e} {:<16.6e} {:<8}",
            load / p_y,
            sigma,
            u_end,
            u_closed,
            regime
        );

        assert!(
            (u_end - u_closed).abs() / u_closed < 1e-8,
            "end displacement {u_end:.6e} vs closed form {u_closed:.6e}"
        );
        // Loading is monotone, so the displacement must be too.
        assert!(u_end > prev_u, "displacement must increase with load");
        prev_u = u_end;

        // The displacement field of a uniform bar is linear in x.
        for (n, &un) in u.iter().enumerate() {
            let x = mesh.node_coords(n)[0];
            let want = u_end * x / length;
            assert!(
                (un - want).abs() <= 1e-8 * u_end.max(1e-30),
                "node {n} at x={x}: u={un:.6e}, expected linear profile {want:.6e}"
            );
        }
    }

    let stiff_ratio = {
        let ep = p.young * h / (p.young + h);
        ep / p.young
    };
    println!(
        "\n  post-yield tangent E_p/E = {stiff_ratio:.4} (E_p = E*H/(E+H) = {:.3e} Pa)",
        p.young * h / (p.young + h)
    );
    println!("\nOK: every load step matches the closed-form uniform-bar response and the");
    println!("    displacement stays linear along the rod.");
}
