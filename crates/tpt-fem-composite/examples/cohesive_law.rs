//! Sweeps a bilinear cohesive (traction-separation) law and verifies that the
//! area under the curve equals the prescribed fracture toughness Gc.
//!
//! Run with:
//! ```text
//! cargo run -p tpt-fem-composite --example cohesive_law
//! ```

use tpt_fem_composite::CohesiveLaw;

fn main() {
    // Bilinear law: peak traction 10 MPa, critical opening 0.1 mm, Gc = 1 N/mm.
    let gc = 1.0; // N/mm
    let peak = 10.0; // MPa = N/mm^2
    let dc = 0.1; // mm
    let law = CohesiveLaw::from_toughness(gc, peak, dc);

    println!("Bilinear cohesive law");
    println!("  peak traction   sigma*  = {:.3} MPa", law.peak_traction);
    println!("  critical open.  delta_c = {:.3} mm", law.critical_opening);
    println!("  final open.     delta_f = {:.3} mm", law.final_opening);
    println!();
    println!("  {:>8} {:>10}", "delta", "traction");

    let n = 20; // 0.2/20 = 0.01 so dc = 0.1 is an exact sample point
    let mut area = 0.0_f64;
    let mut prev_d = 0.0_f64;
    let mut prev_t = 0.0_f64;
    for i in 0..=n {
        let d = law.final_opening * i as f64 / n as f64;
        let t = law.traction(d);
        println!("  {:8.4} {:10.4}", d, t);
        // Trapezoidal integration of the area under the curve.
        area += 0.5 * (t + prev_t) * (d - prev_d);
        prev_d = d;
        prev_t = t;
    }

    println!();
    println!("integrated Gc  = {:.6} N/mm", area);
    println!("prescribed Gc  = {:.6} N/mm", law.toughness());
    assert!((area - gc).abs() / gc < 1e-3, "area must equal Gc");
    assert!((law.toughness() - gc).abs() < 1e-9);
    println!("OK: area under the traction-separation curve recovers Gc.");
}
