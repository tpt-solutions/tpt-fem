use tpt_fem_assembly::apply_neumann;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    // Single triangle with vertices (0,0), (1,0), (0,1).
    let mut b = MeshBuilder::new();
    let a = b.add_node(vec![0.0, 0.0]);
    let c = b.add_node(vec![1.0, 0.0]);
    let d = b.add_node(vec![0.0, 1.0]);
    b.add_element(CellType::Tri, vec![a, c, d]);
    let mesh = b.build();

    // Constant unit flux over the whole boundary. Because the shape functions on
    // a face sum to 1, the total assembled load equals the boundary measure.
    let mut rhs = vec![0.0; 3];
    apply_neumann(&mesh, 1, |_, _| 1.0, &mut rhs);
    let total: f64 = rhs.iter().sum();
    let expected = 2.0 + 2.0_f64.sqrt();
    assert!(
        (total - expected).abs() < 1e-9,
        "got {total} expected {expected}"
    );
    println!("Total Neumann load = {total} (expected {expected})");
}
