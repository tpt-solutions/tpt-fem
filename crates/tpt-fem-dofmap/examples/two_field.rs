//! A genuine two-field map: 3-component velocity + 1-component pressure.
//!
//! Run with: `cargo run -p tpt-fem-dofmap --example two_field`
//!
//! Builds a 4-node mesh, assigns a velocity field with 3 components and a
//! pressure field with 1 component under `Interleaved` (node-major) layout, and
//! prints the resulting global DOF layout. The node/field/component -> global
//! DOF translation is hand-checked against the expected numbering and asserted.

use tpt_fem_dofmap::{FieldSpec, Layout, MultiFieldDofMap};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0]);
    let n3 = b.add_node(vec![1.0, 1.0]);
    b.add_element(CellType::Quad, vec![n0, n1, n2, n3]);
    let mesh = b.build();

    let map = MultiFieldDofMap::new(
        &mesh,
        &[FieldSpec::new("velocity", 3), FieldSpec::new("pressure", 1)],
        Layout::Interleaved,
    );

    println!("two-field map: velocity(3) + pressure(1), Interleaved layout");
    println!(
        "  fields        = {:?}",
        map.fields.iter().map(|f| &f.name).collect::<Vec<_>>()
    );
    println!("  total dofs     = {}", map.ndof);
    println!("  field_ranges   = {:?}", map.field_ranges);

    for n in 0..mesh.node_count() {
        println!("  node {n}: global dofs = {:?}", map.dofs_of(n));
    }

    // Interleaved: node n owns dofs [4n, 4n+3]; components 0..3 are velocity,
    // component 3 is pressure.
    let total = 4;
    assert_eq!(map.ndof, mesh.node_count() * total);
    for n in 0..mesh.node_count() {
        assert_eq!(map.node_field_dof(n, 0, 0), n * total);
        assert_eq!(map.node_field_dof(n, 0, 2), n * total + 2);
        assert_eq!(map.node_field_dof(n, 1, 0), n * total + 3);
        assert_eq!(
            map.dofs_of(n),
            vec![n * total, n * total + 1, n * total + 2, n * total + 3]
        );
    }
    println!("\nverified: (node, field, component) -> global DOF translation is correct");
}
