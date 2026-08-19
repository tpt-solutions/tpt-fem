//! Single-field numbering degenerating to `Mesh::number_dofs`.
//!
//! Run with: `cargo run -p tpt-fem-dofmap --example single_field`
//!
//! Builds a one-field map over a 2-node line and asserts that its per-node DOF
//! lists and total count match `Mesh::number_dofs` exactly (the single-field
//! map is defined to replicate the legacy uniform numbering).

use tpt_fem_dofmap::{FieldSpec, Layout, MultiFieldDofMap};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let a = b.add_node(vec![0.0]);
    let c = b.add_node(vec![1.0]);
    b.add_element(CellType::Line, vec![a, c]);
    let mesh = b.build();

    let components = 3;
    let map = MultiFieldDofMap::new(&mesh, &[FieldSpec::new("u", components)], Layout::Block);
    let legacy = mesh.number_dofs(components);

    println!("single-field map (3 components per node, Block layout)");
    println!("  ndof            = {}", map.ndof);
    println!("  legacy ndof     = {}", legacy.ndof);
    for n in 0..mesh.node_count() {
        println!("  node {n}: dofs = {:?}", map.dofs_of(n));
        assert_eq!(map.dofs_of(n), &legacy.node_dofs[n]);
    }
    assert_eq!(map.ndof, legacy.ndof);
    println!("\nverified: single-field MultiFieldDofMap == Mesh::number_dofs");
}
