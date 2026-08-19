//! `Layout::Block` vs `Layout::Interleaved`, printed side by side.
//!
//! Run with: `cargo run -p tpt-fem-dofmap --example layout_compare`
//!
//! Builds a 3-node mesh with a 2-component displacement field and a 1-component
//! pressure field, then shows how the two layouts assign global DOF indices to
//! the same `(node, field, component)` triples. The difference is purely the
//! ordering: `Block` groups all DOFs of a field contiguously; `Interleaved`
//! keeps a node's fields adjacent.

use tpt_fem_dofmap::{FieldSpec, Layout, MultiFieldDofMap};
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0]);
    b.add_node(vec![1.0]);
    let n2 = b.add_node(vec![2.0]);
    b.add_element(CellType::Line, vec![n0, n2]);
    let mesh = b.build();

    let fields = [FieldSpec::new("u", 2), FieldSpec::new("p", 1)];
    let block = MultiFieldDofMap::new(&mesh, &fields, Layout::Block);
    let inter = MultiFieldDofMap::new(&mesh, &fields, Layout::Interleaved);

    println!("layout comparison: 3 nodes, field u(2) + field p(1)\n");
    println!("field_ranges (Block)  = {:?}", block.field_ranges);
    println!("field_ranges (Interle) = {:?}", inter.field_ranges);
    println!();

    print!("  (node,field,comp)");
    for n in 0..mesh.node_count() {
        for f in 0..fields.len() {
            for c in 0..fields[f].components {
                print!("  n{n} f{f} c{c}");
            }
        }
    }
    println!();

    print!("  Block        ");
    for n in 0..mesh.node_count() {
        for f in 0..fields.len() {
            for c in 0..fields[f].components {
                print!("  {:4}", block.node_field_dof(n, f, c));
            }
        }
    }
    println!();

    print!("  Interleaved  ");
    for n in 0..mesh.node_count() {
        for f in 0..fields.len() {
            for c in 0..fields[f].components {
                print!("  {:4}", inter.node_field_dof(n, f, c));
            }
        }
    }
    println!();

    // Block groups field u (all nodes) then field p; Interleaved keeps each
    // node's u,p,p adjacent.
    assert_eq!(block.field_range(0), (0, 6));
    assert_eq!(block.field_range(1), (6, 3));
    println!("\nverified: Block field u spans [0,6), field p spans [6,9)");
}
