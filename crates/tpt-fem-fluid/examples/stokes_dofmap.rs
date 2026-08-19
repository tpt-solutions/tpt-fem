//! Mixed velocity/pressure DOF layout produced by `stokes_dofmap`.
//!
//! Run with: `cargo run -p tpt-fem-fluid --example stokes_dofmap`
//!
//! `stokes_dofmap` builds a `tpt-fem-dofmap::MultiFieldDofMap` with two fields —
//! `"velocity"` (one component per spatial dimension) followed by `"pressure"`
//! (one component) — using `Layout::Interleaved`, so every node owns
//! `dim + 1` consecutive global DOFs. This example prints that numbering for a
//! 2-D quad mesh and for a single 3-D hex element, and contrasts it with the
//! velocity-only numbering (`node*dim + component`) that `steady_stokes`
//! actually solves and returns.

use tpt_fem_fluid::stokes_dofmap;
use tpt_fem_mesh::{CellType, MeshBuilder};

fn main() {
    // ---------------------------------------------------------------- 2-D quads
    // Two Quad4 cells sharing an edge: 6 nodes, dim = 2 -> 3 DOFs per node.
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0]);
    let n2 = b.add_node(vec![2.0, 0.0]);
    let n3 = b.add_node(vec![0.0, 1.0]);
    let n4 = b.add_node(vec![1.0, 1.0]);
    let n5 = b.add_node(vec![2.0, 1.0]);
    b.add_element(CellType::Quad, vec![n0, n1, n4, n3]);
    b.add_element(CellType::Quad, vec![n1, n2, n5, n4]);
    let mesh = b.build();

    let map = stokes_dofmap(&mesh);
    let dim = map.components(0);
    println!("2-D Quad4 mesh: {} nodes, dim = {dim}", mesh.node_count());
    println!(
        "fields = [{} ({} comp), {} ({} comp)], layout = {:?}, ndof = {}",
        map.fields[0].name,
        map.fields[0].components,
        map.fields[1].name,
        map.fields[1].components,
        map.layout,
        map.ndof
    );
    println!();
    println!("  node      x     y    u_x   u_y      p   |  velocity-only dofs");
    println!("  ----   -----  ----  ----  ----   ----   |  -----------------");
    for n in 0..mesh.node_count() {
        let c = mesh.node_coords(n);
        let ux = map.node_field_dof(n, 0, 0);
        let uy = map.node_field_dof(n, 0, 1);
        let p = map.node_field_dof(n, 1, 0);
        println!(
            "  {n:>4}   {:5.2} {:5.2}  {ux:>4}  {uy:>4}   {p:>4}   |  {:>4}  {:>4}",
            c[0],
            c[1],
            n * dim,
            n * dim + 1
        );
    }

    // Interleaved layout: node n owns [n*(dim+1) .. n*(dim+1)+dim].
    for n in 0..mesh.node_count() {
        assert_eq!(map.node_field_dof(n, 0, 0), n * (dim + 1));
        assert_eq!(map.node_field_dof(n, 0, 1), n * (dim + 1) + 1);
        assert_eq!(map.node_field_dof(n, 1, 0), n * (dim + 1) + dim);
        assert_eq!(map.dofs_of(n).len(), dim + 1);
    }
    assert_eq!(map.ndof, mesh.node_count() * (dim + 1));
    println!("\nverified: node n owns velocity {{3n, 3n+1}} and pressure {{3n+2}}");

    // ----------------------------------------------------------------- 3-D hex
    // One Hex8 cell: dim = 3 -> 4 DOFs per node (u_x, u_y, u_z, p).
    let mut b3 = MeshBuilder::new();
    let mut ids = Vec::new();
    for &(x, y, z) in &[
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ] {
        ids.push(b3.add_node(vec![x, y, z]));
    }
    b3.add_element(CellType::Hex, ids);
    let mesh3 = b3.build();
    let map3 = stokes_dofmap(&mesh3);
    println!(
        "\n3-D Hex8 mesh: {} nodes, dim = {}, ndof = {} ({} per node)",
        mesh3.node_count(),
        map3.components(0),
        map3.ndof,
        map3.ndof / mesh3.node_count()
    );
    println!("  node 3 dofs = {:?}  (u_x, u_y, u_z, p)", map3.dofs_of(3));
    assert_eq!(map3.components(0), 3);
    assert_eq!(map3.components(1), 1);
    assert_eq!(map3.ndof, 8 * 4);
    assert_eq!(map3.dofs_of(3), &[12, 13, 14, 15]);
    println!("\nverified: dim is taken from the first cell type (Hex8 -> 3)");

    // The penalty solver itself is velocity-only: `steady_stokes` returns a
    // velocity vector indexed by `node*dim + component` and a nodal pressure
    // vector with one entry per node.
    println!(
        "\nsteady_stokes velocity length for the 2-D mesh = {} ({} nodes x {dim}),\
         \npressure length = {} (1 per node)",
        mesh.node_count() * dim,
        mesh.node_count(),
        mesh.node_count()
    );
}
