//! Abaqus `.inp` import: read a deck and report its contents.
//!
//! Parses a small Abaqus input deck with [`read_inp`] / [`read_inp_deck`] (from
//! the umbrella prelude), prints the resulting mesh and the captured
//! analysis metadata, then exports to VTK. Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example abaqus_import
//! ```

use tpt_fem::prelude::*;

fn main() {
    let deck = "\
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 0.0, 1.0, 0.0
4, 1.0, 1.0, 0.0
*ELEMENT, TYPE=CPS4
1, 1, 2, 4, 3
*NSET, NSET=fixed
1, 3
*ELSET, ELSET=plate
1
*MATERIAL, NAME=STEEL
*ELASTIC
200.0, 0.3
*BOUNDARY
1, 1, 0.0
1, 2, 0.0
";

    // Geometry-only import.
    let mesh = read_inp(deck).expect("read mesh");
    println!(
        "Mesh: {} nodes, {} elements (cell type {:?})",
        mesh.node_count(),
        mesh.element_count(),
        mesh.elements[0].cell_type,
    );

    // Full-deck import (retains sets / material / boundary).
    let full = read_inp_deck(deck).expect("read deck");
    println!(
        "  nodeset 'fixed' size: {}",
        full.nsets.get("fixed").map_or(0, |v| v.len())
    );
    println!(
        "  elset 'plate'  size: {}",
        full.elsets.get("plate").map_or(0, |v| v.len())
    );
    println!("  material 'STEEL': {:?}", full.materials.get("STEEL"));
    println!("  prescribed boundary conditions: {}", full.boundary.len());

    write_vtk(&mesh, "abaqus_import.vtk").expect("write vtk");
    println!("Wrote abaqus_import.vtk");
}
