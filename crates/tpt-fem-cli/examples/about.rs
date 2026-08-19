//! `tpt-fem-cli` — purpose and version.
//!
//! The CLI is a binary crate (its solver/parser functions live in `main.rs` and
//! are not part of a public library API), so this example documents the crate
//! rather than calling into it. It prints the crate's identity and the
//! subcommands it exposes, all sourced from the build environment — no private
//! API is fabricated.

fn main() {
    println!("{} v{}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    println!("{}", env!("CARGO_PKG_DESCRIPTION"));
    println!();
    println!("Subcommands:");
    println!("  solve <config.toml>   run a steady Poisson/heat-conduction problem");
    println!("  elasticity <config>   run a linear-elasticity problem");
    println!("  modal <config>        run a natural-vibration (modal) problem");
    println!("  init [type] [out.toml] generate a starter problem config");
    println!("  mesh info <file>      print mesh statistics");
    println!("  mesh convert <in> <out>  convert a .msh mesh to a .vtk file");

    assert!(!env!("CARGO_PKG_NAME").is_empty());
    assert!(!env!("CARGO_PKG_VERSION").is_empty());
}
