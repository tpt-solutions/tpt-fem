# tpt-fem-cli

Command-line driver for the [tpt-fem](https://github.com/tpt-solutions/tpt-fem)
core — the mesh-based finite element method library from
[tpt-solutions](https://github.com/tpt-solutions).

The binary is named `tpt-fem` and is built from `crates/tpt-fem-cli`.

## Overview

Subcommands:

- `solve` — run a steady Poisson/heat-conduction problem from a TOML config.
- `elasticity` — run a linear-elasticity statics problem from a TOML config.
- `modal` — extract natural vibration modes from a TOML config.
- `amr` — adaptive h-refinement Poisson solve on the unit square (quadtree,
  ZZ error estimator, Dörfler marking).
- `mesh info` — print summary statistics about a mesh file.
- `mesh convert` — convert a Gmsh `.msh` mesh to a ParaView `.vtk` file.

Error messages reuse the `Display` impls from the core crates, so malformed
input reports a human-readable cause rather than a panic.

## Installation

```toml
[dependencies]
tpt-fem-cli = "0.1"
```

Or install the binary:

```sh
cargo install --path crates/tpt-fem-cli
```

## Usage

```text
# Solve a steady Poisson/heat-conduction problem from a TOML config.
# `solve` takes a single positional config file (no --config/--out flags).
tpt-fem solve problem.toml

# Print mesh statistics.
tpt-fem mesh info mesh.msh

# Convert a Gmsh .msh mesh to a ParaView .vtk file.
tpt-fem mesh convert mesh.msh mesh.vtk
```

## Source layout

| Item | Description |
|------|-------------|
| `solve` | TOML-configured steady Poisson/heat-conduction run. |
| `elasticity` | TOML-configured linear-elasticity statics run. |
| `modal` | TOML-configured natural-vibration mode extraction. |
| `amr` | Adaptive h-refinement Poisson solve (`--max-elements`, `--theta`, `--constant`, `-o/--output`). |
| `mesh info` | Mesh summary statistics. |
| `mesh convert` | Gmsh `.msh` → ParaView `.vtk` conversion. |

## Position in the crate stack

```text
tpt-fem (umbrella) ──► tpt-fem-cli
```

## Examples

| Example | Command | Description |
|---------|---------|-------------|
| `about` | `cargo run -p tpt-fem-cli --example about` | Prints the CLI's name, version, and subcommand summary (the solver logic lives in `main.rs`, not a public library API). |

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
