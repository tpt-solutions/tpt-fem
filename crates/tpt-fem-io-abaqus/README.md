# tpt-fem-io-abaqus

Abaqus `.inp` mesh reader/writer for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

A minimal but practical subset of the Abaqus input deck is supported:

- `*NODE` and `*ELEMENT, TYPE=…` sections, mapping the common linear-element
  types (`T2D2`, `CPS3`/`CPS4`, `C3D4`, `C3D8`) onto the `tpt-fem-mesh`
  `CellType`s.
- All other sections are ignored.

## Installation

```toml
[dependencies]
tpt-fem-io-abaqus = "0.1"
```

## Usage

```rust
use tpt_fem_io_abaqus::{read_inp, write_inp};

// Parse an Abaqus input deck (a subset) into a tpt-fem Mesh.
let mesh = read_inp(deck_text)?;

// Serialise a Mesh back to a minimal Abaqus `.inp` file.
write_inp(&mesh, "mesh.inp")?;
```

Element-type mapping:

| Abaqus type | `CellType` |
|-------------|------------|
| `T2D2` | `Line` |
| `CPS3` | `Tri` |
| `CPS4` / `CPE4` | `Quad` |
| `C3D4` | `Tet` |
| `C3D8` | `Hex` |

## API highlights

| Item | Description |
|------|-------------|
| `read_inp` | Parse an Abaqus deck (subset) into a `Mesh`. |
| `write_inp` | Serialise a `Mesh` to a `.inp` file. |
| `InpError` | Error type for parsing/serialisation failures. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-io-abaqus ──► Abaqus (.inp)
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
