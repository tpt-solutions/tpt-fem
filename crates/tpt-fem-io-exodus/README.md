# tpt-fem-io-exodus

Exodus II mesh reader/writer for
[tpt-fem](https://github.com/tpt-solutions/tpt-fem) — the mesh-based finite
element core from [tpt-solutions](https://github.com/tpt-solutions).

## Overview

Exodus II is a NetCDF-based format with no mature pure-Rust implementation, so
this crate ships a minimal NetCDF-3 *classic* (CDF-1, big-endian) codec and
builds the subset of Exodus variables required to round-trip linear meshes:
`coords`, one `connectN` set per element block, `elem_blk_id`, `eb_status`,
`eb_prop1`, `eb_names`, `elem_num_map`, `node_num_map`, and `time_whole`. Only
the five linear element types are supported.

## Installation

```toml
[dependencies]
tpt-fem-io-exodus = "0.1"
```

## Usage

```rust
use tpt_fem_io_exodus::{read_exodus, write_exodus, bytes_to_mesh, mesh_to_exodus_bytes};

// Read an Exodus II file from disk.
let mesh = read_exodus("mesh.ex2")?;

// Write a Mesh to an Exodus II file.
write_exodus(&mesh, "mesh.ex2")?;

// Or work with raw bytes (useful for in-memory pipelines).
let bytes = mesh_to_exodus_bytes(&mesh)?;
let mesh2 = bytes_to_mesh(&bytes)?;
```

## API highlights

| Item | Description |
|------|-------------|
| `read_exodus` | Read an Exodus II file → `Mesh`. |
| `write_exodus` | Write a `Mesh` → Exodus II file. |
| `bytes_to_mesh` | Decode NetCDF-3 classic bytes → `Mesh`. |
| `mesh_to_exodus_bytes` | Encode a `Mesh` → NetCDF-3 classic bytes. |
| `ExodusError` | Error type for decode/encode failures. |

## Position in the crate stack

```text
tpt-fem-mesh ──► tpt-fem-io-exodus ──► Exodus II / NetCDF
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
