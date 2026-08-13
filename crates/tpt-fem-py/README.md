# tpt-fem-py

Python bindings for the [tpt-fem](https://github.com/tpt-solutions/tpt-fem)
core — the mesh-based finite element method library from
[tpt-solutions](https://github.com/tpt-solutions). Built with
[`maturin`](https://crates.io/crates/maturin) and
[`pyo3`](https://crates.io/crates/pyo3).

> This crate is excluded from the Cargo workspace (`exclude = ["…/tpt-fem-py"]`)
> and is intended for Python-side development this pass.

## Overview

Exposes a `Mesh` class (`load` / `box_mesh` / `coords` / `nodes_on_plane` /
`nodes_in_box` / `write_vtk`) and a `solve_poisson` function accepting either a
constant volumetric source or a Python callable `f(x, y, z)`. Errors from the
core crates are surfaced as Python exceptions via their `Display` impls.

## Installation (development)

```sh
# From the crate directory, build a debug wheel and install into the active env.
maturin develop
```

Or build a release wheel:

```sh
maturin build --release
pip install target/wheels/*.whl
```

## Usage

```python
import tpt_fem_py as fem

# Build a structured box mesh.
mesh = fem.Mesh.box_mesh((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (8, 8, 8))

# Solve -∇·(∇u) = 1 with a Python source callable.
def source(x, y, z):
    return 1.0

u = fem.solve_poisson(mesh, source)

# Export the solution to VTK.
mesh.write_vtk("result.vtk")
```

## API highlights

| Item | Description |
|------|-------------|
| `Mesh.load(path)` | Load a mesh from a file. |
| `Mesh.box_mesh(min, max, n)` | Build a structured box mesh. |
| `Mesh.coords` | Node coordinate array. |
| `Mesh.nodes_on_plane(...)` / `Mesh.nodes_in_box(...)` | Node selection helpers. |
| `Mesh.write_vtk(path)` | Export to VTK. |
| `solve_poisson(mesh, source)` | Solve Poisson with constant or callable source. |

## Position in the crate stack

```text
tpt-fem (umbrella) ──► tpt-fem-py (pyo3 / maturin bindings)
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
