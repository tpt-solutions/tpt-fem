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
`nodes_in_box` / `write_vtk`) and the solver functions `solve_poisson`
(steady heat conduction), `solve_elasticity` (linear statics), and `solve_modal`
(natural-vibration eigenproblem `K φ = ω² M φ`). The Poisson source may be a
constant `float` or a Python callable `f(x, y, z)`. Errors from the core crates
are surfaced as Python exceptions via their `Display` impls.

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
import tpt_fem as fem

# Build a structured box mesh.
mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [8, 8, 8])

# Solve -∇·(k∇u) = 1 with a Python source callable.
def source(x, y, z):
    return 1.0

u = fem.solve_poisson(mesh, 1.0, 2, source, [])

# Linear elasticity (3-D continuum): clamp one face, free elsewhere.
bcs = [(nid, c, 0.0) for nid in mesh.nodes_on_plane(0, 0.0, 1e-9) for c in range(3)]
u = fem.solve_elasticity(mesh, "3d", 200e9, 0.3, 2, bcs)

# Natural vibration: first 4 modes of the K φ = ω² M φ eigenproblem.
modes = fem.solve_modal(mesh, "3d", 200e9, 0.3, 7800.0, 2, 4, bcs)
for w2, shape in modes:
    print("ω² =", w2)

# Export the solution to VTK.
mesh.write_vtk("result.vtk", "u", u)
```

## API highlights

| Item | Description |
|------|-------------|
| `Mesh.load(path)` | Load a mesh from a file. |
| `Mesh.box_mesh(min, max, n)` | Build a structured box mesh. |
| `Mesh.coords` | Node coordinate array. |
| `Mesh.nodes_on_plane(...)` / `Mesh.nodes_in_box(...)` | Node selection helpers. |
| `Mesh.write_vtk(path)` | Export to VTK. |
| `solve_poisson(mesh, conductivity, quad_order, source, bcs)` | Solve Poisson with constant or callable source. |
| `solve_elasticity(mesh, model, young, poisson, quad_order, bcs)` | Linear statics; `model` ∈ `bar`/`plane-stress`/`plane-strain`/`3d`. `bcs` are `(node, component, value)` tuples. |
| `solve_modal(mesh, model, young, poisson, density, quad_order, num_modes, bcs)` | Natural-vibration eigenproblem; returns a list of `(ω², mode_shape)` pairs. |

## Position in the crate stack

```text
tpt-fem (umbrella) ──► tpt-fem-py (pyo3 / maturin bindings)
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
