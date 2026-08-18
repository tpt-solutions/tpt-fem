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

# Dirichlet u = 0 on every boundary face, then solve -∇·(k∇u) = 1.
poisson_bcs = [
    (nid, 0.0)
    for axis in range(3)
    for coord in (0.0, 1.0)
    for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
]
u = fem.solve_poisson(mesh, 1.0, 2, 1.0, poisson_bcs)

# Linear elasticity (3-D continuum) on a slender clamped bar.
bar = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
bcs = [(nid, c, 0.0) for nid in bar.nodes_on_plane(0, 0.0, 1e-9) for c in range(3)]
u = fem.solve_elasticity(bar, "3d", 200e9, 0.3, 2, bcs)

# Natural vibration: first 4 modes of the K φ = ω² M φ eigenproblem.
modes = fem.solve_modal(bar, "3d", 200e9, 0.3, 7800.0, 2, 4, bcs)
for m in modes:
    print("ω² =", m.omega2, "ω =", m.omega)

# Export a solution to VTK.
bar.write_vtk("result.vtk", "u", u.values)
```

## Visualization & Jupyter

Each solver returns a Jupyter-friendly **result object** (not a bare `list`):
`PoissonSolution`, `ElasticitySolution`, and `ModalSolution` (whose elements
are `ModeShape` objects). They carry rich `__repr__` / `_repr_html_` display and
two convenience accessors:

* `to_numpy()` — returns an `np.ndarray` (`(n_nodes,)` for scalar Poisson,
  `(n_nodes, dim)` for the vector fields). Requires `numpy`.
* `to_pyvista()` — returns a `pyvista.UnstructuredGrid` with the field attached
  as point data (`"u"` / `"disp"` / `"mode"`), so results go straight into
  PyVista / ParaView / matplotlib without a manual VTK export round-trip.
  Requires `pyvista` (and `numpy`).

```python
import numpy as np, pyvista as pv

u = fem.solve_poisson(mesh, 1.0, 2, 1.0, poisson_bcs)
field = u.to_numpy()          # np.ndarray, shape (n_nodes,)
grid = u.to_pyvista()         # pyvista.UnstructuredGrid
grid.plot(scalars="u")        # works in a Jupyter notebook

modes = fem.solve_modal(bar, "3d", 200e9, 0.3, 7800.0, 2, 4, bcs)
m0 = modes[0]                 # ModeShape (indexable / iterable)
print(m0.omega, m0.omega2)    # natural frequency and its square
m0.to_pyvista().plot(scalars="mode")
```

## API highlights

| Item | Description |
|------|-------------|
| `Mesh.load(path)` | Load a mesh from a file. |
| `Mesh.box_mesh(min, max, n)` | Build a structured box mesh. |
| `Mesh.coords` | Node coordinate array. |
| `Mesh.nodes_on_plane(...)` / `Mesh.nodes_in_box(...)` | Node selection helpers. |
| `Mesh.write_vtk(path)` | Export to VTK. |
| `solve_poisson(...)` → `PoissonSolution` | Steady heat conduction; constant or callable source. `.values` is `(n_nodes,)`; `.to_numpy()` / `.to_pyvista()` for interop. |
| `solve_elasticity(...)` → `ElasticitySolution` | Linear statics; `model` ∈ `bar`/`plane-stress`/`plane-strain`/`3d`. `.values` is `(n_nodes * dim,)`; `.to_numpy()` is `(n_nodes, dim)`. |
| `solve_modal(...)` → `ModalSolution` | Natural-vibration eigenproblem; indexable / iterable over `ModeShape` objects (`.omega2`, `.omega`, `.to_numpy()`, `.to_pyvista()`). `.omega2s()` / `.frequencies()` list ω² / ω. |

## Position in the crate stack

```text
tpt-fem (umbrella) ──► tpt-fem-py (pyo3 / maturin bindings)
```

## License

Licensed under either of [MIT](../../LICENSE-MIT) or
[Apache-2.0](../../LICENSE-APACHE) at your option.
