"""Smoke tests for the tpt-fem Python bindings.

Run with `maturin develop` then `pytest`, or `maturin pytest`.
"""

import tpt_fem as fem


def test_box_mesh_and_solve(tmp_path):
    # Unit cube, Dirichlet u=0 on every boundary face, source f=1.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 4, 4])
    assert mesh.coords(0) == [0.0, 0.0, 0.0]

    bcs = []
    for axis in range(3):
        for coord in (0.0, 1.0):
            for nid in mesh.nodes_on_plane(axis, coord, 1e-9):
                bcs.append((nid, 0.0))

    u = fem.solve_poisson(mesh, 1.0, 2, 1.0, bcs)
    assert all(0.0 <= v <= 1.0 for v in u)
    out = tmp_path / "py_test.vtk"
    mesh.write_vtk(str(out), "u", u)
    assert out.exists()


def test_solve_elasticity_3d():
    # Slender 3-D bar, clamp one face, zero body load => the trivial
    # displacement field. Verifies the `solve_elasticity` binding (3-D
    # continuum model with per-component `(node, comp, value)` BCs).
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
    bcs = []
    for nid in mesh.nodes_on_plane(0, 0.0, 1e-9):
        for c in range(3):
            bcs.append((nid, c, 0.0))
    u = fem.solve_elasticity(mesh, "3d", 200e9, 0.3, 2, bcs)
    assert len(u) == mesh.node_count() * 3


def test_solve_modal_3d():
    # Same clamped bar: natural-vibration eigenproblem must yield positive
    # squared frequencies and one mode shape per requested mode.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
    bcs = []
    for nid in mesh.nodes_on_plane(0, 0.0, 1e-9):
        for c in range(3):
            bcs.append((nid, c, 0.0))
    modes = fem.solve_modal(mesh, "3d", 200e9, 0.3, 7800.0, 2, 3, bcs)
    assert len(modes) == 3
    for w2, shape in modes:
        assert w2 > 0.0
        assert len(shape) == mesh.node_count() * 3


def test_python_callback_source():
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3, 3, 3])
    bcs = []
    for axis in range(3):
        for coord in (0.0, 1.0):
            for nid in mesh.nodes_on_plane(axis, coord, 1e-9):
                bcs.append((nid, 0.0))

    # Source = x + y + z evaluated at the quadrature point.
    def src(x, y, z):
        return x + y + z

    u = fem.solve_poisson(mesh, 1.0, 2, src, bcs)
    assert len(u) == mesh.node_count()

