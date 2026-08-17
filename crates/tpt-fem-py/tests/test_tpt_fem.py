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
