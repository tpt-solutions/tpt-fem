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
    assert all(0.0 <= v <= 1.0 for v in u.values)
    assert u.mesh is mesh
    out = tmp_path / "py_test.vtk"
    mesh.write_vtk(str(out), "u", u.values)
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
    assert len(u.values) == mesh.node_count() * 3
    assert u.dim == 3


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
    for m in modes:
        assert m.omega2 > 0.0
        assert len(m.shape) == mesh.node_count() * 3
    # Indexing and the omega/frequency accessors work.
    assert modes[0].omega2 == modes.omega2s()[0]
    assert modes[0].omega == modes.frequencies()[0]
    assert modes[0].omega == modes[0].omega2 ** 0.5


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


def test_readme_snippet():
    # Drift guard: the crate README documents this exact usage (Poisson +
    # elasticity + modal on a 3-D bar). If the bound API changes, this fails
    # instead of the docs silently diverging.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 4, 4])
    poisson_bcs = [
        (nid, 0.0)
        for axis in range(3)
        for coord in (0.0, 1.0)
        for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
    ]
    u = fem.solve_poisson(mesh, 1.0, 2, 1.0, poisson_bcs)
    assert len(u) == mesh.node_count()

    bar = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
    bcs = [(nid, c, 0.0) for nid in bar.nodes_on_plane(0, 0.0, 1e-9) for c in range(3)]
    disp = fem.solve_elasticity(bar, "3d", 200e9, 0.3, 2, bcs)
    assert len(disp) == bar.node_count() * 3
    modes = fem.solve_modal(bar, "3d", 200e9, 0.3, 7800.0, 2, 4, bcs)
    assert len(modes) == 4


def test_to_numpy_shapes():
    # Drift guard for the Jupyter-friendly result-object accessors: `to_numpy()`
    # must reshape the flat fields to the expected `(n,)` / `(n, dim)` arrays.
    np = pytest.importorskip("numpy")

    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 4, 4])
    bcs = [
        (nid, 0.0)
        for axis in range(3)
        for coord in (0.0, 1.0)
        for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
    ]
    u = fem.solve_poisson(mesh, 1.0, 2, 1.0, bcs)
    arr = u.to_numpy()
    assert arr.shape == (mesh.node_count(),)
    np.testing.assert_allclose(np.asarray(u.values), arr)

    bar = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
    ebcs = [(nid, c, 0.0) for nid in bar.nodes_on_plane(0, 0.0, 1e-9) for c in range(3)]
    disp = fem.solve_elasticity(bar, "3d", 200e9, 0.3, 2, ebcs)
    darr = disp.to_numpy()
    assert darr.shape == (bar.node_count(), 3)
    np.testing.assert_allclose(np.asarray(disp.values), darr.reshape(-1))

    modes = fem.solve_modal(bar, "3d", 200e9, 0.3, 7800.0, 2, 4, ebcs)
    marr = modes[0].to_numpy()
    assert marr.shape == (bar.node_count(), 3)


def test_to_pyvista_round_trips():
    # Drift guard for the pyvista interop: `to_pyvista()` must return a grid
    # whose point data matches the solution field. Skipped if pyvista isn't
    # installed.
    pytest.importorskip("pyvista")

    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3, 3, 3])
    bcs = [
        (nid, 0.0)
        for axis in range(3)
        for coord in (0.0, 1.0)
        for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
    ]
    u = fem.solve_poisson(mesh, 1.0, 2, 1.0, bcs)
    grid = u.to_pyvista()
    assert grid.n_points == mesh.node_count()
    assert "u" in grid.point_data


