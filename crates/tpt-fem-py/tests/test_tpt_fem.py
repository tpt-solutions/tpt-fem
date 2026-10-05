"""Smoke tests for the tpt-fem Python bindings.

Run with `maturin develop` then `pytest`, or `maturin pytest`.
"""

import pytest

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




def test_topopt_cantilever():
    sol = fem.topopt_cantilever(12, 6, 0.5, max_iter=5)
    assert sol.nx == 12 and sol.ny == 6
    assert len(sol.densities) == 12 * 6
    assert len(sol.compliance) == sol.iterations + 1
    assert sol.compliance[-1] < sol.compliance[0]
    assert abs(sum(sol.densities) / len(sol.densities) - 0.5) < 1e-3


def test_topopt_cantilever_rejects_bad_params():
    with pytest.raises(RuntimeError):
        fem.topopt_cantilever(12, 6, 0.0)
    with pytest.raises(RuntimeError):
        fem.topopt_cantilever(0, 6, 0.5)


def test_transient_heat_cools_toward_boundary_value():
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3, 3, 3])
    bcs = [
        (nid, 0.0)
        for axis in range(3)
        for coord in (0.0, 1.0)
        for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
    ]
    hist = fem.solve_transient_heat(mesh, 1.0, 1.0, 0.01, 20, 1.0, bcs)
    assert len(hist) == 21
    assert hist.times[0] == 0.0 and abs(hist.times[-1] - 0.2) < 1e-12
    interior = [i for i in range(mesh.node_count()) if i not in {b[0] for b in bcs}]
    first, last = hist[0].values, hist[-1].values
    assert all(last[i] < first[i] for i in interior)
    with pytest.raises(RuntimeError):
        fem.solve_transient_heat(mesh, 1.0, 1.0, 0.0, 5, 1.0, bcs)


def test_j2_uniaxial_response_yields_and_hardens():
    e, sy, h = 200e9, 250e6, 20e9
    eps_y = sy / e
    strains = [0.5 * eps_y, eps_y, 2 * eps_y, 4 * eps_y]
    elastic = fem.j2_uniaxial_response(e, 0.3, sy, strains)
    assert abs(elastic[0] - 0.5 * sy) < 1.0
    # Perfect plasticity: stress is capped at the yield stress.
    assert abs(elastic[2] - sy) < 1e3 and abs(elastic[3] - sy) < 1e3
    hard = fem.j2_uniaxial_response(e, 0.3, sy, strains, iso_hardening=h)
    assert hard[3] > hard[2] > sy
    with pytest.raises(RuntimeError):
        fem.j2_uniaxial_response(-1.0, 0.3, sy, strains)


def test_neo_hookean_uniaxial():
    out = fem.neo_hookean_uniaxial(2.0, [1.0, 2.0])
    assert out[0] == 0.0
    assert abs(out[1] - 2.0 * (2.0 - 0.25)) < 1e-12
    with pytest.raises(RuntimeError):
        fem.neo_hookean_uniaxial(2.0, [0.0])


def test_solve_darcy_linear_pressure_drop():
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 2, 2])
    left = mesh.nodes_on_plane(0, 0.0, 1e-9)
    right = mesh.nodes_on_plane(0, 1.0, 1e-9)
    bcs = [(n, 1.0) for n in left] + [(n, 0.0) for n in right]
    p = fem.solve_darcy(mesh, 2.0, bcs)
    for i in range(mesh.node_count()):
        x = mesh.coords(i)[0]
        assert abs(p.values[i] - (1.0 - x)) < 1e-9


def test_solve_stokes_poiseuille_profile():
    # Channel on the unit square driven in x; no-slip walls at y = 0, 1.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4, 4, 4])
    # 3-D duct: no-slip on the y/z walls, driven along x.
    walls = [n for ax in (1, 2) for c in (0.0, 1.0) for n in mesh.nodes_on_plane(ax, c, 1e-9)]
    bcs = [(n, k, 0.0) for n in set(walls) for k in range(3)]
    vel, pres = fem.solve_stokes(mesh, 1.0, [1.0, 0.0, 0.0], bcs)
    assert len(vel.values) == 3 * mesh.node_count()
    assert len(pres.values) == len(pres)
    assert max(vel.values) > 0.0
    with pytest.raises(RuntimeError):
        fem.solve_stokes(mesh, 1.0, [1.0, 0.0], bcs)


def test_laminate_abd_symmetric_has_zero_coupling():
    ply = (140e9, 10e9, 0.3, 5e9, 0.125e-3)
    stack = [ply + (0.0,), ply + (90.0,), ply + (90.0,), ply + (0.0,)]
    abd = fem.laminate_abd(stack)
    assert len(abd) == 6 and all(len(r) == 6 for r in abd)
    for i in range(3):
        for j in range(3):
            assert abs(abd[i][3 + j]) < 1e-6 * abs(abd[i][j] or 1.0)
    assert abd[0][0] > 0 and abd[3][3] > 0
    with pytest.raises(RuntimeError):
        fem.laminate_abd([])


def test_newmark_sdof_matches_closed_form():
    import math
    hist = fem.newmark([[1.0]], [[0.0]], [[4.0]], [1.0], [0.0], [0.0], 0.01, 200)
    t, u = hist[-1]
    assert abs(u[0] - math.cos(2.0 * t)) < 1e-3
    # Callable load + validation.
    hist = fem.newmark([[1.0]], [[0.1]], [[4.0]], [0.0], [0.0], lambda t: [1.0], 0.01, 10)
    assert len(hist) == 11
    with pytest.raises(RuntimeError):
        fem.newmark([[0.0]], [[0.0]], [[1.0]], [1.0], [0.0], [0.0], 0.01, 5)
    with pytest.raises(RuntimeError):
        fem.newmark([[1.0, 0.0]], [[0.0]], [[1.0]], [1.0], [0.0], [0.0], 0.01, 5)


def test_contact_pairs_nearest_and_empty():
    a = [(0, [0.0, 0.0, 0.0]), (1, [10.0, 0.0, 0.0])]
    b = [(7, [1.0, 0.0, 0.0]), (8, [9.0, 0.0, 0.0])]
    pairs = dict(fem.contact_pairs(a, b))
    # The paired value is the *index* into `b`, not its id.
    assert pairs[0][0] == 0 and abs(pairs[0][1] - 1.0) < 1e-12
    assert pairs[1][0] == 1
    assert fem.contact_pairs(a, []) == [(0, None), (1, None)]


def test_thermal_structural_free_expansion():
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2])
    n = mesh.node_count()
    # Kinematic minimum (no rigid-body motion): pin the origin, the x-axis end
    # in y/z and the y-axis end in z; the body then expands freely.
    def at(x, y, z):
        e = 1e-9
        return mesh.nodes_in_box([x - e, y - e, z - e], [x + e, y + e, z + e])[0]

    origin, xend, yend = at(0, 0, 0), at(1, 0, 0), at(0, 1, 0)
    bcs = [(origin, k, 0.0) for k in range(3)]
    bcs += [(xend, 1, 0.0), (xend, 2, 0.0), (yend, 2, 0.0)]
    sol = fem.solve_thermal_structural(mesh, "3d", 1.0, 0.3, 1e-3, [10.0] * n, bcs)
    far = xend
    # Free expansion: u_x(1,0,0) = alpha * dT * L = 1e-2.
    assert abs(sol.values[3 * far] - 1e-2) < 1e-6
    with pytest.raises(RuntimeError):
        fem.solve_thermal_structural(mesh, "3d", 1.0, 0.3, 1e-3, [1.0], bcs)


def test_contact_augmented_lagrangian_holds_the_wall():
    # Spring K=10 pushed toward -x by 4; wall at x >= 0: u -> 0, reaction 4.
    u, lam = fem.contact_augmented_lagrangian([[10.0]], [-4.0], [(0, 0.0)])
    assert abs(u[0]) < 1e-6 and abs(lam[0] - 4.0) < 1e-3
    with pytest.raises(RuntimeError):
        fem.contact_augmented_lagrangian([[10.0]], [-4.0], [(3, 0.0)])
    with pytest.raises(RuntimeError):
        fem.contact_augmented_lagrangian([[10.0]], [-4.0], [(0, 0.0)], penalty=-1.0)


def test_fsi_interface_loads_resultant_matches_pressure_times_area():
    # Unit cube; uniform pressure 5 on the top face (z = 1), area 1 => total
    # outward force 5 in +z, summed over the interface nodes.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2])
    top = mesh.nodes_on_plane(2, 1.0, 1e-9)
    interface = [(n, n) for n in top]
    pressure = [5.0] * mesh.node_count()
    loads = fem.fsi_interface_loads(mesh, interface, pressure)
    assert len(loads) == 3 * mesh.node_count()
    fz = sum(loads[3 * n + 2] for n in range(mesh.node_count()))
    assert abs(abs(fz) - 5.0) < 1e-9
    with pytest.raises(RuntimeError):
        fem.fsi_interface_loads(mesh, [(0, 10_000)], pressure)


def test_gpu_api_is_consistent_with_build():
    if not fem.gpu_enabled():
        with pytest.raises(RuntimeError):
            fem.gpu_solve_cg([(0, 0, 1.0)], [1.0])
        return
    try:
        name = fem.gpu_adapter()
    except RuntimeError:
        pytest.skip("no GPU adapter")
    assert name
    n = 30
    trip = []
    for i in range(n * n):
        trip.append((i, i, 4.0))
        if i % n + 1 < n:
            trip += [(i, i + 1, -1.0), (i + 1, i, -1.0)]
        if i + n < n * n:
            trip += [(i, i + n, -1.0), (i + n, i, -1.0)]
    x, its, res = fem.gpu_solve_cg(trip, [1.0] * (n * n))
    assert res <= 1e-10 and its > 0 and len(x) == n * n
