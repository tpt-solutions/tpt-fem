"""End-to-end example: thermal solve, elasticity, and modal analysis.

Run after `maturin develop`:

    python examples/walkthrough.py
"""

import tpt_fem as fem


def main():
    # 1. Steady Poisson / heat conduction on a unit cube.
    mesh = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [8, 8, 8])
    poisson_bcs = [
        (nid, 0.0)
        for axis in range(3)
        for coord in (0.0, 1.0)
        for nid in mesh.nodes_on_plane(axis, coord, 1e-9)
    ]
    u = fem.solve_poisson(mesh, 1.0, 2, 1.0, poisson_bcs)
    print(f"Poisson: {len(u)} DOFs, u in [{min(u):.4f}, {max(u):.4f}]")
    mesh.write_vtk("walkthrough_poisson.vtk", "u", u)

    # 2. Linear elasticity (3-D continuum) on a slender clamped bar.
    bar = fem.Mesh.box_mesh([0.0, 0.0, 0.0], [1.0, 0.2, 0.2], [8, 2, 2])
    elas_bcs = [
        (nid, c, 0.0)
        for nid in bar.nodes_on_plane(0, 0.0, 1e-9)
        for c in range(3)
    ]
    disp = fem.solve_elasticity(bar, "3d", 200e9, 0.3, 2, elas_bcs)
    print(f"Elasticity: {len(disp)} DOFs (3 per node)")

    # 3. Natural-vibration modes of the same bar.
    modes = fem.solve_modal(bar, "3d", 200e9, 0.3, 7800.0, 2, 4, elas_bcs)
    for i, (w2, shape) in enumerate(modes, 1):
        print(f"  mode {i}: omega^2 = {w2:.4e}")
    bar.write_vtk("walkthrough_modal.vtk", "mode1", modes[0][1])


if __name__ == "__main__":
    main()
