"""Type stubs for the compiled `tpt_fem` extension module.

These annotations describe the public Python surface of the `tpt-fem` core
bindings. They are intentionally conservative: `to_numpy()` and `to_pyvista()`
return `Any` rather than `numpy.ndarray` / `pyvista.UnstructuredGrid` so the
stub stays usable without `numpy`/`pyvista` installed (those are optional
`viz` extras). `source` for `solve_poisson` accepts either a `float` or a
callable `f(x, y, z) -> float`.
"""

from typing import Any, Callable, Sequence

__all__ = [
    "Mesh",
    "PoissonSolution",
    "ElasticitySolution",
    "ModalSolution",
    "ModeShape",
    "TopOptSolution",
    "HeatHistory",
    "solve_poisson",
    "solve_transient_heat",
    "j2_uniaxial_response",
    "solve_darcy",
    "laminate_abd",
    "solve_thermal_structural",
    "contact_pairs",
    "contact_augmented_lagrangian",
    "fsi_interface_loads",
    "gpu_enabled",
    "gpu_adapter",
    "gpu_solve_cg",
    "newmark",
    "solve_stokes",
    "neo_hookean_uniaxial",
    "solve_elasticity",
    "solve_modal",
    "topopt_cantilever",
]

Vector3 = Sequence[float]
IntVector3 = Sequence[int]


class Mesh:
    """A finite-element mesh (nodes + elements)."""

    @staticmethod
    def load(path: str) -> "Mesh":
        """Load a Gmsh `.msh` (v4.1) file into a mesh."""

    @staticmethod
    def box_mesh(min: Vector3, max: Vector3, n: IntVector3) -> "Mesh":
        """Build a structured box mesh of ``[min, max]`` with ``n`` cells per axis."""

    def node_count(self) -> int:
        """Number of nodes in the mesh."""

    def coords(self, i: int) -> list[float]:
        """Coordinates of node ``i``."""

    def nodes_on_plane(self, axis: int, coord: float, tol: float) -> list[int]:
        """Node ids whose ``axis`` coordinate is within ``tol`` of ``coord``."""

    def nodes_in_box(self, min: Vector3, max: Vector3) -> list[int]:
        """Node ids within the axis-aligned box ``[min, max]``."""

    def write_vtk(
        self, path: str, field_name: str = "u", values: Sequence[float] | None = None
    ) -> None:
        """Write the mesh (with an optional per-node scalar field) to a ``.vtk`` file."""


class PoissonSolution:
    """Result of :func:`solve_poisson`: a scalar field on the mesh."""

    @property
    def mesh(self) -> Mesh:
        """The mesh this solution lives on."""

    @property
    def values(self) -> list[float]:
        """Nodal solution values (one per mesh node)."""

    def __len__(self) -> int:
        """Number of nodal values."""

    def __iter__(self) -> Any:
        """Iterate over the nodal values."""

    def to_numpy(self) -> Any:
        """The solution as an ``np.ndarray`` of shape ``(n_nodes,)``."""

    def to_pyvista(self) -> Any:
        """A ``pyvista.UnstructuredGrid`` with the field attached as point data ``"u"``."""


class ElasticitySolution:
    """Result of :func:`solve_elasticity`: a vector displacement field on the mesh."""

    @property
    def mesh(self) -> Mesh:
        """The mesh this solution lives on."""

    @property
    def values(self) -> list[float]:
        """Nodal displacement values, flat ``node_count * dim`` layout."""

    @property
    def dim(self) -> int:
        """Spatial dimension of the displacement field."""

    def __len__(self) -> int:
        """Number of scalar components in the field."""

    def __iter__(self) -> Any:
        """Iterate over the flat displacement values."""

    def to_numpy(self) -> Any:
        """The field as an ``np.ndarray`` of shape ``(n_nodes, dim)``."""

    def to_pyvista(self) -> Any:
        """A ``pyvista.UnstructuredGrid`` with the field attached as point data ``"disp"``."""


class ModeShape:
    """A single eigenmode of :func:`solve_modal`: its squared frequency and shape."""

    @property
    def mesh(self) -> Mesh:
        """The mesh this shape lives on."""

    @property
    def omega2(self) -> float:
        """Squared natural frequency ω²."""

    @property
    def omega(self) -> float:
        """Natural frequency ω = √(ω²)."""

    @property
    def shape(self) -> list[float]:
        """Mode shape, flat ``node_count * dim`` layout."""

    @property
    def dim(self) -> int:
        """Spatial dimension of the mode shape."""

    def to_numpy(self) -> Any:
        """The mode shape as an ``np.ndarray`` of shape ``(n_nodes, dim)``."""

    def to_pyvista(self) -> Any:
        """A ``pyvista.UnstructuredGrid`` with the mode attached as point data ``"mode"``."""


class ModalSolution:
    """Result of :func:`solve_modal`: the natural-vibration eigenproblem."""

    @property
    def mesh(self) -> Mesh:
        """The mesh these modes live on."""

    @property
    def dim(self) -> int:
        """Spatial dimension of the mode shapes."""

    def omega2s(self) -> list[float]:
        """Squared natural frequencies ω², in ascending order."""

    def frequencies(self) -> list[float]:
        """Natural frequencies ω = √(ω²), in ascending order."""

    def __len__(self) -> int:
        """Number of extracted modes."""

    def __getitem__(self, i: int) -> ModeShape:
        """The ``i``-th :class:`ModeShape`."""

    def __iter__(self) -> Any:
        """Iterate over the :class:`ModeShape` objects."""


class HeatHistory:
    """Nodal temperature at every step of :func:`solve_transient_heat`."""

    @property
    def mesh(self) -> Mesh:
        """The mesh this history lives on."""

    @property
    def times(self) -> list[float]:
        """Time of each stored step (``nsteps + 1`` entries, starting at 0)."""

    def __len__(self) -> int: ...
    def __getitem__(self, i: int) -> PoissonSolution:
        """The temperature field at step ``i``."""

    def to_numpy(self) -> Any:
        """``(nsteps + 1, n_nodes)`` array."""


def solve_transient_heat(
    mesh: Mesh,
    conductivity: float,
    rho_c: float,
    dt: float,
    nsteps: int,
    initial: float,
    bcs: Sequence[tuple[int, float]],
    source: float = 0.0,
    theta: float = 1.0,
    quad_order: int = 2,
) -> HeatHistory:
    """Transient heat conduction ``rho_c dT/dt - div(k grad T) = source``."""


def solve_thermal_structural(
    mesh: Mesh,
    model: str,
    young: float,
    poisson: float,
    alpha: float,
    delta_t: Sequence[float],
    bcs: Sequence[tuple[int, int, float]],
) -> ElasticitySolution:
    """Free thermal expansion under a per-node temperature rise ``delta_t``."""


def gpu_enabled() -> bool:
    """``True`` if built with GPU support (``--features gpu``)."""


def gpu_adapter() -> str:
    """Name and backend of the GPU adapter, e.g. ``"NVIDIA ... (Vulkan)"``."""


def gpu_solve_cg(
    triplets: Sequence[tuple[int, int, float]],
    rhs: Sequence[float],
    tol: float = 1e-10,
) -> tuple[list[float], int, float]:
    """GPU Jacobi-PCG for SPD ``A x = b``; returns ``(x, iterations, rel_residual)``."""


def fsi_interface_loads(
    structure: Mesh,
    interface: Sequence[tuple[int, int]],
    fluid_pressure: Sequence[float],
) -> list[float]:
    """Consistent structure load vector from interface pressures.

    ``interface`` pairs ``(structure_node, fluid_node)``; the result has
    ``node * dim + component`` ordering."""


def contact_augmented_lagrangian(
    stiffness: Sequence[Sequence[float]],
    load: Sequence[float],
    constraints: Sequence[tuple[int, float]],
    penalty: float = 1e4,
    max_iter: int = 50,
    tol: float = 1e-9,
) -> tuple[list[float], list[float]]:
    """Solve ``K u = f`` with ``u[dof] >= lower`` constraints; returns ``(u, lambda)``."""


def contact_pairs(
    a: Sequence[tuple[int, Sequence[float]]],
    b: Sequence[tuple[int, Sequence[float]]],
) -> list[tuple[int, tuple[int, float] | None]]:
    """Nearest point of ``b`` for every point of ``a``: ``(a_id, (index_in_b, distance))``
    or ``(a_id, None)`` if ``b`` is empty. Octree-accelerated."""


def laminate_abd(
    plies: Sequence[tuple[float, float, float, float, float, float]],
) -> list[list[float]]:
    """ABD matrix (6x6) from a bottom-to-top stack of
    ``(e1, e2, nu12, g12, thickness, angle_deg)`` plies."""


def newmark(
    mass: Sequence[Sequence[float]],
    damping: Sequence[Sequence[float]],
    stiffness: Sequence[Sequence[float]],
    u0: Sequence[float],
    v0: Sequence[float],
    load: Sequence[float] | Callable[[float], Sequence[float]],
    dt: float,
    nsteps: int,
    beta: float = 0.25,
    gamma: float = 0.5,
) -> list[tuple[float, list[float]]]:
    """Newmark-beta integration of ``M u'' + C u' + K u = f(t)``; returns
    ``(t, u)`` for steps ``0..=nsteps``."""


def solve_darcy(
    mesh: Mesh, permeability: float, bcs: Sequence[tuple[int, float]]
) -> PoissonSolution:
    """Steady Darcy flow ``-div(k grad p) = 0``; returns the nodal pressure."""


def solve_stokes(
    mesh: Mesh,
    viscosity: float,
    body_force: Sequence[float],
    bcs: Sequence[tuple[int, int, float]],
    penalty: float = 1e6,
) -> tuple[ElasticitySolution, PoissonSolution]:
    """Steady Stokes flow (penalty method); returns ``(velocity, pressure)``."""


def j2_uniaxial_response(
    young: float,
    poisson: float,
    yield_stress: float,
    strains: Sequence[float],
    iso_hardening: float = 0.0,
    kin_hardening: float = 0.0,
) -> list[float]:
    """Axial stress of a J2 elastic-plastic material along a monotonic strain path."""


def neo_hookean_uniaxial(mu: float, stretches: Sequence[float]) -> list[float]:
    """Nominal stress ``mu (l - l**-2)`` of an incompressible neo-Hookean solid."""


def solve_poisson(
    mesh: Mesh,
    conductivity: float,
    quad_order: int,
    source: float | Callable[[float, float, float], float],
    bcs: Sequence[tuple[int, float]],
) -> PoissonSolution:
    """Solve the steady Poisson problem ``-∇·(k ∇u) = f`` on ``mesh``."""


def solve_elasticity(
    mesh: Mesh,
    model: str,
    young: float,
    poisson: float,
    quad_order: int,
    bcs: Sequence[tuple[int, int, float]],
) -> ElasticitySolution:
    """Solve a linear-elasticity (static) problem ``K u = 0`` on ``mesh``."""


def solve_modal(
    mesh: Mesh,
    model: str,
    young: float,
    poisson: float,
    density: float,
    quad_order: int,
    num_modes: int,
    bcs: Sequence[tuple[int, int, float]],
) -> ModalSolution:
    """Solve the natural-vibration eigenproblem ``K φ = ω² M φ`` on ``mesh``."""


class TopOptSolution:
    """Result of :func:`topopt_cantilever`: optimized element densities."""

    @property
    def nx(self) -> int:
        """Number of elements along ``x``."""

    @property
    def ny(self) -> int:
        """Number of elements along ``y``."""

    @property
    def densities(self) -> list[float]:
        """Final element densities, row-major (``ny`` rows of ``nx`` columns)."""

    @property
    def compliance(self) -> list[float]:
        """Compliance at each iteration (index 0 is the uniform start)."""

    @property
    def iterations(self) -> int:
        """Number of optimality-criteria iterations performed."""

    def to_numpy(self) -> Any:
        """Densities as a ``(ny, nx)`` ``numpy.ndarray``."""


def topopt_cantilever(
    nx: int,
    ny: int,
    vol_frac: float,
    penal: float = 3.0,
    filter_radius: float = 1.5,
    max_iter: int = 50,
) -> TopOptSolution:
    """SIMP minimum-compliance optimization of a 2-D cantilever (clamped left
    edge, unit downward load at the bottom-right corner)."""
