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
    "solve_poisson",
    "solve_elasticity",
    "solve_modal",
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
