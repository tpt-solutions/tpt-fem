//! Python bindings for the `tpt-fem` core (maturin-based, dev-only this pass).
//!
//! Exposes a `Mesh` class (`load` / `box_mesh` / `coords` / `nodes_on_plane` /
//! `nodes_in_box` / `write_vtk`) and solver functions: `solve_poisson`
//! (steady heat conduction), `solve_elasticity` (linear statics), and
//! `solve_modal` (natural-vibration eigenproblem `K φ = ω² M φ`). The Poisson
//! source may be a constant `float` or a Python callable `f(x, y, z)`; errors
//! from the core crates are surfaced as Python exceptions via their `Display`
//! impls.
//!
//! Each solver returns a Jupyter-friendly result object
//! (`PoissonSolution` / `ElasticitySolution` / `ModalSolution` + `ModeShape`)
//! rather than a bare list: it carries rich `__repr__` / `_repr_html_` display,
//! a `to_numpy()` accessor (returns an `np.ndarray`), and a `to_pyvista()`
//! accessor (returns a `pyvista.UnstructuredGrid`), so results plug straight
//! into PyVista / matplotlib / Jupyter without a manual VTK export-and-reimport
//! round-trip. `numpy` is required for `to_numpy`; `pyvista` for `to_pyvista`.

// pyo3 0.23 deprecated `IntoPy::into_py`; the `&Mesh` -> `Py<Mesh>` path we use
// for result objects is the correct one and has no non-deprecated equivalent
// here, so the single warning is intentionally allowed.
#![allow(deprecated)]

use ::tpt_fem::{
    box_mesh as rs_box_mesh, solve_elasticity as rs_solve_elasticity,
    solve_modal as rs_solve_modal, solve_poisson as rs_solve_poisson, write_vtk_with_data,
    CellType, ElasticModel, Mesh as RsMesh, PointData,
};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::sync::PoisonError;
use pyo3::types::{PyDict, PyList};

#[pyclass]
struct Mesh {
    inner: RsMesh,
}

#[pymethods]
impl Mesh {
    /// Load a Gmsh `.msh` (v4.1) file into a mesh.
    #[staticmethod]
    fn load(path: &str) -> PyResult<Mesh> {
        let bytes = std::fs::read(path)
            .map_err(|e| PyRuntimeError::new_err(format!("read {path}: {e}")))?;
        let inner =
            RsMesh::from_msh_bytes(&bytes).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(Mesh { inner })
    }

    /// Build a structured tetrahedral box mesh of `[min, max]` with `n` cells
    /// per axis.
    #[staticmethod]
    fn box_mesh(min: [f64; 3], max: [f64; 3], n: [usize; 3]) -> Mesh {
        Mesh {
            inner: rs_box_mesh(min, max, n),
        }
    }

    /// Number of nodes in the mesh.
    fn node_count(&self) -> usize {
        self.inner.node_count()
    }

    /// Coordinates of node `i`.
    fn coords(&self, i: usize) -> Vec<f64> {
        self.inner.node_coords(i).to_vec()
    }

    /// Node ids whose `axis` coordinate is within `tol` of `coord`.
    fn nodes_on_plane(&self, axis: usize, coord: f64, tol: f64) -> Vec<usize> {
        self.inner.nodes_on_plane(axis, coord, tol)
    }

    /// Node ids within the axis-aligned box `[min, max]`.
    fn nodes_in_box(&self, min: [f64; 3], max: [f64; 3]) -> Vec<usize> {
        self.inner.nodes_in_box(min, max)
    }

    /// Write the mesh (with an optional per-node scalar field) to a ParaView
    /// `.vtk` file.
    #[pyo3(signature = (path, field_name="u", values=None))]
    fn write_vtk(&self, path: &str, field_name: &str, values: Option<Vec<f64>>) -> PyResult<()> {
        let data = match values {
            Some(v) => vec![PointData::new(field_name, v)],
            None => vec![],
        };
        write_vtk_with_data(&self.inner, &data, path)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }
}

/// Solve the steady Poisson problem `-∇·(k ∇u) = f` on `mesh`.
///
/// * `conductivity` — constant `k`.
/// * `quad_order` — quadrature order.
/// * `source` — either a `float` (constant `f`) or a callable `f(x, y, z)`.
/// * `bcs` — list of `(node_id, value)` Dirichlet conditions.
///
/// Returns a [`PoissonSolution`] (a Jupyter-friendly result object with numpy
/// / pyvista interop and rich display), not a bare `list`.
#[pyfunction]
#[pyo3(signature = (mesh, conductivity, quad_order, source, bcs))]
fn solve_poisson(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    conductivity: f64,
    quad_order: usize,
    source: &Bound<'_, PyAny>,
    bcs: Vec<(usize, f64)>,
) -> PyResult<PoissonSolution> {
    let constant = source.extract::<f64>().ok();
    let callback = if constant.is_none() {
        Some(source.clone().unbind())
    } else {
        None
    };
    // Captures a Python exception raised by the user callback so it can be
    // surfaced as a real Python error after the (GIL-released) solve returns,
    // instead of silently falling back to 0.0. `Arc<Mutex<_>>` (not
    // `Rc<RefCell<_>>`) is required so the closure is `Send` across the
    // `allow_threads` boundary.
    let callback_error: std::sync::Arc<std::sync::Mutex<Option<PyErr>>> =
        std::sync::Arc::new(std::sync::Mutex::new(None));

    // Run the (potentially GIL-unaware) solve with the GIL released; the
    // Python callback re-acquires the GIL per call via `with_gil`.
    let rs_mesh = mesh.borrow().inner.clone();
    let result = py.allow_threads({
        let callback_error = callback_error.clone();
        move || {
            let f = move |x: &[f64]| -> f64 {
                if let Some(c) = constant {
                    return c;
                }
                if let Some(cb) = &callback {
                    Python::with_gil(|py| {
                        let args = (
                            x.first().copied().unwrap_or(0.0),
                            x.get(1).copied().unwrap_or(0.0),
                            x.get(2).copied().unwrap_or(0.0),
                        );
                        match cb.bind(py).call1(args) {
                            Ok(v) => match v.extract::<f64>() {
                                Ok(f) => f,
                                Err(e) => {
                                    *callback_error
                                        .lock()
                                        .unwrap_or_else(PoisonError::into_inner) = Some(e);
                                    0.0
                                }
                            },
                            Err(e) => {
                                *callback_error
                                    .lock()
                                    .unwrap_or_else(PoisonError::into_inner) = Some(e);
                                0.0
                            }
                        }
                    })
                } else {
                    0.0
                }
            };
            rs_solve_poisson(&rs_mesh, conductivity, quad_order, f, &bcs, None, None)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))
        }
    });
    if result.is_ok() {
        if let Some(e) = callback_error
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .take()
        {
            return Err(e);
        }
    }
    let values = result?;
    Ok(PoissonSolution {
        mesh: mesh.unbind(),
        values,
    })
}

/// Reference (spatial) dimension of a mesh's first cell.
fn dim_of(mesh: &RsMesh) -> PyResult<usize> {
    let cell = mesh.elements.first().map(|e| e.cell_type);
    match cell {
        Some(CellType::Line) => Ok(1),
        Some(
            CellType::Tri | CellType::Quad | CellType::Tri6 | CellType::Quad8 | CellType::Quad9,
        ) => Ok(2),
        Some(
            CellType::Tet | CellType::Hex | CellType::Tet10 | CellType::Hex20 | CellType::Hex27,
        ) => Ok(3),
        None => Err(PyRuntimeError::new_err("mesh has no elements")),
    }
}

/// Parse an elasticity-model string into [`ElasticModel`].
fn parse_model(s: &str) -> PyResult<ElasticModel> {
    match s.to_ascii_lowercase().replace(['-', '_'], "").as_str() {
        "bar" | "baraxial" => Ok(ElasticModel::BarAxial),
        "planestress" => Ok(ElasticModel::PlaneStress),
        "planestrain" => Ok(ElasticModel::PlaneStrain),
        "3d" | "continuum" | "continuum3d" => Ok(ElasticModel::Continuum3D),
        other => Err(PyRuntimeError::new_err(format!(
            "unknown elasticity model '{other}' (bar | plane-stress | plane-strain | 3d)"
        ))),
    }
}

/// Solve a linear-elasticity (static) problem `K u = 0` on `mesh`.
///
/// * `model` — `"bar"`, `"plane-stress"`, `"plane-strain"`, or `"3d"`.
/// * `young` / `poisson` — material constants.
/// * `quad_order` — quadrature order.
/// * `bcs` — list of `(node_id, component, value)` Dirichlet conditions (the
///   global DOF is `node_id * dim + component`).
///
/// Returns an [`ElasticitySolution`] (a Jupyter-friendly result object with
/// numpy / pyvista interop and rich display), not a bare `list`.
#[pyfunction]
#[pyo3(signature = (mesh, model, young, poisson, quad_order, bcs))]
fn solve_elasticity(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    model: &str,
    young: f64,
    poisson: f64,
    quad_order: usize,
    bcs: Vec<(usize, usize, f64)>,
) -> PyResult<ElasticitySolution> {
    let model = parse_model(model)?;
    let dim = dim_of(&mesh.borrow().inner)?;
    let rs_mesh = mesh.borrow().inner.clone();
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    let values = py.allow_threads(move || {
        rs_solve_elasticity(
            &rs_mesh,
            model,
            young,
            poisson,
            quad_order,
            |_| vec![0.0; dim],
            &dir,
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    })?;
    Ok(ElasticitySolution {
        mesh: mesh.unbind(),
        values,
        dim,
    })
}

/// Solve the natural-vibration eigenproblem `K φ = ω² M φ` on `mesh`.
///
/// * `model` — as for [`solve_elasticity`].
/// * `young` / `poisson` / `density` — material constants.
/// * `quad_order` — quadrature order.
/// * `num_modes` — number of modes to extract.
/// * `bcs` — list of `(node_id, component, value)` Dirichlet conditions (the
///   constrained DOFs are removed from both `K` and `M`).
///
/// Returns a list of `(ω², φ)` pairs: the squared natural frequency and its
/// mode shape (a `node_count * dim` vector, zero on fixed DOFs).
///
/// Returns a [`ModalSolution`] (a Jupyter-friendly result object with numpy /
/// pyvista interop and rich display) whose elements are [`ModeShape`] objects,
/// not a bare list of `(ω², shape)` tuples.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (mesh, model, young, poisson, density, quad_order, num_modes, bcs))]
fn solve_modal(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    model: &str,
    young: f64,
    poisson: f64,
    density: f64,
    quad_order: usize,
    num_modes: usize,
    bcs: Vec<(usize, usize, f64)>,
) -> PyResult<ModalSolution> {
    let model = parse_model(model)?;
    let dim = dim_of(&mesh.borrow().inner)?;
    let rs_mesh = mesh.borrow().inner.clone();
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    let modes = py.allow_threads(move || {
        rs_solve_modal(
            &rs_mesh, model, young, poisson, density, quad_order, num_modes, &dir,
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    })?;
    let (omega2s, shapes): (Vec<f64>, Vec<Vec<f64>>) = modes.into_iter().unzip();
    Ok(ModalSolution {
        mesh: mesh.unbind(),
        dim,
        omega2s,
        shapes,
    })
}

// ---------------------------------------------------------------------------
// Result objects (Jupyter-friendly: numpy / pyvista interop + rich display).
// ---------------------------------------------------------------------------

/// Borrow the inner `tpt-fem` mesh out of the `Py<Mesh>` stored on a result
/// object (kept as `Py<Mesh>` rather than relying on the feature-gated
/// `Py<T>: Clone`).
fn borrow_rs_mesh(mesh: &Py<Mesh>, py: Python<'_>) -> RsMesh {
    mesh.borrow(py).inner.clone()
}

/// Return `values` as an `np.ndarray` reshaped to `shape`.
fn to_numpy_array(py: Python<'_>, values: &[f64], shape: &[usize]) -> PyResult<PyObject> {
    let np = py.import("numpy").map_err(|_| {
        PyRuntimeError::new_err("numpy is not installed; run `pip install numpy` to use to_numpy()")
    })?;
    let arr = np.call_method1("array", (values.to_vec(),))?;
    Ok(arr.call_method1("reshape", (shape.to_vec(),))?.into())
}

/// Write `mesh` to a temp `.vtk`, read it back with `pyvista`, and attach each
/// `(name, flat_values, ncomp)` field as `ncomp`-wide point data. This reuses
/// the crate's own (well-tested) VTK writer so it works for every cell type the
/// core supports, and avoids fragile hand-rolled pyvista grid construction.
fn to_pyvista_grid(
    py: Python<'_>,
    mesh: &RsMesh,
    fields: &[(&str, Vec<f64>, usize)],
) -> PyResult<PyObject> {
    let pv = py.import("pyvista").map_err(|_| {
        PyRuntimeError::new_err(
            "pyvista is not installed; run `pip install pyvista` to use to_pyvista()",
        )
    })?;
    let _np = py.import("numpy").map_err(|_| {
        PyRuntimeError::new_err(
            "numpy is not installed; run `pip install numpy` to use to_pyvista()",
        )
    })?;
    let tmp = py.import("tempfile")?;
    let kwargs = PyDict::new(py);
    kwargs.set_item("delete", false)?;
    kwargs.set_item("suffix", ".vtk")?;
    let fh = tmp.call_method("NamedTemporaryFile", (), Some(&kwargs))?;
    let path: String = fh.getattr("name")?.extract()?;
    write_vtk_with_data(mesh, &[], &path).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let grid = pv.call_method1("read", (path.clone(),))?;
    let n = mesh.node_count();
    for (name, values, ncomp) in fields {
        let arr = if *ncomp == 1 {
            _np.call_method1("array", (values.to_vec(),))?
        } else {
            let flat = _np.call_method1("array", (values.to_vec(),))?;
            flat.call_method1("reshape", (n, *ncomp))?
        };
        grid.getattr("point_data")?.set_item(*name, arr)?;
    }
    fh.call_method0("close")?;
    let _ = std::fs::remove_file(&path);
    Ok(grid.into())
}

/// `(min, max, mean)` over a scalar field.
fn field_stats(values: &[f64]) -> (f64, f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for &v in values {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    let mean = if values.is_empty() {
        0.0
    } else {
        sum / values.len() as f64
    };
    (min, max, mean)
}

/// Largest L2 magnitude across the nodal vectors of a `n*dim` field.
fn max_vector_magnitude(values: &[f64], dim: usize) -> f64 {
    let mut m = 0.0;
    for c in values.chunks(dim) {
        let mut s = 0.0;
        for &v in c {
            s += v * v;
        }
        let mag = s.sqrt();
        if mag > m {
            m = mag;
        }
    }
    m
}

/// Result of [`solve_poisson`]: a scalar field on the mesh.
///
/// Rich display in Jupyter via `__repr__` / `_repr_html_`; `to_numpy()` returns
/// an `(n_nodes,)` array and `to_pyvista()` a `pyvista.UnstructuredGrid`.
#[pyclass]
struct PoissonSolution {
    mesh: Py<Mesh>,
    values: Vec<f64>,
}

#[pymethods]
impl PoissonSolution {
    /// The mesh this solution lives on.
    #[getter]
    fn mesh(&self, py: Python<'_>) -> Py<Mesh> {
        self.mesh.clone_ref(py)
    }

    /// Nodal solution values (one per mesh node).
    #[getter]
    fn values(&self) -> Vec<f64> {
        self.values.clone()
    }

    /// Number of nodal values.
    fn __len__(&self) -> usize {
        self.values.len()
    }

    /// Iterate over the nodal values.
    fn __iter__(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(PyList::new(py, self.values.clone())?
            .call_method0("__iter__")?
            .into())
    }

    /// The solution as an `np.ndarray` of shape `(n_nodes,)`.
    fn to_numpy(&self, py: Python<'_>) -> PyResult<PyObject> {
        to_numpy_array(py, &self.values, &[self.values.len()])
    }

    /// A `pyvista.UnstructuredGrid` with the field attached as point data `"u"`.
    fn to_pyvista(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = borrow_rs_mesh(&self.mesh, py);
        to_pyvista_grid(py, &m, &[("u", self.values.clone(), 1)])
    }

    fn __repr__(&self) -> String {
        let (min, max, mean) = field_stats(&self.values);
        format!(
            "PoissonSolution(nodes={}, min={:.4e}, max={:.4e}, mean={:.4e})",
            self.values.len(),
            min,
            max,
            mean
        )
    }

    fn _repr_html_(&self) -> String {
        let (min, max, mean) = field_stats(&self.values);
        format!(
            "<table><tr><th colspan=\"2\">PoissonSolution</th></tr>\
             <tr><td>nodes</td><td>{}</td></tr>\
             <tr><td>min</td><td>{:.4e}</td></tr>\
             <tr><td>max</td><td>{:.4e}</td></tr>\
             <tr><td>mean</td><td>{:.4e}</td></tr></table>",
            self.values.len(),
            min,
            max,
            mean
        )
    }
}

/// Result of [`solve_elasticity`]: a vector displacement field on the mesh.
///
/// `values` is the flat `node_count * dim` field; `to_numpy()` reshapes it to
/// `(n_nodes, dim)`, and `to_pyvista()` attaches it as point data `"disp"`.
#[pyclass]
struct ElasticitySolution {
    mesh: Py<Mesh>,
    values: Vec<f64>,
    dim: usize,
}

#[pymethods]
impl ElasticitySolution {
    /// The mesh this solution lives on.
    #[getter]
    fn mesh(&self, py: Python<'_>) -> Py<Mesh> {
        self.mesh.clone_ref(py)
    }

    /// Nodal displacement values, flat `node_count * dim` layout.
    #[getter]
    fn values(&self) -> Vec<f64> {
        self.values.clone()
    }

    /// Spatial dimension of the displacement field.
    #[getter]
    fn dim(&self) -> usize {
        self.dim
    }

    /// Number of scalar components in the field.
    fn __len__(&self) -> usize {
        self.values.len()
    }

    /// Iterate over the flat displacement values.
    fn __iter__(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(PyList::new(py, self.values.clone())?
            .call_method0("__iter__")?
            .into())
    }

    /// The field as an `np.ndarray` of shape `(n_nodes, dim)`.
    fn to_numpy(&self, py: Python<'_>) -> PyResult<PyObject> {
        let n = self.values.len() / self.dim;
        to_numpy_array(py, &self.values, &[n, self.dim])
    }

    /// A `pyvista.UnstructuredGrid` with the field attached as point data
    /// `"disp"` (a vector field).
    fn to_pyvista(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = borrow_rs_mesh(&self.mesh, py);
        to_pyvista_grid(py, &m, &[("disp", self.values.clone(), self.dim)])
    }

    fn __repr__(&self) -> String {
        let mag = max_vector_magnitude(&self.values, self.dim);
        format!(
            "ElasticitySolution(nodes={}, dim={}, max|u|={:.4e})",
            self.values.len() / self.dim,
            self.dim,
            mag
        )
    }

    fn _repr_html_(&self) -> String {
        let mag = max_vector_magnitude(&self.values, self.dim);
        let n = self.values.len() / self.dim;
        format!(
            "<table><tr><th colspan=\"2\">ElasticitySolution</th></tr>\
             <tr><td>nodes</td><td>{}</td></tr>\
             <tr><td>dim</td><td>{}</td></tr>\
             <tr><td>max |u|</td><td>{:.4e}</td></tr></table>",
            n, self.dim, mag
        )
    }
}

/// A single eigenmode of [`solve_modal`]: its squared frequency and shape.
///
/// `shape` is the flat `node_count * dim` mode vector; `to_numpy()` reshapes
/// it to `(n_nodes, dim)`, and `to_pyvista()` attaches it as point data
/// `"mode"` (a vector field).
#[pyclass]
struct ModeShape {
    mesh: Py<Mesh>,
    dim: usize,
    omega2: f64,
    shape: Vec<f64>,
}

#[pymethods]
impl ModeShape {
    /// The mesh this shape lives on.
    #[getter]
    fn mesh(&self, py: Python<'_>) -> Py<Mesh> {
        self.mesh.clone_ref(py)
    }

    /// Squared natural frequency ω².
    #[getter]
    fn omega2(&self) -> f64 {
        self.omega2
    }

    /// Natural frequency ω = √(ω²).
    #[getter]
    fn omega(&self) -> f64 {
        self.omega2.sqrt()
    }

    /// Mode shape, flat `node_count * dim` layout.
    #[getter]
    fn shape(&self) -> Vec<f64> {
        self.shape.clone()
    }

    /// Spatial dimension of the mode shape.
    #[getter]
    fn dim(&self) -> usize {
        self.dim
    }

    /// The mode shape as an `np.ndarray` of shape `(n_nodes, dim)`.
    fn to_numpy(&self, py: Python<'_>) -> PyResult<PyObject> {
        let n = self.shape.len() / self.dim;
        to_numpy_array(py, &self.shape, &[n, self.dim])
    }

    /// A `pyvista.UnstructuredGrid` with the mode attached as point data
    /// `"mode"` (a vector field).
    fn to_pyvista(&self, py: Python<'_>) -> PyResult<PyObject> {
        let m = borrow_rs_mesh(&self.mesh, py);
        to_pyvista_grid(py, &m, &[("mode", self.shape.clone(), self.dim)])
    }

    fn __repr__(&self) -> String {
        let mag = max_vector_magnitude(&self.shape, self.dim);
        format!(
            "ModeShape(ω²={:.4e}, ω={:.4e}, max|φ|={:.4e})",
            self.omega2,
            self.omega2.sqrt(),
            mag
        )
    }

    fn _repr_html_(&self) -> String {
        let mag = max_vector_magnitude(&self.shape, self.dim);
        format!(
            "<table><tr><th colspan=\"2\">ModeShape</th></tr>\
             <tr><td>ω²</td><td>{:.4e}</td></tr>\
             <tr><td>ω</td><td>{:.4e}</td></tr>\
             <tr><td>max |φ|</td><td>{:.4e}</td></tr></table>",
            self.omega2,
            self.omega2.sqrt(),
            mag
        )
    }
}

/// Result of [`solve_modal`]: the natural-vibration eigenproblem.
///
/// Indexable / iterable over [`ModeShape`] objects; `omega2s()` and
/// `frequencies()` return the squared and unsquared frequencies respectively.
#[pyclass]
struct ModalSolution {
    mesh: Py<Mesh>,
    dim: usize,
    omega2s: Vec<f64>,
    shapes: Vec<Vec<f64>>,
}

#[pymethods]
impl ModalSolution {
    /// The mesh these modes live on.
    #[getter]
    fn mesh(&self, py: Python<'_>) -> Py<Mesh> {
        self.mesh.clone_ref(py)
    }

    /// Spatial dimension of the mode shapes.
    #[getter]
    fn dim(&self) -> usize {
        self.dim
    }

    /// Squared natural frequencies ω², in ascending order.
    fn omega2s(&self) -> Vec<f64> {
        self.omega2s.clone()
    }

    /// Natural frequencies ω = √(ω²), in ascending order.
    fn frequencies(&self) -> Vec<f64> {
        self.omega2s.iter().map(|&w| w.sqrt()).collect()
    }

    /// Number of extracted modes.
    fn __len__(&self) -> usize {
        self.omega2s.len()
    }

    /// The `i`-th [`ModeShape`].
    fn __getitem__(&self, py: Python<'_>, i: usize) -> PyResult<ModeShape> {
        if i >= self.omega2s.len() {
            return Err(PyRuntimeError::new_err(format!(
                "mode index {} out of range ({} modes)",
                i,
                self.omega2s.len()
            )));
        }
        Ok(ModeShape {
            mesh: self.mesh.clone_ref(py),
            dim: self.dim,
            omega2: self.omega2s[i],
            shape: self.shapes[i].clone(),
        })
    }

    /// Iterate over the [`ModeShape`] objects.
    fn __iter__(&self, py: Python<'_>) -> PyResult<PyObject> {
        let mut items = Vec::with_capacity(self.omega2s.len());
        for (i, &w2) in self.omega2s.iter().enumerate() {
            items.push(ModeShape {
                mesh: self.mesh.clone_ref(py),
                dim: self.dim,
                omega2: w2,
                shape: self.shapes[i].clone(),
            });
        }
        Ok(PyList::new(py, items)?.call_method0("__iter__")?.into())
    }

    fn __repr__(&self) -> String {
        let fund = self.omega2s.first().map(|&w| w.sqrt()).unwrap_or(0.0);
        format!(
            "ModalSolution(modes={}, dim={}, fundamental ω={:.4e})",
            self.omega2s.len(),
            self.dim,
            fund
        )
    }

    fn _repr_html_(&self) -> String {
        let rows: String = self
            .omega2s
            .iter()
            .enumerate()
            .map(|(i, &w)| {
                format!(
                    "<tr><td>{}</td><td>{:.4e}</td><td>{:.4e}</td></tr>",
                    i,
                    w,
                    w.sqrt()
                )
            })
            .collect();
        format!(
            "<table><tr><th>mode</th><th>ω²</th><th>ω</th></tr>{}</table>",
            rows
        )
    }
}

#[pymodule]
fn tpt_fem(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Mesh>()?;
    m.add_class::<PoissonSolution>()?;
    m.add_class::<ElasticitySolution>()?;
    m.add_class::<ModalSolution>()?;
    m.add_class::<ModeShape>()?;
    m.add_function(wrap_pyfunction!(solve_poisson, py)?)?;
    m.add_function(wrap_pyfunction!(solve_elasticity, py)?)?;
    m.add_function(wrap_pyfunction!(solve_modal, py)?)?;
    Ok(())
}
