//! Python bindings for the `tpt-fem` core (maturin-based, dev-only this pass).
//!
//! Exposes a `Mesh` class (`load` / `box_mesh` / `coords` / `nodes_on_plane` /
//! `nodes_in_box` / `write_vtk`) and solver functions: `solve_poisson`
//! (steady heat conduction), `solve_elasticity` (linear statics), and
//! `solve_modal` (natural-vibration eigenproblem `K φ = ω² M φ`). The Poisson
//! source may be a constant `float` or a Python callable `f(x, y, z)`; errors
//! from the core crates are surfaced as Python exceptions via their `Display`
//! impls.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use ::tpt_fem::{
    box_mesh as rs_box_mesh, solve_elasticity as rs_solve_elasticity,
    solve_modal as rs_solve_modal, solve_poisson as rs_solve_poisson, write_vtk_with_data,
    CellType, ElasticModel, Mesh as RsMesh, PointData,
};

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
        let inner = RsMesh::from_msh_bytes(&bytes)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
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
    fn write_vtk(
        &self,
        path: &str,
        field_name: &str,
        values: Option<Vec<f64>>,
    ) -> PyResult<()> {
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
#[pyfunction]
#[pyo3(signature = (mesh, conductivity, quad_order, source, bcs))]
fn solve_poisson(
    py: Python<'_>,
    mesh: &Mesh,
    conductivity: f64,
    quad_order: usize,
    source: &Bound<'_, PyAny>,
    bcs: Vec<(usize, f64)>,
) -> PyResult<Vec<f64>> {
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
                            x.get(0).copied().unwrap_or(0.0),
                            x.get(1).copied().unwrap_or(0.0),
                            x.get(2).copied().unwrap_or(0.0),
                        );
                        match cb.bind(py).call1(args) {
                            Ok(v) => match v.extract::<f64>() {
                                Ok(f) => f,
                                Err(e) => {
                                    *callback_error.lock().unwrap() = Some(e);
                                    0.0
                                }
                            },
                            Err(e) => {
                                *callback_error.lock().unwrap() = Some(e);
                                0.0
                            }
                        }
                    })
                } else {
                    0.0
                }
            };
            rs_solve_poisson(&mesh.inner, conductivity, quad_order, f, &bcs, None, None)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))
        }
    });
    if result.is_ok() {
        if let Some(e) = callback_error.lock().unwrap().take() {
            return Err(e);
        }
    }
    result
}

/// Reference (spatial) dimension of a mesh's first cell.
fn dim_of(mesh: &RsMesh) -> PyResult<usize> {
    let cell = mesh.elements.first().map(|e| e.cell_type);
    match cell {
        Some(CellType::Line) => Ok(1),
        Some(CellType::Tri | CellType::Quad | CellType::Tri6 | CellType::Quad8 | CellType::Quad9) => {
            Ok(2)
        }
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
#[pyfunction]
#[pyo3(signature = (mesh, model, young, poisson, quad_order, bcs))]
fn solve_elasticity(
    py: Python<'_>,
    mesh: &Mesh,
    model: &str,
    young: f64,
    poisson: f64,
    quad_order: usize,
    bcs: Vec<(usize, usize, f64)>,
) -> PyResult<Vec<f64>> {
    let model = parse_model(model)?;
    let dim = dim_of(&mesh.inner)?;
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    py.allow_threads(move || {
        rs_solve_elasticity(
            &mesh.inner,
            model,
            young,
            poisson,
            quad_order,
            |_| vec![0.0; dim],
            &dir,
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
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
#[pyfunction]
#[pyo3(signature = (mesh, model, young, poisson, density, quad_order, num_modes, bcs))]
fn solve_modal(
    py: Python<'_>,
    mesh: &Mesh,
    model: &str,
    young: f64,
    poisson: f64,
    density: f64,
    quad_order: usize,
    num_modes: usize,
    bcs: Vec<(usize, usize, f64)>,
) -> PyResult<Vec<(f64, Vec<f64>)>> {
    let model = parse_model(model)?;
    let dim = dim_of(&mesh.inner)?;
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    py.allow_threads(move || {
        rs_solve_modal(&mesh.inner, model, young, poisson, density, quad_order, num_modes, &dir)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    })
}

#[pymodule]
fn tpt_fem(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Mesh>()?;
    m.add_function(wrap_pyfunction!(solve_poisson, py)?)?;
    m.add_function(wrap_pyfunction!(solve_elasticity, py)?)?;
    m.add_function(wrap_pyfunction!(solve_modal, py)?)?;
    Ok(())
}
