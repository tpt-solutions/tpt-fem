//! Python bindings for the `tpt-fem` core (maturin-based, dev-only this pass).
//!
//! Exposes a `Mesh` class (`load` / `box_mesh` / `coords` / `nodes_on_plane` /
//! `nodes_in_box` / `write_vtk`) and a `solve_poisson` function accepting either
//! a constant volumetric source or a Python callable `f(x, y, z)`. Errors from
//! the core crates are surfaced as Python exceptions via their `Display` impls.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use tpt_fem::{
    box_mesh as rs_box_mesh, solve_poisson as rs_solve_poisson, write_vtk_with_data, CellType,
    Mesh as RsMesh, MeshBuilder, PointData,
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
        Some(source.into_py(py))
    } else {
        None
    };

    // Run the (potentially GIL-unaware) solve with the GIL released; the
    // Python callback re-acquires the GIL per call via `with_gil`.
    py.allow_threads(move || {
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
                        Ok(v) => v.extract::<f64>().unwrap_or(0.0),
                        Err(_) => 0.0,
                    }
                })
            } else {
                0.0
            }
        };
        rs_solve_poisson(&mesh.inner, conductivity, quad_order, f, &bcs, None, None)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    })
}

#[pymodule]
fn tpt_fem_py(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Mesh>()?;
    m.add_function(wrap_pyfunction!(solve_poisson, py)?)?;
    Ok(())
}
