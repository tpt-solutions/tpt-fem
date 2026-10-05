//! Python bindings for the `tpt-fem` core (maturin-based, dev-only this pass).
//!
//! Exposes a `Mesh` class (`load` / `box_mesh` / `coords` / `nodes_on_plane` /
//! `nodes_in_box` / `write_vtk`) and solver functions: `solve_poisson`
//! (steady heat conduction), `solve_elasticity` (linear statics), and
//! `solve_modal` (natural-vibration eigenproblem `K φ = ω² M φ`), and
//! `topopt_cantilever` (SIMP topology optimization). The Poisson
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
    augmented_lagrangian as rs_augmented_lagrangian, contact_pairs as rs_contact_pairs,
    fsi_interface_loads as rs_fsi_interface_loads, laminate_abd as rs_laminate_abd,
    newmark as rs_newmark, solve_darcy as rs_solve_darcy, steady_stokes as rs_steady_stokes,
    thermal_structural as rs_thermal_structural, ContactConstraint, Coo, NewmarkOptions, Ply,
};
use ::tpt_fem::{
    box_mesh as rs_box_mesh, cantilever_load, solve_elasticity as rs_solve_elasticity,
    solve_modal as rs_solve_modal, solve_poisson as rs_solve_poisson,
    solve_transient_heat as rs_solve_transient_heat, topopt_simp, write_vtk_with_data, CellType,
    ElasticModel, Grid, Mesh as RsMesh, PlasticityParams, PointData, TopOptParams,
    TransientHeatOptions,
};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::sync::PoisonError;

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

/// Transient heat conduction `rho_c dT/dt - div(k grad T) = source` by the
/// theta-method.
///
/// * `mesh` — the mesh.
/// * `conductivity`, `rho_c` — constant `k` and volumetric heat capacity.
/// * `dt`, `nsteps` — time step and number of steps.
/// * `initial` — uniform initial temperature.
/// * `bcs` — list of `(node_id, value)` Dirichlet conditions.
/// * `source` — constant volumetric source (default `0.0`).
/// * `theta` — `1.0` backward Euler (default), `0.5` Crank-Nicolson.
/// * `quad_order` — quadrature order (default `2`).
///
/// Returns a [`HeatHistory`] with one nodal field per step `0..=nsteps`.
#[pyfunction]
#[pyo3(signature = (mesh, conductivity, rho_c, dt, nsteps, initial, bcs, source=0.0, theta=1.0, quad_order=2))]
#[allow(clippy::too_many_arguments)]
fn solve_transient_heat(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    conductivity: f64,
    rho_c: f64,
    dt: f64,
    nsteps: usize,
    initial: f64,
    bcs: Vec<(usize, f64)>,
    source: f64,
    theta: f64,
    quad_order: usize,
) -> PyResult<HeatHistory> {
    let rs_mesh = mesh.borrow().inner.clone();
    let opts = TransientHeatOptions {
        conductivity,
        rho_c,
        quad_order,
        dt,
        nsteps,
        theta,
    };
    let init = vec![initial; rs_mesh.node_count()];
    let hist = py
        .allow_threads(move || {
            rs_solve_transient_heat(&rs_mesh, &opts, &init, move |_, _| source, &bcs)
        })
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let (times, fields): (Vec<f64>, Vec<Vec<f64>>) = hist.into_iter().unzip();
    Ok(HeatHistory {
        mesh: mesh.unbind(),
        times,
        fields,
    })
}

/// Result of [`solve_transient_heat`]: the nodal temperature at every step.
#[pyclass]
struct HeatHistory {
    mesh: Py<Mesh>,
    times: Vec<f64>,
    fields: Vec<Vec<f64>>,
}

#[pymethods]
impl HeatHistory {
    /// The mesh this history lives on.
    #[getter]
    fn mesh(&self, py: Python<'_>) -> Py<Mesh> {
        self.mesh.clone_ref(py)
    }

    /// Time of each stored step (`nsteps + 1` entries, starting at 0).
    #[getter]
    fn times(&self) -> Vec<f64> {
        self.times.clone()
    }

    /// Number of stored steps.
    fn __len__(&self) -> usize {
        self.fields.len()
    }

    /// The nodal temperature at step `i` as a `PoissonSolution`-style field.
    fn __getitem__(&self, py: Python<'_>, i: isize) -> PyResult<PoissonSolution> {
        let n = self.fields.len() as isize;
        let idx = if i < 0 { i + n } else { i };
        if idx < 0 || idx >= n {
            return Err(pyo3::exceptions::PyIndexError::new_err(
                "step index out of range",
            ));
        }
        Ok(PoissonSolution {
            mesh: self.mesh.clone_ref(py),
            values: self.fields[idx as usize].clone(),
        })
    }

    /// The history as an `np.ndarray` of shape `(nsteps + 1, n_nodes)`.
    fn to_numpy(&self, py: Python<'_>) -> PyResult<PyObject> {
        let n = self.fields.first().map_or(0, |f| f.len());
        let flat: Vec<f64> = self.fields.iter().flatten().copied().collect();
        to_numpy_array(py, &flat, &[self.fields.len(), n])
    }

    fn __repr__(&self) -> String {
        let last = self
            .fields
            .last()
            .map(|f| field_stats(f))
            .unwrap_or((0.0, 0.0, 0.0));
        format!(
            "HeatHistory(steps={}, t_end={:.4e}, final min={:.4e}, max={:.4e})",
            self.fields.len(),
            self.times.last().copied().unwrap_or(0.0),
            last.0,
            last.1
        )
    }
}

/// Uniaxial stress-strain response of a J2 (von Mises) elastic-plastic material
/// with isotropic and kinematic hardening, for a *monotonic* strain path.
///
/// `strains` is the sequence of total axial strains (applied in order, starting
/// from the virgin state); returns the axial stress after each one. Material
/// constants: Young's modulus `young`, Poisson's ratio `poisson`, initial yield
/// stress `yield_stress`, hardening moduli `iso_hardening` / `kin_hardening`.
/// A material-point driver, useful for calibrating a model against test data.
#[pyfunction]
#[pyo3(signature = (young, poisson, yield_stress, strains, iso_hardening=0.0, kin_hardening=0.0))]
fn j2_uniaxial_response(
    young: f64,
    poisson: f64,
    yield_stress: f64,
    strains: Vec<f64>,
    iso_hardening: f64,
    kin_hardening: f64,
) -> PyResult<Vec<f64>> {
    if !(young.is_finite() && young > 0.0) || !(yield_stress.is_finite() && yield_stress > 0.0) {
        return Err(PyRuntimeError::new_err(
            "young and yield_stress must be finite and positive",
        ));
    }
    let params = PlasticityParams {
        young,
        poisson,
        yield_stress,
        iso_hardening,
        kin_hardening,
    };
    let mut eps_eq = 0.0;
    let mut out = Vec::with_capacity(strains.len());
    for &eps in &strains {
        let (sigma, new_eq, _) = ::tpt_fem::plastic_1d(&params, eps, eps_eq);
        eps_eq = new_eq;
        out.push(sigma);
    }
    Ok(out)
}

/// Nominal (first Piola-Kirchhoff) stress `P = mu (lambda - lambda^-2)` of an
/// incompressible neo-Hookean solid in uniaxial extension, for each stretch
/// `lambda > 0` in `stretches`.
#[pyfunction]
fn neo_hookean_uniaxial(mu: f64, stretches: Vec<f64>) -> PyResult<Vec<f64>> {
    if let Some(bad) = stretches.iter().find(|l| !(l.is_finite() && **l > 0.0)) {
        return Err(PyRuntimeError::new_err(format!(
            "stretches must be finite and positive (got {bad})"
        )));
    }
    Ok(stretches
        .iter()
        .map(|&l| ::tpt_fem::neo_hookean_1d(l, mu))
        .collect())
}

/// Steady Darcy flow `-div(k grad p) = 0` for the pressure field `p`.
///
/// * `permeability` — constant `k`.
/// * `bcs` — list of `(node_id, pressure)` Dirichlet conditions.
///
/// Returns a [`PoissonSolution`] holding the nodal pressure.
#[pyfunction]
fn solve_darcy(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    permeability: f64,
    bcs: Vec<(usize, f64)>,
) -> PyResult<PoissonSolution> {
    let rs_mesh = mesh.borrow().inner.clone();
    let values = py
        .allow_threads(move || rs_solve_darcy(&rs_mesh, permeability, &[], &bcs))
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok(PoissonSolution {
        mesh: mesh.unbind(),
        values,
    })
}

/// Steady Stokes (creeping) flow by the penalty method.
///
/// * `viscosity` — dynamic viscosity `mu`.
/// * `body_force` — constant body-force vector (length = mesh dimension).
/// * `bcs` — list of `(node_id, component, value)` velocity conditions.
/// * `penalty` — incompressibility penalty (default `1e6`).
///
/// Returns `(velocity, pressure)` as an [`ElasticitySolution`]-style vector
/// field and a [`PoissonSolution`]-style scalar field.
#[pyfunction]
#[pyo3(signature = (mesh, viscosity, body_force, bcs, penalty=1e6))]
fn solve_stokes(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    viscosity: f64,
    body_force: Vec<f64>,
    bcs: Vec<(usize, usize, f64)>,
    penalty: f64,
) -> PyResult<(ElasticitySolution, PoissonSolution)> {
    let dim = dim_of(&mesh.borrow().inner)?;
    if body_force.len() != dim {
        return Err(PyRuntimeError::new_err(format!(
            "body_force must have {dim} components, got {}",
            body_force.len()
        )));
    }
    let rs_mesh = mesh.borrow().inner.clone();
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    let (u, p) = py
        .allow_threads(move || {
            rs_steady_stokes(
                &rs_mesh,
                viscosity,
                move |_| body_force.clone(),
                &dir,
                penalty,
            )
        })
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let handle = mesh.unbind();
    let vel = ElasticitySolution {
        mesh: handle.clone_ref(py),
        values: u,
        dim,
    };
    let pres = PoissonSolution {
        mesh: handle,
        values: p,
    };
    Ok((vel, pres))
}

/// Thermal-structural coupling: free thermal expansion of an elastic body under
/// a nodal temperature rise.
///
/// * `model`, `young`, `poisson` — as for `solve_elasticity`.
/// * `alpha` — coefficient of thermal expansion.
/// * `delta_t` — temperature rise per node (length = node count).
/// * `bcs` — list of `(node_id, component, value)` displacement conditions.
///
/// Returns the displacement field as an [`ElasticitySolution`].
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn solve_thermal_structural(
    py: Python<'_>,
    mesh: Bound<'_, Mesh>,
    model: &str,
    young: f64,
    poisson: f64,
    alpha: f64,
    delta_t: Vec<f64>,
    bcs: Vec<(usize, usize, f64)>,
) -> PyResult<ElasticitySolution> {
    let model = parse_model(model)?;
    let dim = dim_of(&mesh.borrow().inner)?;
    let rs_mesh = mesh.borrow().inner.clone();
    if delta_t.len() != rs_mesh.node_count() {
        return Err(PyRuntimeError::new_err(format!(
            "delta_t has {} entries, mesh has {} nodes",
            delta_t.len(),
            rs_mesh.node_count()
        )));
    }
    let dir: Vec<(usize, f64)> = bcs.iter().map(|(n, c, v)| (n * dim + c, *v)).collect();
    let values = py
        .allow_threads(move || {
            rs_thermal_structural(&rs_mesh, model, young, poisson, alpha, &delta_t, &dir)
        })
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok(ElasticitySolution {
        mesh: mesh.unbind(),
        values,
        dim,
    })
}

/// Work-consistent fluid-structure interface load vector.
///
/// * `structure` — the structure mesh (`Tri`/`Quad` in 2-D, `Tet`/`Hex` in 3-D).
/// * `interface` — list of `(structure_node, fluid_node)` pairs.
/// * `fluid_pressure` — pressure at every fluid node.
///
/// Integrates the traction `p n` over the interface faces of the structure with
/// outward normals and returns the structure load vector (`node * dim +
/// component` ordering).
#[pyfunction]
fn fsi_interface_loads(
    structure: Bound<'_, Mesh>,
    interface: Vec<(usize, usize)>,
    fluid_pressure: Vec<f64>,
) -> PyResult<Vec<f64>> {
    let rs_mesh = structure.borrow().inner.clone();
    rs_fsi_interface_loads(&rs_mesh, &interface, &fluid_pressure)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Augmented-Lagrangian unilateral contact on a small dense system.
///
/// Solves `K u = f` subject to `u[dof] >= lower` for every `(dof, lower)`
/// constraint (non-penetration against a rigid obstacle), enforcing the
/// constraints exactly through multiplier updates. `stiffness` is an `n x n`
/// nested list, `load` has length `n`. Returns `(u, lambda)`: the displacement
/// vector and the contact forces (one multiplier per constraint).
#[pyfunction]
#[pyo3(signature = (stiffness, load, constraints, penalty=1e4, max_iter=50, tol=1e-9))]
fn contact_augmented_lagrangian(
    stiffness: Vec<Vec<f64>>,
    load: Vec<f64>,
    constraints: Vec<(usize, f64)>,
    penalty: f64,
    max_iter: usize,
    tol: f64,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let k = dense_to_coo("stiffness", &stiffness, load.len())?;
    let cons: Vec<ContactConstraint> = constraints
        .iter()
        .map(|&(dof, lower)| ContactConstraint { dof, lower })
        .collect();
    rs_augmented_lagrangian(&k, &load, &cons, penalty, max_iter, tol)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

/// Nearest-node contact pairing between two point sets.
///
/// `a` and `b` are lists of `(id, [x, y, z])`. For every point of `a` returns
/// `(a_id, (b_index, distance))` for its nearest point in `b` (`b_index` is the
/// position in the list `b`, not its id), or `(a_id, None)` when `b` is empty. Octree-accelerated (`O(|a| log |b|)`).
#[pyfunction]
#[allow(clippy::type_complexity)]
fn contact_pairs(
    a: Vec<(usize, Vec<f64>)>,
    b: Vec<(usize, Vec<f64>)>,
) -> Vec<(usize, Option<(usize, f64)>)> {
    rs_contact_pairs(&a, &b)
}

/// Classical lamination theory: the `6x6` ABD stiffness matrix of a laminate.
///
/// `plies` is the bottom-to-top stack of `(e1, e2, nu12, g12, thickness,
/// angle_deg)` tuples. Returns the row-major `[[A, B], [B, D]]` matrix
/// relating in-plane force / moment resultants to mid-surface strain /
/// curvature. A symmetric stack gives `B = 0`.
#[pyfunction]
fn laminate_abd(plies: Vec<(f64, f64, f64, f64, f64, f64)>) -> PyResult<Vec<Vec<f64>>> {
    if plies.is_empty() {
        return Err(PyRuntimeError::new_err("laminate needs at least one ply"));
    }
    let mut stack = Vec::with_capacity(plies.len());
    for (i, &(e1, e2, nu12, g12, thickness, angle_deg)) in plies.iter().enumerate() {
        let positive = [e1, e2, g12, thickness]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0);
        if !positive || !nu12.is_finite() || !angle_deg.is_finite() {
            return Err(PyRuntimeError::new_err(format!(
                "ply {i}: e1, e2, g12 and thickness must be finite and positive"
            )));
        }
        stack.push(Ply {
            e1,
            e2,
            nu12,
            g12,
            thickness,
            angle_deg,
        });
    }
    Ok(rs_laminate_abd(&stack).iter().map(|r| r.to_vec()).collect())
}

fn dense_to_coo(name: &str, m: &[Vec<f64>], n: usize) -> PyResult<Coo> {
    if m.len() != n || m.iter().any(|r| r.len() != n) {
        return Err(PyRuntimeError::new_err(format!(
            "{name} must be a {n}x{n} matrix"
        )));
    }
    let mut c = Coo::new();
    for (i, row) in m.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if v != 0.0 {
                c.push(i, j, v);
            }
        }
    }
    Ok(c)
}

/// Implicit Newmark-beta time integration of `M u'' + C u' + K u = f(t)` for a
/// small dense system.
///
/// `mass`, `damping`, `stiffness` are `n x n` nested lists; `u0`, `v0` the
/// initial displacement / velocity (length `n`); `load` is either a constant
/// force vector or a callable `f(t) -> list[float]`. Returns a list of
/// `(t, u)` pairs for steps `0..=nsteps`. Defaults `beta=0.25`, `gamma=0.5`
/// (average acceleration, unconditionally stable).
#[pyfunction]
#[pyo3(signature = (mass, damping, stiffness, u0, v0, load, dt, nsteps, beta=0.25, gamma=0.5))]
#[allow(clippy::too_many_arguments)]
fn newmark(
    mass: Vec<Vec<f64>>,
    damping: Vec<Vec<f64>>,
    stiffness: Vec<Vec<f64>>,
    u0: Vec<f64>,
    v0: Vec<f64>,
    load: &Bound<'_, PyAny>,
    dt: f64,
    nsteps: usize,
    beta: f64,
    gamma: f64,
) -> PyResult<Vec<(f64, Vec<f64>)>> {
    let n = u0.len();
    if v0.len() != n {
        return Err(PyRuntimeError::new_err(
            "u0 and v0 must have the same length",
        ));
    }
    let m = dense_to_coo("mass", &mass, n)?;
    let c = dense_to_coo("damping", &damping, n)?;
    let k = dense_to_coo("stiffness", &stiffness, n)?;
    let constant = load.extract::<Vec<f64>>().ok();
    let callback_error: std::cell::RefCell<Option<PyErr>> = std::cell::RefCell::new(None);
    let f = |t: f64| -> Vec<f64> {
        if let Some(v) = &constant {
            return v.clone();
        }
        match load.call1((t,)).and_then(|r| r.extract::<Vec<f64>>()) {
            Ok(v) => v,
            Err(e) => {
                *callback_error.borrow_mut() = Some(e);
                vec![0.0; n]
            }
        }
    };
    // Guard the callback result length so a short vector can't index-panic.
    let checked = |t: f64| -> Vec<f64> {
        let mut v = f(t);
        v.resize(n, 0.0);
        v
    };
    let opts = NewmarkOptions { dt, beta, gamma };
    let result = rs_newmark(&m, &c, &k, &u0, &v0, checked, &opts, nsteps)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()));
    if let Some(e) = callback_error.into_inner() {
        return Err(e);
    }
    result
}

/// `True` if this build of the extension was compiled with GPU support
/// (`maturin develop --features gpu`).
#[pyfunction]
fn gpu_enabled() -> bool {
    cfg!(feature = "gpu")
}

/// Name and backend of the GPU adapter used by `gpu_solve_cg`, e.g.
/// `"NVIDIA GeForce RTX 3050 (Vulkan)"`. Raises `RuntimeError` when the
/// extension was built without the `gpu` feature or no adapter exists.
#[pyfunction]
fn gpu_adapter() -> PyResult<String> {
    #[cfg(feature = "gpu")]
    {
        let ctx =
            tpt_fem_gpu::GpuContext::new().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(format!("{} ({})", ctx.adapter_name(), ctx.backend()))
    }
    #[cfg(not(feature = "gpu"))]
    {
        Err(PyRuntimeError::new_err(
            "tpt_fem was built without GPU support (rebuild with --features gpu)",
        ))
    }
}

/// Solve the symmetric positive-definite sparse system `A x = b` with
/// Jacobi-preconditioned conjugate gradients **on the GPU** (single-precision
/// kernels + double-precision iterative refinement).
///
/// * `triplets` — list of `(row, col, value)` entries of `A` (duplicates are
///   summed).
/// * `rhs` — right-hand side `b`.
/// * `tol` — target relative residual `||b - A x|| / ||b||` (default `1e-10`).
///
/// Returns `(x, gpu_iterations, relative_residual)`. Far faster than the CPU
/// solvers beyond roughly 50 000 unknowns. Requires the `gpu` build feature.
#[pyfunction]
#[pyo3(signature = (triplets, rhs, tol=1e-10))]
fn gpu_solve_cg(
    py: Python<'_>,
    triplets: Vec<(usize, usize, f64)>,
    rhs: Vec<f64>,
    tol: f64,
) -> PyResult<(Vec<f64>, usize, f64)> {
    #[cfg(feature = "gpu")]
    {
        let mut coo = Coo::new();
        for (r, c, v) in triplets {
            coo.push(r, c, v);
        }
        let opts = tpt_fem_gpu::GpuCgOptions {
            tol,
            ..tpt_fem_gpu::GpuCgOptions::default()
        };
        let sol = py
            .allow_threads(move || {
                let ctx = tpt_fem_gpu::GpuContext::new()?;
                tpt_fem_gpu::solve_cg_gpu(&ctx, &coo, &rhs, &opts)
            })
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok((sol.x, sol.gpu_iterations, sol.relative_residual))
    }
    #[cfg(not(feature = "gpu"))]
    {
        let _ = (py, triplets, rhs, tol);
        Err(PyRuntimeError::new_err(
            "tpt_fem was built without GPU support (rebuild with --features gpu)",
        ))
    }
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

/// SIMP minimum-compliance topology optimization of a 2-D cantilever.
///
/// The design lives on an `nx × ny` grid of unit square elements: the left
/// edge is clamped and a unit downward point load acts at the bottom-right
/// corner. Returns a [`TopOptSolution`] with the final element densities and
/// the compliance history.
#[pyfunction]
#[pyo3(signature = (nx, ny, vol_frac, penal=3.0, filter_radius=1.5, max_iter=50))]
fn topopt_cantilever(
    py: Python<'_>,
    nx: usize,
    ny: usize,
    vol_frac: f64,
    penal: f64,
    filter_radius: f64,
    max_iter: usize,
) -> PyResult<TopOptSolution> {
    if nx < 1 || ny < 1 {
        return Err(PyRuntimeError::new_err(format!(
            "nx and ny must be >= 1, got {nx} x {ny}"
        )));
    }
    if !(vol_frac > 0.0 && vol_frac <= 1.0) {
        return Err(PyRuntimeError::new_err(format!(
            "vol_frac must be in (0, 1], got {vol_frac}"
        )));
    }
    if penal < 1.0 {
        return Err(PyRuntimeError::new_err(format!(
            "penal must be >= 1, got {penal}"
        )));
    }
    if filter_radius < 0.0 {
        return Err(PyRuntimeError::new_err(format!(
            "filter_radius must be >= 0, got {filter_radius}"
        )));
    }
    let res = py.allow_threads(move || {
        let grid = Grid::new(nx + 1, ny + 1, 1.0);
        let (f, bcs) = cantilever_load(&grid, 1.0);
        let params = TopOptParams {
            grid,
            e0: 1.0,
            nu: 0.3,
            vol_frac,
            penal,
            filter_radius,
            max_iter,
            move_limit: 0.2,
        };
        topopt_simp(&params, &f, &bcs).map_err(|e| PyRuntimeError::new_err(e.to_string()))
    })?;
    Ok(TopOptSolution {
        nx,
        ny,
        densities: res.densities,
        compliance: res.compliance,
        iterations: res.iterations,
    })
}

/// Result of [`topopt_cantilever`]: optimized element densities + history.
#[pyclass]
struct TopOptSolution {
    nx: usize,
    ny: usize,
    densities: Vec<f64>,
    compliance: Vec<f64>,
    iterations: usize,
}

#[pymethods]
impl TopOptSolution {
    /// Number of elements along `x`.
    #[getter]
    fn nx(&self) -> usize {
        self.nx
    }

    /// Number of elements along `y`.
    #[getter]
    fn ny(&self) -> usize {
        self.ny
    }

    /// Final element densities, row-major (`y` rows of `x` columns).
    #[getter]
    fn densities(&self) -> Vec<f64> {
        self.densities.clone()
    }

    /// Compliance at each iteration (index 0 is the uniform start).
    #[getter]
    fn compliance(&self) -> Vec<f64> {
        self.compliance.clone()
    }

    /// Number of optimality-criteria iterations performed.
    #[getter]
    fn iterations(&self) -> usize {
        self.iterations
    }

    /// Densities as a `(ny, nx)` `numpy.ndarray`.
    fn to_numpy(&self, py: Python<'_>) -> PyResult<PyObject> {
        to_numpy_array(py, &self.densities, &[self.ny, self.nx])
    }

    fn __repr__(&self) -> String {
        format!(
            "TopOptSolution({}x{} elements, iterations={}, compliance {:.4e} -> {:.4e})",
            self.nx,
            self.ny,
            self.iterations,
            self.compliance.first().copied().unwrap_or(0.0),
            self.compliance.last().copied().unwrap_or(0.0)
        )
    }
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
    m.add_class::<TopOptSolution>()?;
    m.add_class::<HeatHistory>()?;
    m.add_function(wrap_pyfunction!(solve_transient_heat, py)?)?;
    m.add_function(wrap_pyfunction!(j2_uniaxial_response, py)?)?;
    m.add_function(wrap_pyfunction!(solve_darcy, py)?)?;
    m.add_function(wrap_pyfunction!(laminate_abd, py)?)?;
    m.add_function(wrap_pyfunction!(solve_thermal_structural, py)?)?;
    m.add_function(wrap_pyfunction!(contact_pairs, py)?)?;
    m.add_function(wrap_pyfunction!(contact_augmented_lagrangian, py)?)?;
    m.add_function(wrap_pyfunction!(fsi_interface_loads, py)?)?;
    m.add_function(wrap_pyfunction!(gpu_enabled, py)?)?;
    m.add_function(wrap_pyfunction!(gpu_adapter, py)?)?;
    m.add_function(wrap_pyfunction!(gpu_solve_cg, py)?)?;
    m.add_function(wrap_pyfunction!(newmark, py)?)?;
    m.add_function(wrap_pyfunction!(solve_stokes, py)?)?;
    m.add_function(wrap_pyfunction!(neo_hookean_uniaxial, py)?)?;
    m.add_function(wrap_pyfunction!(solve_poisson, py)?)?;
    m.add_function(wrap_pyfunction!(solve_elasticity, py)?)?;
    m.add_function(wrap_pyfunction!(solve_modal, py)?)?;
    m.add_function(wrap_pyfunction!(topopt_cantilever, py)?)?;
    Ok(())
}
