//! WebAssembly bindings for the `tpt-fem` SIMP topology-optimization demo.
//!
//! Exposes a single entry point, [`optimize_cantilever`], that runs the same
//! solver as `tpt-fem topopt` / `tpt_fem.topopt_cantilever` in the browser and
//! returns the result for the `www/` page to draw.

use tpt_fem_topopt::{cantilever_load, topopt_simp, Grid, TopOptParams};
use wasm_bindgen::prelude::*;

/// Largest grid accepted; the in-house dense LU solve makes bigger grids too
/// slow for the main thread.
const MAX_ELEMENTS_X: usize = 40;
const MAX_ELEMENTS_Y: usize = 16;

/// Result of [`optimize_cantilever`].
#[wasm_bindgen]
pub struct TopOptResult {
    nx: usize,
    ny: usize,
    densities: Vec<f64>,
    compliance: Vec<f64>,
}

#[wasm_bindgen]
impl TopOptResult {
    /// Elements along `x`.
    #[wasm_bindgen(getter)]
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Elements along `y`.
    #[wasm_bindgen(getter)]
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Final element densities, row-major from the bottom row (`y = 0`) up.
    #[wasm_bindgen(getter)]
    pub fn densities(&self) -> Vec<f64> {
        self.densities.clone()
    }

    /// Compliance per iteration (index 0 is the uniform start).
    #[wasm_bindgen(getter)]
    pub fn compliance(&self) -> Vec<f64> {
        self.compliance.clone()
    }
}

/// Optimize a 2-D cantilever (clamped left edge, unit downward load at the
/// bottom-right corner) on an `nx × ny` element grid.
#[wasm_bindgen]
pub fn optimize_cantilever(
    nx: usize,
    ny: usize,
    vol_frac: f64,
    penal: f64,
    filter_radius: f64,
    max_iter: usize,
) -> Result<TopOptResult, JsError> {
    if !(1..=MAX_ELEMENTS_X).contains(&nx) || !(1..=MAX_ELEMENTS_Y).contains(&ny) {
        return Err(JsError::new(&format!(
            "grid must be 1..={MAX_ELEMENTS_X} x 1..={MAX_ELEMENTS_Y} elements, got {nx} x {ny}"
        )));
    }
    if !(vol_frac > 0.0 && vol_frac <= 1.0) {
        return Err(JsError::new("vol_frac must be in (0, 1]"));
    }
    if penal < 1.0 || filter_radius < 0.0 {
        return Err(JsError::new("penal must be >= 1 and filter_radius >= 0"));
    }
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
    let res = topopt_simp(&params, &f, &bcs).map_err(|e| JsError::new(&e.to_string()))?;
    Ok(TopOptResult {
        nx,
        ny,
        densities: res.densities,
        compliance: res.compliance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optimizes_small_cantilever() {
        let r = optimize_cantilever(12, 6, 0.5, 3.0, 1.5, 5).ok().unwrap();
        assert_eq!(r.densities.len(), 72);
        assert!(r.compliance.last().unwrap() < &r.compliance[0]);
    }
}
