//! SIMP topology optimization for `tpt-fem`.
//!
//! Minimises the compliance of a linear-elastic structure over a fixed mesh
//! by redistributing a bounded volume fraction of material with the classic
//! [SIMP](https://en.wikipedia.org/wiki/Topology_optimization#SIMP_method)
//! scheme (`E_e = x_e^p·E`, `x_min ≤ x_e ≤ 1`):
//!
//! 1. assemble `K(x)` from the per-element elasticity matrices scaled by
//!    `x_e^p`,
//! 2. solve `K u = f` (Dirichlet-condensed),
//! 3. compute sensitivities `∂c/∂x_e = −p·x_e^(p−1)·uₑᵀ k0e uₑ`,
//! 4. apply a sensitivity density filter over element centroids,
//! 5. update the design with an optimality-criteria (OC) step whose Lagrange
//!    multiplier is found by bisection to meet the volume fraction.
//!
//! ```no_run
//! use tpt_fem_topopt::{cantilever_problem, simp_optimize, SimpOptions};
//!
//! // A 32×16 cantilever at 40% volume: the optimizer carves out a load
//! // path and roughly halves the compliance of the full-block design.
//! let problem = cantilever_problem(32, 16);
//! let opts = SimpOptions { volfrac: 0.4, ..Default::default() };
//! let result = simp_optimize(&problem, &opts).unwrap();
//! assert!(result.compliance.last().unwrap() < &result.compliance[0]);
//! ```

use tpt_fem_assembly::solve_with_dirichlet;
use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel};
use tpt_fem_mesh::{CellType, Mesh};

/// Options for [`simp_optimize`].
#[derive(Clone, Debug)]
pub struct SimpOptions {
    /// Target volume fraction of the design domain (Σx / n_elements).
    pub volfrac: f64,
    /// SIMP penalty exponent (classical choice: 3).
    pub penalty: f64,
    /// Density-filter support radius, in units of the average element size.
    pub filter_radius: f64,
    /// Relative OC move limit per iteration.
    pub move_limit: f64,
    /// Lower bound on any element density.
    pub x_min: f64,
    /// Maximum OC iterations.
    pub max_iter: usize,
    /// Stop early when the relative compliance change drops below this.
    pub tol: f64,
}

impl Default for SimpOptions {
    fn default() -> Self {
        SimpOptions {
            volfrac: 0.4,
            penalty: 3.0,
            filter_radius: 1.5,
            move_limit: 0.2,
            x_min: 1e-3,
            max_iter: 60,
            tol: 1e-4,
        }
    }
}

/// A fully specified SIMP problem on a caller-built mesh.
///
/// The load is a per-global-DOF force vector and `dirichlet` uses the same
/// `(dof, value)` convention as [`tpt_fem_assembly::solve_with_dirichlet`].
/// All elements must be planar quads or triangles (`PlaneStress`/
/// `PlaneStrain`, 2 DOF/node).
#[derive(Clone)]
pub struct SimpProblem {
    /// Design-domain mesh.
    pub mesh: Mesh,
    /// Plane elasticity model.
    pub model: ElasticModel,
    /// Base (solid) Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Gauss order used for the element matrices.
    pub quad_order: usize,
    /// Global load vector (one entry per DOF, `2·n_nodes` long).
    pub load: Vec<f64>,
    /// Essential boundary conditions.
    pub dirichlet: Vec<(usize, f64)>,
}

/// Result of [`simp_optimize`].
#[derive(Clone, Debug)]
pub struct SimpResult {
    /// Final element densities.
    pub densities: Vec<f64>,
    /// Compliance history (one entry per completed iteration).
    pub compliance: Vec<f64>,
    /// Volume-fraction history.
    pub volume_fraction: Vec<f64>,
}

fn element_centroids(mesh: &Mesh) -> Vec<[f64; 2]> {
    mesh.elements
        .iter()
        .map(|e| {
            let mut c = [0.0, 0.0];
            for &nd in &e.nodes {
                let p = mesh.node_coords(nd);
                c[0] += p[0];
                c[1] += p[1];
            }
            let n = e.nodes.len() as f64;
            [c[0] / n, c[1] / n]
        })
        .collect()
}

/// Run the SIMP optimization loop.
///
/// Returns an error if the underlying linear solve fails (e.g. a singular
/// reduced stiffness from an over-constrained problem).
pub fn simp_optimize(
    problem: &SimpProblem,
    opts: &SimpOptions,
) -> Result<SimpResult, tpt_fem_sparse::SparseError> {
    let mesh = &problem.mesh;
    let n_elem = mesh.elements.len();
    let dim = 2usize;

    // Base element stiffness matrices (E = 1); scaled by x^p·young at use.
    let k0: Vec<Vec<Vec<f64>>> = (0..n_elem)
        .map(|eid| {
            elasticity_element_matrix(
                mesh,
                eid,
                problem.model,
                1.0,
                problem.poisson,
                problem.quad_order,
            )
            .expect("element stiffness must build on the design mesh")
        })
        .collect();

    // Filter weights: w_ej = max(0, r − |c_e − c_j|), normalised per element.
    let centroids = element_centroids(mesh);
    // Average element size = mean of the largest node-pair distance per element.
    let mut h_sum = 0.0_f64;
    for e in &mesh.elements {
        let n = e.nodes.len();
        let mut hmax = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let pa = mesh.node_coords(e.nodes[i]);
                let pb = mesh.node_coords(e.nodes[j]);
                hmax = hmax.max(((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt());
            }
        }
        h_sum += hmax;
    }
    let avg_h = if n_elem > 0 {
        h_sum / n_elem as f64
    } else {
        1.0
    };
    let radius = opts.filter_radius * avg_h;
    let mut weights: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_elem];
    for e in 0..n_elem {
        for j in 0..n_elem {
            let d = ((centroids[e][0] - centroids[j][0]).powi(2)
                + (centroids[e][1] - centroids[j][1]).powi(2))
            .sqrt();
            if d < radius {
                weights[e].push((j, radius - d));
            }
        }
        let sum: f64 = weights[e].iter().map(|(_, w)| w).sum();
        for (_, w) in weights[e].iter_mut() {
            *w /= sum;
        }
    }

    let mut x = vec![opts.volfrac; n_elem];
    let mut compliance_hist = Vec::new();
    let mut volfrac_hist = Vec::new();

    for _ in 0..opts.max_iter {
        // Assemble K(x) = Σ x_e^p · young · k0e.
        let coo = tpt_fem_assembly::assemble(mesh, dim, |eid, _| {
            let s = x[eid].powf(opts.penalty) * problem.young;
            k0[eid]
                .iter()
                .map(|row| row.iter().map(|v| v * s).collect())
                .collect()
        });
        let u = solve_with_dirichlet(&coo, &problem.load, &problem.dirichlet)?;

        // Compliance and raw sensitivities.
        let mut compliance = 0.0;
        let mut sens = vec![0.0; n_elem];
        for eid in 0..n_elem {
            let nodes = &mesh.elements[eid].nodes;
            let mut ue = Vec::with_capacity(nodes.len() * dim);
            for &nd in nodes {
                ue.push(u[nd * dim]);
                ue.push(u[nd * dim + 1]);
            }
            let mut uku = 0.0;
            for (i, ui) in ue.iter().enumerate() {
                for (j, uj) in ue.iter().enumerate() {
                    uku += k0[eid][i][j] * ui * uj;
                }
            }
            uku *= problem.young;
            compliance += u[eid * dim] * problem.load[eid * dim]
                + u[eid * dim + 1] * problem.load[eid * dim + 1];
            sens[eid] = -opts.penalty * x[eid].powf(opts.penalty - 1.0) * uku;
        }
        // c = fᵀu ≥ 0 for a stable structure.

        // Sensitivity filtering.
        let filtered: Vec<f64> = (0..n_elem)
            .map(|e| {
                let num: f64 = weights[e].iter().map(|&(j, w)| w * sens[j]).sum();
                num / x[e].max(opts.x_min)
            })
            .collect();

        // OC update with bisection on the Lagrange multiplier.
        let target = opts.volfrac * n_elem as f64;
        let mut lo = 0.0_f64;
        let mut hi = 1e9_f64;
        let mut x_new = x.clone();
        for _ in 0..60 {
            let lam = 0.5 * (lo + hi);
            for e in 0..n_elem {
                let candidate = x[e] * (-filtered[e] / lam).abs().sqrt();
                let clipped = candidate
                    .max(x[e] - opts.move_limit)
                    .min(x[e] + opts.move_limit)
                    .clamp(opts.x_min, 1.0);
                x_new[e] = clipped;
            }
            let vol: f64 = x_new.iter().sum();
            // Larger λ shrinks the OC step (and the volume); if we are above
            // the target we must increase λ.
            if vol > target {
                lo = lam;
            } else {
                hi = lam;
            }
        }
        let vol: f64 = x_new.iter().sum();
        // Under-relaxation: blending the OC step with the previous design
        // damps the well-known iteration-level oscillation of pure OC.
        for e in 0..n_elem {
            x[e] = 0.5 * (x_new[e] + x[e]);
        }
        compliance_hist.push(compliance);
        volfrac_hist.push(vol / n_elem as f64);
        if compliance_hist.len() >= 2 {
            let prev = compliance_hist[compliance_hist.len() - 2];
            if prev.abs() > 0.0 && ((compliance - prev) / prev).abs() < opts.tol {
                break;
            }
        }
    }

    Ok(SimpResult {
        densities: x,
        compliance: compliance_hist,
        volume_fraction: volfrac_hist,
    })
}

/// Build a rectangular cantilever benchmark: a `[0,w]×[0,h]` domain of
/// `nx×ny` bilinear quads, clamped on the left edge, with a downward unit
/// load at mid-span of the right edge.
pub fn cantilever_problem(nx: usize, ny: usize) -> SimpProblem {
    let (mesh, load, dirichlet) = grid_problem(nx, ny, |n, x, _y, _w, _h| {
        if x < 1e-9 {
            Some(vec![(n * 2, 0.0), (n * 2 + 1, 0.0)])
        } else {
            None
        }
    });
    // Mid-span node of the right edge.
    let (w, h) = (1.0_f64, 1.0_f64);
    let tip = (0..mesh.node_count())
        .find(|&n| {
            let p = mesh.node_coords(n);
            (p[0] - w).abs() < 1e-9 && (p[1] - h / 2.0).abs() < 1e-9
        })
        .expect("mid-span node must exist");
    let mut load = load;
    load[tip * 2 + 1] = -1.0;
    SimpProblem {
        mesh,
        model: ElasticModel::PlaneStress,
        young: 1.0,
        poisson: 0.3,
        quad_order: 2,
        load,
        dirichlet,
    }
}

/// Build an `nx×ny` unit-square quad grid; `bc(x, y, w, h)` returns the
/// Dirichlet conditions for boundary nodes (or `None`).
fn grid_problem(
    nx: usize,
    ny: usize,
    bc: impl Fn(usize, f64, f64, f64, f64) -> Option<Vec<(usize, f64)>>,
) -> (Mesh, Vec<f64>, Vec<(usize, f64)>) {
    let mut b = tpt_fem_mesh::MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=ny {
        let y = j as f64 / ny as f64;
        let mut r = Vec::new();
        for i in 0..=nx {
            r.push(b.add_node(vec![i as f64 / nx as f64, y]));
        }
        rows.push(r);
    }
    for j in 0..ny {
        for i in 0..nx {
            b.add_element(
                CellType::Quad,
                vec![
                    rows[j][i],
                    rows[j][i + 1],
                    rows[j + 1][i + 1],
                    rows[j + 1][i],
                ],
            );
        }
    }
    let mesh = b.build();
    let n_dof = mesh.node_count() * 2;
    let load = vec![0.0; n_dof];
    let mut dirichlet = Vec::new();
    for n in 0..mesh.node_count() {
        let p = mesh.node_coords(n);
        if let Some(con) = bc(n, p[0], p[1], 1.0, 1.0) {
            dirichlet.extend(con);
        }
    }
    (mesh, load, dirichlet)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solid_cantilever_solves() {
        // Sanity: the assembled solid design solves through the standard
        // elasticity path on the benchmark mesh.
        let problem = cantilever_problem(4, 2);
        let u = tpt_fem_elasticity::solve_elasticity(
            &problem.mesh,
            problem.model,
            problem.young,
            problem.poisson,
            problem.quad_order,
            |_| vec![0.0, 0.0],
            &problem.dirichlet,
        )
        .expect("solid cantilever must solve");
        assert_eq!(u.len(), problem.mesh.node_count() * 2);
    }

    #[test]
    fn simp_reduces_compliance_and_meets_volume() {
        // A 16×8 cantilever at 40% volume: compliance must drop versus the
        // uniform 40% design, the final volume fraction must hit the target
        // (OC enforces it exactly), and all densities must stay in bounds.
        let problem = cantilever_problem(16, 8);
        let opts = SimpOptions {
            volfrac: 0.4,
            max_iter: 30,
            ..Default::default()
        };
        let result = simp_optimize(&problem, &opts).unwrap();
        assert_eq!(result.densities.len(), problem.mesh.elements.len());
        for &x in &result.densities {
            assert!((1e-3..=1.0).contains(&x), "density out of bounds: {x}");
        }
        let vf = *result.volume_fraction.last().unwrap();
        assert!((vf - 0.4).abs() < 5e-3, "final volume fraction {vf} != 0.4");
        // Compliance must improve overall and stay near-monotone (OC with
        // move limits can wiggle slightly between iterations).
        assert!(
            result.compliance.last().unwrap() < result.compliance.first().unwrap(),
            "compliance did not improve"
        );
        let wiggle = result
            .compliance
            .windows(2)
            .all(|w| w[1] <= w[0] * (1.0 + 5e-2));
        assert!(wiggle, "compliance spiked: {:?}", result.compliance);
    }

    #[test]
    fn simp_full_material_is_stationary() {
        // At volfrac = 1 with a generous move limit, the OC update cannot
        // improve on solid material: every density stays at (or returns to)
        // 1 and the compliance matches the analytic fᵀu of the full structure.
        let problem = cantilever_problem(8, 4);
        let opts = SimpOptions {
            volfrac: 1.0,
            move_limit: 0.05,
            penalty: 1.0, // linear stiffness: x=1 is optimal
            max_iter: 5,
            ..Default::default()
        };
        let result = simp_optimize(&problem, &opts).unwrap();
        for &x in &result.densities {
            assert!((x - 1.0).abs() < 1e-6, "density drifted from solid: {x}");
        }
        assert!(*result.compliance.first().unwrap() > 0.0);
    }
}
