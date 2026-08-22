//! Porous-media models for `tpt-fem`: steady Darcy flow and transient Biot
//! consolidation.
//!
//! Both are built from a scalar (pressure / excess-pore-pressure) field using
//! the same isoparametric gradient machinery as the displacement elements:
//!
//! * [`solve_darcy`] assembles `K = ∫ k ∇Nᵀ∇N dΩ` (steady-state
//!   `−∇·(k∇p) = f`) and solves it with `tpt-fem-sparse`.
//! * [`terzaghi_consolidation`] integrates the 1-D Terzaghi consolidation
//!   equation `∂u/∂t = cᵥ ∂²u/∂z²` (a diffusion problem) with a stable
//!   backward-Euler step and recovers the settlement history, checked against
//!   the closed-form final (drained) settlement.
//!
//! ```
//! use tpt_fem_porous::solve_darcy;
//! use tpt_fem_mesh::{CellType, MeshBuilder};
//!
//! // 1-D column of unit length, permeability k=1, drained at both ends (p=0),
//! // zero source: the solution is the trivial p≡0.
//! let mut b = MeshBuilder::new();
//! let a = b.add_node(vec![0.0]);
//! let m = b.add_node(vec![0.5]);
//! let c = b.add_node(vec![1.0]);
//! b.add_element(CellType::Line, vec![a, m]);
//! b.add_element(CellType::Line, vec![m, c]);
//! let mesh = b.build();
//! let src: Vec<Vec<f64>> = Vec::new();
//! let p = solve_darcy(&mesh, 1.0, &src, &[(a, 0.0), (c, 0.0)]).unwrap();
//! assert!(p.iter().all(|&x| x.abs() < 1e-9));
//! ```

use tpt_fem_assembly::{assemble, solve_with_dirichlet};
use tpt_fem_element::{Hex8, Line2, Map, Quad4, ReferenceElement, Tet4, Tri3};
use tpt_fem_mesh::{CellType, Mesh};
use tpt_fem_quadrature::{
    gauss_legendre, tensor_cube, tensor_square, tetrahedron, triangle, TetrahedronRule,
    TriangleRule,
};

/// Cell-type dispatch for a scalar element's reference gradients and quadrature.
fn cell_scalar(cell: CellType, order: usize) -> (Vec<Vec<f64>>, Vec<f64>, usize) {
    match cell {
        CellType::Line => {
            let r = gauss_legendre(order);
            let qpts: Vec<Vec<f64>> = r.points.iter().map(|x| vec![*x]).collect();
            (qpts, r.weights, Line2::DIM)
        }
        CellType::Tri => {
            let r = triangle(TriangleRule::Degree2);
            let qpts = r.points.iter().map(|p| vec![p[0], p[1]]).collect();
            (qpts, r.weights, Tri3::DIM)
        }
        CellType::Quad => {
            let r = tensor_square(&gauss_legendre(order));
            let qpts = r.points.iter().map(|p| vec![p[0], p[1]]).collect();
            (qpts, r.weights, Quad4::DIM)
        }
        CellType::Tet => {
            let r = tetrahedron(TetrahedronRule::Degree2);
            let qpts = r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect();
            (qpts, r.weights, Tet4::DIM)
        }
        CellType::Hex => {
            let r = tensor_cube(&gauss_legendre(order));
            let qpts = r.points.iter().map(|p| vec![p[0], p[1], p[2]]).collect();
            (qpts, r.weights, Hex8::DIM)
        }
        other => panic!("porous: unsupported cell {other:?} for scalar field"),
    }
}

fn ref_grad_scalar(cell: CellType, xi: &[f64]) -> Vec<Vec<f64>> {
    match cell {
        CellType::Line => Line2::grad(xi),
        CellType::Tri => Tri3::grad(xi),
        CellType::Quad => Quad4::grad(xi),
        CellType::Tet => Tet4::grad(xi),
        CellType::Hex => Hex8::grad(xi),
        other => panic!("porous: unsupported cell {other:?}"),
    }
}

fn ref_shape_scalar(cell: CellType, xi: &[f64]) -> Vec<f64> {
    match cell {
        CellType::Line => Line2::shape(xi),
        CellType::Tri => Tri3::shape(xi),
        CellType::Quad => Quad4::shape(xi),
        CellType::Tet => Tet4::shape(xi),
        CellType::Hex => Hex8::shape(xi),
        other => panic!("porous: unsupported cell {other:?}"),
    }
}

/// Element conductivity `K_e = k ∫ ∇Nᵀ∇N dΩ` and mass `M_e = ∫ NᵀN dΩ` for a
/// scalar field, returned in node order (1 DOF per node).
fn scalar_elem(mesh: &Mesh, eid: usize, k: f64, order: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let elem = &mesh.elements[eid];
    let cell = elem.cell_type;
    let phys: Vec<Vec<f64>> = elem
        .nodes
        .iter()
        .map(|&n| mesh.node_coords(n).to_vec())
        .collect();
    let n = elem.nodes.len();
    let (qpts, qw, _dim) = cell_scalar(cell, order);
    let mut ke = vec![vec![0.0; n]; n];
    let mut me = vec![vec![0.0; n]; n];
    for (qp, w) in qpts.iter().zip(&qw) {
        let local = ref_grad_scalar(cell, qp);
        let map = Map::from_nodes_and_grad(&phys, &local);
        let det = map.determinant.abs();
        let g: Vec<Vec<f64>> = local.iter().map(|g| map.physical_grad(g)).collect();
        for i in 0..n {
            for j in 0..n {
                let mut d = 0.0;
                for d_ in 0..g[0].len() {
                    d += g[i][d_] * g[j][d_];
                }
                ke[i][j] += w * det * k * d;
            }
        }
        let ns = ref_shape_scalar(cell, qp);
        for i in 0..n {
            for j in 0..n {
                me[i][j] += w * det * ns[i] * ns[j];
            }
        }
    }
    (ke, me)
}

/// Steady Darcy flow: solve `−∇·(k∇p) = f` for the pressure field `p`, with
/// `dirichlet` fixing `(global_dof, value)` pressure nodes. `source(eid) -> f_e`
/// (node-order vector) supplies the element body source; pass `&[]` for none.
pub fn solve_darcy(
    mesh: &Mesh,
    permeability: f64,
    source: &[Vec<f64>],
    dirichlet: &[(usize, f64)],
) -> Result<Vec<f64>, tpt_fem_sparse::SparseError> {
    let order = 2;
    let coo = assemble(mesh, 1, |eid, m| scalar_elem(m, eid, permeability, order).0);
    let mut rhs = vec![0.0; mesh.node_count()];
    for (eid, elem) in mesh.elements.iter().enumerate() {
        if eid < source.len() {
            for (i, &node) in elem.nodes.iter().enumerate() {
                rhs[node] += source[eid][i];
            }
        }
    }
    solve_with_dirichlet(&coo, &rhs, dirichlet)
}

/// Terzaghi 1-D consolidation of a saturated column of height `H` (a `Line2`
/// mesh spanning `[0, H]`) under an initial uniform excess pore pressure `q0`.
///
/// The top (`z = H`) is drained (`u = 0`); the base (`z = 0`) is impermeable.
/// Returns the settlement history `(t, settlement)` together with the maximum
/// excess pore pressure at each step. The final (drained) settlement equals the
/// closed-form `q0·H / E_v`.
///
/// * `cv` — consolidation coefficient `cᵥ = k·E_v / γ_w` (length²/time).
/// * `ev` — 1-D drained modulus `E_v`.
/// * `dt` — time step. The backward-Euler step is unconditionally stable, so
///   there is no Courant-type restriction on `dt`; larger steps simply reduce
///   the temporal accuracy.
pub fn terzaghi_consolidation(
    mesh: &Mesh,
    q0: f64,
    cv: f64,
    ev: f64,
    total_time: f64,
    dt: f64,
) -> Vec<(f64, f64, f64)> {
    // Uniform element length from the first element.
    let e0 = &mesh.elements[0];
    let _dz = (mesh.node_coords(e0.nodes[1])[0] - mesh.node_coords(e0.nodes[0])[0]).abs();

    // Diffusion matrices: K = ∫ k ∇Nᵀ∇N, M = ∫ NᵀN. The Terzaghi equation is
    // ∂u/∂t = cᵥ ∂²u/∂z², whose FE form is M u̇ + cᵥ·K' u = 0 (K' = ∫∇Nᵀ∇N), so
    // the scalar conductivity fed to the element matrix is k = cᵥ.
    let k = cv;
    let coo_k = assemble(mesh, 1, |eid, m| scalar_elem(m, eid, k, 2).0);
    let coo_m = assemble(mesh, 1, |eid, m| scalar_elem(m, eid, 1.0, 2).1);

    // Backward-Euler: (M + dt K) u_{n+1} = M u_n.
    let mut lhs = coo_m.clone();
    // Add dt*K to lhs by pushing scaled entries.
    for i in 0..coo_k.rows.len() {
        lhs.push(coo_k.rows[i], coo_k.cols[i], dt * coo_k.vals[i]);
    }
    let m_csr = coo_m.to_csr();

    let h = mesh.node_count();
    let mut u = vec![q0; h];
    u[0] = q0; // impermeable base: no-flux handled by K structure (no top BC at z=0)
               // Drained top is the last node (largest coordinate).
    let top = (0..h)
        .max_by(|&a, &b| {
            mesh.node_coords(a)[0]
                .partial_cmp(&mesh.node_coords(b)[0])
                .unwrap()
        })
        .unwrap();
    u[top] = 0.0;

    let dirichlet = [(top, 0.0)];
    let mut out = Vec::new();
    let mut t = 0.0;
    let nsteps = (total_time / dt).round() as usize;
    for _ in 0..nsteps {
        // r = M u_n
        let mut r = vec![0.0; h];
        for row in 0..h {
            let mut s = 0.0;
            for c in m_csr.row_ptrs[row]..m_csr.row_ptrs[row + 1] {
                s += m_csr.values[c] * u[m_csr.col_ind[c]];
            }
            r[row] = s;
        }
        let u_next = solve_with_dirichlet(&lhs, &r, &dirichlet).expect("terzaghi solve");
        u = u_next;
        u[top] = 0.0;
        t += dt;
        // Settlement = (1/E_v) ∫ (q0 - u) dz, approximated by trapezoidal rule.
        let mut integral = 0.0;
        for row in 0..h {
            let mut s = 0.0;
            for c in m_csr.row_ptrs[row]..m_csr.row_ptrs[row + 1] {
                s += m_csr.values[c] * (q0 - u[m_csr.col_ind[c]]);
            }
            integral += s;
        }
        let settlement = integral / ev;
        let max_u = u.iter().cloned().fold(0.0_f64, f64::max);
        out.push((t, settlement, max_u));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    fn line_mesh(n: usize, length: f64) -> Mesh {
        let mut b = MeshBuilder::new();
        let mut prev = b.add_node(vec![0.0]);
        for i in 1..=n {
            let node = b.add_node(vec![length * i as f64 / n as f64]);
            b.add_element(CellType::Line, vec![prev, node]);
            prev = node;
        }
        b.build()
    }

    #[test]
    fn darcy_zero_source_zero_pressure() {
        let mesh = line_mesh(4, 1.0);
        let p = solve_darcy(&mesh, 1.0, &[], &[(0, 0.0), (4, 0.0)]).unwrap();
        assert!(p.iter().all(|&x| x.abs() < 1e-9));
    }

    #[test]
    fn darcy_linear_pressure_drop() {
        // Constant permeability, zero source, p=1 at x=0 and p=0 at x=1:
        // steady solution is linear, p(x)=1-x.
        let mesh = line_mesh(4, 1.0);
        let p = solve_darcy(&mesh, 2.0, &[], &[(0, 1.0), (4, 0.0)]).unwrap();
        // sample a couple of nodes by coordinate
        let at = |x: f64| -> f64 {
            let idx = (x * 4.0).round() as usize;
            p[idx]
        };
        assert!((at(0.0) - 1.0).abs() < 1e-9);
        assert!((at(1.0) - 0.0).abs() < 1e-9);
        assert!((at(0.5) - 0.5).abs() < 1e-2);
    }

    #[test]
    fn terzaghi_final_settlement_matches_closed_form() {
        let mesh = line_mesh(10, 1.0); // H = 1.0
        let q0 = 100.0;
        let ev = 1e6;
        let cv = 0.1;
        let dt = 0.01; // stable: Δz²/(2 cᵥ) = 0.01/0.2 = 0.05
        let total = 30.0; // Tᵥ = cᵥ t / H² = 3 -> >95% consolidated
        let hist = terzaghi_consolidation(&mesh, q0, cv, ev, total, dt);
        let final_settlement = hist.last().unwrap().1;
        let want = q0 * 1.0 / ev; // closed-form drained settlement
        assert!(
            (final_settlement - want).abs() / want < 0.05,
            "settlement {} vs {}",
            final_settlement,
            want
        );
        // Excess pore pressure should have decayed substantially.
        assert!(hist.last().unwrap().2 < q0 * 0.2);
        // Settlement must be monotonically non-decreasing.
        for w in 1..hist.len() {
            assert!(hist[w].1 >= hist[w - 1].1 - 1e-12);
        }
    }
}
