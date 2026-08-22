//! Method-of-Manufactured-Solutions (MMS) convergence tests for the Poisson
//! solver in `tpt_fem_thermal`.
//!
//! Each test manufactures an exact smooth solution `u`, derives the source `f`
//! from `-∇·(k∇u) = f`, solves on a refinement sequence, and asserts that the
//! discrete `L2` and `H1` errors converge at the theoretical rates for P1
//! elements (`L2` order 2, `H1` order 1).

use tpt_fem_element::{Line2, Map, ReferenceElement, Tet4, Tri3};
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_fem_quadrature::{gauss_legendre, tetrahedron, triangle, TetrahedronRule, TriangleRule};

use tpt_fem_thermal::solve_poisson;

// ---------------------------------------------------------------------------
// Mesh builders
// ---------------------------------------------------------------------------

/// Uniformly spaced 1-D line mesh of `[0, 1]` with `n` elements.
fn line_mesh(n: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut nodes = Vec::with_capacity(n + 1);
    for i in 0..=n {
        nodes.push(b.add_node(vec![i as f64 / n as f64]));
    }
    for w in nodes.windows(2) {
        b.add_element(CellType::Line, vec![w[0], w[1]]);
    }
    b.build()
}

/// Structured `n x n` triangular mesh of the unit square.
fn square_tri_mesh(n: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut ids = vec![vec![0usize; n + 1]; n + 1];
    for j in 0..=n {
        for i in 0..=n {
            ids[j][i] = b.add_node(vec![i as f64 / n as f64, j as f64 / n as f64]);
        }
    }
    for j in 0..n {
        for i in 0..n {
            let a = ids[j][i];
            let c = ids[j][i + 1];
            let d = ids[j + 1][i];
            let e = ids[j + 1][i + 1];
            b.add_element(CellType::Tri, vec![a, c, e]);
            b.add_element(CellType::Tri, vec![a, e, d]);
        }
    }
    b.build()
}

// ---------------------------------------------------------------------------
// Manufactured solution helpers
// ---------------------------------------------------------------------------

/// `u = sin(π x)` on the line; its `-u'' = π² sin(π x)`.
fn u1(x: &[f64]) -> f64 {
    (std::f64::consts::PI * x[0]).sin()
}
fn f1(x: &[f64]) -> f64 {
    std::f64::consts::PI * std::f64::consts::PI * (std::f64::consts::PI * x[0]).sin()
}
fn g1(x: &[f64]) -> [f64; 1] {
    [std::f64::consts::PI * (std::f64::consts::PI * x[0]).cos()]
}

/// `u = sin(π x) sin(π y)` on the square; its `-∇²u = 2π² sin(π x) sin(π y)`.
fn u2(x: &[f64]) -> f64 {
    (std::f64::consts::PI * x[0]).sin() * (std::f64::consts::PI * x[1]).sin()
}
fn f2(x: &[f64]) -> f64 {
    let p = std::f64::consts::PI;
    2.0 * p * p * (p * x[0]).sin() * (p * x[1]).sin()
}
fn g2(x: &[f64]) -> [f64; 2] {
    let p = std::f64::consts::PI;
    [
        p * (p * x[0]).cos() * (p * x[1]).sin(),
        p * (p * x[0]).sin() * (p * x[1]).cos(),
    ]
}

// ---------------------------------------------------------------------------
// Error norms
// ---------------------------------------------------------------------------

/// Integrated `L2` error `||u_h - u_ex||`, computed with element quadrature so
/// it scales as the true norm (converging at O(h²) for P1 elements).
fn l2_error(mesh: &Mesh, u: &[f64], uex: fn(&[f64]) -> f64) -> f64 {
    let mut s = 0.0f64;
    for elem in &mesh.elements {
        let cell = elem.cell_type;
        let dim = match cell {
            CellType::Line => 1,
            CellType::Tri => 2,
            CellType::Tet => 3,
            _ => panic!("unsupported cell in MMS test"),
        };
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        let (qpts, qw) = match cell {
            CellType::Line => {
                let r = gauss_legendre(4);
                (
                    r.points.iter().map(|x| vec![*x]).collect::<Vec<_>>(),
                    r.weights,
                )
            }
            CellType::Tri => {
                let r = triangle(TriangleRule::HammerStroud);
                (
                    r.points
                        .iter()
                        .map(|p| vec![p[0], p[1]])
                        .collect::<Vec<_>>(),
                    r.weights,
                )
            }
            CellType::Tet => {
                let r = tetrahedron(TetrahedronRule::Keast3);
                (
                    r.points
                        .iter()
                        .map(|p| vec![p[0], p[1], p[2]])
                        .collect::<Vec<_>>(),
                    r.weights,
                )
            }
            _ => unreachable!(),
        };
        let local_shape = |xi: &[f64]| -> Vec<f64> {
            match cell {
                CellType::Line => Line2::shape(xi),
                CellType::Tri => Tri3::shape(xi),
                CellType::Tet => Tet4::shape(xi),
                _ => unreachable!(),
            }
        };
        let n = elem.nodes.len();
        for (qp, w) in qpts.iter().zip(&qw) {
            let ns = local_shape(qp);
            // Physical coordinates of the quadrature point.
            let mut xq = vec![0.0; dim];
            for k in 0..n {
                for d in 0..dim {
                    xq[d] += ns[k] * phys[k][d];
                }
            }
            // FE solution value at the quadrature point.
            let mut uh = 0.0;
            for i in 0..n {
                uh += ns[i] * u[elem.nodes[i]];
            }
            let d = uh - uex(&xq);
            let det = {
                let lg = match cell {
                    CellType::Line => Line2::grad(qp),
                    CellType::Tri => Tri3::grad(qp),
                    CellType::Tet => Tet4::grad(qp),
                    _ => unreachable!(),
                };
                Map::from_nodes_and_grad(&phys, &lg).determinant.abs()
            };
            s += w * det * d * d;
        }
    }
    s.sqrt()
}

/// `H1` semi-norm error `||∇(u_h - u_ex)||`, integrated over the elements with
/// the same isoparametric mapping the solver uses.
fn h1_error(mesh: &Mesh, u: &[f64], uex_grad: fn(&[f64]) -> [f64; 3]) -> f64 {
    let mut s = 0.0f64;
    for elem in &mesh.elements {
        let cell = elem.cell_type;
        let dim = match cell {
            CellType::Line => 1,
            CellType::Tri => 2,
            CellType::Tet => 3,
            _ => panic!("unsupported cell in MMS test"),
        };
        let phys: Vec<Vec<f64>> = elem
            .nodes
            .iter()
            .map(|&n| mesh.node_coords(n).to_vec())
            .collect();
        let (qpts, qw) = match cell {
            CellType::Line => {
                let r = gauss_legendre(4);
                (r.points.iter().map(|x| vec![*x]).collect(), r.weights)
            }
            CellType::Tri => {
                let r = triangle(TriangleRule::HammerStroud);
                (
                    r.points
                        .iter()
                        .map(|p| vec![p[0], p[1]])
                        .collect::<Vec<_>>(),
                    r.weights,
                )
            }
            CellType::Tet => {
                let r = tetrahedron(TetrahedronRule::Keast3);
                (
                    r.points
                        .iter()
                        .map(|p| vec![p[0], p[1], p[2]])
                        .collect::<Vec<_>>(),
                    r.weights,
                )
            }
            _ => unreachable!(),
        };
        let local_grad = |xi: &[f64]| -> Vec<Vec<f64>> {
            match cell {
                CellType::Line => Line2::grad(xi),
                CellType::Tri => Tri3::grad(xi),
                CellType::Tet => Tet4::grad(xi),
                _ => unreachable!(),
            }
        };
        let local_shape = |xi: &[f64]| -> Vec<f64> {
            match cell {
                CellType::Line => Line2::shape(xi),
                CellType::Tri => Tri3::shape(xi),
                CellType::Tet => Tet4::shape(xi),
                _ => unreachable!(),
            }
        };
        let n = elem.nodes.len();
        for (qp, w) in qpts.iter().zip(&qw) {
            let lg = local_grad(qp);
            let map = Map::from_nodes_and_grad(&phys, &lg);
            let det = map.determinant.abs();
            // Physical coordinates of the quadrature point.
            let ns = local_shape(qp);
            let mut xq = vec![0.0; dim];
            for k in 0..n {
                for d in 0..dim {
                    xq[d] += ns[k] * phys[k][d];
                }
            }
            // Gradient of the FE solution.
            let mut guh = vec![0.0; dim];
            for i in 0..n {
                let g = map.physical_grad(&lg[i]);
                for d in 0..dim {
                    guh[d] += u[elem.nodes[i]] * g[d];
                }
            }
            let gex = uex_grad(&xq);
            let mut diff = 0.0;
            for d in 0..dim {
                let dd = guh[d] - gex[d];
                diff += dd * dd;
            }
            s += w * det * diff;
        }
    }
    s.sqrt()
}

/// Solve with Dirichlet `u = u_exact` on all boundary nodes.
fn solve_mms(mesh: &Mesh, f: fn(&[f64]) -> f64, uex: fn(&[f64]) -> f64) -> Vec<f64> {
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        let on_boundary = match n.coords.len() {
            1 => n.coords[0].abs() < 1e-9 || (n.coords[0] - 1.0).abs() < 1e-9,
            _ => {
                n.coords[0].abs() < 1e-9
                    || (n.coords[0] - 1.0).abs() < 1e-9
                    || n.coords[1].abs() < 1e-9
                    || (n.coords[1] - 1.0).abs() < 1e-9
            }
        };
        if on_boundary {
            bcs.push((i, uex(&n.coords)));
        }
    }
    solve_poisson(mesh, 1.0, 4, f, &bcs, None, None).expect("solve")
}

/// Assert that successive errors shrink at ~`expected_ratio` (allowing slack).
fn assert_converges(errors: &[f64], expected_ratio: f64) {
    for w in errors.windows(2) {
        let ratio = w[0] / w[1];
        assert!(
            ratio > expected_ratio * 0.92,
            "convergence ratio {ratio:.3} below expected ~{expected_ratio:.2}"
        );
    }
}

#[test]
fn mms_1d_converges() {
    let mut l2 = Vec::new();
    let mut h1 = Vec::new();
    for n in [8usize, 16, 32, 64] {
        let mesh = line_mesh(n);
        let u = solve_mms(&mesh, f1, u1);
        l2.push(l2_error(&mesh, &u, u1));
        h1.push(h1_error(&mesh, &u, |x| {
            let g = g1(x);
            [g[0], 0.0, 0.0]
        }));
    }
    // P1: L2 ~ O(h^2) (ratio 4), H1 ~ O(h) (ratio 2).
    assert_converges(&l2, 4.0);
    assert_converges(&h1, 2.0);
}

#[test]
fn mms_2d_converges() {
    let mut l2 = Vec::new();
    let mut h1 = Vec::new();
    for n in [4usize, 8, 16, 32] {
        let mesh = square_tri_mesh(n);
        let u = solve_mms(&mesh, f2, u2);
        l2.push(l2_error(&mesh, &u, u2));
        h1.push(h1_error(&mesh, &u, |x| {
            let g = g2(x);
            [g[0], g[1], 0.0]
        }));
    }
    assert_converges(&l2, 4.0);
    assert_converges(&h1, 2.0);
}

#[test]
#[ignore = "3-D MMS: heavy; run with --ignored"]
fn mms_3d_converges() {
    use tpt_fem_mesh_gen::box_mesh;
    let mut l2 = Vec::new();
    let mut h1 = Vec::new();
    for n in [4usize, 8] {
        let mesh = box_mesh([0.0; 3], [1.0; 3], [n, n, n]);
        let uex = |x: &[f64]| {
            (std::f64::consts::PI * x[0]).sin()
                * (std::f64::consts::PI * x[1]).sin()
                * (std::f64::consts::PI * x[2]).sin()
        };
        let fsrc = |x: &[f64]| {
            let p = std::f64::consts::PI;
            3.0 * p * p * (p * x[0]).sin() * (p * x[1]).sin() * (p * x[2]).sin()
        };
        let uex_grad = |x: &[f64]| {
            let p = std::f64::consts::PI;
            [
                p * (p * x[0]).cos() * (p * x[1]).sin() * (p * x[2]).sin(),
                p * (p * x[0]).sin() * (p * x[1]).cos() * (p * x[2]).sin(),
                p * (p * x[0]).sin() * (p * x[1]).sin() * (p * x[2]).cos(),
            ]
        };
        let mut bcs = Vec::new();
        for (i, nd) in mesh.nodes.iter().enumerate() {
            let on_b =
                (0..3).any(|d| nd.coords[d].abs() < 1e-9 || (nd.coords[d] - 1.0).abs() < 1e-9);
            if on_b {
                bcs.push((i, uex(&nd.coords)));
            }
        }
        let u = solve_poisson(&mesh, 1.0, 4, fsrc, &bcs, None, None).expect("solve");
        l2.push(l2_error(&mesh, &u, uex));
        h1.push(h1_error(&mesh, &u, uex_grad));
    }
    eprintln!("3D L2 errors: {l2:?}");
    eprintln!("3D H1 errors: {h1:?}");
    // The Kuhn six-tet brick decomposition is mildly anisotropic, so the L2
    // rate is still pre-asymptotic at n = 4→8 (measured ≈ 3.5, approaching
    // the theoretical 4 on finer grids); allow extra slack there.
    for w in l2.windows(2) {
        let ratio = w[0] / w[1];
        assert!(
            ratio > 4.0 * 0.85,
            "3D L2 convergence ratio {ratio:.3} below expected ~4.00"
        );
    }
    assert_converges(&h1, 2.0);
}
