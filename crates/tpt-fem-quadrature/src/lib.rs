//! Gauss quadrature rules for finite-element reference elements.
#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_range_loop)]
//!
//! This crate provides exact fixed-order quadrature rules on the standard
//! reference elements used by `tpt-fem-element`:
//!
//! * 1-D Gauss–Legendre rules of orders 1–5 on `[-1, 1]` (and on `[0, 1]`),
//! * tensor-product rules on the square `[-1, 1]²` and cube `[-1, 1]³`,
//! * fixed low-order rules on the reference triangle `(0,0),(1,0),(0,1)`,
//! * fixed low-order rules on the reference tetrahedron
//!   `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`.
//!
//! Every rule is exact for polynomials up to a stated degree and is verified
//! by the unit tests, which integrate monomials against their closed-form
//! values.
//!
//! # Example
//!
//! ```
//! use tpt_fem_quadrature::gauss_legendre;
//!
//! // Order-2 Gauss-Legendre is exact for cubics on [-1, 1].
//! let rule = gauss_legendre(2);
//! let approx: f64 = rule.weights.iter().zip(&rule.points)
//!     .map(|(w, x)| w * x * x * x)
//!     .sum();
//! assert!((approx - 0.0).abs() < 1e-12);
//! ```

/// A 1-D quadrature rule: evaluation points and matching weights.
#[derive(Clone, Debug, PartialEq)]
pub struct Quad1D {
    /// Points on the integration interval.
    pub points: Vec<f64>,
    /// Quadrature weights, one per point.
    pub weights: Vec<f64>,
}

/// A 2-D quadrature rule on a reference element.
#[derive(Clone, Debug, PartialEq)]
pub struct Quad2D {
    /// Points `(x, y)` on the reference element.
    pub points: Vec<[f64; 2]>,
    /// Quadrature weights, one per point.
    pub weights: Vec<f64>,
}

/// A 3-D quadrature rule on a reference element.
#[derive(Clone, Debug, PartialEq)]
pub struct Quad3D {
    /// Points `(x, y, z)` on the reference element.
    pub points: Vec<[f64; 3]>,
    /// Quadrature weights, one per point.
    pub weights: Vec<f64>,
}

impl Quad1D {
    /// Number of quadrature points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if the rule has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Sum of the weights (should equal the reference interval length, `2`).
    pub fn weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }
}

impl Quad2D {
    /// Number of quadrature points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if the rule has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Sum of the weights (equals the reference-element area, `0.5` on the
    /// standard triangle, `4` on the standard square).
    pub fn weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }
}

impl Quad3D {
    /// Number of quadrature points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if the rule has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Sum of the weights (equals the reference-element volume, `1/6` on the
    /// standard tetrahedron, `8` on the standard cube).
    pub fn weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }
}

/// Gauss–Legendre quadrature of the given order (1–5) on `[-1, 1]`.
///
/// The rule is exact for polynomials of degree `2*order - 1`.
///
/// # Panics
///
/// Panics if `order` is not in `1..=5`.
pub fn gauss_legendre(order: usize) -> Quad1D {
    let (points, weights): (&[f64], &[f64]) = match order {
        1 => (&[0.0], &[2.0]),
        2 => (&[-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()], &[1.0, 1.0]),
        3 => (
            &[0.0, -(3.0_f64 / 5.0).sqrt(), (3.0_f64 / 5.0).sqrt()],
            &[8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0],
        ),
        4 => (
            &[
                -0.3399810435848563,
                0.3399810435848563,
                -0.8611363115940526,
                0.8611363115940526,
            ],
            &[
                0.6521451548625461,
                0.6521451548625461,
                0.3478548451374539,
                0.3478548451374539,
            ],
        ),
        5 => (
            &[
                0.0,
                -0.5384693101056831,
                0.5384693101056831,
                -0.9061798459386640,
                0.9061798459386640,
            ],
            &[
                0.5688888888888889,
                0.4786286704993665,
                0.4786286704993665,
                0.2369268850561891,
                0.2369268850561891,
            ],
        ),
        _ => panic!("gauss_legendre: order must be in 1..=5, got {order}"),
    };
    Quad1D {
        points: points.to_vec(),
        weights: weights.to_vec(),
    }
}

/// Gauss–Legendre quadrature of the given order (1–5) mapped to `[0, 1]`.
///
/// The rule is exact for polynomials of degree `2*order - 1` on `[0, 1]`.
///
/// # Panics
///
/// Panics if `order` is not in `1..=5`.
pub fn gauss_legendre_unit(order: usize) -> Quad1D {
    let r = gauss_legendre(order);
    Quad1D {
        points: r.points.iter().map(|x| (x + 1.0) * 0.5).collect(),
        weights: r.weights.iter().map(|w| w * 0.5).collect(),
    }
}

/// Tensor-product quadrature on the square `[-1, 1]²` from a 1-D rule.
pub fn tensor_square(rule: &Quad1D) -> Quad2D {
    let n = rule.points.len();
    let mut points = Vec::with_capacity(n * n);
    let mut weights = Vec::with_capacity(n * n);
    for (xi, wi) in rule.points.iter().zip(&rule.weights) {
        for (xj, wj) in rule.points.iter().zip(&rule.weights) {
            points.push([*xi, *xj]);
            weights.push(wi * wj);
        }
    }
    Quad2D { points, weights }
}

/// Tensor-product quadrature on the cube `[-1, 1]³` from a 1-D rule.
pub fn tensor_cube(rule: &Quad1D) -> Quad3D {
    let n = rule.points.len();
    let mut points = Vec::with_capacity(n * n * n);
    let mut weights = Vec::with_capacity(n * n * n);
    for (xi, wi) in rule.points.iter().zip(&rule.weights) {
        for (xj, wj) in rule.points.iter().zip(&rule.weights) {
            for (xk, wk) in rule.points.iter().zip(&rule.weights) {
                points.push([*xi, *xj, *xk]);
                weights.push(wi * wj * wk);
            }
        }
    }
    Quad3D { points, weights }
}

/// Reference-triangle rule selector.
///
/// The reference triangle has vertices `(0,0)`, `(1,0)`, `(0,1)` (area `0.5`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriangleRule {
    /// 1 point, degree 1 (centroid).
    Degree1,
    /// 3 points, degree 2.
    Degree2,
    /// 7 points, Hammer–Stroud, degree 5.
    HammerStroud,
}

/// Quadrature rule on the reference triangle `(0,0),(1,0),(0,1)`.
pub fn triangle(rule: TriangleRule) -> Quad2D {
    // Points are `(x, y)`; the reference triangle is `x >= 0, y >= 0, x + y <= 1`.
    let (pts, wts): (&[[f64; 2]], &[f64]) = match rule {
        TriangleRule::Degree1 => (&[[1.0 / 3.0, 1.0 / 3.0]], &[0.5]),
        TriangleRule::Degree2 => (
            &[
                [2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0],
                [1.0 / 6.0, 1.0 / 6.0],
            ],
            &[1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
        ),
        TriangleRule::HammerStroud => (
            &[
                [3.3333333333333331e-01, 3.3333333333333331e-01],
                [5.9716185317142512e-02, 4.7014190734142874e-01],
                [4.7014190734142874e-01, 4.7014190734142874e-01],
                [4.7014190734142874e-01, 5.9716185317142512e-02],
                [7.9742706745232927e-01, 1.0128646627383536e-01],
                [1.0128646627383536e-01, 1.0128646627383536e-01],
                [1.0128646627383536e-01, 7.9742706745232927e-01],
            ],
            &[
                1.1249952634765750e-01,
                6.6197272786762434e-02,
                6.6197272786762434e-02,
                6.6197272786762434e-02,
                6.2969551762997400e-02,
                6.2969551762997400e-02,
                6.2969551762997400e-02,
            ],
        ),
    };
    Quad2D {
        points: pts.to_vec(),
        weights: wts.to_vec(),
    }
}

/// Reference-tetrahedron rule selector.
///
/// The reference tetrahedron has vertices `(0,0,0)`, `(1,0,0)`, `(0,1,0)`,
/// `(0,0,1)` (volume `1/6`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TetrahedronRule {
    /// 1 point, degree 1 (centroid).
    Degree1,
    /// 4 points, degree 2.
    Degree2,
    /// 5 points, Keast degree 3.
    Keast3,
    /// 11 points, Keast degree 4.
    Keast4,
}

/// Quadrature rule on the reference tetrahedron.
pub fn tetrahedron(rule: TetrahedronRule) -> Quad3D {
    let (pts, wts): (&[[f64; 3]], &[f64]) = match rule {
        TetrahedronRule::Degree1 => (&[[1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0]], &[1.0 / 6.0]),
        TetrahedronRule::Degree2 => {
            let a = 0.5854101966249685;
            let b = 0.1381966011250105;
            (
                &[[a, b, b], [b, a, b], [b, b, a], [b, b, b]],
                &[1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0],
            )
        }
        TetrahedronRule::Keast3 => {
            let w_c = -2.0 / 15.0;
            let w_o = 3.0 / 40.0;
            (
                &[
                    [1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0],
                    [1.0 / 2.0, 1.0 / 6.0, 1.0 / 6.0],
                    [1.0 / 6.0, 1.0 / 2.0, 1.0 / 6.0],
                    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 2.0],
                    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
                ],
                &[w_c, w_o, w_o, w_o, w_o],
            )
        }
        TetrahedronRule::Keast4 => (
            &[
                [
                    3.4755987341592781e-01,
                    8.9833494257077623e-02,
                    2.7444361775607534e-02,
                ],
                [
                    5.6324396794297564e-02,
                    1.0824401145734477e-01,
                    1.8455109205301029e-01,
                ],
                [
                    2.1308877184004529e-01,
                    2.7235826324187240e-01,
                    3.1429349534962059e-01,
                ],
                [
                    3.9737822358754359e-02,
                    1.0367439369541243e-01,
                    7.0201200748180603e-01,
                ],
                [
                    7.4704100235526549e-01,
                    3.4504439235972663e-02,
                    1.1583086727771186e-01,
                ],
                [
                    3.6107462289211784e-01,
                    1.9351669025520701e-01,
                    7.0607393788550299e-01,
                ],
                [
                    5.7831079909913763e-02,
                    5.0453520347486192e-01,
                    3.2171662411990731e-01,
                ],
                [
                    3.2121308072437715e-01,
                    4.4103399120113559e-02,
                    4.2394778054055604e-01,
                ],
                [
                    4.8682903810900369e-01,
                    3.0856376319184292e-01,
                    9.6688241005899650e-02,
                ],
                [
                    9.1155031509020074e-02,
                    4.9741562468248268e-01,
                    3.6520406819370660e-02,
                ],
                [
                    1.3824099494800238e-01,
                    7.8071364095426199e-01,
                    5.5402341306998303e-02,
                ],
            ],
            &[
                1.3369217678467679e-02,
                1.8226860983725623e-02,
                2.7290508173583462e-02,
                1.2871564902353203e-02,
                9.2296479881365471e-03,
                1.8129933227631694e-03,
                1.6616322747826342e-02,
                2.2159660825794909e-02,
                2.4310916171224342e-02,
                1.4220267110531069e-02,
                6.5587066189644319e-03,
            ],
        ),
    };
    Quad3D {
        points: pts.to_vec(),
        weights: wts.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Closed-form integrals over the reference triangle x,y>=0, x+y<=1.
    fn tri_int(i: usize, j: usize) -> f64 {
        // int x^i y^j = i! j! / (i+j+2)!
        let f = |n: usize| (1..=n).product::<usize>() as f64;
        f(i) * f(j) / f(i + j + 2)
    }

    // Closed-form integrals over the reference tetrahedron.
    fn tet_int(i: usize, j: usize, k: usize) -> f64 {
        let f = |n: usize| (1..=n).product::<usize>() as f64;
        f(i) * f(j) * f(k) / f(i + j + k + 3)
    }

    fn approx_1d(r: &Quad1D, i: usize) -> f64 {
        r.weights
            .iter()
            .zip(&r.points)
            .map(|(w, x)| w * x.powi(i as i32))
            .sum()
    }

    fn approx_tri(r: &Quad2D, i: usize, j: usize) -> f64 {
        r.weights
            .iter()
            .zip(&r.points)
            .map(|(w, p)| w * p[0].powi(i as i32) * p[1].powi(j as i32))
            .sum()
    }

    fn approx_tet(r: &Quad3D, i: usize, j: usize, k: usize) -> f64 {
        r.weights
            .iter()
            .zip(&r.points)
            .map(|(w, p)| w * p[0].powi(i as i32) * p[1].powi(j as i32) * p[2].powi(k as i32))
            .sum()
    }

    #[test]
    fn gauss_legendre_degrees() {
        for order in 1..=5 {
            let r = gauss_legendre(order);
            assert_eq!(r.weight_sum(), 2.0);
            let deg = 2 * order - 1;
            for m in 0..=deg {
                let exact = if m % 2 == 1 {
                    0.0
                } else {
                    2.0 / (m as f64 + 1.0)
                };
                assert!(
                    (approx_1d(&r, m) - exact).abs() < 1e-9,
                    "order {order} monomial x^{m}"
                );
            }
            // One degree beyond exactness should generally fail.
            assert!((approx_1d(&r, deg + 1) - 2.0 / (deg as f64 + 2.0)).abs() > 1e-7);
        }
    }

    #[test]
    fn gauss_legendre_unit_maps() {
        let r = gauss_legendre_unit(3);
        assert!((r.weight_sum() - 1.0).abs() < 1e-12);
        // int_0^1 x^5 dx = 1/6, order 3 is exact up to degree 5.
        assert!((approx_1d(&r, 5) - 1.0 / 6.0).abs() < 1e-9);
    }

    #[test]
    fn tensor_square_is_exact() {
        let r = tensor_square(&gauss_legendre(2));
        assert!((r.weight_sum() - 4.0).abs() < 1e-9);
        // Exact integral of x^i y^j over [-1,1]^2: 0 when either power is odd,
        // else 2/(i+1) * 2/(j+1).
        let exact = |i: usize, j: usize| -> f64 {
            let one = |n: usize| {
                if n % 2 == 1 {
                    0.0
                } else {
                    2.0 / (n as f64 + 1.0)
                }
            };
            one(i) * one(j)
        };
        for i in 0..=3 {
            for j in 0..=3 {
                let approx = approx_tri(&r, i, j);
                assert!((approx - exact(i, j)).abs() < 1e-9, "x^{i} y^{j}");
            }
        }
    }

    #[test]
    fn tensor_cube_is_exact() {
        let r = tensor_cube(&gauss_legendre(2));
        assert!((r.weight_sum() - 8.0).abs() < 1e-9);
        let exact = |i: usize, j: usize, k: usize| -> f64 {
            let one = |n: usize| {
                if n % 2 == 1 {
                    0.0
                } else {
                    2.0 / (n as f64 + 1.0)
                }
            };
            one(i) * one(j) * one(k)
        };
        for i in 0..=3 {
            for j in 0..=3 {
                for k in 0..=3 {
                    let approx = approx_tet(&r, i, j, k);
                    assert!((approx - exact(i, j, k)).abs() < 1e-9, "x^{i} y^{j} z^{k}");
                }
            }
        }
    }

    #[test]
    fn triangle_rules() {
        let t1 = triangle(TriangleRule::Degree1);
        assert!((t1.weight_sum() - 0.5).abs() < 1e-12);
        assert!((approx_tri(&t1, 1, 0) - tri_int(1, 0)).abs() < 1e-9);

        let t2 = triangle(TriangleRule::Degree2);
        assert!((t2.weight_sum() - 0.5).abs() < 1e-12);
        for i in 0..=2 {
            for j in 0..=(2 - i) {
                assert!(
                    (approx_tri(&t2, i, j) - tri_int(i, j)).abs() < 1e-9,
                    "tri deg2 x^{i} y^{j}"
                );
            }
        }

        let hs = triangle(TriangleRule::HammerStroud);
        assert!((hs.weight_sum() - 0.5).abs() < 1e-10);
        for i in 0..=5 {
            for j in 0..=(5 - i) {
                assert!(
                    (approx_tri(&hs, i, j) - tri_int(i, j)).abs() < 1e-9,
                    "HS x^{i} y^{j}"
                );
            }
        }
    }

    #[test]
    fn tetrahedron_rules() {
        let t1 = tetrahedron(TetrahedronRule::Degree1);
        assert!((t1.weight_sum() - 1.0 / 6.0).abs() < 1e-12);
        assert!((approx_tet(&t1, 1, 0, 0) - tet_int(1, 0, 0)).abs() < 1e-9);

        let t2 = tetrahedron(TetrahedronRule::Degree2);
        assert!((t2.weight_sum() - 1.0 / 6.0).abs() < 1e-9);
        for i in 0..=2 {
            for j in 0..=(2 - i) {
                for k in 0..=(2 - i - j) {
                    assert!(
                        (approx_tet(&t2, i, j, k) - tet_int(i, j, k)).abs() < 1e-9,
                        "tet deg2 x^{i} y^{j} z^{k}"
                    );
                }
            }
        }

        let k3 = tetrahedron(TetrahedronRule::Keast3);
        assert!((k3.weight_sum() - 1.0 / 6.0).abs() < 1e-9);
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                for k in 0..=(3 - i - j) {
                    assert!(
                        (approx_tet(&k3, i, j, k) - tet_int(i, j, k)).abs() < 1e-9,
                        "keast3 x^{i} y^{j} z^{k}"
                    );
                }
            }
        }

        let k4 = tetrahedron(TetrahedronRule::Keast4);
        assert!((k4.weight_sum() - 1.0 / 6.0).abs() < 1e-9);
        for i in 0..=4 {
            for j in 0..=(4 - i) {
                for k in 0..=(4 - i - j) {
                    assert!(
                        (approx_tet(&k4, i, j, k) - tet_int(i, j, k)).abs() < 1e-7,
                        "keast4 x^{i} y^{j} z^{k}"
                    );
                }
            }
        }
    }
}
