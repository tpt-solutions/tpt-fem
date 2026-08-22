//! Contact mechanics for `tpt-fem`.
//!
//! Provides unilateral (non-penetration) contact enforcement written from
//! scratch: a **penalty** method and an **augmented-Lagrangian** iteration, both
//! operating on a caller-supplied base stiffness [`Coo`] plus a set of
//! [`ContactConstraint`]s (each a global DOF constrained to stay at or above a
//! lower bound — e.g. a node's normal coordinate relative to a rigid obstacle).
//!
//! A minimal from-scratch [`contact_pairs`] helper performs brute-force
//! nearest-node pairing between two surfaces (a full BVH/octree is a future
//! optimisation; the constraint enforcement below is independent of the search).
//!
//! ```
//! use tpt_fem_contact::{augmented_lagrangian, ContactConstraint};
//! use tpt_fem_sparse::Coo;
//!
//! // One DOF on a spring (k=10) to ground, pushed into a wall at x=0 by F=-4
//! // (force toward −x). Hard contact should hold x≈0 with reaction λ≈4.
//! let base = Coo { rows: vec![0], cols: vec![0], vals: vec![10.0] };
//! let load = vec![-4.0];
//! let con = ContactConstraint { dof: 0, lower: 0.0 };
//! let (u, lambda) = augmented_lagrangian(&base, &load, &[con], 1e4, 50, 1e-9);
//! assert!(u[0].abs() < 1e-6);
//! assert!((lambda[0] - 4.0).abs() < 1e-3);
//! ```

use tpt_fem_sparse::{solve, Coo};

/// A unilateral constraint `x_dof ≥ lower` (non-penetration against a rigid
/// obstacle located at `lower` along the constrained DOF's axis).
#[derive(Clone, Copy, Debug)]
pub struct ContactConstraint {
    /// Global DOF the constraint acts on.
    pub dof: usize,
    /// Lower bound the DOF must satisfy (`x ≥ lower`).
    pub lower: f64,
}

/// Add a penalty contact contribution to the global system.
///
/// For each constraint, the energy penalty `½ κ (min(0, lower − x))²` contributes
/// `κ` to the diagonal stiffness at `dof` and a load `κ·lower` (so the
/// equilibrium shifts toward the obstacle). Returns the augmented `(stiffness,
/// load)`.
pub fn penalty_contact(
    base: &Coo,
    load: &[f64],
    constraints: &[ContactConstraint],
    penalty: f64,
) -> (Coo, Vec<f64>) {
    let mut k = base.clone();
    let mut f = load.to_vec();
    for c in constraints {
        k.push(c.dof, c.dof, penalty);
        f[c.dof] += penalty * c.lower;
    }
    (k, f)
}

/// Augmented-Lagrangian contact iteration.
///
/// Each step solves the penalty-augmented system with the *current* Lagrange
/// multipliers folded into the load, then updates
/// `λ ← max(0, λ + κ·(lower − x))` (the standard ALM multiplier update for a
/// unilateral constraint). Converges to the hard-contact solution (exact
/// non-penetration) for a sufficiently large `penalty`.
pub fn augmented_lagrangian(
    base: &Coo,
    load: &[f64],
    constraints: &[ContactConstraint],
    penalty: f64,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut lambda = vec![0.0; constraints.len()];
    let mut u = vec![0.0; load.len()];
    for _ in 0..max_iter {
        // Fold multipliers into the load: f_eff = f + Σ λ n  (n = +1 for x ≥ lower).
        let mut f_eff = load.to_vec();
        for (i, c) in constraints.iter().enumerate() {
            f_eff[c.dof] += lambda[i];
        }
        let (k_aug, f_aug) = penalty_contact(base, &f_eff, constraints, penalty);
        let u_new = solve(&k_aug, &f_aug).expect("contact augmented system must solve");
        // Multiplier update on the violation (lower − x).
        let mut max_viol = 0.0_f64;
        for (i, c) in constraints.iter().enumerate() {
            let viol = c.lower - u_new[c.dof];
            let new_lambda = (lambda[i] + penalty * viol).max(0.0);
            lambda[i] = new_lambda;
            max_viol = max_viol.max(viol.abs());
        }
        u = u_new;
        if max_viol < tol {
            break;
        }
    }
    (u, lambda)
}

/// Brute-force nearest-node pairing between two surfaces `a` and `b` (each a
/// list of `(node_id, coords)`).
///
/// For every node in `a`, returns `(node_id, Some((index_in_b, gap)))` for its
/// closest node in `b`, or `(node_id, None)` if `b` is empty (there is nothing
/// to pair with — callers must handle the `None` rather than indexing a
/// sentinel). The scan is O(|a|·|b|) and dependency-free; a spatial acceleration
/// structure (BVH/octree) is a future optimisation and the constraint
/// enforcement below is independent of the search.
pub fn contact_pairs(
    a: &[(usize, Vec<f64>)],
    b: &[(usize, Vec<f64>)],
) -> Vec<(usize, Option<(usize, f64)>)> {
    let mut out = Vec::with_capacity(a.len());
    for (na, ca) in a.iter() {
        let mut best: Option<(usize, f64)> = None;
        for (ib, (_nb, cb)) in b.iter().enumerate() {
            let d: f64 = ca
                .iter()
                .zip(cb)
                .map(|(x, y)| (x - y) * (x - y))
                .sum::<f64>()
                .sqrt();
            if best.map_or(true, |(_, bd)| d < bd) {
                best = Some((ib, d));
            }
        }
        out.push((*na, best));
    }
    out
}

/// Uniform octree over a fixed point set, for fast nearest-neighbour queries.
///
/// The root cell is the smallest axis-aligned cube enclosing all points;
/// cells are subdivided into 8 children until they hold at most `max_leaf`
/// points or reach `max_depth`. A nearest-neighbour query descends cells in
/// order of their minimum possible distance to the query point, pruning any
/// subtree whose bounding cube cannot beat the current best — O(log n) per
/// query on evenly distributed points versus O(n) for [`contact_pairs`].
#[derive(Debug)]
pub struct Octree {
    /// Point coordinates, normalised internally to 3-D (2-D points are
    /// zero-padded).
    pts: Vec<[f64; 3]>,
    /// Flattened cell storage; `children` indices refer back into this vec.
    cells: Vec<OctCell>,
}

#[derive(Debug)]
struct OctCell {
    center: [f64; 3],
    half: f64,
    /// `Some([child; 8])` for interior cells (child order: bit 0 = −x/+x,
    /// bit 1 = −y/+y, bit 2 = −z/+z), `None` for leaves holding `points`.
    children: Option<[usize; 8]>,
    points: Vec<usize>,
}

impl Octree {
    fn min_dist2(cell: &OctCell, q: &[f64; 3]) -> f64 {
        let mut d2 = 0.0;
        for a in 0..3 {
            let delta = (q[a] - cell.center[a]).abs() - cell.half;
            if delta > 0.0 {
                d2 += delta * delta;
            }
        }
        d2
    }

    /// Build an octree over `points` (any dimensionality ≤ 3).
    pub fn new(points: &[Vec<f64>], max_leaf: usize, max_depth: u32) -> Octree {
        let pts: Vec<[f64; 3]> = points
            .iter()
            .map(|p| {
                [
                    p.first().copied().unwrap_or(0.0),
                    p.get(1).copied().unwrap_or(0.0),
                    p.get(2).copied().unwrap_or(0.0),
                ]
            })
            .collect();
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [-f64::INFINITY; 3];
        for p in &pts {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        if pts.is_empty() {
            return Octree {
                pts,
                cells: Vec::new(),
            };
        }
        let center = [
            (lo[0] + hi[0]) / 2.0,
            (lo[1] + hi[1]) / 2.0,
            (lo[2] + hi[2]) / 2.0,
        ];
        let half = hi
            .iter()
            .zip(&lo)
            .map(|(h, l)| (h - l) / 2.0)
            .fold(0.0_f64, f64::max)
            .max(1e-12);
        let mut tree = Octree {
            pts,
            cells: Vec::new(),
        };
        let all: Vec<usize> = (0..tree.pts.len()).collect();
        tree.subdivide(center, half, all, 0, max_leaf.max(1), max_depth);
        tree
    }

    fn subdivide(
        &mut self,
        center: [f64; 3],
        half: f64,
        points: Vec<usize>,
        depth: u32,
        max_leaf: usize,
        max_depth: u32,
    ) -> usize {
        if points.len() <= max_leaf || depth >= max_depth {
            self.cells.push(OctCell {
                center,
                half,
                children: None,
                points,
            });
            return self.cells.len() - 1;
        }
        let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); 8];
        for &pi in &points {
            let p = self.pts[pi];
            let mut code = 0usize;
            for a in 0..3 {
                if p[a] >= center[a] {
                    code |= 1 << a;
                }
            }
            buckets[code].push(pi);
        }
        // Reserve this slot first so children can reference siblings.
        self.cells.push(OctCell {
            center,
            half,
            children: None,
            points: Vec::new(),
        });
        let my_idx = self.cells.len() - 1;
        let child_half = half / 2.0;
        let mut children = [my_idx; 8];
        for (oct, bucket) in buckets.into_iter().enumerate() {
            if bucket.is_empty() {
                continue;
            }
            let mut cc = center;
            for a in 0..3 {
                cc[a] += if oct & (1 << a) != 0 {
                    child_half
                } else {
                    -child_half
                };
            }
            let ci = self.subdivide(cc, child_half, bucket, depth + 1, max_leaf, max_depth);
            children[oct] = ci;
        }
        self.cells[my_idx].children = Some(children);
        my_idx
    }

    /// Index of the nearest point to `q`, with its squared distance
    /// (`None` only for an empty tree).
    pub fn nearest(&self, q: &[f64]) -> Option<(usize, f64)> {
        if self.cells.is_empty() {
            return None;
        }
        let qq = [
            q.first().copied().unwrap_or(0.0),
            q.get(1).copied().unwrap_or(0.0),
            q.get(2).copied().unwrap_or(0.0),
        ];
        let mut best: Option<(usize, f64)> = None;
        // DFS with bound pruning: every cell whose cube cannot hold a point
        // closer than the current best is skipped.
        let mut stack = vec![0usize];
        while let Some(ci) = stack.pop() {
            let cell = &self.cells[ci];
            if let Some((_, bd)) = best {
                if Self::min_dist2(cell, &qq) > bd {
                    continue;
                }
            }
            match &cell.children {
                None => {
                    for &pi in &cell.points {
                        let d2: f64 = (0..3)
                            .map(|a| {
                                let e = qq[a] - self.pts[pi][a];
                                e * e
                            })
                            .sum();
                        if best.map_or(true, |(_, bd)| d2 < bd) {
                            best = Some((pi, d2));
                        }
                    }
                }
                Some(children) => {
                    // Visit nearer children first.
                    let mut order: Vec<usize> =
                        children.iter().copied().filter(|&c| c != ci).collect();
                    order.sort_by(|&x, &y| {
                        let dx = Self::min_dist2(&self.cells[x], &qq);
                        let dy = Self::min_dist2(&self.cells[y], &qq);
                        dx.partial_cmp(&dy).unwrap_or(std::cmp::Ordering::Equal)
                    });
                    stack.extend(order);
                }
            }
        }
        best
    }

    /// Number of points indexed.
    pub fn len(&self) -> usize {
        self.pts.len()
    }

    /// Whether the tree indexes no points.
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }
}

/// Octree-accelerated nearest-node pairing between two surfaces — same
/// contract as [`contact_pairs`], but builds an [`Octree`] over surface `b`
/// once and answers each query in O(log |b|) instead of O(|b|).
///
/// Results are identical to [`contact_pairs`] up to tie-breaking between
/// exactly equidistant candidates.
pub fn contact_pairs_octree(
    a: &[(usize, Vec<f64>)],
    b: &[(usize, Vec<f64>)],
) -> Vec<(usize, Option<(usize, f64)>)> {
    let coords: Vec<Vec<f64>> = b.iter().map(|(_, c)| c.clone()).collect();
    let tree = Octree::new(&coords, 8, 24);
    a.iter()
        .map(|(na, ca)| (*na, tree.nearest(ca).map(|(i, d2)| (i, d2.sqrt()))))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn penalty_penetrates_wall() {
        // Soft penalty lets the node sink slightly below the wall.
        let base = Coo {
            rows: vec![0],
            cols: vec![0],
            vals: vec![10.0],
        };
        let load = vec![-4.0];
        let con = ContactConstraint { dof: 0, lower: 0.0 };
        let (k, f) = penalty_contact(&base, &load, &[con], 1e3);
        let u = solve(&k, &f).unwrap();
        // x = (F + κ·lower)/(k+κ) < 0 (penetration).
        assert!(u[0] < 0.0);
        assert!(u[0] > -0.1);
    }

    #[test]
    fn augmented_lagrangian_hard_contact() {
        let base = Coo {
            rows: vec![0],
            cols: vec![0],
            vals: vec![10.0],
        };
        let load = vec![-4.0];
        let con = ContactConstraint { dof: 0, lower: 0.0 };
        let (u, lambda) = augmented_lagrangian(&base, &load, &[con], 1e4, 100, 1e-9);
        assert!(
            u[0].abs() < 1e-6,
            "x should be held at the wall, got {}",
            u[0]
        );
        assert!(
            (lambda[0] - 4.0).abs() < 1e-2,
            "reaction should balance load, got {}",
            lambda[0]
        );
    }

    #[test]
    fn augmented_lagrangian_two_dof() {
        // Two springs in series to ground (k=5 each), middle node pushed by
        // F=-3 into a wall at x=0. Hard contact: middle node x≈0, reaction λ≈3.
        let base = Coo {
            rows: vec![0, 0, 1, 1],
            cols: vec![0, 1, 0, 1],
            vals: vec![5.0, -5.0, -5.0, 5.0],
        };
        let load = vec![0.0, -3.0];
        let con = ContactConstraint { dof: 1, lower: 0.0 };
        let (u, lambda) = augmented_lagrangian(&base, &load, &[con], 1e4, 200, 1e-8);
        assert!(u[1].abs() < 1e-6);
        assert!((lambda[0] - 3.0).abs() < 1e-2);
    }

    #[test]
    fn contact_pairs_finds_nearest() {
        let a = vec![(0usize, vec![0.0, 0.0])];
        let b = vec![(1usize, vec![1.0, 0.0]), (2usize, vec![0.1, 0.0])];
        let pairs = contact_pairs(&a, &b);
        assert_eq!(pairs.len(), 1);
        let (na, best) = &pairs[0];
        assert_eq!(*na, 0);
        let (ib, d) = best.expect("surface b is non-empty");
        assert_eq!(ib, 1); // node 2 (index 1) is closest
        assert!((d - 0.1).abs() < 1e-12);
    }

    #[test]
    fn contact_pairs_empty_surface_is_none() {
        let a = vec![(0usize, vec![0.0, 0.0])];
        let b: Vec<(usize, Vec<f64>)> = Vec::new();
        let pairs = contact_pairs(&a, &b);
        assert_eq!(pairs.len(), 1);
        assert!(pairs[0].1.is_none());
    }

    /// Deterministic LCG pseudo-random coordinates in [0,1)^dim.
    fn lcg_points(n: usize, dim: usize) -> Vec<Vec<f64>> {
        let mut s: u64 = 0x9E3779B97F4A7C15;
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let mut p = Vec::with_capacity(dim);
            for _ in 0..dim {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                p.push((s >> 11) as f64 / (1u64 << 53) as f64);
            }
            out.push(p);
        }
        out
    }

    #[test]
    fn octree_matches_brute_force_nearest() {
        // 400 points in a unit cube; every query point must find the same
        // nearest neighbour (and distance) as the brute-force scan.
        let pts = lcg_points(400, 3);
        let queries = lcg_points(60, 3);
        let tree = Octree::new(&pts, 8, 24);
        assert_eq!(tree.len(), 400);
        for q in &queries {
            let brute = contact_pairs(
                &[(0usize, q.clone())],
                &(0..pts.len())
                    .map(|i| (i, pts[i].clone()))
                    .collect::<Vec<_>>(),
            )[0]
            .1
            .unwrap();
            let fast = tree.nearest(q).unwrap();
            assert_eq!(fast.0, brute.0, "different node for {q:?}");
            assert!((fast.1.sqrt() - brute.1).abs() < 1e-12);
        }
    }

    #[test]
    fn octree_2d_and_edge_cases() {
        // 2-D points (zero-padded internally) and duplicate points must work;
        // an empty tree reports no neighbour.
        let pts = vec![vec![0.0, 0.0], vec![2.0, 2.0], vec![2.0, 2.0]];
        let tree = Octree::new(&pts, 1, 16);
        let (i, d2) = tree.nearest(&[1.9, 1.9]).unwrap();
        assert!(i >= 1 && d2.sqrt() < 0.15, "got ({i}, {})", d2.sqrt());
        assert!(tree.nearest(&[10.0, 10.0]).unwrap().0 >= 1);
        let empty = Octree::new(&[], 4, 8);
        assert!(empty.is_empty());
        assert!(empty.nearest(&[0.0, 0.0]).is_none());
    }

    #[test]
    fn contact_pairs_octree_agrees_with_brute_force() {
        let pa = lcg_points(50, 3);
        let pb = lcg_points(500, 3);
        let a: Vec<(usize, Vec<f64>)> = pa
            .iter()
            .enumerate()
            .map(|(i, c)| (100 + i, c.clone()))
            .collect();
        let b: Vec<(usize, Vec<f64>)> = pb
            .iter()
            .enumerate()
            .map(|(i, c)| (2000 + i, c.clone()))
            .collect();
        let brute = contact_pairs(&a, &b);
        let fast = contact_pairs_octree(&a, &b);
        assert_eq!(brute.len(), fast.len());
        for ((_, x), (_, y)) in brute.iter().zip(&fast) {
            match (x, y) {
                (Some((bi, bd)), Some((fi, fd))) => {
                    assert_eq!(bi, fi, "index mismatch");
                    assert!((bd - fd).abs() < 1e-12);
                }
                (None, None) => {}
                other => panic!("option mismatch: {other:?}"),
            }
        }
    }
}
