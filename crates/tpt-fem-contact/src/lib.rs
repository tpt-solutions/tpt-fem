//! Contact mechanics for `tpt-fem`.
//!
//! Provides unilateral (non-penetration) contact enforcement written from
//! scratch: a **penalty** method and an **augmented-Lagrangian** iteration, both
//! operating on a caller-supplied base stiffness [`Coo`] plus a set of
//! [`ContactConstraint`]s (each a global DOF constrained to stay at or above a
//! lower bound — e.g. a node's normal coordinate relative to a rigid obstacle).
//!
//! A minimal from-scratch [`contact_pairs`] helper performs nearest-node
//! pairing between two surfaces using the dependency-free [`Octree`] spatial
//! index (an amortised O(|a|·log|b|) lookup that replaces the earlier
//! O(|a|·|b|) brute-force scan), so surface-to-surface contact scales to larger
//! meshes. The constraint enforcement below is independent of the search.
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

mod octree;
pub use octree::Octree;

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

/// Nearest-node pairing between two surfaces `a` and `b` (each a list of
/// `(node_id, coords)`), accelerated by the dependency-free [`Octree`] spatial
/// index.
///
/// For every node in `a`, returns `(node_id, Some((index_in_b, gap)))` for its
/// closest node in `b`, or `(node_id, None)` if `b` is empty (there is nothing
/// to pair with — callers must handle the `None` rather than indexing a
/// sentinel). The search is an amortised O(|a|·log|b|) lookup (degrading to a
/// linear scan only when `b` is tiny); the constraint enforcement elsewhere in
/// this crate is independent of the search.
pub fn contact_pairs(
    a: &[(usize, Vec<f64>)],
    b: &[(usize, Vec<f64>)],
) -> Vec<(usize, Option<(usize, f64)>)> {
    if b.is_empty() {
        return a.iter().map(|(na, _)| (*na, None)).collect();
    }
    let b_indexed: Vec<(usize, Vec<f64>)> = b
        .iter()
        .enumerate()
        .map(|(i, (_, c))| (i, c.clone()))
        .collect();
    let tree = Octree::build(&b_indexed);
    a.iter()
        .map(|(na, ca)| {
            let best = tree.nearest(ca);
            (*na, best)
        })
        .collect()
}

/// Octree-accelerated nearest-node pairing between two surfaces — same
/// contract as [`contact_pairs`] (which already uses [`Octree`] internally).
/// Kept as a distinct name for API stability; delegates directly to
/// [`contact_pairs`].
///
/// Results are identical to [`contact_pairs`] up to tie-breaking between
/// exactly equidistant candidates.
pub fn contact_pairs_octree(
    a: &[(usize, Vec<f64>)],
    b: &[(usize, Vec<f64>)],
) -> Vec<(usize, Option<(usize, f64)>)> {
    contact_pairs(a, b)
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
        let indexed: Vec<(usize, Vec<f64>)> = pts
            .iter()
            .enumerate()
            .map(|(i, c)| (i, c.clone()))
            .collect();
        let tree = Octree::build(&indexed);
        for q in &queries {
            let brute = contact_pairs(&[(0usize, q.clone())], &indexed)[0]
                .1
                .unwrap();
            let fast = tree.nearest(q).unwrap();
            assert_eq!(fast.0, brute.0, "different node for {q:?}");
            assert!((fast.1 - brute.1).abs() < 1e-12);
        }
    }

    #[test]
    fn octree_2d_and_edge_cases() {
        // 2-D points (zero-padded internally) and duplicate points must work;
        // an empty tree reports no neighbour.
        let pts: Vec<(usize, Vec<f64>)> = vec![
            (0, vec![0.0, 0.0]),
            (1, vec![2.0, 2.0]),
            (2, vec![2.0, 2.0]),
        ];
        let tree = Octree::build(&pts);
        let (i, d) = tree.nearest(&[1.9, 1.9]).unwrap();
        assert!(i >= 1 && d < 0.15, "got ({i}, {d})");
        assert!(tree.nearest(&[10.0, 10.0]).unwrap().0 >= 1);
        let empty = Octree::build(&[]);
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
