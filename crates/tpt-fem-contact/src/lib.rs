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
}
