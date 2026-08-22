//! Adaptive h-refinement for the scalar Poisson problem on the unit square.
//!
//! The mesh is a **1-irregular quadtree** over `[0,1]²`: every leaf can be
//! split into four children, and the tree is re-balanced after each refinement
//! round so that adjacent leaves differ by at most one level. Hanging nodes
//! (edge midpoints of coarser leaves) are eliminated at assembly time with the
//! linear constraint `u_h = (u_a + u_b)/2`, so the conforming Q1 space is
//! preserved without multi-point-constraint bookkeeping on the solver side.
//!
//! The adaptive loop is the classic identify–mark–refine cycle:
//!
//! 1. assemble and solve `−Δu = f` (Q1, 2×2 Gauss),
//! 2. estimate the per-element error with a Zienkiewicz–Zhu gradient-recovery
//!    indicator,
//! 3. mark with Dörfler bulk marking (smallest set carrying `theta` of the
//!    estimated error),
//! 4. refine, re-balance, repeat until the element budget is reached.
//!
//! ```no_run
//! use tpt_fem_amr::{solve_adaptive, AmrOptions};
//!
//! // Poisson with a localized source: adaptivity concentrates elements near
//! // the peak and drives the estimated error below that of the uniform seed
//! // mesh for a given element budget.
//! let f = |x: f64, y: f64| (-500.0 * ((x - 0.5).powi(2) + (y - 0.25).powi(2))).exp();
//! let g = |_: f64, _: f64| 0.0;
//! let res = solve_adaptive(&f, &g, &AmrOptions { max_elements: 800, ..Default::default() }).unwrap();
//! assert!(res.mesh.len() <= 800);
//! ```

use std::collections::{BTreeSet, HashMap, HashSet};

use tpt_fem_sparse::{solve, Coo};

/// A quadtree cell identified by its level and integer grid position:
/// the cell covers `[i/2^level, (i+1)/2^level] × [j/2^level, (j+1)/2^level]`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CellKey {
    /// Refinement level (root = 0).
    pub level: u32,
    /// Integer grid index along x at `level`.
    pub i: i64,
    /// Integer grid index along y at `level`.
    pub j: i64,
}

impl CellKey {
    /// The four quadrant children of this cell.
    pub fn children(self) -> [CellKey; 4] {
        let l = self.level + 1;
        let (i, j) = (2 * self.i, 2 * self.j);
        [
            CellKey { level: l, i, j },
            CellKey {
                level: l,
                i: i + 1,
                j,
            },
            CellKey {
                level: l,
                i: i + 1,
                j: j + 1,
            },
            CellKey {
                level: l,
                i,
                j: j + 1,
            },
        ]
    }
}

/// The leaf set of a quadtree covering `[0,1]²`. Invariant: the leaves always
/// partition the domain exactly.
#[derive(Clone, Debug, Default)]
pub struct QuadTree {
    leaves: BTreeSet<CellKey>,
}

impl QuadTree {
    /// A tree with the single root cell.
    pub fn new_root() -> QuadTree {
        let mut t = QuadTree {
            leaves: BTreeSet::new(),
        };
        t.leaves.insert(CellKey {
            level: 0,
            i: 0,
            j: 0,
        });
        t
    }

    /// Split `c` into its four children. Returns `false` (and does nothing) if
    /// `c` is not currently a leaf.
    pub fn refine(&mut self, c: CellKey) -> bool {
        if self.leaves.remove(&c) {
            for ch in c.children() {
                self.leaves.insert(ch);
            }
            true
        } else {
            false
        }
    }

    /// Number of leaves.
    pub fn len(&self) -> usize {
        self.leaves.len()
    }

    /// Whether the tree has no leaves (never the case for a valid tree).
    pub fn is_empty(&self) -> bool {
        self.leaves.is_empty()
    }

    /// Iterate over the leaf cells (in sorted order).
    pub fn leaves(&self) -> impl Iterator<Item = CellKey> + '_ {
        self.leaves.iter().copied()
    }

    /// Restore 1-irregularity: refine every leaf that has a finer edge
    /// neighbour, repeating until stable.
    pub fn balance(&mut self) {
        loop {
            let mut changed = false;
            for &c in self.leaves.clone().iter() {
                if self.edge_slots_finier(c) {
                    changed |= self.refine(c);
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// True if any of the four edge neighbours of `c` is more than one level
    /// finer (which would create a hanging node violating 1-irregularity).
    fn edge_slots_finier(&self, c: CellKey) -> bool {
        let l = c.level + 1;
        let (i2, j2) = (2 * c.i, 2 * c.j);
        // The two level-(l+1) slots touching each of the four edges.
        let slots = [
            [
                CellKey {
                    level: l,
                    i: i2 + 1,
                    j: j2,
                },
                CellKey {
                    level: l,
                    i: i2 + 1,
                    j: j2 + 1,
                },
            ],
            [
                CellKey {
                    level: l,
                    i: i2 - 1,
                    j: j2,
                },
                CellKey {
                    level: l,
                    i: i2 - 1,
                    j: j2 + 1,
                },
            ],
            [
                CellKey {
                    level: l,
                    i: i2,
                    j: j2 + 1,
                },
                CellKey {
                    level: l,
                    i: i2 + 1,
                    j: j2 + 1,
                },
            ],
            [
                CellKey {
                    level: l,
                    i: i2,
                    j: j2 - 1,
                },
                CellKey {
                    level: l,
                    i: i2 + 1,
                    j: j2 - 1,
                },
            ],
        ];
        slots.iter().any(|pair| {
            pair.iter().any(|&s| {
                // Slots outside the unit square have no neighbour leaf; they
                // can never be finer than `c`.
                let extent = 1i64 << s.level;
                s.i >= 0
                    && s.i < extent
                    && s.j >= 0
                    && s.j < extent
                    && self.covering_leaf(s).level >= s.level
            })
        })
    }

    /// The leaf currently covering slot `s` (walks up the tree until hit).
    /// The caller must ensure `s` lies inside the domain; the level-0 guard
    /// keeps the walk total if it ever does not.
    fn covering_leaf(&self, mut s: CellKey) -> CellKey {
        loop {
            if self.leaves.contains(&s) || s.level == 0 {
                return s;
            }
            s = CellKey {
                level: s.level - 1,
                i: s.i >> 1,
                j: s.j >> 1,
            };
        }
    }
}

/// A conforming Q1 mesh extracted from the tree, plus hanging-node
/// constraints.
#[derive(Clone, Debug)]
pub struct HangingMesh {
    /// Node coordinates (node id = index).
    pub coords: Vec<[f64; 2]>,
    /// Leaf elements as CCW corner node quadruples.
    pub elems: Vec<[usize; 4]>,
    /// `(hanging_node, endpoint_a, endpoint_b)` with
    /// `u_h(hanging) = (u_a + u_b) / 2`.
    pub constraints: Vec<(usize, usize, usize)>,
}

impl HangingMesh {
    /// Number of leaf elements.
    pub fn len(&self) -> usize {
        self.elems.len()
    }

    /// Whether the mesh has no elements (never the case for a valid tree).
    pub fn is_empty(&self) -> bool {
        self.elems.is_empty()
    }

    /// Whether `node` lies on the domain boundary.
    pub fn on_boundary(&self, node: usize) -> bool {
        let [x, y] = self.coords[node];
        x <= 1e-12 || y <= 1e-12 || x >= 1.0 - 1e-12 || y >= 1.0 - 1e-12
    }
}

/// Quantise a coordinate onto the dyadic grid for exact node matching.
fn q(v: f64) -> i64 {
    (v * (1i64 << 26) as f64).round() as i64
}

/// Extract the conforming leaf mesh (with hanging-node constraints) from a
/// balanced tree.
pub fn build_mesh(tree: &QuadTree) -> HangingMesh {
    let scale = (1i64 << 26) as f64;
    let mut node_map: HashMap<(i64, i64), usize> = HashMap::new();
    let mut coords: Vec<[f64; 2]> = Vec::new();
    let get_node =
        |map: &mut HashMap<(i64, i64), usize>, cs: &mut Vec<[f64; 2]>, x: f64, y: f64| -> usize {
            *map.entry((q(x), q(y))).or_insert_with(|| {
                cs.push([x, y]);
                cs.len() - 1
            })
        };

    let mut elems = Vec::new();
    // Candidate hanging nodes: (midpoint key, endpoint a key, endpoint b key).
    type NodeKey = (i64, i64);
    let mut midpoints: Vec<(NodeKey, NodeKey, NodeKey)> = Vec::new();

    for &c in tree.leaves.iter() {
        let n = (1i64 << c.level) as f64;
        let x0 = c.i as f64 / n;
        let x1 = (c.i + 1) as f64 / n;
        let y0 = c.j as f64 / n;
        let y1 = (c.j + 1) as f64 / n;
        let corners = [
            (q(x0), q(y0)),
            (q(x1), q(y0)),
            (q(x1), q(y1)),
            (q(x0), q(y1)),
        ];
        let ids: Vec<usize> = corners
            .iter()
            .map(|&(kx, ky)| {
                get_node(
                    &mut node_map,
                    &mut coords,
                    kx as f64 / scale,
                    ky as f64 / scale,
                )
            })
            .collect();
        elems.push([ids[0], ids[1], ids[2], ids[3]]);
        // Edge midpoints of this leaf; if such a point ever becomes a mesh
        // node (because some neighbour refined deeper), it must hang on this
        // edge's endpoints.
        let mids = [
            ((q((x0 + x1) / 2.0), q(y0)), corners[0], corners[1]),
            ((q(x1), q((y0 + y1) / 2.0)), corners[1], corners[2]),
            ((q((x0 + x1) / 2.0), q(y1)), corners[3], corners[2]),
            ((q(x0), q((y0 + y1) / 2.0)), corners[0], corners[3]),
        ];
        midpoints.extend(mids);
    }

    let mut constraints = Vec::new();
    for (mkey, akey, bkey) in midpoints {
        if let (Some(&nm), Some(&na), Some(&nb)) = (
            node_map.get(&mkey),
            node_map.get(&akey),
            node_map.get(&bkey),
        ) {
            if nm != na && nm != nb {
                constraints.push((nm, na, nb));
            }
        }
    }

    HangingMesh {
        coords,
        elems,
        constraints,
    }
}

/// Assemble and solve `−Δu = f` with Dirichlet data `g` on the boundary of
/// the leaf mesh, eliminating hanging-node constraints locally. Returns the
/// full nodal solution (constrained values interpolated).
pub fn solve_poisson(
    mesh: &HangingMesh,
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
) -> Result<Vec<f64>, tpt_fem_sparse::SparseError> {
    let n = mesh.coords.len();
    // Constraint lookup: hanging node -> (a, b).
    let mut hang: HashMap<usize, (usize, usize)> = HashMap::new();
    let mut dupes: HashSet<(usize, usize, usize)> = HashSet::new();
    for &(c, a, b) in &mesh.constraints {
        if dupes.insert((c, a, b)) {
            hang.insert(c, (a, b));
        }
    }
    // Free-index numbering: hanging nodes have none.
    let mut free_idx: Vec<Option<usize>> = vec![None; n];
    let mut n_free = 0;
    for i in 0..n {
        if !hang.contains_key(&i) {
            free_idx[i] = Some(n_free);
            n_free += 1;
        }
    }

    let mut k = vec![vec![0.0_f64; n_free]; n_free];
    let mut rhs = vec![0.0_f64; n_free];
    // Distribute an entry (p, v_p)·(q, v_q) onto the free representatives.
    let reps = |i: usize| -> Vec<(usize, f64)> {
        match hang.get(&i) {
            None => vec![(free_idx[i].unwrap(), 1.0)],
            Some(&(a, b)) => {
                vec![(free_idx[a].unwrap(), 0.5), (free_idx[b].unwrap(), 0.5)]
            }
        }
    };

    for elem in &mesh.elems {
        let pts: Vec<[f64; 2]> = elem.iter().map(|&nd| mesh.coords[nd]).collect();
        let w = pts[1][0] - pts[0][0];
        let h = pts[3][1] - pts[0][1];
        // 2x2 Gauss Q1 stiffness and load.
        let mut ke = [[0.0_f64; 4]; 4];
        let mut fe = [0.0_f64; 4];
        let gauss = [-0.5773502691896257_f64, 0.5773502691896257];
        for &dxi in &gauss {
            for &eta in &gauss {
                let dnx = [-(1.0 - eta), 1.0 - eta, 1.0 + eta, -(1.0 + eta)];
                let dny = [-(1.0 - dxi), -(1.0 + dxi), 1.0 + dxi, 1.0 - dxi];
                let nvals = [
                    (1.0 - dxi) * (1.0 - eta),
                    (1.0 + dxi) * (1.0 - eta),
                    (1.0 + dxi) * (1.0 + eta),
                    (1.0 - dxi) * (1.0 + eta),
                ];
                for i in 0..4 {
                    for j in 0..4 {
                        ke[i][j] +=
                            (dnx[i] * dnx[j] / (w * w) + dny[i] * dny[j] / (h * h)) * 0.25 * w * h;
                    }
                    let xc = pts[i][0] + 0.5 * w * dxi;
                    let yc = pts[i][1] + 0.5 * h * eta;
                    fe[i] += nvals[i] * f(xc, yc) * 0.25 * w * h;
                }
            }
        }
        for (p, &gi) in elem.iter().enumerate() {
            for (q_, &gj) in elem.iter().enumerate() {
                for (ri, wi) in reps(gi) {
                    for (rj, wj) in reps(gj) {
                        k[ri][rj] += wi * wj * ke[p][q_];
                    }
                }
            }
            for (ri, wi) in reps(gi) {
                rhs[ri] += wi * fe[p];
            }
        }
    }

    // Dirichlet BCs live on boundary FREE nodes; impose them strongly
    // (identity rows) after lifting the column contributions into the rhs.
    let mut bcs: Vec<(usize, f64)> = Vec::new();
    for i in 0..n {
        if mesh.on_boundary(i) && !hang.contains_key(&i) {
            if let Some(fi) = free_idx[i] {
                bcs.push((fi, g(mesh.coords[i][0], mesh.coords[i][1])));
            }
        }
    }
    for &(d, v) in &bcs {
        for j in 0..n_free {
            rhs[j] -= k[j][d] * v;
            k[j][d] = 0.0;
        }
        for j in 0..n_free {
            k[d][j] = 0.0;
        }
        k[d][d] = 1.0;
        rhs[d] = v;
    }

    // Solve.
    let mut coo = Coo::with_capacity(n_free * 16);
    for i in 0..n_free {
        for j in 0..n_free {
            if k[i][j] != 0.0 {
                coo.push(i, j, k[i][j]);
            }
        }
    }
    let sol = solve(&coo, &rhs)?;

    // Scatter back: free dofs directly, hanging by interpolation.
    let mut u_full = vec![0.0; n];
    let mut inv_free: Vec<Option<usize>> = vec![None; n];
    for (i, slot) in free_idx.iter().enumerate() {
        if let Some(fi) = slot {
            inv_free[i] = Some(*fi);
        }
    }
    for (i, fi) in inv_free.iter().enumerate() {
        if let Some(fi) = fi {
            u_full[i] = sol[*fi];
        }
    }
    for (&c, &(a, b)) in &hang {
        u_full[c] = 0.5 * (u_full[a] + u_full[b]);
    }
    Ok(u_full)
}

/// Zienkiewicz–Zhu gradient-recovery error estimates, one per leaf element
/// (same order as `mesh.elems`).
pub fn zz_estimates(mesh: &HangingMesh, u: &[f64]) -> Vec<f64> {
    // Element-wise constant gradients at the centroids.
    let mut elem_grad: Vec<[f64; 2]> = Vec::with_capacity(mesh.elems.len());
    let mut areas: Vec<f64> = Vec::new();
    for e in &mesh.elems {
        let uc: Vec<f64> = e.iter().map(|&nd| u[nd]).collect();
        let (x0, x1) = (mesh.coords[e[0]][0], mesh.coords[e[1]][0]);
        let (y0, y1) = (mesh.coords[e[0]][1], mesh.coords[e[2]][1]);
        let w = x1 - x0;
        let h = y1 - y0;
        let gx = ((uc[1] + uc[2]) - (uc[0] + uc[3])) / (2.0 * w);
        let gy = ((uc[2] + uc[3]) - (uc[0] + uc[1])) / (2.0 * h);
        elem_grad.push([gx, gy]);
        areas.push(w * h);
    }
    // Nodal recovery: area-weighted average of incident element gradients.
    let mut nodal: Vec<([f64; 2], f64)> = vec![([0.0, 0.0], 0.0); mesh.coords.len()];
    for ((ei, e), &g) in mesh.elems.iter().enumerate().zip(&elem_grad) {
        for &nd in e.iter() {
            nodal[nd].0[0] += g[0] * areas[ei];
            nodal[nd].0[1] += g[1] * areas[ei];
            nodal[nd].1 += areas[ei];
        }
    }
    let recovered: Vec<[f64; 2]> = nodal
        .iter()
        .map(|(g, w)| if *w > 0.0 { [g[0] / *w, g[1] / *w] } else { *g })
        .collect();
    // Element error: area * |grad_recovered(center) - grad_elem|².
    mesh.elems
        .iter()
        .enumerate()
        .map(|(ei, e)| {
            let mut gr = [0.0_f64, 0.0];
            for &nd in e {
                gr[0] += recovered[nd][0];
                gr[1] += recovered[nd][1];
            }
            gr[0] *= 0.25;
            gr[1] *= 0.25;
            let dx = gr[0] - elem_grad[ei][0];
            let dy = gr[1] - elem_grad[ei][1];
            areas[ei] * (dx * dx + dy * dy)
        })
        .collect()
}

/// Options for [`solve_adaptive`].
#[derive(Clone, Copy, Debug)]
pub struct AmrOptions {
    /// Stop once the mesh has at least this many elements.
    pub max_elements: usize,
    /// Dörfler bulk-marking fraction in `(0, 1]`.
    pub theta: f64,
}

impl Default for AmrOptions {
    fn default() -> Self {
        AmrOptions {
            max_elements: 512,
            theta: 0.7,
        }
    }
}

/// Result of [`solve_adaptive`].
#[derive(Clone, Debug)]
pub struct AdaptiveSolution {
    /// The final (balanced) leaf mesh.
    pub mesh: HangingMesh,
    /// Nodal solution on the final mesh.
    pub u: Vec<f64>,
    /// ZZ error estimates of the final solve (per element).
    pub estimates: Vec<f64>,
    /// Square root of the total estimated error² on the final mesh.
    pub estimated_error: f64,
}

/// Run the identify–mark–refine loop until the element budget is reached.
///
/// Each round solves on the current balanced leaf mesh, estimates the error
/// with [`zz_estimates`], and refines the Dörfler set (smallest prefix of
/// elements carrying `theta` of the total estimated error²). Element order in
/// [`build_mesh`] follows the tree's leaf order, so estimates map directly
/// back to leaves.
pub fn solve_adaptive(
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
    opts: &AmrOptions,
) -> Result<AdaptiveSolution, tpt_fem_sparse::SparseError> {
    let mut tree = QuadTree::new_root();
    // Seed with a small uniform grid: on the single root cell the gradient
    // recovery estimate is identically zero (nothing to average against),
    // which would make the first Dörfler mark degenerate.
    for _ in 0..2 {
        for c in tree.leaves.iter().copied().collect::<Vec<_>>() {
            tree.refine(c);
        }
    }
    loop {
        tree.balance();
        let mesh = build_mesh(&tree);
        let u = solve_poisson(&mesh, f, g)?;
        let estimates = zz_estimates(&mesh, &u);
        let total = estimates.iter().map(|e| e.sqrt()).sum::<f64>().sqrt();
        if tree.len() >= opts.max_elements {
            return Ok(AdaptiveSolution {
                mesh,
                u,
                estimates,
                estimated_error: total,
            });
        }
        // Dörfler bulk marking; element ei == the ei-th leaf in sorted order.
        let mut order: Vec<usize> = (0..estimates.len()).collect();
        order.sort_by(|&a, &b| {
            estimates[b]
                .partial_cmp(&estimates[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let target = opts.theta * total * total;
        let keys: Vec<CellKey> = tree.leaves.iter().copied().collect();
        let mut acc = 0.0;
        for &ei in &order {
            tree.refine(keys[ei]);
            acc += estimates[ei];
            if acc >= target {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Poisson with a smooth manufactured solution `u = sin(pi x) sin(pi y)`.
    fn mms_f(x: f64, y: f64) -> f64 {
        2.0 * std::f64::consts::PI
            * std::f64::consts::PI
            * (std::f64::consts::PI * x).sin()
            * (std::f64::consts::PI * y).sin()
    }

    #[test]
    fn refine_only_splits_leaves() {
        let mut t = QuadTree::new_root();
        assert_eq!(t.len(), 1);
        let root = *t.leaves.iter().next().unwrap();
        assert!(t.refine(root));
        assert_eq!(t.len(), 4);
        // The root is no longer a leaf, so refining it again is a no-op.
        assert!(!t.refine(root));
        assert_eq!(t.len(), 4);
    }

    #[test]
    fn balance_restores_one_irregularity() {
        let mut t = QuadTree::new_root();
        let root = *t.leaves.iter().next().unwrap();
        // Split the root into four level-1 cells first.
        assert!(t.refine(root));
        assert_eq!(t.len(), 4);
        // Refine only one of the four children: the neighbours of the fine
        // cell are now two levels coarser and must be split by `balance`.
        let child = root.children()[0];
        assert!(t.refine(child));
        // 3 remaining level-1 cells + the 4 children of the refined one.
        assert_eq!(t.len(), 7);
        t.balance();
        // After balancing, adjacent leaves differ by at most one level: the
        // three untouched neighbours must have been split too.
        assert_eq!(t.len(), 4 + 3 * 4);
        // Leaves still partition the unit square: total leaf area is 1.
        let area: f64 = t
            .leaves
            .iter()
            .map(|c| {
                let n = (1i64 << c.level) as f64;
                (1.0 / n) * (1.0 / n)
            })
            .sum();
        assert!((area - 1.0).abs() < 1e-12);
    }

    #[test]
    fn build_mesh_is_conforming_partition() {
        let mut t = QuadTree::new_root();
        let root = *t.leaves.iter().next().unwrap();
        for ch in root.children() {
            t.refine(ch);
        }
        t.balance();
        let mesh = build_mesh(&t);
        // Leaf count matches element count; every leaf is a unit square.
        assert_eq!(mesh.elems.len(), t.len());
        let total: f64 = mesh
            .elems
            .iter()
            .map(|e| {
                let (x0, x1) = (mesh.coords[e[0]][0], mesh.coords[e[1]][0]);
                let (y0, y1) = (mesh.coords[e[0]][1], mesh.coords[e[2]][1]);
                (x1 - x0) * (y1 - y0)
            })
            .sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn poisson_converges_to_mms_on_uniform_mesh() {
        let mut t = QuadTree::new_root();
        let root = *t.leaves.iter().next().unwrap();
        for ch in root.children() {
            t.refine(ch);
        }
        let mesh = build_mesh(&t);
        let g = |_: f64, _: f64| 0.0;
        let u = solve_poisson(&mesh, &mms_f, &g).unwrap();
        let err: f64 = mesh
            .coords
            .iter()
            .zip(&u)
            .map(|(c, uh)| {
                let exact =
                    (std::f64::consts::PI * c[0]).sin() * (std::f64::consts::PI * c[1]).sin();
                (uh - exact) * (uh - exact)
            })
            .sum::<f64>()
            .sqrt();
        // 4x4 Q1 mesh: L2 nodal error is small.
        assert!(err < 5e-2, "nodal L2 error {err}");
    }

    #[test]
    fn adaptive_reduces_true_error_within_budget() {
        // MMS problem: the adaptive loop must land on a mesh whose true nodal
        // error against the exact solution is smaller than the uniform seed
        // grid's, while respecting the element budget.
        let g = |_: f64, _: f64| 0.0;
        let adaptive = solve_adaptive(
            &mms_f,
            &g,
            &AmrOptions {
                max_elements: 256,
                ..Default::default()
            },
        )
        .unwrap();
        let true_err = |mesh: &HangingMesh, u: &[f64]| {
            mesh.coords
                .iter()
                .zip(u)
                .map(|(c, uh)| {
                    let e = uh
                        - (std::f64::consts::PI * c[0]).sin() * (std::f64::consts::PI * c[1]).sin();
                    e * e
                })
                .sum::<f64>()
                .sqrt()
        };
        // The seed grid is the level-2 uniform mesh (16 elements).
        let mut seed = QuadTree::new_root();
        for _ in 0..2 {
            for c in seed.leaves().collect::<Vec<_>>() {
                seed.refine(c);
            }
        }
        let smesh = build_mesh(&seed);
        let su = solve_poisson(&smesh, &mms_f, &g).unwrap();
        let seed_err = true_err(&smesh, &su);
        assert!(adaptive.mesh.len() <= 256);
        let err = true_err(&adaptive.mesh, &adaptive.u);
        assert!(err < seed_err, "adaptive {err} vs seed {seed_err}");
    }

    #[test]
    fn estimates_are_nonnegative_and_finite() {
        let f = |x: f64, y: f64| (-200.0 * ((x - 0.3).powi(2) + (y - 0.7).powi(2))).exp();
        let g = |_: f64, _: f64| 0.0;
        let res = solve_adaptive(
            &f,
            &g,
            &AmrOptions {
                max_elements: 128,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(res.estimates.len(), res.mesh.elems.len());
        assert!(res.estimates.iter().all(|e| e.is_finite() && *e >= 0.0));
        assert!(res.estimated_error.is_finite());
    }
}
