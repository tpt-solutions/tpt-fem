//! Dependency-free octree spatial index for nearest-neighbour queries.
//!
//! Generalised to 1-D/2-D/3-D by the dimension of the first point supplied to
//! [`Octree::build`]. Build once, query many times. Each stored point carries an
//! opaque `usize` payload (e.g. an index into the caller's point list) that is
//! returned by [`Octree::nearest`].
//!
//! This backs [`crate::contact_pairs`], replacing the former O(|a|·|b|)
//! brute-force scan with an amortised O(|a|·log|b|) lookup so that
//! surface-to-surface contact can scale to larger meshes.

/// Maximum number of points kept in a leaf before the node subdivides.
const BUCKET: usize = 16;
/// Maximum tree depth; guards against degenerate (e.g. identical) point clouds
/// that would otherwise subdivide forever.
const MAX_DEPTH: usize = 24;

/// An octree over a fixed set of points, supporting nearest-neighbour queries.
///
/// ```
/// use tpt_fem_contact::Octree;
/// let pts = vec![(0usize, vec![0.0, 0.0]), (1, vec![1.0, 0.0]), (2, vec![0.1, 0.0])];
/// let tree = Octree::build(&pts);
/// let (idx, d) = tree.nearest(&[0.12, 0.0]).unwrap();
/// assert_eq!(idx, 2);
/// assert!((d - 0.02).abs() < 1e-12);
/// ```
pub struct Octree {
    dim: usize,
    root: OctNode,
}

struct OctNode {
    lo: Vec<f64>,
    hi: Vec<f64>,
    /// `None` ⇒ this is a leaf whose `points` hold the data.
    /// `Some` ⇒ an internal node whose `children` hold the sub-trees.
    children: Option<Vec<OctNode>>,
    points: Vec<(usize, Vec<f64>)>,
}

impl Octree {
    /// Build an octree over `points`, where each entry is `(payload, coords)`.
    ///
    /// Returns an empty tree (whose [`Octree::nearest`] always yields `None`)
    /// when `points` is empty.
    pub fn build(points: &[(usize, Vec<f64>)]) -> Octree {
        if points.is_empty() {
            return Octree {
                dim: 0,
                root: OctNode {
                    lo: Vec::new(),
                    hi: Vec::new(),
                    children: Some(Vec::new()),
                    points: Vec::new(),
                },
            };
        }
        let dim = points[0].1.len();
        let (lo, hi) = bounding_box(points, dim);
        let root = build_node(points.to_vec(), lo, hi, dim, 0);
        Octree { dim, root }
    }

    /// Nearest point to `q`, returning `(payload, distance)`.
    ///
    /// Returns `None` if the tree is empty.
    pub fn nearest(&self, q: &[f64]) -> Option<(usize, f64)> {
        if self.dim == 0 {
            return None;
        }
        let mut best: Option<(usize, f64)> = None;
        nearest_in(&self.root, q, &mut best);
        best
    }
}

fn build_node(
    mut pts: Vec<(usize, Vec<f64>)>,
    lo: Vec<f64>,
    hi: Vec<f64>,
    dim: usize,
    depth: usize,
) -> OctNode {
    if pts.len() <= BUCKET || depth >= MAX_DEPTH {
        return OctNode {
            lo,
            hi,
            children: None,
            points: pts,
        };
    }
    let mid: Vec<f64> = lo.iter().zip(&hi).map(|(a, b)| 0.5 * (a + b)).collect();
    let n_children = 1usize << dim;
    let mut buckets: Vec<Vec<(usize, Vec<f64>)>> = vec![Vec::new(); n_children];
    for (payload, p) in pts.drain(..) {
        let mut idx = 0usize;
        for d in 0..dim {
            if p[d] >= mid[d] {
                idx |= 1 << d;
            }
        }
        buckets[idx].push((payload, p));
    }
    // Degenerate split (all points in one child, e.g. coincident points): stop
    // and store as a leaf to avoid an infinite subdivision loop.
    if buckets.iter().filter(|b| !b.is_empty()).count() <= 1 {
        return OctNode {
            lo,
            hi,
            children: None,
            points: buckets.into_iter().flatten().collect(),
        };
    }
    let children: Vec<OctNode> = buckets
        .into_iter()
        .enumerate()
        .map(|(ci, cpts)| {
            let mut clo = lo.clone();
            let mut chi = hi.clone();
            for d in 0..dim {
                if ci & (1 << d) != 0 {
                    clo[d] = mid[d];
                } else {
                    chi[d] = mid[d];
                }
            }
            build_node(cpts, clo, chi, dim, depth + 1)
        })
        .collect();
    OctNode {
        lo,
        hi,
        children: Some(children),
        points: Vec::new(),
    }
}

fn nearest_in(node: &OctNode, q: &[f64], best: &mut Option<(usize, f64)>) {
    match &node.children {
        None => {
            for (payload, p) in &node.points {
                let d = dist(q, p);
                if best.map_or(true, |(_, bd)| d < bd) {
                    *best = Some((*payload, d));
                }
            }
        }
        Some(children) => {
            let mut order: Vec<usize> = (0..children.len()).collect();
            order.sort_by(|&a, &b| {
                bbox_min_dist(q, &children[a])
                    .partial_cmp(&bbox_min_dist(q, &children[b]))
                    .unwrap()
            });
            for ci in order {
                let c = &children[ci];
                if c.points.is_empty() && c.children.is_none() {
                    continue;
                }
                let mind = bbox_min_dist(q, c);
                if best.is_some_and(|(_, bd)| mind > bd) {
                    continue;
                }
                nearest_in(c, q, best);
            }
        }
    }
}

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

fn bbox_min_dist(q: &[f64], node: &OctNode) -> f64 {
    let mut s = 0.0_f64;
    for d in 0..q.len() {
        let v = if q[d] < node.lo[d] {
            node.lo[d] - q[d]
        } else if q[d] > node.hi[d] {
            q[d] - node.hi[d]
        } else {
            0.0
        };
        s += v * v;
    }
    s.sqrt()
}

fn bounding_box(points: &[(usize, Vec<f64>)], dim: usize) -> (Vec<f64>, Vec<f64>) {
    let mut lo = vec![f64::INFINITY; dim];
    let mut hi = vec![f64::NEG_INFINITY; dim];
    for (_, p) in points {
        for d in 0..dim {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    (lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octree_matches_bruteforce() {
        // Deterministic LCG so the test is reproducible without `rand`.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut rng = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / (1u64 << 31) as f64 - 1.0
        };
        let mut b = Vec::new();
        for i in 0..500usize {
            b.push((i, vec![rng() * 10.0, rng() * 10.0, rng() * 10.0]));
        }
        let tree = Octree::build(&b);
        for _ in 0..200 {
            let q = vec![rng() * 10.0, rng() * 10.0, rng() * 10.0];
            let got = tree.nearest(&q);
            let mut exp: Option<(usize, f64)> = None;
            for (i, (_, p)) in b.iter().enumerate() {
                let d = dist(&q, p);
                if exp.map_or(true, |(_, bd)| d < bd) {
                    exp = Some((i, d));
                }
            }
            assert_eq!(
                got.unwrap().0,
                exp.unwrap().0,
                "nearest index mismatch for q={q:?}"
            );
            assert!((got.unwrap().1 - exp.unwrap().1).abs() < 1e-9);
        }
    }

    #[test]
    fn octree_empty_is_none() {
        let tree = Octree::build(&[]);
        assert!(tree.nearest(&[0.0, 0.0]).is_none());
    }
}
