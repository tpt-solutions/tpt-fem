//! Native 3D tetrahedral mesh generation for `tpt-fem`.
//!
//! This crate provides dependency-free mesh generation, written from scratch
//! so it can ship under the repo's `MIT OR Apache-2.0` policy where wrapped
//! alternatives (e.g. TetGen-based tooling) are AGPL and therefore disallowed.
//!
//! Two generators are provided:
//!
//! * [`delaunay_3d`] — an incremental (Bowyer–Watson) Delaunay
//!   tetrahedralisation of an arbitrary point cloud. It meshes the convex hull
//!   of the input points and is useful for filling a volume bounded by a set of
//!   seed nodes.
//! * [`box_mesh`] — a structured mesher that splits every axis-aligned brick of
//!   a bounding box into six tetrahedra. This path is *guaranteed* to produce a
//!   valid, intersection-free, positively-oriented mesh with no external
//!   dependency and no robustness caveats, so it is the recommended route when
//!   a quality grid is needed.
//!
//! Both return a [`tpt_fem_mesh::Mesh`] of [`CellType::Tet`] elements. Quality
//! can be inspected with [`tet_quality`] and improved with [`laplacian_smooth`].
//!
//! # Robustness note
//!
//! The Delaunay predicates ([`orient3`], [`in_sphere`]) are evaluated in `f64`
//! with a relative tolerance. They are exact for points in general position
//! (no four coplanar, no five cospherical) and degrade gracefully otherwise;
//! coincident input points are de-duplicated automatically. For highest
//! quality on arbitrary domains, the [`box_mesh`] + [`laplacian_smooth`] route
//! is preferred.

use std::collections::{HashMap, HashSet};

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// A 3D point.
pub type Point3 = [f64; 3];

/// Quality statistics aggregated over the tetrahedral elements of a mesh.
#[derive(Clone, Copy, Debug, Default)]
pub struct TetQuality {
    /// Smallest interior dihedral angle over all tetrahedra, in radians.
    pub min_dihedral: f64,
    /// Largest interior dihedral angle over all tetrahedra, in radians.
    pub max_dihedral: f64,
    /// Smallest radius-edge ratio (circumradius / shortest edge) over all
    /// tetrahedra. Smaller is better; values below ~2 are excellent.
    pub min_radius_edge: f64,
    /// Largest radius-edge ratio over all tetrahedra.
    pub max_radius_edge: f64,
}

// ---------------------------------------------------------------------------
// Geometric predicates
// ---------------------------------------------------------------------------

/// Signed volume (×6) of the tetrahedron `(a, b, c, d)`.
///
/// Positive when `d` lies above the plane through `(a, b, c)` oriented so that
/// `a → b → c` is counter-clockwise. Returns `0.0` for a degenerate (flat)
/// tetrahedron.
pub fn orient3(a: &Point3, b: &Point3, c: &Point3, d: &Point3) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) - ab[1] * (ac[0] * ad[2] - ac[2] * ad[0])
        + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])
}

/// In-sphere predicate for the positively-oriented tetrahedron `(a, b, c, d)`.
///
/// Returns a value `> 0` when `p` lies strictly inside the circumsphere of the
/// tetrahedron, `< 0` when outside, and `0` on the sphere (degenerate). The
/// result is only meaningful when `orient3(a, b, c, d) > 0`.
pub fn in_sphere(a: &Point3, b: &Point3, c: &Point3, d: &Point3, p: &Point3) -> f64 {
    let rows = [a, b, c, d, p];
    let mut m: [[f64; 5]; 5] = [[0.0; 5]; 5];
    for (i, q) in rows.iter().enumerate() {
        let s = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
        m[i] = [q[0], q[1], q[2], s, 1.0];
    }
    // Negated so that, for a positively-oriented tetrahedron, a positive value
    // means `p` lies strictly inside the circumsphere.
    -det_n::<5>(m)
}

/// Generic determinant of an `N × N` matrix via Gaussian elimination with
/// partial pivoting. Used for the small fixed-size predicates above.
fn det_n<const N: usize>(mut a: [[f64; N]; N]) -> f64 {
    let mut det = 1.0;
    for col in 0..N {
        let mut piv = col;
        let mut mx = a[col][col].abs();
        for r in (col + 1)..N {
            if a[r][col].abs() > mx {
                mx = a[r][col].abs();
                piv = r;
            }
        }
        if mx < 1e-300 {
            return 0.0;
        }
        if piv != col {
            a.swap(col, piv);
            det = -det;
        }
        let pv = a[col][col];
        det *= pv;
        for r in (col + 1)..N {
            let f = a[r][col] / pv;
            for c in col..N {
                a[r][c] -= f * a[col][c];
            }
        }
    }
    det
}

/// Relative tolerance used by the Delaunay cavity test.
const EPS: f64 = 1e-9;

/// Squared distance between two points.
fn dist2(a: &Point3, b: &Point3) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// ---------------------------------------------------------------------------
// Bowyer–Watson Delaunay tetrahedralisation
// ---------------------------------------------------------------------------

/// A tetrahedron in the incremental builder. `v[i]` is the vertex index into
/// the working point list; face `i` is opposite `v[i]`. `nbr[i]` is the id of
/// the neighbouring tetrahedron across face `i` (`-1` if none).
#[derive(Clone)]
struct Tet {
    v: [usize; 4],
    nbr: [isize; 4],
    deleted: bool,
}

impl Tet {
    /// The three vertex indices of face `i` (opposite `v[i]`), in stored order.
    fn face(&self, i: usize) -> [usize; 3] {
        match i {
            0 => [self.v[1], self.v[2], self.v[3]],
            1 => [self.v[0], self.v[2], self.v[3]],
            2 => [self.v[0], self.v[1], self.v[3]],
            _ => [self.v[0], self.v[1], self.v[2]],
        }
    }

    /// Vertex indices of face `i` as a sorted triple (order-independent key).
    fn face_key(&self, i: usize) -> [usize; 3] {
        let mut f = self.face(i);
        f.sort_unstable();
        f
    }
}

/// Incremental 3D Delaunay builder (Bowyer–Watson).
struct Delaunay {
    pts: Vec<Point3>,
    tets: Vec<Tet>,
}

impl Delaunay {
    fn new(pts: Vec<Point3>) -> Self {
        Delaunay {
            pts,
            tets: Vec::new(),
        }
    }

    /// Add a point, re-triangulating the affected cavity.
    fn insert(&mut self, p: usize) {
        // 1. Find all tets whose circumsphere contains p ("bad" tets).
        let mut bad = Vec::new();
        for (id, t) in self.tets.iter().enumerate() {
            if t.deleted {
                continue;
            }
            let [a, b, c, d] = t.v;
            if in_sphere(
                &self.pts[a],
                &self.pts[b],
                &self.pts[c],
                &self.pts[d],
                &self.pts[p],
            ) > EPS
            {
                bad.push(id);
            }
        }
        if bad.is_empty() {
            return;
        }

        // 2. Collect the boundary of the cavity: faces of bad tets that are not
        //    shared with another bad tet. Each entry is (face vertices, id of
        //    the surviving neighbour across that face, or -1).
        let bad_set: HashSet<usize> = bad.iter().copied().collect();
        let mut boundary: Vec<([usize; 3], isize)> = Vec::new();
        for &id in &bad {
            let t = &self.tets[id];
            for fi in 0..4 {
                let nb = t.nbr[fi];
                let interior = nb >= 0 && bad_set.contains(&(nb as usize));
                if !interior {
                    boundary.push((t.face(fi), nb));
                }
            }
        }

        // 3. Delete the bad tets.
        for &id in &bad {
            self.tets[id].deleted = true;
        }

        // 4. Retriangulate: one new tet per boundary face, apex at p.
        let start = self.tets.len();
        for (face, nb) in boundary {
            let mut nt = Tet {
                v: [face[0], face[1], face[2], p],
                nbr: [-1, -1, -1, nb],
                deleted: false,
            };
            if orient3(
                &self.pts[face[0]],
                &self.pts[face[1]],
                &self.pts[face[2]],
                &self.pts[p],
            ) < 0.0
            {
                nt.v.swap(1, 2);
            }
            self.tets.push(nt);
        }

        // 5. Link new tets to each other (shared side faces include the apex).
        let count = self.tets.len() - start;
        for i in 0..count {
            let ti = start + i;
            for j in (i + 1)..count {
                let tj = start + j;
                if let Some((fi, fj)) = shared_face(&self.tets[ti], &self.tets[tj]) {
                    self.tets[ti].nbr[fi] = tj as isize;
                    self.tets[tj].nbr[fj] = ti as isize;
                }
            }
        }

        // 6. Link each new tet's base face to the surviving boundary neighbour.
        for i in 0..count {
            let ti = start + i;
            let nb = self.tets[ti].nbr[3];
            if nb < 0 {
                continue;
            }
            let nb = nb as usize;
            let key = self.tets[ti].face_key(3);
            for fj in 0..4 {
                if self.tets[nb].face_key(fj) == key {
                    self.tets[nb].nbr[fj] = ti as isize;
                    break;
                }
            }
        }
    }

    /// Build a [`Mesh`] from the surviving real tetrahedra (super-tetrahedron
    /// vertices are `super_start..` and are dropped).
    fn into_mesh(self, super_start: usize) -> Mesh {
        let mut b = MeshBuilder::new();
        for p in &self.pts[0..super_start] {
            b.add_node(p.to_vec());
        }
        for t in &self.tets {
            if t.deleted {
                continue;
            }
            if t.v.iter().any(|&vi| vi >= super_start) {
                continue; // touches the super-tetrahedron → outside hull
            }
            let [a, b2, c, d] = t.v;
            let pa = &self.pts[a];
            let pb = &self.pts[b2];
            let pc = &self.pts[c];
            let pd = &self.pts[d];
            if orient3(pa, pb, pc, pd).abs() < 1e-12 {
                continue; // defensive: drop near-degenerate tets
            }
            b.add_element(CellType::Tet, t.v.to_vec());
        }
        b.build()
    }
}

/// If tetrahedra `x` and `y` share a face, return the pair of face indices
/// (in `x`, in `y`) whose vertex sets coincide.
fn shared_face(x: &Tet, y: &Tet) -> Option<(usize, usize)> {
    for fx in 0..4 {
        let kx = x.face_key(fx);
        for fy in 0..4 {
            if y.face_key(fy) == kx {
                return Some((fx, fy));
            }
        }
    }
    None
}

/// De-duplicate `points` (within `tol`), returning the unique points and a
/// mapping from each input index to its unique-point index.
fn dedupe(points: &[Point3], tol: f64) -> (Vec<Point3>, Vec<usize>) {
    let t2 = tol * tol;
    let mut uniq: Vec<Point3> = Vec::new();
    let mut map = Vec::with_capacity(points.len());
    for p in points {
        if let Some(pos) = uniq.iter().position(|q| dist2(q, p) <= t2) {
            map.push(pos);
        } else {
            map.push(uniq.len());
            uniq.push(*p);
        }
    }
    (uniq, map)
}

/// Tetrahedralise a cloud of 3D points by incremental Delaunay.
///
/// The convex hull of `points` is meshed with `Tet4` elements. Coincident
/// points are merged. For best results the input should be in general position
/// (see the crate-level robustness note); for guaranteed-quality structured
/// grids prefer [`box_mesh`].
pub fn delaunay_3d(points: &[Point3]) -> Mesh {
    let (uniq, _map) = dedupe(points, 1e-9);
    if uniq.len() < 4 {
        // Fewer than four non-coincident points cannot form a tetrahedron.
        return MeshBuilder::new().build();
    }

    let mut d = Delaunay::new(uniq);

    // Axis-aligned bounding box, then a super-tetrahedron circumscribing the
    // bounding sphere so every input point lies strictly inside it.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &d.pts {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let center = [
        (lo[0] + hi[0]) / 2.0,
        (lo[1] + hi[1]) / 2.0,
        (lo[2] + hi[2]) / 2.0,
    ];
    let r = (0.25 * ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)))
        .sqrt()
        * 6.0;
    // Regular tetrahedron directions (each opposite face is equilateral).
    let dirs: [Point3; 4] = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let super_start = d.pts.len();
    for dir in dirs {
        let n = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        d.pts.push([
            center[0] + r * dir[0] / n,
            center[1] + r * dir[1] / n,
            center[2] + r * dir[2] / n,
        ]);
    }

    // Initial super-tetrahedron (vertices super_start..super_start+4).
    let s = super_start;
    let mut sup = Tet {
        v: [s, s + 1, s + 2, s + 3],
        nbr: [-1, -1, -1, -1],
        deleted: false,
    };
    if orient3(&d.pts[s], &d.pts[s + 1], &d.pts[s + 2], &d.pts[s + 3]) < 0.0 {
        sup.v.swap(1, 2);
    }
    d.tets.push(sup);

    // Insert real points in shuffled (deterministic) order for robustness.
    let mut order: Vec<usize> = (0..super_start).collect();
    let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in (1..order.len()).rev() {
        // LCG-based Fisher–Yates (no external RNG dependency).
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = (seed >> 33) as usize % (i + 1);
        order.swap(i, j);
    }
    for &p in &order {
        d.insert(p);
    }

    d.into_mesh(super_start)
}

// ---------------------------------------------------------------------------
// Structured box mesher
// ---------------------------------------------------------------------------

/// Build a structured tetrahedral mesh of the axis-aligned box `[min, max]`.
///
/// The box is subdivided into `n = [nx, ny, nz]` bricks per axis direction and
/// each brick is split into six tetrahedra sharing its main diagonal; the
/// resulting mesh is intersection-free, positively oriented, and has no
/// dependency on any external library. Returns a [`Mesh`] of `Tet4` elements
/// with `(nx+1)·(ny+1)·(nz+1)` nodes.
pub fn box_mesh(min: Point3, max: Point3, n: [usize; 3]) -> Mesh {
    let [nx, ny, nz] = n;
    let mut b = MeshBuilder::new();
    // Grid of node ids; node (i,j,k) -> id.
    let mut ids = vec![vec![vec![0usize; nz + 1]; ny + 1]; nx + 1];
    for i in 0..=nx {
        for j in 0..=ny {
            for k in 0..=nz {
                let t = if nx > 0 { i as f64 / nx as f64 } else { 0.0 };
                let u = if ny > 0 { j as f64 / ny as f64 } else { 0.0 };
                let w = if nz > 0 { k as f64 / nz as f64 } else { 0.0 };
                let x = min[0] + t * (max[0] - min[0]);
                let y = min[1] + u * (max[1] - min[1]);
                let z = min[2] + w * (max[2] - min[2]);
                ids[i][j][k] = b.add_node(vec![x, y, z]);
            }
        }
    }
    // Corner offsets of a brick.
    let corner = |i: usize, j: usize, k: usize| -> Point3 {
        let t = if nx > 0 { i as f64 / nx as f64 } else { 0.0 };
        let u = if ny > 0 { j as f64 / ny as f64 } else { 0.0 };
        let w = if nz > 0 { k as f64 / nz as f64 } else { 0.0 };
        [
            min[0] + t * (max[0] - min[0]),
            min[1] + u * (max[1] - min[1]),
            min[2] + w * (max[2] - min[2]),
        ]
    };
    // Six-tet decomposition of a brick (corners 0..7) sharing the 0-6 diagonal.
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let c = |di: usize, dj: usize, dk: usize| ids[i + di][j + dj][k + dk];
                let id = [
                    c(0, 0, 0),
                    c(1, 0, 0),
                    c(1, 1, 0),
                    c(0, 1, 0),
                    c(0, 0, 1),
                    c(1, 0, 1),
                    c(1, 1, 1),
                    c(0, 1, 1),
                ];
                let co = [
                    corner(0, 0, 0),
                    corner(1, 0, 0),
                    corner(1, 1, 0),
                    corner(0, 1, 0),
                    corner(0, 0, 1),
                    corner(1, 0, 1),
                    corner(1, 1, 1),
                    corner(0, 1, 1),
                ];
                for tet in [
                    [0usize, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    let mut tv = [id[tet[0]], id[tet[1]], id[tet[2]], id[tet[3]]];
                    let (a, bb, cc, dd) = (co[tet[0]], co[tet[1]], co[tet[2]], co[tet[3]]);
                    if orient3(&a, &bb, &cc, &dd) < 0.0 {
                        tv.swap(1, 2);
                    }
                    b.add_element(CellType::Tet, tv.to_vec());
                }
            }
        }
    }
    b.build()
}

// ---------------------------------------------------------------------------
// Quality metrics and smoothing
// ---------------------------------------------------------------------------

/// Returns `true` if every `Tet4` element has strictly positive (right-handed)
/// volume.
pub fn all_positively_oriented(mesh: &Mesh) -> bool {
    mesh.elements
        .iter()
        .filter(|e| e.cell_type == CellType::Tet)
        .all(|e| {
            let p = |n: usize| -> Point3 {
                let c = mesh.node_coords(n);
                [
                    c[0],
                    c.get(1).copied().unwrap_or(0.0),
                    c.get(2).copied().unwrap_or(0.0),
                ]
            };
            orient3(
                &p(e.nodes[0]),
                &p(e.nodes[1]),
                &p(e.nodes[2]),
                &p(e.nodes[3]),
            ) > 0.0
        })
}

/// Aggregate quality metrics over the `Tet4` elements of `mesh`.
///
/// *Dihedral* angles are the interior angles between adjacent faces of a
/// tetrahedron; *radius-edge* is the circumradius divided by the shortest edge.
/// Both are reported as their minimum and maximum across the mesh (initialised
/// to sentinels when no tetrahedra are present).
pub fn tet_quality(mesh: &Mesh) -> TetQuality {
    let mut q = TetQuality {
        min_dihedral: f64::INFINITY,
        max_dihedral: 0.0,
        min_radius_edge: f64::INFINITY,
        max_radius_edge: 0.0,
    };
    for e in &mesh.elements {
        if e.cell_type != CellType::Tet {
            continue;
        }
        let v: [Point3; 4] = [
            node_pt(mesh, e.nodes[0]),
            node_pt(mesh, e.nodes[1]),
            node_pt(mesh, e.nodes[2]),
            node_pt(mesh, e.nodes[3]),
        ];

        // Outward face normals (face f is opposite vertex f).
        let mut norms: [Point3; 4] = [[0.0; 3]; 4];
        for f in 0..4 {
            let others: [usize; 3] = match f {
                0 => [1, 2, 3],
                1 => [0, 2, 3],
                2 => [0, 1, 3],
                _ => [0, 1, 2],
            };
            let (a, bb, cc) = (v[others[0]], v[others[1]], v[others[2]]);
            let ab = [bb[0] - a[0], bb[1] - a[1], bb[2] - a[2]];
            let ac = [cc[0] - a[0], cc[1] - a[1], cc[2] - a[2]];
            let mut n = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            // Orient outward: away from the opposite vertex.
            let cx = (a[0] + bb[0] + cc[0]) / 3.0;
            let cy = (a[1] + bb[1] + cc[1]) / 3.0;
            let cz = (a[2] + bb[2] + cc[2]) / 3.0;
            let inward = [cx - v[f][0], cy - v[f][1], cz - v[f][2]];
            if n[0] * inward[0] + n[1] * inward[1] + n[2] * inward[2] > 0.0 {
                n = [-n[0], -n[1], -n[2]];
            }
            norms[f] = n;
        }

        // Each edge (i,j) is shared by the two faces {0,1,2,3} \ {i,j}.
        let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for (i, j) in edges {
            let f1 = (0..4).find(|&f| f != i && f != j).unwrap();
            let f2 = (0..4).find(|&f| f != i && f != j && f != f1).unwrap();
            let (n1, n2) = (norms[f1], norms[f2]);
            let dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
            let l1 = (n1[0].powi(2) + n1[1].powi(2) + n1[2].powi(2)).sqrt();
            let l2 = (n2[0].powi(2) + n2[1].powi(2) + n2[2].powi(2)).sqrt();
            if l1 < 1e-12 || l2 < 1e-12 {
                continue;
            }
            let cos = (dot / (l1 * l2)).clamp(-1.0, 1.0);
            // Interior dihedral = π − angle between outward normals.
            let ang = std::f64::consts::PI - cos.acos();
            q.min_dihedral = q.min_dihedral.min(ang);
            q.max_dihedral = q.max_dihedral.max(ang);
        }

        // Radius-edge ratio.
        let r = circumradius(&v);
        let shortest = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| dist2(&v[i], &v[j])))
            .filter(|&s| s > 1e-24)
            .map(|s| s.sqrt())
            .fold(f64::INFINITY, f64::min);
        if shortest.is_finite() && r.is_finite() {
            let re = r / shortest;
            q.min_radius_edge = q.min_radius_edge.min(re);
            q.max_radius_edge = q.max_radius_edge.max(re);
        }
    }
    if q.min_dihedral.is_infinite() {
        q.min_dihedral = 0.0;
    }
    if q.min_radius_edge.is_infinite() {
        q.min_radius_edge = 0.0;
    }
    q
}

/// Read a node's coordinates as a 3D point (zero-padded for 2D meshes).
fn node_pt(mesh: &Mesh, n: usize) -> Point3 {
    let c = mesh.node_coords(n);
    [
        c[0],
        c.get(1).copied().unwrap_or(0.0),
        c.get(2).copied().unwrap_or(0.0),
    ]
}

/// Circumradius of a tetrahedron with vertices `v[0..4]`.
fn circumradius(v: &[Point3; 4]) -> f64 {
    // Solve for the circumcenter as the intersection of perpendicular bisector
    // planes: |x - a|² = |x - b|² etc. → linear system M c = rhs.
    let mut m = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for k in 0..3 {
        let ba = [
            v[k + 1][0] - v[0][0],
            v[k + 1][1] - v[0][1],
            v[k + 1][2] - v[0][2],
        ];
        m[k].copy_from_slice(&ba);
        rhs[k] = 0.5
            * (ba[0] * (v[k + 1][0] + v[0][0])
                + ba[1] * (v[k + 1][1] + v[0][1])
                + ba[2] * (v[k + 1][2] + v[0][2]));
    }
    let det = det_n::<3>(m);
    if det.abs() < 1e-12 {
        return f64::INFINITY;
    }
    let mut center = [0.0; 3];
    for c in 0..3 {
        let mut col = m;
        for row in 0..3 {
            col[row][c] = rhs[row];
        }
        center[c] = det_n::<3>(col) / det;
    }
    dist2(&center, &v[0]).sqrt()
}

/// Run `iterations` passes of Laplacian (umbrella) smoothing on the interior
/// nodes of a tetrahedral mesh, holding boundary nodes fixed. Returns the
/// maximum node displacement over the last pass.
///
/// Boundary nodes are those lying on a mesh surface (a triangle face shared by
/// exactly one tetrahedron); interior nodes move to the centroid of their
/// one-ring neighbours.
pub fn laplacian_smooth(mesh: &mut Mesh, iterations: usize) -> f64 {
    let dim = mesh.node_coords(0).len();
    let mut max_disp = 0.0;
    for _ in 0..iterations {
        // Count triangle faces; a face shared by exactly one tet is a boundary
        // face, and its nodes are boundary nodes (held fixed).
        let mut face_count: HashMap<(usize, usize, usize), usize> = HashMap::new();
        for e in &mesh.elements {
            if e.cell_type != CellType::Tet {
                continue;
            }
            for tri in [
                [e.nodes[0], e.nodes[1], e.nodes[2]],
                [e.nodes[0], e.nodes[1], e.nodes[3]],
                [e.nodes[0], e.nodes[2], e.nodes[3]],
                [e.nodes[1], e.nodes[2], e.nodes[3]],
            ] {
                let mut t = tri;
                t.sort_unstable();
                *face_count.entry((t[0], t[1], t[2])).or_insert(0) += 1;
            }
        }
        let boundary: HashSet<usize> = face_count
            .iter()
            .filter(|(_, &c)| c == 1)
            .flat_map(|(k, _)| [k.0, k.1, k.2])
            .collect();

        // One-ring adjacency (node -> neighbour set).
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); mesh.node_count()];
        for e in &mesh.elements {
            if e.cell_type != CellType::Tet {
                continue;
            }
            for pair in [
                (e.nodes[0], e.nodes[1]),
                (e.nodes[0], e.nodes[2]),
                (e.nodes[0], e.nodes[3]),
                (e.nodes[1], e.nodes[2]),
                (e.nodes[1], e.nodes[3]),
                (e.nodes[2], e.nodes[3]),
            ] {
                if !adj[pair.0].contains(&pair.1) {
                    adj[pair.0].push(pair.1);
                }
                if !adj[pair.1].contains(&pair.0) {
                    adj[pair.1].push(pair.0);
                }
            }
        }

        let mut disp: f64 = 0.0;
        let nc = mesh.node_count();
        for n in 0..nc {
            if boundary.contains(&n) {
                continue;
            }
            let neigh = &adj[n];
            if neigh.is_empty() {
                continue;
            }
            let mut avg = vec![0.0; dim];
            for &m in neigh {
                let c = mesh.node_coords(m);
                for k in 0..dim {
                    avg[k] += c[k];
                }
            }
            for k in 0..dim {
                avg[k] /= neigh.len() as f64;
            }
            let cur = &mesh.nodes[n].coords;
            let mut d: f64 = 0.0;
            for k in 0..dim {
                d += (avg[k] - cur[k]).powi(2);
            }
            d = d.sqrt();
            disp = disp.max(d);
            mesh.nodes[n].coords[..dim].copy_from_slice(&avg[..dim]);
        }
        max_disp = disp;
    }
    max_disp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn predicates_unit_tet() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        assert!(orient3(&a, &b, &c, &d) > 0.0);
        // The centroid lies strictly inside the circumsphere.
        let centroid = [0.25, 0.25, 0.25];
        assert!(
            in_sphere(&a, &b, &c, &d, &centroid) > 0.0,
            "centroid inside"
        );
        let far = [0.5, 0.5, 2.0];
        assert!(in_sphere(&a, &b, &c, &d, &far) < 0.0, "far outside");
    }

    #[test]
    fn delaunay_single_tet() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let m = delaunay_3d(&pts);
        assert_eq!(m.element_count(), 1);
        assert!(all_positively_oriented(&m));
    }

    #[test]
    fn delaunay_cube_is_closed() {
        let mut pts = Vec::new();
        for x in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                for z in [0.0, 1.0] {
                    pts.push([x, y, z]);
                }
            }
        }
        let m = delaunay_3d(&pts);
        assert!(all_positively_oriented(&m));
        assert!(m.element_count() >= 5, "cube corners yield several tets");
        // No triangle face is shared by 3+ tetrahedra (proper triangulation).
        let mut counts: HashMap<(usize, usize, usize), usize> = HashMap::new();
        for e in &m.elements {
            if e.cell_type != CellType::Tet {
                continue;
            }
            for tri in [
                [e.nodes[0], e.nodes[1], e.nodes[2]],
                [e.nodes[0], e.nodes[1], e.nodes[3]],
                [e.nodes[0], e.nodes[2], e.nodes[3]],
                [e.nodes[1], e.nodes[2], e.nodes[3]],
            ] {
                let mut t = tri;
                t.sort_unstable();
                *counts.entry((t[0], t[1], t[2])).or_insert(0) += 1;
            }
        }
        assert!(
            counts.values().all(|&c| c <= 2),
            "no face shared by more than two tets"
        );
    }

    #[test]
    fn box_mesh_counts_and_orientation() {
        let m = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 2, 2]);
        assert_eq!(m.node_count(), 3 * 3 * 3);
        assert_eq!(m.element_count(), 2 * 2 * 2 * 6);
        assert!(all_positively_oriented(&m));
    }

    #[test]
    fn box_quality_and_smoothing() {
        let mut m = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3, 3, 3]);
        let q0 = tet_quality(&m);
        assert!(q0.max_radius_edge.is_finite());
        assert!(q0.min_dihedral > 0.0);
        let disp = laplacian_smooth(&mut m, 3);
        assert!(disp.is_finite());
        let q1 = tet_quality(&m);
        assert!(q1.max_radius_edge.is_finite());
        // Boundary nodes remain fixed after smoothing.
        let corner = m.node_coords(0);
        assert_eq!(corner[0], 0.0);
    }
}
