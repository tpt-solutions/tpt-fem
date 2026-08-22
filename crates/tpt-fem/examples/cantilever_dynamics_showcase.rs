//! Structural-dynamics pipeline: assembly -> modal analysis -> time history.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-fem --example cantilever_dynamics_showcase --features modal,dynamic
//! ```
//!
//! A thin plane-stress cantilever (Quad4) is modelled end to end:
//!
//! 1. **Assembly** — the stiffness comes from `tpt-fem-elasticity` +
//!    `tpt-fem-assembly`; the mass is lumped consistently with the element
//!    areas (each Quad4 shares its `rho*t*A` equally between its four nodes).
//!    The clamped DOFs are eliminated by zeroing rows/columns with unit
//!    diagonals, which keeps the eigenproblem symmetric.
//!
//! 2. **Modal analysis** — `modal_analysis` extracts the lowest modes via a
//!    generalized Lanczos solve. The first bending frequency is compared with
//!    Euler-Bernoulli beam theory,
//!    `omega_1 = (1.875)^2 / L^2 * sqrt(E*I / (rho*A))`; a coarse Quad4 mesh is
//!    slightly stiff in bending (shear locking of bilinear quads), so FEM sits
//!    a few percent above theory.
//!
//! 3. **Transient response** — a step load at the tip is integrated two ways:
//!    full-system implicit Newmark and truncated modal superposition. With
//!    enough modes retained the two histories must agree closely, which
//!    validates both paths against each other on the same FE model.

use tpt_fem_assembly::assemble;
use tpt_fem_dynamic::{newmark, NewmarkOptions};
use tpt_fem_elasticity::{elasticity_element_matrix, ElasticModel};
use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_fem_modal::modal_analysis;
use tpt_fem_sparse::Coo;

/// Structured `nx` x `ny` Quad4 strip of length `l` and thickness `h`.
fn strip(l: f64, h: f64, nx: usize, ny: usize) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut rows = Vec::new();
    for j in 0..=ny {
        let y = h * j as f64 / ny as f64;
        let mut r = Vec::new();
        for i in 0..=nx {
            r.push(b.add_node(vec![l * i as f64 / nx as f64, y]));
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
    b.build()
}

/// Drop all triplets touching a `fixed` DOF and renumber the remaining ones
/// compactly, returning the reduced matrix and the `old -> new` index map
/// (`usize::MAX` for eliminated DOFs).
fn reduce(coo: &Coo, fixed: &[usize], n: usize) -> (Coo, Vec<usize>) {
    let mut map = vec![usize::MAX; n];
    let mut next = 0usize;
    for d in 0..n {
        if !fixed.contains(&d) {
            map[d] = next;
            next += 1;
        }
    }
    let mut out = Coo {
        rows: Vec::new(),
        cols: Vec::new(),
        vals: Vec::new(),
    };
    for i in 0..coo.rows.len() {
        let (r, c) = (map[coo.rows[i]], map[coo.cols[i]]);
        if r != usize::MAX && c != usize::MAX {
            out.rows.push(r);
            out.cols.push(c);
            out.vals.push(coo.vals[i]);
        }
    }
    (out, map)
}

fn main() {
    let l = 1.0;
    let h = 0.05;
    let (nx, ny) = (20usize, 2usize);
    let (e, nu, rho) = (70.0e9, 0.33, 2700.0); // aluminium

    // NOTE: the plane-stress operators are per unit thickness, so the whole
    // model is unit-thick; frequencies are thickness-independent anyway since
    // both bending stiffness and mass scale linearly with thickness.
    let mesh = strip(l, h, nx, ny);
    let ndof = mesh.node_count() * 2;

    // Stiffness: plane-stress Quad4 assembly.
    let k_full: Coo = assemble(&mesh, 2, |eid, _m| {
        elasticity_element_matrix(&mesh, eid, ElasticModel::PlaneStress, e, nu, 2)
            .expect("element stiffness")
    });

    // Lumped diagonal mass: element area shared by its 4 nodes (unit thickness).
    let dx = l / nx as f64;
    let dy = h / ny as f64;
    let m_elem = rho * dx * dy;
    let mut mdiag = vec![0.0_f64; ndof];
    for elem in &mesh.elements {
        for &n in &elem.nodes {
            mdiag[n * 2] += m_elem / 4.0;
            mdiag[n * 2 + 1] += m_elem / 4.0;
        }
    }
    let m_full = Coo {
        rows: (0..ndof).collect(),
        cols: (0..ndof).collect(),
        vals: mdiag,
    };

    // Clamp x = 0 in both components and reduce to the free DOFs.
    let fixed: Vec<usize> = (0..mesh.node_count())
        .filter(|&n| mesh.node_coords(n)[0] < 1e-9)
        .flat_map(|n| [n * 2, n * 2 + 1])
        .collect();
    let (k, map) = reduce(&k_full, &fixed, ndof);
    let (m, _) = reduce(&m_full, &fixed, ndof);
    let n_red = k
        .rows
        .iter()
        .cloned()
        .chain(std::iter::once(0))
        .max()
        .unwrap()
        + 1;

    // ---- Modal analysis vs beam theory ---------------------------------
    let num_modes = 8;
    let data = modal_analysis(&k, &m, 0.0, num_modes, 2 * num_modes, 0.0).expect("modal analysis");
    println!(
        "cantilever L = {l}, h = {h}, {nx} x {ny} Quad4, {} DOFs",
        ndof - fixed.len()
    );
    for (i, w) in data.omega.iter().enumerate() {
        println!(
            "  f{:2} = {:10.3} Hz",
            i + 1,
            w / (2.0 * std::f64::consts::PI)
        );
    }
    let beam_omega1 = 1.875_f64.powi(2) / l / l * (e * h * h / (12.0 * rho)).sqrt();
    let ratio = data.omega[0] / beam_omega1;
    println!(
        "\n  omega_1 FEM  = {:.4e} rad/s\n  omega_1 beam = {:.4e} rad/s\n  ratio        = {:.4}",
        data.omega[0], beam_omega1, ratio
    );
    assert!(
        (0.95..=1.25).contains(&ratio),
        "first bending mode must be near Euler-Bernoulli theory (slightly stiff)"
    );

    // ---- Transient response: Newmark vs modal superposition ------------
    let tip_node = (0..mesh.node_count())
        .filter(|&n| (mesh.node_coords(n)[1] - h / 2.0).abs() < 1e-12)
        .max_by(|&a, &b| {
            mesh.node_coords(a)[0]
                .partial_cmp(&mesh.node_coords(b)[0])
                .unwrap()
        })
        .expect("mid-surface tip node");
    let tip_dof = tip_node * 2 + 1;
    let tip_red = map[tip_dof];
    assert!(tip_red != usize::MAX, "tip dof must be free");

    let f0 = 1000.0; // N step load at the tip
    let mut force = vec![0.0_f64; n_red];
    force[tip_red] = f0;
    #[allow(unused_mut)]
    let load = move |_t: f64| force.clone();

    let opts = NewmarkOptions {
        dt: 2.0e-5,
        beta: 0.25,
        gamma: 0.5,
    };
    let nsteps = 500;
    let zero: Vec<f64> = vec![0.0; n_red];
    let no_damping = Coo {
        rows: vec![],
        cols: vec![],
        vals: vec![],
    };
    let hist_full = newmark(
        &m,
        &no_damping,
        &k,
        &zero,
        &zero,
        load.clone(),
        &opts,
        nsteps,
    );
    let hist_modal = data.modal_superposition(&zero, &zero, &load, &opts, nsteps);

    println!("\ntip displacement history (step load {f0:.0} N):");
    println!(
        "  {:>8}  {:>13}  {:>13}  {:>8}",
        "t [ms]", "Newmark", "modal", "diff"
    );
    let mut max_rel = 0.0_f64;
    for step in [50usize, 150, 250, 350, 500] {
        let ((_, uf), (_, um)) = (&hist_full[step], &hist_modal[step]);
        let d = (uf[tip_red] - um[tip_red]).abs();
        max_rel = max_rel.max(d / uf[tip_red].abs().max(1e-30));
        println!(
            "  {:8.2}  {:13.6e}  {:13.6e}  {:8.1e}",
            step as f64 * opts.dt * 1.0e3,
            uf[tip_red],
            um[tip_red],
            d
        );
    }
    assert!(
        max_rel < 0.02,
        "modal superposition must track the full Newmark solve (max rel diff {max_rel:.3e})"
    );
    println!(
        "\nOK: modes match beam theory within {:.1}% and both integrators agree.",
        100.0 * (ratio - 1.0)
    );
}
