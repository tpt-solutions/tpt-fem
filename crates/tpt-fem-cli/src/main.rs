//! `tpt-fem` command-line driver.
//!
//! Subcommands:
//! * `solve` — run a steady Poisson/heat-conduction problem from a TOML config.
//! * `elasticity` / `modal` — elasticity and modal problems from a TOML config.
//! * `amr` — adaptive h-refinement Poisson solve on the unit square
//!   (`tpt-fem-amr::solve_adaptive`): quadtree refinement driven by a
//!   Zienkiewicz–Zhu error estimator with Dörfler marking.
//! * `mesh info` — print summary statistics about a mesh file.
//! * `mesh convert` — convert a Gmsh `.msh` mesh to a ParaView `.vtk` file.
//!
//! Error messages reuse the `Display` impls from the core crates, so malformed
//! input reports a human-readable cause rather than a panic.

use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, Subcommand};
use serde::Deserialize;
use tpt_fem::{
    boundary_faces, box_mesh, read_exodus, read_inp, read_vtk, solve_adaptive, solve_elasticity,
    solve_modal, solve_poisson, write_vtk, write_vtk_with_data, AmrOptions, CellType, ElasticModel,
    Error, Mesh, MeshBuilder, PointData,
};

type Err = Error;

// ---------------------------------------------------------------------------
// CLI surface
// ---------------------------------------------------------------------------

#[derive(Parser)]
#[command(name = "tpt-fem", about = "Finite-element core driver (solve / mesh).")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Solve a problem defined by a TOML config file.
    Solve {
        /// Path to the TOML problem description.
        config: PathBuf,
    },
    /// Solve a linear-elasticity problem (TOML config, `problem.type = "elasticity"`).
    Elasticity {
        /// Path to the TOML problem description.
        config: PathBuf,
    },
    /// Solve a natural-vibration (modal) problem (TOML config, `problem.type = "modal"`).
    Modal {
        /// Path to the TOML problem description.
        config: PathBuf,
    },
    /// Adaptive h-refinement Poisson solve on `[0,1]²`: refine where the
    /// ZZ error estimator says it matters until the element budget is spent.
    Amr {
        /// Element budget for the adaptive loop.
        #[arg(long, default_value_t = 512)]
        max_elements: usize,
        /// Dörfler bulk-marking fraction in `(0, 1]`.
        #[arg(long, default_value_t = 0.7)]
        theta: f64,
        /// Constant volumetric source `f(x) = constant`.
        #[arg(long, default_value_t = 1.0)]
        constant: f64,
        /// Path of the exported ParaView `.vtk` result.
        #[arg(short, long, default_value = "amr.vtk")]
        output: PathBuf,
    },
    /// Generate a starter `problem.toml` for a chosen problem type.
    Init {
        /// Problem type: `poisson`, `elasticity`, or `modal`.
        #[arg(default_value = "poisson")]
        problem: String,
        /// Output config path.
        #[arg(default_value = "problem.toml")]
        output: PathBuf,
    },
    /// Mesh inspection and conversion utilities.
    Mesh {
        #[command(subcommand)]
        action: MeshAction,
    },
}

#[derive(Subcommand)]
enum MeshAction {
    /// Print node/element statistics for a mesh (`.msh` or `.vtk`).
    Info { file: PathBuf },
    /// Convert a Gmsh `.msh` mesh into a ParaView `.vtk` file.
    Convert { input: PathBuf, output: PathBuf },
}

// ---------------------------------------------------------------------------
// Config schema
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct Config {
    #[serde(default)]
    problem: Problem,
    mesh: MeshSpec,
    #[serde(default)]
    material: Material,
    #[serde(default)]
    source: Source,
    #[serde(default)]
    bc: Vec<Bc>,
    #[serde(default)]
    output: Output,
}

#[derive(Deserialize, Default)]
struct Problem {
    /// Problem type: `"poisson"` (steady heat conduction), `"elasticity"`
    /// (linear statics), or `"modal"` (natural vibration modes).
    #[serde(default = "default_problem")]
    r#type: String,
    /// Elasticity model for `problem.type = "elasticity"` / `"modal"`:
    /// `"bar"`, `"plane-stress"`, `"plane-strain"`, or `"3d"`.
    #[serde(default = "default_model")]
    model: String,
    /// Number of vibration modes to extract for `problem.type = "modal"`.
    #[serde(default = "four")]
    num_modes: usize,
}

fn default_problem() -> String {
    "poisson".into()
}

fn default_model() -> String {
    "plane-stress".into()
}

fn four() -> usize {
    4
}

#[derive(Deserialize, Default)]
struct Material {
    /// Constant conductivity `k` in `-∇·(k∇u) = f`.
    #[serde(default = "one")]
    conductivity: f64,
    /// Young's modulus `E` for elasticity / modal problems.
    #[serde(default = "one")]
    young: f64,
    /// Poisson's ratio `ν` for elasticity / modal problems.
    #[serde(default = "default_poisson")]
    poisson: f64,
    /// Mass density `ρ` for modal problems.
    #[serde(default = "one")]
    density: f64,
}

fn one() -> f64 {
    1.0
}

fn default_poisson() -> f64 {
    0.3
}

#[derive(Deserialize, Default)]
struct Source {
    /// Constant volumetric source `f(x) = constant`.
    #[serde(default)]
    constant: f64,
}

#[derive(Deserialize, Default)]
struct Output {
    /// Path of the exported ParaView `.vtk` result.
    #[serde(default = "default_output")]
    vtk: PathBuf,
}

fn default_output() -> PathBuf {
    PathBuf::from("solution.vtk")
}

#[derive(Deserialize)]
struct MeshSpec {
    /// Import an existing Gmsh `.msh` file instead of generating a box.
    #[serde(default)]
    file: Option<PathBuf>,
    /// Spatial dimension of the generated box (2 or 3).
    #[serde(default = "two")]
    dim: usize,
    /// Lower corner of the box.
    #[serde(default = "origin")]
    min: Vec<f64>,
    /// Upper corner of the box.
    #[serde(default = "unit")]
    max: Vec<f64>,
    /// Number of elements per axis.
    #[serde(default = "ten")]
    n: Vec<usize>,
}

fn two() -> usize {
    2
}
fn origin() -> Vec<f64> {
    vec![0.0, 0.0, 0.0]
}
fn unit() -> Vec<f64> {
    vec![1.0, 1.0, 1.0]
}
fn ten() -> Vec<usize> {
    vec![10, 10, 10]
}

#[derive(Deserialize, Default)]
struct PlaneSel {
    axis: usize,
    coord: f64,
    #[serde(default = "default_tol")]
    tol: f64,
}
fn default_tol() -> f64 {
    1e-9
}

#[derive(Deserialize, Default)]
struct BoxSel {
    min: Vec<f64>,
    max: Vec<f64>,
}

#[derive(Deserialize)]
struct Bc {
    /// Prescribed value (Dirichlet).
    value: f64,
    #[serde(default)]
    nodes: Option<Vec<usize>>,
    #[serde(default)]
    plane: Option<PlaneSel>,
    #[serde(default)]
    box_sel: Option<BoxSel>,
    #[serde(default)]
    region: Option<i32>,
    #[serde(default)]
    boundary: bool,
    /// Which displacement/DOF components to constrain (0-based). Defaults to
    /// *all* components of the problem's spatial dimension. This is ignored for
    /// scalar (Poisson) problems, which have a single DOF per node. For vector
    /// problems it lets you pin only, say, the x-component (`[0]`) of a node
    /// while leaving the others free.
    #[serde(default)]
    dofs: Option<Vec<usize>>,
}

/// A resolved boundary-condition hit on a single node, retaining the optional
/// per-component DOF mask so vector solvers can constrain selected components
/// rather than the whole node.
struct BcHit {
    node: usize,
    value: f64,
    dofs: Option<Vec<usize>>,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Err> {
    let cli = Cli::parse();
    match cli.command {
        Command::Solve { config } => solve_config(&config, None),
        Command::Elasticity { config } => solve_config(&config, Some("elasticity")),
        Command::Modal { config } => solve_config(&config, Some("modal")),
        Command::Amr {
            max_elements,
            theta,
            constant,
            output,
        } => run_amr(max_elements, theta, constant, &output),
        Command::Init { problem, output } => init_config(&problem, &output),
        Command::Mesh { action } => match action {
            MeshAction::Info { file } => mesh_info(&file),
            MeshAction::Convert { input, output } => mesh_convert(&input, &output),
        },
    }
}

// ---------------------------------------------------------------------------
// Mesh generation / loading helpers
// ---------------------------------------------------------------------------

/// Build or import the mesh described by `MeshSpec`.
fn build_mesh(spec: &MeshSpec) -> Result<Mesh, Err> {
    if let Some(file) = &spec.file {
        let bytes = std::fs::read(file)?;
        return Ok(Mesh::from_msh_bytes(&bytes)?);
    }
    let dim = if spec.dim == 3 { 3 } else { 2 };
    let min = pad(&spec.min, dim, 0.0);
    let max = pad(&spec.max, dim, 1.0);
    let n = pad_usize(&spec.n, dim, 10);
    if dim == 3 {
        Ok(box_mesh(
            [min[0], min[1], min[2]],
            [max[0], max[1], max[2]],
            [n[0], n[1], n[2]],
        ))
    } else {
        Ok(tri_box_2d([min[0], min[1]], [max[0], max[1]], [n[0], n[1]]))
    }
}

/// Structured triangulated `[min,max]²` box.
fn tri_box_2d(min: [f64; 2], max: [f64; 2], n: [usize; 2]) -> Mesh {
    let mut b = MeshBuilder::new();
    let mut ids = vec![vec![0usize; n[1] + 1]; n[0] + 1];
    for j in 0..=n[1] {
        for i in 0..=n[0] {
            let x = min[0] + (max[0] - min[0]) * i as f64 / n[0] as f64;
            let y = min[1] + (max[1] - min[1]) * j as f64 / n[1] as f64;
            ids[j][i] = b.add_node(vec![x, y]);
        }
    }
    for j in 0..n[1] {
        for i in 0..n[0] {
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

fn pad(v: &[f64], dim: usize, fill: f64) -> Vec<f64> {
    let mut out = vec![fill; dim];
    let n = dim.min(v.len());
    out[..n].copy_from_slice(&v[..n]);
    out
}

fn pad_usize(v: &[usize], dim: usize, fill: usize) -> Vec<usize> {
    let mut out = vec![fill; dim];
    let n = dim.min(v.len());
    out[..n].copy_from_slice(&v[..n]);
    out
}

/// Local node indices of face `fi` for a cell type (mirrors the reference
/// element face definitions used by `tpt_fem_assembly::boundary_faces`).
fn cell_face_nodes(cell: CellType, fi: usize) -> &'static [usize] {
    match cell {
        CellType::Line => &[[0], [1]][fi],
        CellType::Tri => &[[1, 2], [2, 0], [0, 1]][fi],
        CellType::Quad => &[[0, 1], [1, 2], [2, 3], [3, 0]][fi],
        CellType::Tet => &[[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]][fi],
        CellType::Hex => &[
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [3, 2, 6, 7],
            [0, 4, 7, 3],
            [1, 5, 6, 2],
        ][fi],
        // Quadratic (P2) faces are keyed by corner nodes, matching
        // `tpt-fem-assembly`'s `faces_of` definitions.
        CellType::Tri6 => &[[1, 2], [2, 0], [0, 1]][fi],
        CellType::Quad8 => &[[0, 1], [1, 2], [2, 3], [3, 0]][fi],
        CellType::Quad9 => &[[0, 6], [6, 8], [8, 2], [2, 0]][fi],
        CellType::Tet10 => &[[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]][fi],
        CellType::Hex20 => &[
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [3, 2, 6, 7],
            [0, 4, 7, 3],
            [1, 5, 6, 2],
        ][fi],
        CellType::Hex27 => &[
            [0, 18, 24, 6],
            [2, 20, 26, 8],
            [0, 18, 20, 2],
            [6, 24, 26, 8],
            [0, 2, 8, 6],
            [18, 20, 26, 24],
        ][fi],
    }
}

/// Expand a boundary condition into a list of per-node [`BcHit`]s (carrying the
/// optional per-component DOF mask from `bc.dofs`).
fn bc_nodes(mesh: &Mesh, bc: &Bc) -> Vec<BcHit> {
    let mut ids: Vec<usize> = Vec::new();
    if let Some(list) = &bc.nodes {
        ids.extend(list.iter().copied());
    }
    if let Some(p) = &bc.plane {
        ids.extend(mesh.nodes_on_plane(p.axis, p.coord, p.tol));
    }
    if let Some(s) = &bc.box_sel {
        let mn = pad(&s.min, 3, f64::NEG_INFINITY);
        let mx = pad(&s.max, 3, f64::INFINITY);
        let arr = [mn[0], mn[1], mn[2]];
        let arr2 = [mx[0], mx[1], mx[2]];
        ids.extend(mesh.nodes_in_box(arr, arr2));
    }
    if let Some(tag) = bc.region {
        for n in &mesh.nodes {
            if n.region == Some(tag) {
                ids.push(n.id);
            }
        }
    }
    if bc.boundary {
        for (eid, fi) in boundary_faces(mesh) {
            let elem = &mesh.elements[eid];
            for &local in cell_face_nodes(elem.cell_type, fi) {
                ids.push(elem.nodes[local]);
            }
        }
    }
    ids.into_iter()
        .map(|id| BcHit {
            node: id,
            value: bc.value,
            dofs: bc.dofs.clone(),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Subcommands
// ---------------------------------------------------------------------------

/// Starter `problem.toml` templates, keyed by problem type.
const POISSON_TEMPLATE: &str = "\
[problem]
type = \"poisson\"

[mesh]
dim = 2
min = [0.0, 0.0]
max = [1.0, 1.0]
n   = [20, 20]

[material]
conductivity = 1.0

[source]
constant = 1.0

[[bc]]
value = 0.0
boundary = true

[output]
vtk = \"solution.vtk\"
";

const ELASTICITY_TEMPLATE: &str = "\
[problem]
type = \"elasticity\"
model = \"plane-stress\"

[mesh]
dim = 2
min = [0.0, 0.0]
max = [1.0, 1.0]
n   = [20, 20]

[material]
young   = 1.0
poisson = 0.3

# Pin selected DOFs with `dofs = [0]` etc.; omit to fix every component.
[[bc]]
value = 0.0
boundary = true

[output]
vtk = \"displacement.vtk\"
";

const MODAL_TEMPLATE: &str = "\
[problem]
type = \"modal\"
model = \"plane-stress\"
num_modes = 4

[mesh]
dim = 2
min = [0.0, 0.0]
max = [1.0, 1.0]
n   = [20, 20]

[material]
young   = 1.0
poisson = 0.3
density = 1.0

[[bc]]
value = 0.0
boundary = true

[output]
vtk = \"mode1.vtk\"
";

/// Generate a starter problem config for `problem` and write it to `output`.
fn init_config(problem: &str, output: &PathBuf) -> Result<(), Err> {
    let body = match problem.to_ascii_lowercase().as_str() {
        "elasticity" => ELASTICITY_TEMPLATE,
        "modal" => MODAL_TEMPLATE,
        "poisson" => POISSON_TEMPLATE,
        other => {
            return Err(Error::Msg(format!(
                "unknown problem type '{other}' (supported: poisson, elasticity, modal)"
            )));
        }
    };
    std::fs::write(output, body)?;
    println!("Wrote starter '{problem}' config to {}", output.display());
    Ok(())
}

fn solve_config(path: &PathBuf, expected: Option<&str>) -> Result<(), Err> {
    let text = std::fs::read_to_string(path)?;
    let cfg: Config = toml::from_str(&text)?;

    if let Some(expected) = expected {
        if cfg.problem.r#type != expected {
            return Err(Error::Msg(format!(
                "subcommand '{}' requires problem.type = \"{}\", but config declares \"{}\"",
                expected, expected, cfg.problem.r#type
            )));
        }
    }

    let mesh = build_mesh(&cfg.mesh)?;
    println!(
        "Mesh: {} nodes, {} elements",
        mesh.node_count(),
        mesh.element_count()
    );

    let mut hits = Vec::new();
    for bc in &cfg.bc {
        hits.extend(bc_nodes(&mesh, bc));
    }
    println!("Applied {} Dirichlet conditions", hits.len());

    let vtk = &cfg.output.vtk;
    match cfg.problem.r#type.as_str() {
        "poisson" => {
            let f = cfg.source.constant;
            let ndof = mesh.node_count();
            // Scalar problem: one DOF per node, so the per-component `dofs`
            // mask is irrelevant here.
            let bcs: Vec<(usize, f64)> = hits.iter().map(|h| (h.node, h.value)).collect();
            let t0 = Instant::now();
            let u = solve_poisson(
                &mesh,
                cfg.material.conductivity,
                4,
                move |_| f,
                &bcs,
                None,
                None,
            )?;
            let elapsed = t0.elapsed();
            let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
            let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            // Brief report summary so results can be sanity-checked without
            // opening ParaView (DOF count, solve time, solution range).
            println!("DOFs:      {}", ndof);
            println!("Solve time: {:.3?}", elapsed);
            println!("Solution u in [{:.6e}, {:.6e}]", umin, umax);
            write_vtk_with_data(&mesh, &[PointData::new("u", u)], vtk)?;
        }
        "elasticity" => {
            let model = parse_model(&cfg.problem.model)?;
            let dim = cell_dim(mesh.elements[0].cell_type);
            let dir = expand_dirichlet(&hits, dim);
            let ndof = mesh.node_count() * dim;
            let t0 = Instant::now();
            let u = solve_elasticity(
                &mesh,
                model,
                cfg.material.young,
                cfg.material.poisson,
                4,
                |_| vec![0.0; dim],
                &dir,
            )?;
            let elapsed = t0.elapsed();
            let disp = disp_magnitude(&u, dim);
            let dmax = disp.iter().cloned().fold(0.0_f64, f64::max);
            println!("DOFs:      {}", ndof);
            println!("Solve time: {:.3?}", elapsed);
            println!("Max |displacement|: {:.6e}", dmax);
            write_vtk_with_data(&mesh, &[PointData::new("displacement", disp)], vtk)?;
        }
        "modal" => {
            let model = parse_model(&cfg.problem.model)?;
            let dim = cell_dim(mesh.elements[0].cell_type);
            let dir = expand_dirichlet(&hits, dim);
            let ndof = mesh.node_count() * dim;
            let t0 = Instant::now();
            let modes = solve_modal(
                &mesh,
                model,
                cfg.material.young,
                cfg.material.poisson,
                cfg.material.density,
                4,
                cfg.problem.num_modes,
                &dir,
            )?;
            let elapsed = t0.elapsed();
            let (lam, phi) = &modes[0];
            let shape = disp_magnitude(phi, dim);
            println!("DOFs:      {}", ndof);
            println!("Solve time: {:.3?}", elapsed);
            println!(
                "Fundamental frequency ω = {:.6e} (ω² = {:.6e})",
                lam.sqrt(),
                lam
            );
            write_vtk_with_data(&mesh, &[PointData::new("mode1", shape)], vtk)?;
        }
        other => {
            return Err(Error::Msg(format!(
                "unsupported problem type '{other}' (supported: poisson, elasticity, modal)"
            )));
        }
    }
    println!("Wrote {}", vtk.display());
    Ok(())
}

/// Adaptive h-refinement Poisson driver (`amr` subcommand).
///
/// Runs [`tpt_fem::solve_adaptive`] on `[0,1]²` with a constant source and
/// homogeneous Dirichlet boundary data, then exports the final quadtree leaf
/// mesh and solution as a ParaView `.vtk` file.
fn run_amr(max_elements: usize, theta: f64, constant: f64, output: &PathBuf) -> Result<(), Err> {
    if theta <= 0.0 || theta > 1.0 {
        return Err(Error::Msg(format!(
            "--theta must be in (0, 1], got {theta}"
        )));
    }
    let f = move |_: f64, _: f64| constant;
    let g = |_: f64, _: f64| 0.0;
    let t0 = Instant::now();
    let res = solve_adaptive(
        &f,
        &g,
        &AmrOptions {
            max_elements,
            theta,
        },
    )?;
    let elapsed = t0.elapsed();

    // Convert the quadtree leaf mesh into the workspace `Mesh` model so the
    // standard VTK export path can be reused.
    let mut b = MeshBuilder::new();
    for c in &res.mesh.coords {
        b.add_node(vec![c[0], c[1]]);
    }
    for e in &res.mesh.elems {
        b.add_element(CellType::Quad, e.to_vec());
    }
    let mesh = b.build();

    let umin = res.u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = res.u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // Brief report summary so results can be sanity-checked without opening
    // ParaView (element count, estimated error, wall-clock time).
    println!("Elements:            {}", res.mesh.len());
    println!("Nodes:               {}", res.mesh.coords.len());
    println!("DOFs:                {}", res.mesh.coords.len());
    println!("Solve time:          {:.3?}", elapsed);
    println!("Estimated error:     {:.6e}", res.estimated_error);
    println!("Solution u in [{:.6e}, {:.6e}]", umin, umax);
    write_vtk_with_data(&mesh, &[PointData::new("u", res.u)], output)?;
    println!("Wrote {}", output.display());
    Ok(())
}

/// Reference (spatial) dimension of a cell type.
fn cell_dim(cell: CellType) -> usize {
    match cell {
        CellType::Line => 1,
        CellType::Tri | CellType::Quad | CellType::Tri6 | CellType::Quad8 | CellType::Quad9 => 2,
        CellType::Tet | CellType::Hex | CellType::Tet10 | CellType::Hex20 | CellType::Hex27 => 3,
    }
}

/// Parse an elasticity/model string into [`ElasticModel`].
fn parse_model(s: &str) -> Result<ElasticModel, Err> {
    match s.to_ascii_lowercase().replace(['-', '_'], "").as_str() {
        "bar" | "baraxial" => Ok(ElasticModel::BarAxial),
        "planestress" => Ok(ElasticModel::PlaneStress),
        "planestrain" => Ok(ElasticModel::PlaneStrain),
        "3d" | "continuum" | "continuum3d" => Ok(ElasticModel::Continuum3D),
        other => Err(Error::Msg(format!("unknown elasticity model '{other}"))),
    }
}

/// Expand resolved BC hits into per-DOF `(global_dof, value)` pairs for a
/// `dim`-component vector problem. A hit without an explicit `dofs` mask
/// constrains every component of the node; otherwise only the listed
/// components (clamped to `[0, dim)`) are fixed.
fn expand_dirichlet(hits: &[BcHit], dim: usize) -> Vec<(usize, f64)> {
    let mut dir = Vec::new();
    for h in hits {
        let comps: Vec<usize> = match &h.dofs {
            Some(v) => v.iter().copied().filter(|&c| c < dim).collect(),
            None => (0..dim).collect(),
        };
        for c in comps {
            dir.push((h.node * dim + c, h.value));
        }
    }
    dir
}

/// Per-node displacement magnitude from a `node_count * dim` solution vector.
fn disp_magnitude(u: &[f64], dim: usize) -> Vec<f64> {
    let n = u.len() / dim;
    (0..n)
        .map(|i| {
            let mut s = 0.0;
            for c in 0..dim {
                s += u[i * dim + c] * u[i * dim + c];
            }
            s.sqrt()
        })
        .collect()
}

fn mesh_info(path: &PathBuf) -> Result<(), Err> {
    let mesh = load_mesh(path)?;
    println!("nodes:    {}", mesh.node_count());
    println!("elements: {}", mesh.element_count());
    let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for e in &mesh.elements {
        *counts.entry(e.cell_type.name().to_string()).or_insert(0) += 1;
    }
    for (cell, c) in counts {
        println!("  {:>6} : {}", cell, c);
    }
    if mesh.node_count() > 0 {
        let dim = mesh.node_coords(0).len();
        let mut lo = vec![f64::INFINITY; dim];
        let mut hi = vec![f64::NEG_INFINITY; dim];
        for n in &mesh.nodes {
            for d in 0..dim {
                lo[d] = lo[d].min(n.coords[d]);
                hi[d] = hi[d].max(n.coords[d]);
            }
        }
        print!("bounds:   ");
        for d in 0..dim {
            print!("[{}..{}] ", lo[d], hi[d]);
        }
        println!();
    }
    Ok(())
}

fn mesh_convert(input: &PathBuf, output: &PathBuf) -> Result<(), Err> {
    let mesh = load_mesh(input)?;
    write_vtk(&mesh, output)?;
    println!(
        "Converted {} nodes / {} elements -> {}",
        mesh.node_count(),
        mesh.element_count(),
        output.display()
    );
    Ok(())
}

/// Load a mesh: `.msh` (Gmsh), `.inp` (Abaqus), `.ex`/`.ex2`/`.e` (Exodus), or
/// `.vtk` (re-imported through `tpt-fem-io-vtk`'s reader).
fn load_mesh(path: &PathBuf) -> Result<Mesh, Err> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase();
    match ext.as_str() {
        "msh" => {
            let bytes = std::fs::read(path)?;
            Ok(Mesh::from_msh_bytes(&bytes)?)
        }
        "inp" => {
            let text = std::fs::read_to_string(path)?;
            Ok(read_inp(&text).map_err(|e| Error::Msg(format!("abaqus import: {e}")))?)
        }
        // Exodus II (NetCDF-3) mesh.
        "ex" | "ex2" | "e" => {
            Ok(read_exodus(path).map_err(|e| Error::Msg(format!("exodus import: {e}")))?)
        }
        // VTK (re-imported through `tpt-fem-io-vtk`'s crate-level reader,
        // promoted out of the CLI in Phase 10c).
        "vtk" => Ok(read_vtk(path).map_err(|e| Error::Msg(format!("vtk import: {e}")))?),
        other => Err(Error::Msg(format!(
            "unsupported mesh extension '.{other}' (use .msh, .inp, .ex, or .vtk)"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    /// A minimal valid Gmsh v4.1 `.msh` (two triangles on the unit square).
    const TRI_MSH: &str = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
2 1 0 4
1 2 3 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
1.0 1.0 0.0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
1 1 2 3
2 2 4 3
$EndElements
";

    fn write_temp(name: &str, contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, contents).expect("write temp file");
        path
    }

    #[test]
    fn solve_runs_and_writes_vtk() {
        // Box mesh, homogeneous Dirichlet on the boundary, zero source => the
        // trivial solution u = 0. The test just verifies the end-to-end path
        // (config parse -> build -> assemble -> solve -> VTK export) succeeds
        // and produces a non-empty output file.
        let out = std::env::temp_dir().join("tpt_fem_cli_solve_test.vtk");
        // Windows temp paths contain backslashes, which TOML would parse as
        // escape sequences; use forward slashes inside the generated config.
        let out_toml = out.to_string_lossy().replace('\\', "/");
        let cfg = format!(
            "[problem]\ntype = \"poisson\"\n\
             [mesh]\ndim = 2\nmin = [0.0, 0.0]\nmax = [1.0, 1.0]\nn = [4, 4]\n\
             [material]\nconductivity = 1.0\n[source]\nconstant = 0.0\n\
             [[bc]]\nvalue = 0.0\nboundary = true\n\
             [output]\nvtk = \"{}\"\n",
            out_toml
        );
        let cfg_path = write_temp("tpt_fem_cli_solve_test.toml", &cfg);
        solve_config(&cfg_path, None).expect("solve_config");
        let meta = std::fs::metadata(&out).expect("output written");
        assert!(meta.len() > 0, "VTK output should be non-empty");
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(&cfg_path);
    }

    #[test]
    fn mesh_info_reports_counts() {
        let path = write_temp("tpt_fem_cli_info_test.msh", TRI_MSH);
        // Should not error and should report 4 nodes / 2 elements.
        mesh_info(&path).expect("mesh_info");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn mesh_convert_produces_vtk() {
        let msh = write_temp("tpt_fem_cli_convert_in.msh", TRI_MSH);
        let out = std::env::temp_dir().join("tpt_fem_cli_convert_out.vtk");
        mesh_convert(&msh, &out).expect("mesh_convert");
        let meta = std::fs::metadata(&out).expect("converted file written");
        assert!(meta.len() > 0);
        let _ = std::fs::remove_file(&msh);
        let _ = std::fs::remove_file(&out);
    }

    #[test]
    fn load_mesh_handles_inp_and_exodus() {
        // Both the Abaqus and Exodus readers must round-trip the same triangle.
        let msh = write_temp("tpt_fem_cli_load.msh", TRI_MSH);
        let mesh0 = load_mesh(&msh).expect("load .msh");
        let n = mesh0.node_count();
        let e = mesh0.element_count();

        // Abaqus: export via the writer, re-import through the CLI loader.
        let inp_path = std::env::temp_dir().join("tpt_fem_cli_load.inp");
        tpt_fem::write_inp(&mesh0, &inp_path).expect("write .inp");
        let mesh_inp = load_mesh(&inp_path).expect("load .inp");
        assert_eq!(mesh_inp.node_count(), n);
        assert_eq!(mesh_inp.element_count(), e);

        // Exodus: export via the writer, re-import through the CLI loader.
        let ex_path = std::env::temp_dir().join("tpt_fem_cli_load.ex");
        tpt_fem::write_exodus(&mesh0, &ex_path).expect("write .ex");
        let mesh_ex = load_mesh(&ex_path).expect("load .ex");
        assert_eq!(mesh_ex.node_count(), n);
        assert_eq!(mesh_ex.element_count(), e);

        let _ = std::fs::remove_file(&msh);
        let _ = std::fs::remove_file(&inp_path);
        let _ = std::fs::remove_file(&ex_path);
    }

    #[test]
    fn init_writes_starter_config() {
        // `init` must produce a parseable starter config for each problem
        // type; the generated Poisson config then drives a full solve.
        for ptype in ["poisson", "elasticity", "modal"] {
            let cfg_path = std::env::temp_dir().join(format!("tpt_fem_cli_init_{ptype}.toml"));
            init_config(ptype, &cfg_path).unwrap_or_else(|e| panic!("{ptype}: {e}"));
            let text = std::fs::read_to_string(&cfg_path).expect("read generated config");
            assert!(
                text.contains("[problem]"),
                "{ptype} config missing [problem]"
            );
            assert!(text.contains("type = "), "{ptype} config missing type");
            let _ = std::fs::remove_file(&cfg_path);
        }
        // And an unknown problem type is rejected, not silently written.
        let bad = std::env::temp_dir().join("tpt_fem_cli_init_bad.toml");
        assert!(init_config("nope", &bad).is_err());
        let _ = std::fs::remove_file(&bad);
    }

    #[test]
    fn cli_usage_matches_readme() {
        // Drift guard: the README documents `cargo run -p tpt-fem-cli -- solve
        // problem.toml`, i.e. `solve` takes a single *positional* config file.
        // If the clap definition changes (e.g. becomes `--config`/flags), this
        // fails loudly instead of the docs silently diverging.
        let mut cmd = Cli::command();
        for sub in ["solve", "elasticity", "modal"] {
            let usage = cmd
                .find_subcommand_mut(sub)
                .expect("subcommand exists")
                .render_usage()
                .to_string();
            assert!(
                usage.contains(format!("{sub} <CONFIG>").as_str()),
                "{sub} usage drifted from README: {usage}"
            );
        }
        // And the bare-positional form must actually parse.
        Cli::command().debug_assert()
    }

    #[test]
    fn expand_dirichlet_honors_dof_mask() {
        // Per-component Dirichlet: node 0 pinned in x only (dofs=[0]); node 1
        // pinned in every component (default). The expansion must constrain
        // exactly the listed global DOFs.
        let hits = vec![
            BcHit {
                node: 0,
                value: 0.0,
                dofs: Some(vec![0]),
            },
            BcHit {
                node: 1,
                value: 0.0,
                dofs: None,
            },
        ];
        let dir = expand_dirichlet(&hits, 2);
        assert!(dir.contains(&(0, 0.0))); // node 0, dof x
        assert!(!dir.contains(&(1, 0.0))); // node 0, dof y must NOT be fixed
        assert!(dir.contains(&(2, 0.0))); // node 1, dof x
        assert!(dir.contains(&(3, 0.0))); // node 1, dof y
                                          // Out-of-range components are clamped to [0, dim).
        let hits2 = vec![BcHit {
            node: 2,
            value: 0.0,
            dofs: Some(vec![0, 5]),
        }];
        assert_eq!(expand_dirichlet(&hits2, 2), vec![(4, 0.0)]);
    }

    #[test]
    fn elasticity_runs_and_writes_vtk() {
        // Clamped 2-D box, zero body load => the trivial displacement field.
        // Verifies the `elasticity` subcommand path (config -> mesh -> assemble
        // -> solve_elasticity -> VTK) succeeds and writes a non-empty file.
        let out = std::env::temp_dir().join("tpt_fem_cli_elasticity_test.vtk");
        let out_toml = out.to_string_lossy().replace('\\', "/");
        let cfg = format!(
            "[problem]\ntype = \"elasticity\"\nmodel = \"plane-stress\"\n\
             [mesh]\ndim = 2\nmin = [0.0, 0.0]\nmax = [1.0, 1.0]\nn = [4, 4]\n\
             [material]\nyoung = 1.0\npoisson = 0.3\n\
             [[bc]]\nvalue = 0.0\nboundary = true\n\
             [output]\nvtk = \"{}\"\n",
            out_toml
        );
        let cfg_path = write_temp("tpt_fem_cli_elasticity_test.toml", &cfg);
        solve_config(&cfg_path, Some("elasticity")).expect("elasticity solve_config");
        let meta = std::fs::metadata(&out).expect("output written");
        assert!(meta.len() > 0, "VTK output should be non-empty");
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(&cfg_path);
    }

    #[test]
    fn modal_runs_and_writes_vtk() {
        // Clamped 2-D box modal problem: must extract the fundamental mode
        // without error and write a non-empty file.
        let out = std::env::temp_dir().join("tpt_fem_cli_modal_test.vtk");
        let out_toml = out.to_string_lossy().replace('\\', "/");
        let cfg = format!(
            "[problem]\ntype = \"modal\"\nmodel = \"plane-stress\"\nnum_modes = 2\n\
             [mesh]\ndim = 2\nmin = [0.0, 0.0]\nmax = [1.0, 1.0]\nn = [6, 6]\n\
             [material]\nyoung = 1.0\npoisson = 0.3\ndensity = 1.0\n\
             [[bc]]\nvalue = 0.0\nboundary = true\n\
             [output]\nvtk = \"{}\"\n",
            out_toml
        );
        let cfg_path = write_temp("tpt_fem_cli_modal_test.toml", &cfg);
        solve_config(&cfg_path, Some("modal")).expect("modal solve_config");
        let meta = std::fs::metadata(&out).expect("output written");
        assert!(meta.len() > 0, "VTK output should be non-empty");
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(&cfg_path);
    }

    #[test]
    fn amr_runs_and_writes_vtk() {
        // Small element budget so the adaptive loop finishes quickly: must
        // produce a non-empty VTK with the refined quadtree mesh + solution.
        let out = std::env::temp_dir().join("tpt_fem_cli_amr_test.vtk");
        run_amr(64, 0.7, 1.0, &out).expect("amr run");
        let meta = std::fs::metadata(&out).expect("output written");
        assert!(meta.len() > 0, "VTK output should be non-empty");
        let _ = std::fs::remove_file(&out);
    }

    #[test]
    fn amr_rejects_invalid_theta() {
        let out = std::env::temp_dir().join("tpt_fem_cli_amr_bad_theta.vtk");
        assert!(run_amr(64, 0.0, 1.0, &out).is_err());
        assert!(run_amr(64, 1.5, 1.0, &out).is_err());
    }

    #[test]
    fn elasticity_subcommand_rejects_non_elastic_config() {
        // Regression for the cosmetic-subcommand bug: invoking `elasticity`
        // on a Poisson config must error rather than silently solving it.
        let out = std::env::temp_dir().join("tpt_fem_cli_typecheck_test.vtk");
        let out_toml = out.to_string_lossy().replace('\\', "/");
        let cfg = format!(
            "[problem]\ntype = \"poisson\"\n\
             [mesh]\ndim = 2\nmin = [0.0, 0.0]\nmax = [1.0, 1.0]\nn = [4, 4]\n\
             [material]\nconductivity = 1.0\n[source]\nconstant = 0.0\n\
             [[bc]]\nvalue = 0.0\nboundary = true\n\
             [output]\nvtk = \"{}\"\n",
            out_toml
        );
        let cfg_path = write_temp("tpt_fem_cli_typecheck_test.toml", &cfg);
        let result = solve_config(&cfg_path, Some("elasticity"));
        assert!(
            result.is_err(),
            "elasticity subcommand must reject poisson config"
        );
        let _ = std::fs::remove_file(&cfg_path);
    }
}
