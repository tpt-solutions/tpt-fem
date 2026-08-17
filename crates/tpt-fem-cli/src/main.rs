//! `tpt-fem` command-line driver.
//!
//! Subcommands:
//! * `solve` — run a steady Poisson/heat-conduction problem from a TOML config.
//! * `mesh info` — print summary statistics about a mesh file.
//! * `mesh convert` — convert a Gmsh `.msh` mesh to a ParaView `.vtk` file.
//!
//! Error messages reuse the `Display` impls from the core crates, so malformed
//! input reports a human-readable cause rather than a panic.

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use serde::Deserialize;
use tpt_fem::{
    boundary_faces, box_mesh, solve_poisson, write_vtk, write_vtk_with_data, CellType, Mesh,
    MeshBuilder, MeshError, PointData,
};

type Err = Box<dyn std::error::Error>;

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
    /// Problem type; only `"poisson"` (heat conduction) is supported this pass.
    #[serde(default = "default_problem")]
    r#type: String,
}

fn default_problem() -> String {
    "poisson".into()
}

#[derive(Deserialize, Default)]
struct Material {
    /// Constant conductivity `k` in `-∇·(k∇u) = f`.
    #[serde(default = "one")]
    conductivity: f64,
}

fn one() -> f64 {
    1.0
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
    /// Prescribed value (Dirichlet, this pass).
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
        Command::Solve { config } => solve_config(&config),
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

/// Expand a boundary condition into a `(node, value)` list.
fn bc_nodes(mesh: &Mesh, bc: &Bc) -> Vec<(usize, f64)> {
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
    ids.into_iter().map(|id| (id, bc.value)).collect()
}

// ---------------------------------------------------------------------------
// Subcommands
// ---------------------------------------------------------------------------

fn solve_config(path: &PathBuf) -> Result<(), Err> {
    let text = std::fs::read_to_string(path)?;
    let cfg: Config = toml::from_str(&text)?;

    if cfg.problem.r#type != "poisson" {
        return Err(format!(
            "unsupported problem type '{}' (only 'poisson' is supported this pass)",
            cfg.problem.r#type
        )
        .into());
    }

    let mesh = build_mesh(&cfg.mesh)?;
    println!(
        "Mesh: {} nodes, {} elements",
        mesh.node_count(),
        mesh.element_count()
    );

    let mut bcs = Vec::new();
    for bc in &cfg.bc {
        let n = bc_nodes(&mesh, bc);
        bcs.extend(n);
    }
    println!("Applied {} Dirichlet conditions", bcs.len());

    let f = cfg.source.constant;
    let u = solve_poisson(
        &mesh,
        cfg.material.conductivity,
        4,
        move |_| f,
        &bcs,
        None,
        None,
    )?;

    let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("Solution u in [{:.6e}, {:.6e}]", umin, umax);

    let vtk = &cfg.output.vtk;
    write_vtk_with_data(&mesh, &[PointData::new("u", u)], vtk)?;
    println!("Wrote {}", vtk.display());
    Ok(())
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

/// Load a mesh: `.msh` (Gmsh) or `.vtk` (re-imported through vtkio).
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
        "vtk" => {
            let vtk = vtkio::model::Vtk::import(path)
                .map_err(|e| MeshError::Parse(format!("vtk import: {e}")))?;
            mesh_from_vtk(&vtk)
        }
        other => Err(format!("unsupported mesh extension '.{other}' (use .msh or .vtk)").into()),
    }
}

/// Minimal VTK → Mesh reader for the linear cell types this crate writes.
fn mesh_from_vtk(vtk: &vtkio::model::Vtk) -> Result<Mesh, Err> {
    use vtkio::model::{DataSet, Piece};
    let ds = &vtk.data;
    let (points, cells) = match ds {
        DataSet::UnstructuredGrid { pieces, .. } => match &pieces[0] {
            Piece::Inline(p) => {
                let pts = match &p.points {
                    vtkio::model::IOBuffer::F64(v) => v.clone(),
                    _ => return Err("unsupported point coordinate type".into()),
                };
                (pts, p.cells.clone())
            }
            _ => return Err("expected inline VTK piece".into()),
        },
        _ => return Err("expected unstructured grid".into()),
    };
    let np = points.len() / 3;
    let mut b = MeshBuilder::new();
    for i in 0..np {
        b.add_node(vec![points[3 * i], points[3 * i + 1], points[3 * i + 2]]);
    }
    let verts = match &cells.cell_verts {
        vtkio::model::VertexNumbers::Legacy { vertices, .. } => vertices.clone(),
        _ => return Err("unsupported cell numbering".into()),
    };
    let mut i = 0;
    while i < verts.len() {
        let cnt = verts[i] as usize;
        let cell = match cnt {
            2 => CellType::Line,
            3 => CellType::Tri,
            4 => CellType::Quad,
            6 => CellType::Tri6,
            8 => CellType::Hex,
            9 => CellType::Quad9,
            10 => CellType::Tet10,
            20 => CellType::Hex20,
            27 => CellType::Hex27,
            _ => return Err(format!("unsupported cell with {cnt} nodes").into()),
        };
        let nodes = verts[i + 1..i + 1 + cnt]
            .iter()
            .map(|&v| v as usize)
            .collect();
        b.add_element(cell, nodes);
        i += 1 + cnt;
    }
    Ok(b.build())
}
