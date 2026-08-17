//! Exodus II mesh reader/writer for `tpt-fem`.
//!
//! Exodus II is a NetCDF-based format with no mature pure-Rust implementation,
//! so this crate ships a minimal NetCDF-3 *classic* (CDF-1, big-endian) codec
//! and builds the subset of Exodus variables required to round-trip linear
//! meshes: `coords`, one `connectN` set per element block, `elem_blk_id`,
//! `eb_status`, `eb_prop1`, `eb_names`, `elem_num_map`, `node_num_map`, and
//! `time_whole`. Only the five linear element types are supported.

use std::collections::HashMap;
use std::path::Path;

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Errors returned by the Exodus reader/writer.
#[derive(Debug)]
pub enum ExodusError {
    /// A malformed NetCDF/Exodus structure.
    Parse(String),
    /// An unsupported construct (e.g. record variables, non-3D coordinates).
    Unsupported(String),
    /// The mesh built from the file failed validation (e.g. a node index out
    /// of range).
    Mesh(tpt_fem_mesh::MeshError),
    /// An I/O error occurred while reading or writing the file.
    Io(std::io::Error),
}

impl std::fmt::Display for ExodusError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExodusError::Parse(m) => write!(f, "malformed Exodus/NetCDF file: {m}"),
            ExodusError::Unsupported(m) => write!(f, "unsupported Exodus construct: {m}"),
            ExodusError::Mesh(e) => write!(f, "invalid mesh in Exodus file: {e}"),
            ExodusError::Io(e) => write!(f, "I/O error reading/writing Exodus file: {e}"),
        }
    }
}

impl std::error::Error for ExodusError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ExodusError::Mesh(e) => Some(e),
            ExodusError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<tpt_fem_mesh::MeshError> for ExodusError {
    fn from(e: tpt_fem_mesh::MeshError) -> Self {
        ExodusError::Mesh(e)
    }
}

impl From<std::io::Error> for ExodusError {
    fn from(e: std::io::Error) -> Self {
        ExodusError::Io(e)
    }
}

// ---------------------------------------------------------------------------
// Minimal NetCDF-3 classic (CDF-1) codec
// ---------------------------------------------------------------------------

const NC_DIMENSION: u32 = 10;
const NC_VARIABLE: u32 = 11;
const NC_FLOAT: u32 = 5;
const NC_INT: u32 = 4;
const NC_CHAR: u32 = 2;
const LEN_NAME: u32 = 33;

struct NcDim {
    name: String,
    size: u32,
}

struct NcVar {
    name: String,
    dtype: u32,
    dim_ids: Vec<u32>,
    /// Pre-encoded big-endian payload (length must equal `vsize`).
    data: Vec<u8>,
}

fn to_u32(v: u32) -> [u8; 4] {
    v.to_be_bytes()
}

/// Encode a NetCDF name (1 length byte + chars, padded to a 4-byte boundary).
fn enc_name(name: &str) -> Vec<u8> {
    let bytes = name.as_bytes();
    let mut out = Vec::with_capacity(4);
    out.push(bytes.len().min(255) as u8);
    out.extend_from_slice(bytes);
    while out.len() % 4 != 0 {
        out.push(0);
    }
    out
}

fn dim_size(dims: &[NcDim], id: u32) -> Result<u32, ExodusError> {
    let d = usize::try_from(id)
        .map_err(|_| ExodusError::Parse(format!("dimension id {id} overflows")))?;
    dims.get(d)
        .map(|dim| dim.size)
        .ok_or_else(|| ExodusError::Parse(format!("dimension id {id} out of range")))
}

fn vsize_of(dims: &[NcDim], v: &NcVar) -> Result<u32, ExodusError> {
    let mut s: u32 = match v.dtype {
        NC_CHAR => 1,
        NC_INT | NC_FLOAT => 4,
        _ => 4,
    };
    for &d in &v.dim_ids {
        s = s
            .checked_mul(dim_size(dims, d)?)
            .ok_or_else(|| ExodusError::Parse("variable size overflows u32".into()))?;
    }
    Ok(s)
}

/// Reject a count read from an untrusted file before using it to size an
/// allocation: a corrupted `numrecs`/dimension/variable count could otherwise
/// ask for an absurd `Vec::with_capacity`, which aborts the process.
fn sane_count(count: usize, bytes: &[u8]) -> Result<usize, ExodusError> {
    // Every entry consumes at least 4 bytes (a name length byte plus padding,
    // or a u32), so the count can never legitimately exceed this bound.
    if count > bytes.len() / 4 + 8 {
        return Err(ExodusError::Parse(format!("implausible count {count}")));
    }
    Ok(count)
}

/// Encode a file from dimensions and variables into NetCDF-3 classic bytes.
///
/// Returns an error (rather than panicking) if a variable's dimension product
/// overflows a `u32`; this is unreachable for meshes we encode in practice, but
/// surfacing it as [`ExodusError`] keeps the writer panic-free.
fn encode_nc3(dims: &[NcDim], vars: &[NcVar]) -> Result<Vec<u8>, ExodusError> {
    // First pass: header with placeholder begins (we recompute afterwards).
    let header = build_header(dims, vars, &vec![0u32; vars.len()])?;
    let mut data_off = header.len();
    data_off = (data_off + 3) & !3;
    let mut begins = Vec::with_capacity(vars.len());
    let mut cursor = data_off;
    for v in vars {
        begins.push(cursor as u32);
        let size = vsize_of(dims, v)? as usize;
        cursor += size;
        cursor = (cursor + 3) & !3;
    }
    let header = build_header(dims, vars, &begins)?;

    let mut out = header;
    // Pad header to 4-byte boundary (already aligned by construction of vars
    // sizes, but ensure).
    while out.len() % 4 != 0 {
        out.push(0);
    }
    for (v, _) in vars.iter().zip(&begins) {
        out.extend_from_slice(&v.data);
        while out.len() % 4 != 0 {
            out.push(0);
        }
    }
    Ok(out)
}

fn build_header(dims: &[NcDim], vars: &[NcVar], begins: &[u32]) -> Result<Vec<u8>, ExodusError> {
    let mut h = Vec::new();
    h.extend_from_slice(&[b'C', b'D', b'F', 1]);
    h.extend_from_slice(&to_u32(0)); // numrecs (no record variables)
    h.extend_from_slice(&to_u32(NC_DIMENSION));
    h.extend_from_slice(&to_u32(dims.len() as u32));
    for d in dims {
        h.extend_from_slice(&enc_name(&d.name));
        h.extend_from_slice(&to_u32(d.size));
    }
    h.extend_from_slice(&to_u32(NC_VARIABLE));
    h.extend_from_slice(&to_u32(vars.len() as u32));
    for (v, &begin) in vars.iter().zip(begins) {
        h.extend_from_slice(&enc_name(&v.name));
        h.extend_from_slice(&to_u32(v.dim_ids.len() as u32));
        for &d in &v.dim_ids {
            h.extend_from_slice(&to_u32(d));
        }
        h.extend_from_slice(&to_u32(0)); // number of attributes
        h.extend_from_slice(&to_u32(v.dtype));
        // `dims`/`vars` describe a mesh this crate just assembled, so `vsize_of`
        // is effectively unreachable; surface it as an error regardless.
        h.extend_from_slice(&to_u32(vsize_of(dims, v)?));
        h.extend_from_slice(&to_u32(begin));
    }
    h.extend_from_slice(&to_u32(0)); // NC_END
    Ok(h)
}

/// Decode dimensions and variables from NetCDF-3 classic bytes.
#[allow(unused_assignments, clippy::type_complexity)]
fn decode_nc3(bytes: &[u8]) -> Result<(Vec<NcDim>, Vec<(NcVar, usize)>), ExodusError> {
    let mut pos = 0usize;
    if bytes.len() < 4 || &bytes[0..4] != b"CDF\x01" {
        return Err(ExodusError::Parse("not a CDF-1 file".into()));
    }
    pos = 4;
    let _numrecs = read_u32(bytes, &mut pos)?;
    let tag = read_u32(bytes, &mut pos)?;
    if tag != NC_DIMENSION {
        return Err(ExodusError::Parse("expected dimension tag".into()));
    }
    let ndims = read_u32(bytes, &mut pos)? as usize;
    let ndims = sane_count(ndims, bytes)?;
    let mut dims = Vec::with_capacity(ndims);
    for _ in 0..ndims {
        let name = read_name(bytes, &mut pos)?;
        let size = read_u32(bytes, &mut pos)?;
        dims.push(NcDim { name, size });
    }
    let tag = read_u32(bytes, &mut pos)?;
    if tag != NC_VARIABLE {
        return Err(ExodusError::Parse("expected variable tag".into()));
    }
    let nvars = read_u32(bytes, &mut pos)? as usize;
    let nvars = sane_count(nvars, bytes)?;
    let mut vars = Vec::with_capacity(nvars);
    for _ in 0..nvars {
        let name = read_name(bytes, &mut pos)?;
        let nd = read_u32(bytes, &mut pos)? as usize;
        let nd = sane_count(nd, bytes)?;
        let mut dim_ids = Vec::with_capacity(nd);
        for _ in 0..nd {
            dim_ids.push(read_u32(bytes, &mut pos)?);
        }
        let _natts = read_u32(bytes, &mut pos)?;
        let dtype = read_u32(bytes, &mut pos)?;
        let _vsize = read_u32(bytes, &mut pos)? as usize; // declared vsize; recomputed below
        let vsize = usize::try_from(vsize_of(
            &dims,
            &NcVar {
                name: name.clone(),
                dtype,
                dim_ids: dim_ids.clone(),
                data: Vec::new(),
            },
        )?)
        .map_err(|_| ExodusError::Parse("variable size too large".into()))?;
        let begin = read_u32(bytes, &mut pos)? as usize;
        let data = bytes
            .get(begin..begin + vsize)
            .ok_or_else(|| ExodusError::Parse("variable data out of range".into()))?
            .to_vec();
        vars.push((
            NcVar {
                name,
                dtype,
                dim_ids,
                data,
            },
            begin,
        ));
    }
    Ok((dims, vars))
}

fn read_u32(bytes: &[u8], pos: &mut usize) -> Result<u32, ExodusError> {
    let s = bytes
        .get(*pos..*pos + 4)
        .ok_or_else(|| ExodusError::Parse("unexpected EOF".into()))?;
    *pos += 4;
    Ok(u32::from_be_bytes([s[0], s[1], s[2], s[3]]))
}

fn read_name(bytes: &[u8], pos: &mut usize) -> Result<String, ExodusError> {
    let len = *bytes
        .get(*pos)
        .ok_or_else(|| ExodusError::Parse("name length truncated".into()))? as usize;
    *pos += 1;
    let end = *pos + len;
    if end > bytes.len() {
        return Err(ExodusError::Parse("name truncated".into()));
    }
    let name = String::from_utf8_lossy(&bytes[*pos..end]).into_owned();
    *pos = end;
    while *pos < bytes.len() && *pos % 4 != 0 {
        *pos += 1;
    }
    Ok(name)
}

fn encode_floats(vals: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(vals.len() * 4);
    for v in vals {
        out.extend_from_slice(&(f32::to_bits(*v as f32)).to_be_bytes());
    }
    out
}

fn encode_ints(vals: &[i64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(vals.len() * 4);
    for v in vals {
        out.extend_from_slice(&(*v as u32).to_be_bytes());
    }
    out
}

fn encode_chars(cols: &[String], ncol: usize) -> Vec<u8> {
    // Each column is padded to LEN_NAME (33) chars.
    let mut out = Vec::with_capacity(cols.len() * LEN_NAME as usize);
    for c in cols {
        let mut s = c.clone();
        while s.len() < LEN_NAME as usize {
            s.push('\0');
        }
        s.truncate(LEN_NAME as usize);
        out.extend_from_slice(s.as_bytes());
    }
    let _ = ncol;
    out
}

fn decode_floats(data: &[u8]) -> Vec<f64> {
    data.chunks_exact(4)
        .map(|c| f32::from_bits(u32::from_be_bytes([c[0], c[1], c[2], c[3]])) as f64)
        .collect()
}

fn decode_ints(data: &[u8]) -> Vec<i64> {
    data.chunks_exact(4)
        .map(|c| u32::from_be_bytes([c[0], c[1], c[2], c[3]]) as i64)
        .collect()
}

fn decode_chars(data: &[u8], ncol: usize) -> Vec<String> {
    let mut out = Vec::with_capacity(ncol);
    for chunk in data.chunks(LEN_NAME as usize) {
        let s = String::from_utf8_lossy(chunk)
            .trim_end_matches('\0')
            .to_string();
        out.push(s);
    }
    out
}

// ---------------------------------------------------------------------------
// Exodus layer
// ---------------------------------------------------------------------------

fn cell_exodus_name(cell: CellType) -> &'static str {
    match cell {
        CellType::Line => "LINE2",
        CellType::Tri => "TRI3",
        CellType::Quad => "QUAD4",
        CellType::Tet => "TET4",
        CellType::Hex => "HEX8",
        CellType::Tri6 => "TRI6",
        CellType::Quad8 => "QUAD8",
        CellType::Quad9 => "QUAD9",
        CellType::Tet10 => "TET10",
        CellType::Hex20 => "HEX20",
        CellType::Hex27 => "HEX27",
    }
}

fn exodus_name_to_cell(name: &str) -> Option<CellType> {
    match name.trim().to_uppercase().as_str() {
        "LINE2" | "L2" | "T2D2" => Some(CellType::Line),
        "TRI3" | "T3" => Some(CellType::Tri),
        "TRI6" => Some(CellType::Tri6),
        "QUAD4" | "Q4" => Some(CellType::Quad),
        "QUAD8" => Some(CellType::Quad8),
        "QUAD9" => Some(CellType::Quad9),
        "TET4" | "T4" => Some(CellType::Tet),
        "TET10" => Some(CellType::Tet10),
        "HEX8" | "H8" => Some(CellType::Hex),
        "HEX20" => Some(CellType::Hex20),
        "HEX27" => Some(CellType::Hex27),
        _ => None,
    }
}

/// Write a `Mesh` to an Exodus II file at `path`.
pub fn write_exodus(mesh: &Mesh, path: impl AsRef<Path>) -> std::io::Result<()> {
    let bytes = mesh_to_exodus_bytes(mesh)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
    std::fs::write(path, bytes)
}

/// Serialise a `Mesh` to NetCDF-3 classic (Exodus II) bytes.
///
/// Returns an error instead of panicking if a variable's dimension product
/// overflows a `u32` (unreachable for the meshes this crate builds, but kept
/// panic-free for completeness).
pub fn mesh_to_exodus_bytes(mesh: &Mesh) -> Result<Vec<u8>, ExodusError> {
    // Group elements by cell type into blocks.
    let mut block_order: Vec<CellType> = Vec::new();
    let mut block_elems: Vec<Vec<usize>> = Vec::new();
    for (ei, e) in mesh.elements.iter().enumerate() {
        if let Some(pos) = block_order.iter().position(|c| *c == e.cell_type) {
            block_elems[pos].push(ei);
        } else {
            block_order.push(e.cell_type);
            block_elems.push(vec![ei]);
        }
    }
    let num_nodes = mesh.node_count() as u32;
    let num_elem = mesh.element_count() as u32;
    let num_blk = block_order.len() as u32;
    let num_dim: u32 = 3;

    let mut dims = vec![
        NcDim {
            name: "num_nodes".into(),
            size: num_nodes,
        },
        NcDim {
            name: "num_elem".into(),
            size: num_elem,
        },
        NcDim {
            name: "num_elem_blk".into(),
            size: num_blk,
        },
        NcDim {
            name: "num_node_sets".into(),
            size: 0,
        },
        NcDim {
            name: "num_side_sets".into(),
            size: 0,
        },
        NcDim {
            name: "len_name".into(),
            size: LEN_NAME,
        },
        NcDim {
            name: "four".into(),
            size: 4,
        },
        NcDim {
            name: "time_step".into(),
            size: 1,
        },
        NcDim {
            name: "num_dim".into(),
            size: num_dim,
        },
    ];
    // Per-block dimensions for the connectN arrays.
    let base = dims.len() as u32;
    let mut d0: Vec<u32> = Vec::new();
    let mut d1: Vec<u32> = Vec::new();
    for (bi, cell) in block_order.iter().enumerate() {
        let n_in_block = block_elems[bi].len() as u32;
        let np = cell.node_count() as u32;
        dims.push(NcDim {
            name: format!("eb{}_n", bi),
            size: n_in_block,
        });
        dims.push(NcDim {
            name: format!("eb{}_np", bi),
            size: np,
        });
        d0.push(base + 2 * bi as u32);
        d1.push(base + 2 * bi as u32 + 1);
    }
    // dim indices:
    let di_nodes = 0;
    let di_elem = 1;
    let di_blk = 2;
    let di_ln = 5;
    let di_ts = 7;
    let di_dim = 8;

    let mut vars: Vec<NcVar> = Vec::new();

    // time_whole
    vars.push(NcVar {
        name: "time_whole".into(),
        dtype: NC_FLOAT,
        dim_ids: vec![di_ts],
        data: encode_floats(&[0.0]),
    });

    // coords
    let mut coords = Vec::with_capacity(mesh.node_count() * 3);
    for n in &mesh.nodes {
        coords.push(n.coords.first().copied().unwrap_or(0.0));
        coords.push(n.coords.get(1).copied().unwrap_or(0.0));
        coords.push(n.coords.get(2).copied().unwrap_or(0.0));
    }
    vars.push(NcVar {
        name: "coords".into(),
        dtype: NC_FLOAT,
        dim_ids: vec![di_nodes, di_dim],
        data: encode_floats(&coords),
    });

    // elem_blk_id, eb_status, eb_prop1, eb_names
    let blk_ids: Vec<i64> = (1..=num_blk as i64).collect();
    vars.push(NcVar {
        name: "elem_blk_id".into(),
        dtype: NC_INT,
        dim_ids: vec![di_blk],
        data: encode_ints(&blk_ids),
    });
    vars.push(NcVar {
        name: "eb_status".into(),
        dtype: NC_INT,
        dim_ids: vec![di_blk],
        data: encode_ints(&vec![1i64; num_blk as usize]),
    });
    vars.push(NcVar {
        name: "eb_prop1".into(),
        dtype: NC_INT,
        dim_ids: vec![di_blk],
        data: encode_ints(&blk_ids),
    });
    let names: Vec<String> = block_order
        .iter()
        .map(|c| cell_exodus_name(*c).into())
        .collect();
    vars.push(NcVar {
        name: "eb_names".into(),
        dtype: NC_CHAR,
        dim_ids: vec![di_ln, di_blk],
        data: encode_chars(&names, num_blk as usize),
    });

    // connect1..k
    for (bi, _cell) in block_order.iter().enumerate() {
        let mut conn: Vec<i64> = Vec::new();
        for &ei in &block_elems[bi] {
            for &n in &mesh.elements[ei].nodes {
                conn.push((n + 1) as i64); // Exodus uses 1-based node ids
            }
        }
        vars.push(NcVar {
            name: format!("connect{}", bi + 1),
            dtype: NC_INT,
            dim_ids: vec![d0[bi], d1[bi]],
            data: encode_ints(&conn),
        });
    }

    // elem_num_map, node_num_map
    let elem_map: Vec<i64> = (1..=num_elem as i64).collect();
    vars.push(NcVar {
        name: "elem_num_map".into(),
        dtype: NC_INT,
        dim_ids: vec![di_elem],
        data: encode_ints(&elem_map),
    });
    let node_map: Vec<i64> = (1..=num_nodes as i64).collect();
    vars.push(NcVar {
        name: "node_num_map".into(),
        dtype: NC_INT,
        dim_ids: vec![di_nodes],
        data: encode_ints(&node_map),
    });

    encode_nc3(&dims, &vars)
}

/// Read an Exodus II mesh from `path`.
pub fn read_exodus(path: impl AsRef<Path>) -> Result<Mesh, ExodusError> {
    let bytes = std::fs::read(path).map_err(|e| ExodusError::Parse(format!("{e}")))?;
    bytes_to_mesh(&bytes)
}

/// Parse NetCDF-3 classic (Exodus II) bytes into a `Mesh`.
pub fn bytes_to_mesh(bytes: &[u8]) -> Result<Mesh, ExodusError> {
    let (dims, vars) = decode_nc3(bytes)?;
    let mut varmap: HashMap<String, (u32, Vec<u8>)> = HashMap::new();
    for (v, _) in &vars {
        varmap.insert(v.name.clone(), (v.dtype, v.data.clone()));
    }

    let coords = varmap
        .get("coords")
        .ok_or_else(|| ExodusError::Parse("missing coords".into()))?;
    if coords.0 != NC_FLOAT {
        return Err(ExodusError::Parse(
            "coords variable must be NC_FLOAT".into(),
        ));
    }
    let coords = decode_floats(&coords.1);
    if coords.len() % 3 != 0 {
        return Err(ExodusError::Parse(
            "coords length is not a multiple of 3".into(),
        ));
    }
    let num_nodes = coords.len() / 3;
    if num_nodes == 0 {
        return Err(ExodusError::Parse("coords contains no nodes".into()));
    }

    let eb_names = varmap
        .get("eb_names")
        .map(|(dt, d)| {
            if *dt != NC_CHAR {
                Err(ExodusError::Parse(
                    "eb_names variable must be NC_CHAR".into(),
                ))
            } else {
                Ok(decode_chars(d, d.len() / LEN_NAME as usize))
            }
        })
        .transpose()?
        .unwrap_or_default();

    let mut builder = MeshBuilder::new();
    for i in 0..num_nodes {
        let x = coords[i * 3];
        let y = coords[i * 3 + 1];
        let z = coords[i * 3 + 2];
        builder.add_node(vec![x, y, z]);
    }

    // For each connectN variable, build a block. Retain the on-disk `begin`
    // offset so tests can locate (and corrupt) the connectivity payload.
    let mut connect_vars: Vec<(usize, &NcVar, usize)> = vars
        .iter()
        .filter(|(v, _)| v.name.starts_with("connect"))
        .map(|(v, begin)| {
            let idx: String = v.name.trim_start_matches("connect").to_string();
            (idx.parse::<usize>().unwrap_or(0), v, *begin)
        })
        .collect();
    connect_vars.sort_by_key(|(i, _, _)| *i);

    for (_, v, _begin) in &connect_vars {
        if v.dtype != NC_INT {
            return Err(ExodusError::Parse(format!(
                "connect variable {} must be NC_INT",
                v.name
            )));
        }
        // dim_ids = [num_elem_in_block, nodes_per_elem] (dimension indices)
        let last_dim = *v
            .dim_ids
            .last()
            .ok_or_else(|| ExodusError::Parse(format!("connect {} has no dims", v.name)))?;
        let npe = dim_size(&dims, last_dim)? as usize;
        let cell = if !eb_names.is_empty() {
            // Determine block index from variable number.
            let idx: usize = v.name.trim_start_matches("connect").parse().unwrap_or(1);
            if idx == 0 {
                return Err(ExodusError::Parse(format!(
                    "connect variable {} has no numeric suffix",
                    v.name
                )));
            }
            let name = eb_names.get(idx - 1).map(|s| s.as_str()).unwrap_or("");
            exodus_name_to_cell(name)
                .ok_or_else(|| ExodusError::Parse(format!("unknown block {name}")))?
        } else {
            cell_from_npe(npe)
        };
        if cell.node_count() != npe {
            return Err(ExodusError::Parse(format!(
                "block {}: cell {:?} expects {} nodes, got {}",
                v.name,
                cell,
                cell.node_count(),
                npe
            )));
        }
        let raw = decode_ints(&v.data);
        if raw.len() % npe != 0 {
            return Err(ExodusError::Parse(format!(
                "connect {} has {} values, not a multiple of {npe}",
                v.name,
                raw.len()
            )));
        }
        let n_in_block = raw.len() / npe;
        for e in 0..n_in_block {
            // 1-based Exodus node ids -> 0-based, validating each against the
            // mesh's node count. A `0` entry (or any id out of range) corrupts
            // the connectivity and would otherwise be accepted silently.
            let mut nodes: Vec<usize> = Vec::with_capacity(npe);
            for k in 0..npe {
                let c = raw[e * npe + k];
                if c < 1 {
                    return Err(ExodusError::Mesh(
                        tpt_fem_mesh::MeshError::NodeIndexOutOfRange {
                            element: e,
                            node: c as usize,
                            node_count: num_nodes,
                        },
                    ));
                }
                let z = (c - 1) as usize;
                if z >= num_nodes {
                    return Err(ExodusError::Mesh(
                        tpt_fem_mesh::MeshError::NodeIndexOutOfRange {
                            element: e,
                            node: z,
                            node_count: num_nodes,
                        },
                    ));
                }
                nodes.push(z);
            }
            builder.try_add_element(cell, nodes)?;
        }
    }

    Ok(builder.try_build()?)
}

/// Infer a cell type from nodes-per-element when block names are absent.
fn cell_from_npe(npe: usize) -> CellType {
    match npe {
        2 => CellType::Line,
        3 => CellType::Tri,
        4 => CellType::Quad, // ambiguity resolved by eb_names when present
        6 => CellType::Tri6,
        8 => CellType::Hex, // ambiguity with Quad8 resolved by eb_names when present
        9 => CellType::Quad9,
        10 => CellType::Tet10,
        20 => CellType::Hex20,
        27 => CellType::Hex27,
        _ => CellType::Tri,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    #[test]
    fn round_trips_tri_mesh() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        let mesh = b.build();

        let bytes = mesh_to_exodus_bytes(&mesh).unwrap();
        assert_eq!(&bytes[0..4], b"CDF\x01");
        let parsed = bytes_to_mesh(&bytes).unwrap();
        assert_eq!(parsed.node_count(), 3);
        assert_eq!(parsed.element_count(), 1);
        assert_eq!(parsed.elements[0].cell_type, CellType::Tri);
        assert_eq!(parsed.elements[0].nodes, vec![0, 1, 2]);
        assert!((parsed.node_coords(1)[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn round_trips_multi_block() {
        // Two element types => two Exodus blocks.
        let mut b = MeshBuilder::new();
        let a = b.add_node(vec![0.0, 0.0]);
        let c = b.add_node(vec![1.0, 0.0]);
        let d = b.add_node(vec![0.0, 1.0]);
        let e = b.add_node(vec![1.0, 1.0]);
        let f = b.add_node(vec![2.0, 0.0]);
        b.add_element(CellType::Tri, vec![a, c, d]);
        b.add_element(CellType::Quad, vec![c, e, d, a]);
        b.add_element(CellType::Tri, vec![c, f, e]);
        let mesh = b.build();

        let bytes = mesh_to_exodus_bytes(&mesh).unwrap();
        let parsed = bytes_to_mesh(&bytes).unwrap();
        assert_eq!(parsed.element_count(), 3);
        let tri = parsed
            .elements
            .iter()
            .filter(|x| x.cell_type == CellType::Tri)
            .count();
        let quad = parsed
            .elements
            .iter()
            .filter(|x| x.cell_type == CellType::Quad)
            .count();
        assert_eq!(tri, 2);
        assert_eq!(quad, 1);
    }

    #[test]
    fn rejects_non_cdf1_magic() {
        let err = bytes_to_mesh(b"NOPE").unwrap_err();
        assert!(matches!(err, ExodusError::Parse(_)));
    }

    #[test]
    fn rejects_truncated_header() {
        let err = bytes_to_mesh(b"CDF\x01").unwrap_err();
        assert!(matches!(err, ExodusError::Parse(_)));
    }

    #[test]
    fn rejects_implausible_dimension_count() {
        let mut b: Vec<u8> = vec![b'C', b'D', b'F', 1];
        b.extend_from_slice(&0u32.to_be_bytes());
        b.extend_from_slice(&NC_DIMENSION.to_be_bytes());
        b.extend_from_slice(&u32::MAX.to_be_bytes());
        let err = bytes_to_mesh(&b).unwrap_err();
        assert!(matches!(err, ExodusError::Parse(_)));
    }

    #[test]
    fn rejects_truncated_name() {
        let mut b: Vec<u8> = vec![b'C', b'D', b'F', 1];
        b.extend_from_slice(&0u32.to_be_bytes());
        b.extend_from_slice(&NC_DIMENSION.to_be_bytes());
        b.extend_from_slice(&1u32.to_be_bytes());
        b.push(200u8);
        let err = bytes_to_mesh(&b).unwrap_err();
        assert!(matches!(err, ExodusError::Parse(_)));
    }

    #[test]
    fn rejects_connect_zero_node_id() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        let bytes = mesh_to_exodus_bytes(&b.build()).unwrap();
        // Locate the connect1 payload on disk and overwrite its first node id
        // (1-based) with 0, which underflows to usize::MAX when subtracted.
        let (_, vars) = decode_nc3(&bytes).unwrap();
        let begin = vars
            .iter()
            .find(|(v, _)| v.name == "connect1")
            .map(|(_, b)| *b)
            .expect("connect1 present");
        let mut corrupted = bytes.clone();
        corrupted[begin..begin + 4].copy_from_slice(&0u32.to_be_bytes());
        let err = bytes_to_mesh(&corrupted).unwrap_err();
        assert!(matches!(err, ExodusError::Mesh(_)));
    }

    #[test]
    fn rejects_connect_out_of_range_node_id() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        let bytes = mesh_to_exodus_bytes(&b.build()).unwrap();
        // Overwrite the third (last) node id of the only element with 9, which
        // maps to 0-based node 8 and is out of range for a 3-node mesh.
        let (_, vars) = decode_nc3(&bytes).unwrap();
        let begin = vars
            .iter()
            .find(|(v, _)| v.name == "connect1")
            .map(|(_, b)| *b)
            .expect("connect1 present");
        let mut corrupted = bytes.clone();
        corrupted[begin + 8..begin + 12].copy_from_slice(&9u32.to_be_bytes());
        let err = bytes_to_mesh(&corrupted).unwrap_err();
        assert!(matches!(err, ExodusError::Mesh(_)));
    }

    #[test]
    fn rejects_coords_wrong_dtype() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        let mut bytes = mesh_to_exodus_bytes(&b.build()).unwrap();
        let name = [6u8, b'c', b'o', b'o', b'r', b'd', b's', 0];
        if let Some(p) = bytes.windows(name.len()).position(|w| w == name) {
            let dtype_pos = p + 8 + 4 + 8 + 4;
            bytes[dtype_pos] = NC_INT as u8;
            let err = bytes_to_mesh(&bytes).unwrap_err();
            assert!(matches!(err, ExodusError::Parse(_)));
        }
    }
}
