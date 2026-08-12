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

fn dim_size(dims: &[NcDim], id: u32) -> u32 {
    dims[id as usize].size
}

fn vsize_of(dims: &[NcDim], v: &NcVar) -> u32 {
    let mut s: u32 = match v.dtype {
        NC_CHAR => 1,
        NC_INT | NC_FLOAT => 4,
        _ => 4,
    };
    for &d in &v.dim_ids {
        s *= dim_size(dims, d);
    }
    s
}

/// Encode a file from dimensions and variables into NetCDF-3 classic bytes.
fn encode_nc3(dims: &[NcDim], vars: &[NcVar]) -> Vec<u8> {
    // First pass: header with placeholder begins (we recompute afterwards).
    let header = build_header(dims, vars, &vec![0u32; vars.len()]);
    let mut data_off = header.len();
    data_off = (data_off + 3) & !3;
    let mut begins = Vec::with_capacity(vars.len());
    let mut cursor = data_off;
    for v in vars {
        begins.push(cursor as u32);
        let size = vsize_of(dims, v) as usize;
        cursor += size;
        cursor = (cursor + 3) & !3;
    }
    let header = build_header(dims, vars, &begins);

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
    out
}

fn build_header(dims: &[NcDim], vars: &[NcVar], begins: &[u32]) -> Vec<u8> {
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
        h.extend_from_slice(&to_u32(vsize_of(dims, v)));
        h.extend_from_slice(&to_u32(begin));
    }
    h.extend_from_slice(&to_u32(0)); // NC_END
    h
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
    let mut vars = Vec::with_capacity(nvars);
    for _ in 0..nvars {
        let name = read_name(bytes, &mut pos)?;
        let nd = read_u32(bytes, &mut pos)? as usize;
        let mut dim_ids = Vec::with_capacity(nd);
        for _ in 0..nd {
            dim_ids.push(read_u32(bytes, &mut pos)?);
        }
        let _natts = read_u32(bytes, &mut pos)?;
        let dtype = read_u32(bytes, &mut pos)?;
        let vsize = read_u32(bytes, &mut pos)? as usize;
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
        .ok_or_else(|| ExodusError::Parse("name length".into()))? as usize;
    *pos += 1;
    let end = *pos + len;
    let name = String::from_utf8_lossy(&bytes[*pos..end]).into_owned();
    *pos = end;
    while *pos % 4 != 0 {
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
    }
}

fn exodus_name_to_cell(name: &str) -> Option<CellType> {
    match name.trim().to_uppercase().as_str() {
        "LINE2" | "L2" | "T2D2" => Some(CellType::Line),
        "TRI3" | "T3" => Some(CellType::Tri),
        "QUAD4" | "Q4" => Some(CellType::Quad),
        "TET4" | "T4" => Some(CellType::Tet),
        "HEX8" | "H8" => Some(CellType::Hex),
        _ => None,
    }
}

/// Write a `Mesh` to an Exodus II file at `path`.
pub fn write_exodus(mesh: &Mesh, path: impl AsRef<Path>) -> std::io::Result<()> {
    let bytes = mesh_to_exodus_bytes(mesh);
    std::fs::write(path, bytes)
}

/// Serialise a `Mesh` to NetCDF-3 classic (Exodus II) bytes.
pub fn mesh_to_exodus_bytes(mesh: &Mesh) -> Vec<u8> {
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
    let coords = decode_floats(&coords.1);
    let num_nodes = coords.len() / 3;

    let eb_names = varmap
        .get("eb_names")
        .map(|(_, d)| decode_chars(d, d.len() / LEN_NAME as usize))
        .unwrap_or_default();

    let mut builder = MeshBuilder::new();
    for i in 0..num_nodes {
        let x = coords[i * 3];
        let y = coords[i * 3 + 1];
        let z = coords[i * 3 + 2];
        builder.add_node(vec![x, y, z]);
    }

    // For each connectN variable, build a block.
    let mut connect_vars: Vec<(usize, &NcVar)> = vars
        .iter()
        .filter(|(v, _)| v.name.starts_with("connect"))
        .map(|(v, _)| {
            let idx: String = v.name.trim_start_matches("connect").to_string();
            (idx.parse::<usize>().unwrap_or(0), v)
        })
        .collect();
    connect_vars.sort_by_key(|(i, _)| *i);

    for (_, v) in &connect_vars {
        // dim_ids = [num_elem_in_block, nodes_per_elem] (dimension indices)
        let npe = dims[*v.dim_ids.last().unwrap() as usize].size as usize;
        let cell = if !eb_names.is_empty() {
            // Determine block index from variable number.
            let idx: usize = v.name.trim_start_matches("connect").parse().unwrap_or(1);
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
        let n_in_block = raw.len() / npe;
        for e in 0..n_in_block {
            let nodes: Vec<usize> = (0..npe)
                .map(|k| (raw[e * npe + k] - 1) as usize) // 1-based -> 0-based
                .collect();
            builder.add_element(cell, nodes);
        }
    }

    Ok(builder.build())
}

/// Infer a cell type from nodes-per-element when block names are absent.
fn cell_from_npe(npe: usize) -> CellType {
    match npe {
        2 => CellType::Line,
        3 => CellType::Tri,
        4 => CellType::Quad, // ambiguity resolved by eb_names when present
        8 => CellType::Hex,
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

        let bytes = mesh_to_exodus_bytes(&mesh);
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

        let bytes = mesh_to_exodus_bytes(&mesh);
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
}
