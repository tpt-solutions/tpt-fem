//! Abaqus `.inp` mesh import/export for `tpt-fem`.
//!
//! A minimal but practical subset of the Abaqus input deck is supported:
//! `*NODE` and `*ELEMENT, TYPE=…` sections, mapping the common linear-element
//! types (`T2D2`, `CPS3`/`CPS4`, `C3D4`, `C3D8`) onto the `tpt-fem-mesh`
//! `CellType`s. All other sections are ignored.

use std::collections::HashMap;
use std::path::Path;

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};

/// Errors returned while parsing an Abaqus input deck.
#[derive(Debug)]
pub enum InpError {
    /// A `*NODE` or `*ELEMENT` record was malformed.
    Parse(String),
    /// An unknown element type was encountered.
    UnknownElementType(String),
    /// The mesh built from the deck failed validation (e.g. a node index out
    /// of range).
    Mesh(tpt_fem_mesh::MeshError),
    /// An I/O error occurred while reading or writing the file.
    Io(std::io::Error),
}

impl std::fmt::Display for InpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InpError::Parse(m) => write!(f, "malformed Abaqus .inp record: {m}"),
            InpError::UnknownElementType(t) => write!(f, "unknown Abaqus element type: {t}"),
            InpError::Mesh(e) => write!(f, "invalid mesh in Abaqus .inp file: {e}"),
            InpError::Io(e) => write!(f, "I/O error reading/writing Abaqus .inp file: {e}"),
        }
    }
}

impl std::error::Error for InpError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            InpError::Mesh(e) => Some(e),
            InpError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<tpt_fem_mesh::MeshError> for InpError {
    fn from(e: tpt_fem_mesh::MeshError) -> Self {
        InpError::Mesh(e)
    }
}

impl From<std::io::Error> for InpError {
    fn from(e: std::io::Error) -> Self {
        InpError::Io(e)
    }
}

/// Map an Abaqus element type name to a `tpt-fem` cell type.
fn abaqus_type_to_cell(ty: &str) -> Option<CellType> {
    match ty.trim().to_uppercase().as_str() {
        "T2D2" | "C1D2" | "B31" => Some(CellType::Line),
        "CPS3" | "CPE3" | "S3" | "STRI3" | "M3D3" | "R3D3" => Some(CellType::Tri),
        "CPS4" | "CPE4" | "S4" | "M3D4" | "R3D4" => Some(CellType::Quad),
        "C3D4" | "S4R" => Some(CellType::Tet),
        "C3D8" | "C3D8R" | "C3D8H" | "S8R" => Some(CellType::Hex),
        _ => None,
    }
}

/// Map a `tpt-fem` cell type to its Abaqus element type name (writer side).
fn cell_to_abaqus_type(cell: CellType) -> &'static str {
    match cell {
        CellType::Line => "T2D2",
        CellType::Tri => "CPS3",
        CellType::Quad => "CPS4",
        CellType::Tet => "C3D4",
        CellType::Hex => "C3D8",
    }
}

/// Parse an Abaqus `.inp` text buffer into a `Mesh`.
pub fn read_inp(text: &str) -> Result<Mesh, InpError> {
    let mut builder = MeshBuilder::new();
    let mut id_map: HashMap<u64, usize> = HashMap::new();
    let mut current: Option<Section> = None;
    let mut element_type: Option<CellType> = None;

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with("**") {
            continue;
        }
        if line.starts_with('*') {
            // Section header.
            let header = line.trim_start_matches('*');
            if header.to_uppercase().starts_with("NODE") {
                current = Some(Section::Node);
                element_type = None;
            } else if header.to_uppercase().starts_with("ELEMENT") {
                current = Some(Section::Element);
                element_type = parse_element_type(header);
                if element_type.is_none() {
                    return Err(InpError::UnknownElementType(header.to_string()));
                }
            } else {
                current = None;
                element_type = None;
            }
            continue;
        }

        match current {
            Some(Section::Node) => {
                let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
                if parts.len() < 2 {
                    return Err(InpError::Parse(format!("node record: {line}")));
                }
                let aid = parts[0]
                    .parse::<u64>()
                    .map_err(|_| InpError::Parse(format!("node id: {line}")))?;
                let mut coords = vec![0.0; 3];
                for (i, c) in parts.iter().skip(1).take(3).enumerate() {
                    coords[i] = c
                        .parse::<f64>()
                        .map_err(|_| InpError::Parse(format!("node coord: {line}")))?;
                }
                let id = builder.add_node(coords);
                id_map.insert(aid, id);
            }
            Some(Section::Element) => {
                let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
                if parts.len() < 2 {
                    return Err(InpError::Parse(format!("element record: {line}")));
                }
                let cell = element_type.unwrap();
                let nodes: Vec<usize> = parts
                    .iter()
                    .skip(1)
                    .map(|p| {
                        p.parse::<u64>()
                            .ok()
                            .and_then(|aid| id_map.get(&aid).copied())
                            .ok_or_else(|| InpError::Parse(format!("element node: {line}")))
                    })
                    .collect::<Result<_, _>>()?;
                if nodes.len() != cell.node_count() {
                    return Err(InpError::Parse(format!(
                        "element node count {} != expected {} for {}",
                        nodes.len(),
                        cell.node_count(),
                        line
                    )));
                }
                builder.try_add_element(cell, nodes)?;
            }
            None => {}
        }
    }

    Ok(builder.try_build()?)
}

#[derive(Clone, Copy)]
enum Section {
    Node,
    Element,
}

/// Extract the `TYPE=xxx` from an `*ELEMENT` header.
fn parse_element_type(header: &str) -> Option<CellType> {
    for tok in header.split(',') {
        let tok = tok.trim();
        if let Some(ty) = tok
            .strip_prefix("TYPE=")
            .or_else(|| tok.strip_prefix("type="))
        {
            return abaqus_type_to_cell(ty);
        }
    }
    None
}

/// Write a `Mesh` to an Abaqus `.inp` file (nodes + a single element section).
pub fn write_inp(mesh: &Mesh, path: impl AsRef<Path>) -> Result<(), InpError> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "*NODE")?;
    for n in &mesh.nodes {
        writeln!(
            f,
            "{}, {}, {}, {}",
            n.id + 1,
            n.coords[0],
            n.coords.get(1).copied().unwrap_or(0.0),
            n.coords.get(2).copied().unwrap_or(0.0),
        )?;
    }
    // Group elements by cell type; emit one *ELEMENT section per type.
    let mut order: Vec<CellType> = Vec::new();
    for e in &mesh.elements {
        if !order.contains(&e.cell_type) {
            order.push(e.cell_type);
        }
    }
    for cell in order {
        writeln!(f, "*ELEMENT, TYPE={}", cell_to_abaqus_type(cell))?;
        for e in &mesh.elements {
            if e.cell_type != cell {
                continue;
            }
            let mut line = format!("{}", e.id + 1);
            for &n in &e.nodes {
                line.push_str(&format!(", {}", n + 1));
            }
            writeln!(f, "{line}")?;
        }
    }
    Ok(())
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

        // Reuse write_inp by writing to a temp path then reading the string.
        let dir = std::env::temp_dir();
        let path = dir.join("tpt_fem_io_abaqus_test.inp");
        write_inp(&mesh, &path).unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        let parsed = read_inp(&s).unwrap();
        assert_eq!(parsed.node_count(), 3);
        assert_eq!(parsed.element_count(), 1);
        assert_eq!(parsed.elements[0].cell_type, CellType::Tri);
        // Connectivity should map back to the same relative ordering.
        assert_eq!(parsed.elements[0].nodes, vec![0, 1, 2]);
        // Coordinates preserved.
        assert!((parsed.node_coords(0)[0] - 0.0).abs() < 1e-12);
        assert!((parsed.node_coords(2)[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn parses_unknown_type_as_error() {
        let text = "*ELEMENT, TYPE=UNKNOWN\n1, 1, 2, 3\n";
        assert!(matches!(
            read_inp(text),
            Err(InpError::UnknownElementType(_))
        ));
    }

    #[test]
    fn ignores_other_sections() {
        let text = "*HEADING\nSome title\n*NODE\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n*ELEMENT, TYPE=T2D2\n1, 1, 2\n*MATERIAL, NAME=STEEL\n*ELASTIC\n200.0, 0.3\n";
        let mesh = read_inp(text).unwrap();
        assert_eq!(mesh.node_count(), 2);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].cell_type, CellType::Line);
    }
}
