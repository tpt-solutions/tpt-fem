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
        "CPS4" | "CPE4" | "S4" | "S4R" | "M3D4" | "R3D4" => Some(CellType::Quad),
        "C3D4" => Some(CellType::Tet),
        "C3D8" | "C3D8R" | "C3D8H" | "S8R" => Some(CellType::Hex),
        // Quadratic (P2) families.
        "CPS6" | "CPE6" | "STRI65" => Some(CellType::Tri6),
        "CPS8" | "CPE8" | "S8" => Some(CellType::Quad8),
        "CPS9" | "CPE9" => Some(CellType::Quad9),
        "C3D10" => Some(CellType::Tet10),
        "C3D20" | "C3D20R" | "C3D20H" => Some(CellType::Hex20),
        "C3D27" | "C3D27R" => Some(CellType::Hex27),
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
        CellType::Tri6 => "CPS6",
        CellType::Quad8 => "CPS8",
        CellType::Quad9 => "CPS9",
        CellType::Tet10 => "C3D10",
        CellType::Hex20 => "C3D20",
        CellType::Hex27 => "C3D27",
    }
}

/// A fully parsed Abaqus input deck: the geometry [`Mesh`] plus the named
/// node/element sets, material definitions, and prescribed-boundary entries
/// that the geometry-only [`read_inp`] would otherwise discard.
#[derive(Debug)]
pub struct InpDeck {
    /// The mesh assembled from `*NODE` / `*ELEMENT` sections.
    pub mesh: Mesh,
    /// `*NSET` definitions: set name -> internal node ids.
    pub nsets: HashMap<String, Vec<usize>>,
    /// `*ELSET` definitions: set name -> internal element ids.
    pub elsets: HashMap<String, Vec<usize>>,
    /// `*MATERIAL`/`*ELASTIC` definitions: material name -> `(young, poisson)`.
    pub materials: HashMap<String, (f64, f64)>,
    /// `*BOUNDARY` prescriptions: `(internal node id, dof, value)`.
    pub boundary: Vec<(usize, usize, f64)>,
}

/// Parse an Abaqus `.inp` text buffer into a geometry [`Mesh`].
///
/// Convenience wrapper around [`read_inp_deck`] that discards the named sets,
/// materials, and boundary conditions.
pub fn read_inp(text: &str) -> Result<Mesh, InpError> {
    Ok(read_inp_deck(text)?.mesh)
}

/// Parse an Abaqus `.inp` text buffer into an [`InpDeck`], retaining the
/// `*NSET` / `*ELSET` / `*MATERIAL` / `*BOUNDARY` sections (not just geometry)
/// so that real analysis decks can be imported, not only their topology.
pub fn read_inp_deck(text: &str) -> Result<InpDeck, InpError> {
    let mut builder = MeshBuilder::new();
    let mut id_map: HashMap<u64, usize> = HashMap::new();
    let mut elem_tag_map: HashMap<u64, usize> = HashMap::new();
    let mut current: Option<Section> = None;
    let mut element_type: Option<CellType> = None;
    let mut nsets: HashMap<String, Vec<usize>> = HashMap::new();
    let mut elsets: HashMap<String, Vec<usize>> = HashMap::new();
    let mut materials: HashMap<String, (f64, f64)> = HashMap::new();
    let mut boundary: Vec<(usize, usize, f64)> = Vec::new();
    let mut cur_nset: Option<String> = None;
    let mut cur_elset: Option<String> = None;
    let mut cur_material: Option<String> = None;

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with("**") {
            continue;
        }
        if line.starts_with('*') {
            // Section header.
            let header = line.trim_start_matches('*');
            let h = header.to_uppercase();
            if h.starts_with("NODE") {
                current = Some(Section::Node);
                element_type = None;
                cur_nset = None;
                cur_elset = None;
                cur_material = None;
            } else if h.starts_with("ELEMENT") {
                current = Some(Section::Element);
                element_type = parse_element_type(header);
                if element_type.is_none() {
                    return Err(InpError::UnknownElementType(header.to_string()));
                }
                cur_nset = None;
                cur_elset = None;
                cur_material = None;
            } else if h.starts_with("NSET") {
                current = Some(Section::Nset);
                cur_nset = set_name(header);
                cur_elset = None;
            } else if h.starts_with("ELSET") {
                current = Some(Section::Elset);
                cur_elset = set_name(header);
                cur_nset = None;
            } else if h.starts_with("MATERIAL") {
                current = Some(Section::Material);
                cur_material = set_name(header);
                cur_nset = None;
                cur_elset = None;
            } else if h.starts_with("ELASTIC") {
                // Associated with the current material; `(young, poisson)`
                // follows on the next line.
                current = Some(Section::Elastic);
                cur_nset = None;
                cur_elset = None;
            } else if h.starts_with("BOUNDARY") {
                current = Some(Section::Boundary);
                cur_nset = None;
                cur_elset = None;
            } else {
                current = None;
                element_type = None;
                cur_nset = None;
                cur_elset = None;
                cur_material = None;
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
                let aid = parts[0]
                    .parse::<u64>()
                    .map_err(|_| InpError::Parse(format!("element id: {line}")))?;
                let eid = builder.try_add_element(cell, nodes)?;
                elem_tag_map.insert(aid, eid);
            }
            Some(Section::Nset) => {
                let name = cur_nset
                    .clone()
                    .ok_or_else(|| InpError::Parse("NSET without a name".into()))?;
                for tok in line
                    .split([',', ' '])
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                {
                    let aid = tok
                        .parse::<u64>()
                        .map_err(|_| InpError::Parse(format!("nset id: {line}")))?;
                    let id = *id_map.get(&aid).ok_or_else(|| {
                        InpError::Parse(format!("nset references unknown node {aid}"))
                    })?;
                    nsets.entry(name.clone()).or_default().push(id);
                }
            }
            Some(Section::Elset) => {
                let name = cur_elset
                    .clone()
                    .ok_or_else(|| InpError::Parse("ELSET without a name".into()))?;
                for tok in line
                    .split([',', ' '])
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                {
                    let aid = tok
                        .parse::<u64>()
                        .map_err(|_| InpError::Parse(format!("elset id: {line}")))?;
                    let id = *elem_tag_map.get(&aid).ok_or_else(|| {
                        InpError::Parse(format!("elset references unknown element {aid}"))
                    })?;
                    elsets.entry(name.clone()).or_default().push(id);
                }
            }
            Some(Section::Material) => {
                // Name already captured via `cur_material`; no data on this line.
            }
            Some(Section::Elastic) => {
                // `young, poisson[, ...]` follows the `*ELASTIC` header.
                let parts: Vec<&str> = line
                    .split(',')
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                    .collect();
                if parts.len() < 2 {
                    return Err(InpError::Parse(format!("elastic record: {line}")));
                }
                let young = parts[0]
                    .parse::<f64>()
                    .map_err(|_| InpError::Parse(format!("young: {line}")))?;
                let nu = parts[1]
                    .parse::<f64>()
                    .map_err(|_| InpError::Parse(format!("poisson: {line}")))?;
                if let Some(m) = &cur_material {
                    materials.insert(m.clone(), (young, nu));
                }
            }
            Some(Section::Boundary) => {
                // `node, first_dof[, last_dof], value`. A single dof when the
                // middle field is absent; otherwise all dofs in `[first, last]`.
                let parts: Vec<&str> = line
                    .split(',')
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                    .collect();
                if parts.len() < 3 {
                    return Err(InpError::Parse(format!("boundary record: {line}")));
                }
                let aid = parts[0]
                    .parse::<u64>()
                    .map_err(|_| InpError::Parse(format!("boundary node: {line}")))?;
                let id = *id_map.get(&aid).ok_or_else(|| {
                    InpError::Parse(format!("boundary references unknown node {aid}"))
                })?;
                let first = parts[1]
                    .parse::<usize>()
                    .map_err(|_| InpError::Parse(format!("boundary dof: {line}")))?;
                let value = parts[parts.len() - 1]
                    .parse::<f64>()
                    .map_err(|_| InpError::Parse(format!("boundary value: {line}")))?;
                if parts.len() == 3 {
                    boundary.push((id, first, value));
                } else {
                    let last = parts[2]
                        .parse::<usize>()
                        .map_err(|_| InpError::Parse(format!("boundary dof: {line}")))?;
                    for d in first..=last {
                        boundary.push((id, d, value));
                    }
                }
            }
            None => {}
        }
    }

    let mesh = builder.try_build()?;
    Ok(InpDeck {
        mesh,
        nsets,
        elsets,
        materials,
        boundary,
    })
}

/// Extract the `NSET=`/`ELSET=`/`NAME=` value from a header token list.
fn set_name(header: &str) -> Option<String> {
    for tok in header.split(',') {
        let tok = tok.trim();
        let v = tok
            .strip_prefix("NSET=")
            .or_else(|| tok.strip_prefix("nset="))
            .or_else(|| tok.strip_prefix("ELSET="))
            .or_else(|| tok.strip_prefix("elset="))
            .or_else(|| tok.strip_prefix("NAME="))
            .or_else(|| tok.strip_prefix("name="));
        if let Some(name) = v {
            return Some(name.trim().to_string());
        }
    }
    None
}

#[derive(Clone, Copy)]
enum Section {
    Node,
    Element,
    Nset,
    Elset,
    Material,
    Elastic,
    Boundary,
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
    fn s4r_shell_quad_maps_to_quad() {
        // S4R is a 4-node shell quadrilateral, not a tetrahedron.
        let text = "*NODE\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n3, 1.0, 1.0, 0.0\n4, 0.0, 1.0, 0.0\n*ELEMENT, TYPE=S4R\n1, 1, 2, 3, 4\n";
        let mesh = read_inp(text).unwrap();
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].cell_type, CellType::Quad);
    }

    #[test]
    fn ignores_other_sections() {
        let text = "*HEADING\nSome title\n*NODE\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n*ELEMENT, TYPE=T2D2\n1, 1, 2\n*MATERIAL, NAME=STEEL\n*ELASTIC\n200.0, 0.3\n";
        let mesh = read_inp(text).unwrap();
        assert_eq!(mesh.node_count(), 2);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].cell_type, CellType::Line);
    }

    #[test]
    fn deck_captures_sets_material_boundary() {
        // A small analysis deck: geometry plus *NSET / *ELSET / *MATERIAL /
        // *ELASTIC / *BOUNDARY. The deck parser must retain all of them, not
        // just the topology.
        let text = "\
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 0.0, 1.0, 0.0
4, 1.0, 1.0, 0.0
*ELEMENT, TYPE=CPS4
1, 1, 2, 4, 3
*NSET, NSET=fixed
1, 3
*ELSET, ELSET=plate
1
*MATERIAL, NAME=STEEL
*ELASTIC
200.0, 0.3
*BOUNDARY
1, 1, 0.0
1, 2, 0.0
";
        let deck = read_inp_deck(text).expect("parse deck");
        assert_eq!(deck.mesh.node_count(), 4);
        assert_eq!(deck.mesh.element_count(), 1);
        assert_eq!(deck.nsets.get("fixed").map(|v| v.len()), Some(2));
        assert_eq!(deck.elsets.get("plate").map(|v| v.len()), Some(1));
        assert_eq!(deck.materials.get("STEEL"), Some(&(200.0, 0.3)));
        // Node 1 (internal id 0) pinned in dofs 1 and 2.
        assert_eq!(deck.boundary.len(), 2);
        assert!(deck.boundary.contains(&(0, 1, 0.0)));
        assert!(deck.boundary.contains(&(0, 2, 0.0)));
    }
}
