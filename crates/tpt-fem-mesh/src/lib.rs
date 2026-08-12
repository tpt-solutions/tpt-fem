//! Mesh data structures, DOF numbering, and Gmsh `.msh` v4.1 import.
//!
//! This crate provides the mesh containers used by the rest of `tpt-fem`:
//!
//! * [`Node`] / [`Element`] / [`Mesh`] — the in-memory mesh model,
//! * [`MeshBuilder`] — a manual mesh-builder API for constructing meshes in code,
//! * [`Mesh::number_dofs`] — configurable (dofs-per-node) degree-of-freedom
//!   numbering,
//! * [`Mesh::from_msh_bytes`] — import of Gmsh `.msh` version 4.1 ASCII files
//!   via the [`mshio`] crate.
//!
//! Only the five linear element types `Line2`, `Tri3`, `Quad4`, `Tet4`, and
//! `Hex8` are imported; other Gmsh element types produce
//! [`MeshError::UnsupportedElementType`]. Higher-order element types are a
//! tracked follow-up.
//!
//! # Example
//!
//! ```
//! use tpt_fem_mesh::{MeshBuilder, CellType};
//!
//! let mut b = MeshBuilder::new();
//! let n0 = b.add_node(vec![0.0, 0.0]);
//! let n1 = b.add_node(vec![1.0, 0.0]);
//! let n2 = b.add_node(vec![0.0, 1.0]);
//! b.add_element(CellType::Tri, vec![n0, n1, n2]);
//! let mesh = b.build();
//! assert_eq!(mesh.node_count(), 3);
//! assert_eq!(mesh.element_count(), 1);
//! ```

use std::collections::HashMap;

/// A node identifier (index into [`Mesh::nodes`]).
pub type NodeId = usize;
/// An element identifier (index into [`Mesh::elements`]).
pub type ElementId = usize;

/// The supported linear cell types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CellType {
    /// 2-node line.
    Line,
    /// 3-node triangle.
    Tri,
    /// 4-node quadrilateral.
    Quad,
    /// 4-node tetrahedron.
    Tet,
    /// 8-node hexahedron.
    Hex,
}

impl CellType {
    /// Number of nodes for this cell type.
    pub fn node_count(self) -> usize {
        match self {
            CellType::Line => 2,
            CellType::Tri => 3,
            CellType::Quad => 4,
            CellType::Tet => 4,
            CellType::Hex => 8,
        }
    }

    /// Reference-element name, used in diagnostics.
    pub fn name(self) -> &'static str {
        match self {
            CellType::Line => "Line2",
            CellType::Tri => "Tri3",
            CellType::Quad => "Quad4",
            CellType::Tet => "Tet4",
            CellType::Hex => "Hex8",
        }
    }
}

/// A mesh node: an identifier and its spatial coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    /// Node identifier.
    pub id: NodeId,
    /// Coordinates (length 3; for 2-D meshes `z` is `0.0`).
    pub coords: Vec<f64>,
    /// Gmsh physical-group tag of the geometric entity this node was
    /// imported from, if any (see [`Mesh::from_msh_bytes`]). `None` for
    /// manually built meshes and meshes imported from formats that don't
    /// carry region tags.
    pub region: Option<i32>,
}

/// A mesh element: an identifier, a cell type, and its node list.
#[derive(Clone, Debug, PartialEq)]
pub struct Element {
    /// Element identifier.
    pub id: ElementId,
    /// Cell type.
    pub cell_type: CellType,
    /// Node identifiers (in the order required by the reference element).
    pub nodes: Vec<NodeId>,
    /// Gmsh physical-group tag of the geometric entity this element was
    /// imported from, if any. See [`Node::region`].
    pub region: Option<i32>,
}

/// An in-memory finite-element mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    /// Nodes, indexed by [`NodeId`].
    pub nodes: Vec<Node>,
    /// Elements, indexed by [`ElementId`].
    pub elements: Vec<Element>,
}

/// Degree-of-freedom numbering for a mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct DofMap {
    /// For each node, the list of DOF indices it owns.
    pub node_dofs: Vec<Vec<usize>>,
    /// Total number of degrees of freedom.
    pub ndof: usize,
}

/// Errors produced while building or importing a mesh.
#[derive(Debug)]
pub enum MeshError {
    /// The Gmsh file could not be parsed.
    Parse(String),
    /// A Gmsh element type that this crate does not support was encountered.
    UnsupportedElementType(u8),
    /// The file contained no `$Nodes` section.
    MissingNodes,
    /// The file contained no `$Elements` section.
    MissingElements,
    /// An `$Elements` entry referenced a node tag absent from `$Nodes`.
    DanglingNodeTag(u64),
    /// An element's node count does not match its cell type.
    NodeCountMismatch {
        /// The element's cell type.
        cell: CellType,
        /// The number of nodes `cell` requires.
        expected: usize,
        /// The number of nodes actually supplied.
        got: usize,
    },
    /// An element referenced a node index beyond the mesh's node count.
    NodeIndexOutOfRange {
        /// The offending element's id.
        element: ElementId,
        /// The out-of-range node id it referenced.
        node: NodeId,
        /// The mesh's actual node count.
        node_count: usize,
    },
}

impl std::fmt::Display for MeshError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MeshError::Parse(m) => write!(f, "failed to parse mesh: {m}"),
            MeshError::UnsupportedElementType(t) => {
                write!(f, "unsupported element type (Gmsh type id {t})")
            }
            MeshError::MissingNodes => write!(f, "mesh file has no $Nodes section"),
            MeshError::MissingElements => write!(f, "mesh file has no $Elements section"),
            MeshError::DanglingNodeTag(t) => write!(
                f,
                "element references node tag {t}, which is not present in the node section"
            ),
            MeshError::NodeCountMismatch {
                cell,
                expected,
                got,
            } => write!(f, "{} expects {expected} nodes, got {got}", cell.name()),
            MeshError::NodeIndexOutOfRange {
                element,
                node,
                node_count,
            } => write!(
                f,
                "element {element} references node {node}, but the mesh has only {node_count} nodes"
            ),
        }
    }
}

impl std::error::Error for MeshError {}

impl Mesh {
    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }

    /// Coordinates of node `id`.
    pub fn node_coords(&self, id: NodeId) -> &[f64] {
        &self.nodes[id].coords
    }

    /// Assign degrees of freedom, `dofs_per_node` per node, numbered
    /// contiguously by node id.
    pub fn number_dofs(&self, dofs_per_node: usize) -> DofMap {
        let node_dofs = (0..self.nodes.len())
            .map(|n| (0..dofs_per_node).map(|k| n * dofs_per_node + k).collect())
            .collect();
        DofMap {
            node_dofs,
            ndof: self.nodes.len() * dofs_per_node,
        }
    }

    /// Parse a Gmsh `.msh` version 4.1 (ASCII) byte buffer into a mesh.
    ///
    /// Only the five linear cell types are supported; see [`MeshError`]. If
    /// the file has an `$Entities` section, each node/element's
    /// [`Node::region`]/[`Element::region`] is set to the first Gmsh
    /// physical-group tag of the geometric entity it belongs to (`None` if
    /// the entity has no physical tag, or if the file has no `$Entities`
    /// section at all).
    pub fn from_msh_bytes(bytes: &[u8]) -> Result<Mesh, MeshError> {
        let parsed =
            mshio::parse_msh_bytes(bytes).map_err(|e| MeshError::Parse(format!("{e:?}")))?;
        let data = &parsed.data;
        let nodes_section = data.nodes.as_ref().ok_or(MeshError::MissingNodes)?;
        let elements_section = data.elements.as_ref().ok_or(MeshError::MissingElements)?;

        // Map (entity_dim, entity_tag) -> first physical-group tag, from the
        // optional $Entities section (points = dim 0, curves = dim 1,
        // surfaces = dim 2, volumes = dim 3).
        let mut region_map: HashMap<(i32, i32), i32> = HashMap::new();
        if let Some(entities) = &data.entities {
            for p in &entities.points {
                if let Some(&tag) = p.physical_tags.first() {
                    region_map.insert((0, p.tag), tag);
                }
            }
            for c in &entities.curves {
                if let Some(&tag) = c.physical_tags.first() {
                    region_map.insert((1, c.tag), tag);
                }
            }
            for s in &entities.surfaces {
                if let Some(&tag) = s.physical_tags.first() {
                    region_map.insert((2, s.tag), tag);
                }
            }
            for v in &entities.volumes {
                if let Some(&tag) = v.physical_tags.first() {
                    region_map.insert((3, v.tag), tag);
                }
            }
        }

        // Build a global node-tag -> NodeId map.
        let mut nodes: Vec<Node> = Vec::new();
        let mut tag_to_id: HashMap<u64, NodeId> = HashMap::new();
        let min_tag = nodes_section.min_node_tag;
        for block in &nodes_section.node_blocks {
            let block_start = nodes.len();
            let region = region_map.get(&(block.entity_dim, block.entity_tag)).copied();
            match &block.node_tags {
                // Sparse tags: the map gives tag -> local index within the block.
                Some(map) => {
                    for (&tag, &local_idx) in map {
                        let node = &block.nodes[local_idx];
                        let id = block_start + local_idx;
                        tag_to_id.insert(tag, id);
                        nodes.push(Node {
                            id,
                            coords: vec![node.x, node.y, node.z],
                            region,
                        });
                    }
                }
                // Sequential tags: node at global index `i` has tag `min_tag + i`.
                None => {
                    for (local_idx, node) in block.nodes.iter().enumerate() {
                        let id = block_start + local_idx;
                        let tag = min_tag + id as u64;
                        tag_to_id.insert(tag, id);
                        nodes.push(Node {
                            id,
                            coords: vec![node.x, node.y, node.z],
                            region,
                        });
                    }
                }
            }
        }

        // Build elements.
        let mut elements: Vec<Element> = Vec::new();
        for block in &elements_section.element_blocks {
            let cell = match block.element_type {
                mshio::mshfile::ElementType::Lin2 => CellType::Line,
                mshio::mshfile::ElementType::Tri3 => CellType::Tri,
                mshio::mshfile::ElementType::Qua4 => CellType::Quad,
                mshio::mshfile::ElementType::Tet4 => CellType::Tet,
                mshio::mshfile::ElementType::Hex8 => CellType::Hex,
                other => return Err(MeshError::UnsupportedElementType(other as u8)),
            };
            let expected = cell.node_count();
            let region = region_map.get(&(block.entity_dim, block.entity_tag)).copied();
            for element in &block.elements {
                if element.nodes.len() != expected {
                    return Err(MeshError::UnsupportedElementType(block.element_type as u8));
                }
                let node_ids: Vec<NodeId> = element
                    .nodes
                    .iter()
                    .map(|t| tag_to_id.get(t).copied().ok_or(MeshError::DanglingNodeTag(*t)))
                    .collect::<Result<_, _>>()?;
                let id = elements.len();
                elements.push(Element {
                    id,
                    cell_type: cell,
                    nodes: node_ids,
                    region,
                });
            }
        }

        Ok(Mesh { nodes, elements })
    }

    /// Validate that every element's node count matches its cell type and
    /// every node id it references is within range.
    ///
    /// [`MeshBuilder::build`]/[`MeshBuilder::add_element`] already guarantee
    /// the first invariant by construction (they panic instead); this method
    /// is for meshes assembled from untrusted external data (e.g. a file
    /// importer) where a `Result` is preferable to a panic. See
    /// [`MeshBuilder::try_build`].
    pub fn validate(&self) -> Result<(), MeshError> {
        let node_count = self.nodes.len();
        for e in &self.elements {
            let expected = e.cell_type.node_count();
            if e.nodes.len() != expected {
                return Err(MeshError::NodeCountMismatch {
                    cell: e.cell_type,
                    expected,
                    got: e.nodes.len(),
                });
            }
            for &n in &e.nodes {
                if n >= node_count {
                    return Err(MeshError::NodeIndexOutOfRange {
                        element: e.id,
                        node: n,
                        node_count,
                    });
                }
            }
        }
        Ok(())
    }

    /// Node ids whose coordinate along `axis` (`0` = x, `1` = y, `2` = z) is
    /// within `tol` of `coord`. Useful for selecting a boundary conditions
    /// on an axis-aligned face of a box-like mesh.
    pub fn nodes_on_plane(&self, axis: usize, coord: f64, tol: f64) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| (n.coords.get(axis).copied().unwrap_or(0.0) - coord).abs() <= tol)
            .map(|n| n.id)
            .collect()
    }

    /// Node ids within the axis-aligned box `[min, max]` (inclusive, with a
    /// small tolerance to absorb floating-point rounding at the boundary).
    pub fn nodes_in_box(&self, min: [f64; 3], max: [f64; 3]) -> Vec<NodeId> {
        const TOL: f64 = 1e-12;
        self.nodes
            .iter()
            .filter(|n| {
                (0..3).all(|d| {
                    let c = n.coords.get(d).copied().unwrap_or(0.0);
                    c >= min[d] - TOL && c <= max[d] + TOL
                })
            })
            .map(|n| n.id)
            .collect()
    }
}

/// A builder for constructing meshes manually in code.
#[derive(Clone, Debug, Default)]
pub struct MeshBuilder {
    nodes: Vec<Node>,
    elements: Vec<Element>,
}

impl MeshBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        MeshBuilder::default()
    }

    /// Add a node with the given coordinates; returns its [`NodeId`].
    pub fn add_node(&mut self, coords: Vec<f64>) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(Node {
            id,
            coords,
            region: None,
        });
        id
    }

    /// Add an element of the given cell type and node list; returns its
    /// [`ElementId`].
    ///
    /// # Panics
    ///
    /// Panics if `nodes.len()` does not match `cell.node_count()`.
    pub fn add_element(&mut self, cell: CellType, nodes: Vec<NodeId>) -> ElementId {
        assert_eq!(
            nodes.len(),
            cell.node_count(),
            "{} expects {} nodes, got {}",
            cell.name(),
            cell.node_count(),
            nodes.len()
        );
        let id = self.elements.len();
        self.elements.push(Element {
            id,
            cell_type: cell,
            nodes,
            region: None,
        });
        id
    }

    /// Like [`MeshBuilder::add_element`], but tagged with a region (Gmsh
    /// physical-group-style) tag instead of leaving it `None`.
    pub fn add_element_with_region(
        &mut self,
        cell: CellType,
        nodes: Vec<NodeId>,
        region: i32,
    ) -> ElementId {
        let id = self.add_element(cell, nodes);
        self.elements[id].region = Some(region);
        id
    }

    /// Like [`MeshBuilder::add_element`], but returns a [`MeshError`]
    /// instead of panicking if `nodes.len()` does not match
    /// `cell.node_count()`. Prefer this over `add_element` when the node
    /// list originates from untrusted external data (e.g. a file importer).
    pub fn try_add_element(
        &mut self,
        cell: CellType,
        nodes: Vec<NodeId>,
    ) -> Result<ElementId, MeshError> {
        let expected = cell.node_count();
        if nodes.len() != expected {
            return Err(MeshError::NodeCountMismatch {
                cell,
                expected,
                got: nodes.len(),
            });
        }
        let id = self.elements.len();
        self.elements.push(Element {
            id,
            cell_type: cell,
            nodes,
            region: None,
        });
        Ok(id)
    }

    /// Finalize into a [`Mesh`].
    pub fn build(self) -> Mesh {
        Mesh {
            nodes: self.nodes,
            elements: self.elements,
        }
    }

    /// Finalize into a [`Mesh`], validating (via [`Mesh::validate`]) that
    /// every element's node ids are in range. Prefer this over `build` when
    /// node ids originate from untrusted external data.
    pub fn try_build(self) -> Result<Mesh, MeshError> {
        let mesh = Mesh {
            nodes: self.nodes,
            elements: self.elements,
        };
        mesh.validate()?;
        Ok(mesh)
    }
}

impl DofMap {
    /// The DOF index of `component` (0-based) at `node`.
    pub fn dof(&self, node: NodeId, component: usize) -> usize {
        self.node_dofs[node][component]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builder_and_dofs() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        let e = b.add_element(CellType::Tri, vec![n0, n1, n2]);
        let mesh = b.build();

        assert_eq!(mesh.node_count(), 3);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[e].cell_type, CellType::Tri);
        assert_eq!(mesh.node_coords(n1), &[1.0, 0.0]);

        let dof1 = mesh.number_dofs(1);
        assert_eq!(dof1.ndof, 3);
        assert_eq!(dof1.dof(n2, 0), 2);

        let dof3 = mesh.number_dofs(3);
        assert_eq!(dof3.ndof, 9);
        assert_eq!(dof3.dof(n1, 0), 3);
        assert_eq!(dof3.dof(n1, 2), 5);
    }

    #[test]
    fn import_gmsh_tri_mesh() {
        let msh = "\
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
        let mesh = Mesh::from_msh_bytes(msh.as_bytes()).expect("parse");
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 2);
        // Node with Gmsh tag 1 is our NodeId 0 at the origin.
        assert_eq!(mesh.node_coords(0), &[0.0, 0.0, 0.0]);
        assert_eq!(mesh.elements[0].cell_type, CellType::Tri);
        assert_eq!(mesh.elements[0].nodes, vec![0, 1, 2]);
        assert_eq!(mesh.elements[1].nodes, vec![1, 3, 2]);

        let dofs = mesh.number_dofs(2);
        assert_eq!(dofs.ndof, 8);
    }

    #[test]
    fn import_rejects_unsupported() {
        // A mesh containing a single 3-node line (Lin3, type 8) is unsupported.
        let msh = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
1 1 0 3
1 2 3
0.0 0.0 0.0
1.0 0.0 0.0
2.0 0.0 0.0
$EndNodes
$Elements
1 1 1 1
1 1 8 1
1 1 2 3
$EndElements
";
        let err = Mesh::from_msh_bytes(msh.as_bytes()).unwrap_err();
        assert!(matches!(err, MeshError::UnsupportedElementType(8)));
    }

    #[test]
    fn import_rejects_dangling_node_tag() {
        // $Elements references node tag 99, which $Nodes never defines.
        let msh = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
1 1 0 3
1 2 3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2 99
$EndElements
";
        let err = Mesh::from_msh_bytes(msh.as_bytes()).unwrap_err();
        assert!(matches!(err, MeshError::DanglingNodeTag(99)));
    }

    #[test]
    fn import_reads_physical_group_regions() {
        // A single triangle whose surface entity (dim=2, tag=1) carries
        // physical-group tag 42.
        let msh = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 1 0
1 0 0 0 1 1 0 1 42 0
$EndEntities
$Nodes
1 3 1 3
2 1 0 3
1 2 3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2 3
$EndElements
";
        let mesh = Mesh::from_msh_bytes(msh.as_bytes()).expect("parse");
        assert_eq!(mesh.nodes[0].region, Some(42));
        assert_eq!(mesh.elements[0].region, Some(42));
    }

    #[test]
    fn validate_rejects_out_of_range_node() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        // Manually corrupt connectivity to reference a nonexistent node 7.
        let mut mesh = b.build();
        mesh.elements[0].nodes[2] = 7;
        let err = mesh.validate().unwrap_err();
        assert!(matches!(
            err,
            MeshError::NodeIndexOutOfRange {
                node: 7,
                node_count: 3,
                ..
            }
        ));
    }

    #[test]
    fn try_add_element_rejects_wrong_count() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let err = b.try_add_element(CellType::Tri, vec![n0, n1]).unwrap_err();
        assert!(matches!(
            err,
            MeshError::NodeCountMismatch {
                cell: CellType::Tri,
                expected: 3,
                got: 2,
            }
        ));
    }

    #[test]
    fn try_build_accepts_valid_mesh() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.try_add_element(CellType::Tri, vec![n0, n1, n2]).unwrap();
        let mesh = b.try_build().expect("valid mesh");
        assert_eq!(mesh.element_count(), 1);
    }

    #[test]
    fn selectors_find_nodes_on_plane_and_in_box() {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0, 0.0]);
        let n3 = b.add_node(vec![1.0, 1.0, 1.0]);
        let mesh = b.build();

        let mut on_x0 = mesh.nodes_on_plane(0, 0.0, 1e-9);
        on_x0.sort();
        assert_eq!(on_x0, vec![n0, n2]);

        let mut in_box = mesh.nodes_in_box([0.0, 0.0, 0.0], [1.0, 1.0, 0.5]);
        in_box.sort();
        assert_eq!(in_box, vec![n0, n1, n2]);
        assert!(!in_box.contains(&n3));
    }

    #[test]
    fn error_display_is_human_readable() {
        let e = MeshError::DanglingNodeTag(7);
        assert!(e.to_string().contains("node tag 7"));
        let e = MeshError::NodeIndexOutOfRange {
            element: 2,
            node: 9,
            node_count: 4,
        };
        assert!(e.to_string().contains("element 2"));
    }
}
