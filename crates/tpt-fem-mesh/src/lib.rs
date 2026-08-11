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
}

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
    /// Only the five linear cell types are supported; see [`MeshError`].
    pub fn from_msh_bytes(bytes: &[u8]) -> Result<Mesh, MeshError> {
        let parsed =
            mshio::parse_msh_bytes(bytes).map_err(|e| MeshError::Parse(format!("{e:?}")))?;
        let data = &parsed.data;
        let nodes_section = data.nodes.as_ref().ok_or(MeshError::MissingNodes)?;
        let elements_section = data.elements.as_ref().ok_or(MeshError::MissingElements)?;

        // Build a global node-tag -> NodeId map.
        let mut nodes: Vec<Node> = Vec::new();
        let mut tag_to_id: HashMap<u64, NodeId> = HashMap::new();
        let min_tag = nodes_section.min_node_tag;
        for block in &nodes_section.node_blocks {
            let block_start = nodes.len();
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
            for element in &block.elements {
                if element.nodes.len() != expected {
                    return Err(MeshError::UnsupportedElementType(block.element_type as u8));
                }
                let node_ids = element
                    .nodes
                    .iter()
                    .map(|t| *tag_to_id.get(t).expect("node tag present in mesh"))
                    .collect();
                let id = elements.len();
                elements.push(Element {
                    id,
                    cell_type: cell,
                    nodes: node_ids,
                });
            }
        }

        Ok(Mesh { nodes, elements })
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
        self.nodes.push(Node { id, coords });
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
        });
        id
    }

    /// Finalize into a [`Mesh`].
    pub fn build(self) -> Mesh {
        Mesh {
            nodes: self.nodes,
            elements: self.elements,
        }
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
}
