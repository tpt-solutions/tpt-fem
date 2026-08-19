//! Multi-field degree-of-freedom numbering for `tpt-fem`.
//!
//! The rest of the core (`tpt-fem-mesh::Mesh::number_dofs`) hard-codes a single
//! *uniform* `dofs_per_node` across the whole mesh — fine for a single physical
//! field (e.g. temperature, or a single displacement vector) but insufficient
//! for coupled problems where different fields carry different component counts
//! (e.g. displacement with 3 components and pressure with 1). [`MultiFieldDofMap`]
//! generalises that to N named fields, each with its own per-node component
//! count, numbered on one [`Mesh`].
//!
//! ```
//! use tpt_fem_dofmap::{FieldSpec, Layout, MultiFieldDofMap};
//! use tpt_fem_mesh::{CellType, MeshBuilder};
//!
//! let mut b = MeshBuilder::new();
//! let a = b.add_node(vec![0.0, 0.0]);
//! let c = b.add_node(vec![1.0, 0.0]);
//! b.add_element(CellType::Line, vec![a, c]);
//! let mesh = b.build();
//!
//! // One displacement field (2 components) + one pressure field (1 component).
//! let map = MultiFieldDofMap::new(
//!     &mesh,
//!     &[FieldSpec::new("u", 2), FieldSpec::new("p", 1)],
//!     Layout::Interleaved,
//! );
//! // Node 0 owns dofs 0,1 (u) and 2 (p); node 1 owns 3,4 and 5.
//! assert_eq!(map.node_field_dof(0, 0, 1), 1);
//! assert_eq!(map.node_field_dof(0, 1, 0), 2);
//! assert_eq!(map.node_field_dof(1, 1, 0), 5);
//! ```

use tpt_fem_mesh::Mesh;

/// A named physical field and its per-node component count.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldSpec {
    /// Field name (e.g. `"u"` for displacement, `"p"` for pressure).
    pub name: String,
    /// Number of scalar components carried per node by this field.
    pub components: usize,
}

impl FieldSpec {
    /// Construct a field spec. `components` must be non-zero.
    pub fn new(name: impl Into<String>, components: usize) -> Self {
        assert!(components > 0, "FieldSpec::components must be non-zero");
        FieldSpec {
            name: name.into(),
            components,
        }
    }
}

/// How the per-field DOFs are laid out in the global numbering.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// Field-major: all DOFs of field 0 for every node come first, then field 1,
    /// and so on. `node n`'s field-`f` dofs are
    /// `[offset(f) + n * components(f), …]`.
    Block,
    /// Node-major: for each node, the component blocks of every field are
    /// contiguous. `node n`'s field-`f` component `c` dof is
    /// `[n * total_components + prefix(f) + c]`.
    Interleaved,
}

/// A multi-field degree-of-freedom map over a [`Mesh`].
///
/// Indexing helpers translate `(node, field, component)` triples into the single
/// global DOF index used by the assembled system matrices. When a map has a
/// single field it degenerates to exactly `Mesh::number_dofs`'s numbering.
#[derive(Clone, Debug)]
pub struct MultiFieldDofMap {
    /// The field descriptors, in map order.
    pub fields: Vec<FieldSpec>,
    /// The layout used to assign global DOFs.
    pub layout: Layout,
    /// Per-node list of *all* DOFs owned by that node, in global order.
    pub node_dofs: Vec<Vec<usize>>,
    /// For each field, the `(start, count)` of its contiguous global DOF block.
    /// Meaningful for both layouts (it is the block for that field, though for
    /// `Interleaved` the field is not globally contiguous — see [`Self::field_range`]
    /// for the block-major view used by `Block`).
    pub field_ranges: Vec<(usize, usize)>,
    /// Total number of degrees of freedom.
    pub ndof: usize,
}

impl MultiFieldDofMap {
    /// Build the map for `mesh` with the given fields and layout.
    pub fn new(mesh: &Mesh, fields: &[FieldSpec], layout: Layout) -> Self {
        let nnodes = mesh.node_count();
        let total_components: usize = fields.iter().map(|f| f.components).sum();
        let mut prefixes = vec![0usize; fields.len()];
        let mut running = 0;
        for (i, f) in fields.iter().enumerate() {
            prefixes[i] = running;
            running += f.components;
        }
        let ndof = nnodes * total_components;

        let mut node_dofs = vec![Vec::new(); nnodes];
        for n in 0..nnodes {
            let mut dofs = Vec::with_capacity(total_components);
            match layout {
                Layout::Block => {
                    for (fi, f) in fields.iter().enumerate() {
                        let base = prefixes[fi] * nnodes + n * f.components;
                        for c in 0..f.components {
                            dofs.push(base + c);
                        }
                    }
                }
                Layout::Interleaved => {
                    for c in 0..total_components {
                        dofs.push(n * total_components + c);
                    }
                }
            }
            node_dofs[n] = dofs;
        }

        let field_ranges = fields
            .iter()
            .enumerate()
            .map(|(fi, f)| {
                let start = prefixes[fi] * nnodes;
                (start, f.components * nnodes)
            })
            .collect();

        MultiFieldDofMap {
            fields: fields.to_vec(),
            layout,
            node_dofs,
            field_ranges,
            ndof,
        }
    }

    /// Global DOF for `(node, field_index, component)`.
    pub fn node_field_dof(&self, node: usize, field: usize, component: usize) -> usize {
        match self.layout {
            Layout::Block => {
                let prefix: usize = self.fields[..field].iter().map(|f| f.components).sum();
                prefix * self.node_dofs.len() + node * self.fields[field].components + component
            }
            Layout::Interleaved => {
                let total: usize = self.fields.iter().map(|f| f.components).sum();
                let prefix: usize = self.fields[..field].iter().map(|f| f.components).sum();
                node * total + prefix + component
            }
        }
    }

    /// Number of components of a given field.
    pub fn components(&self, field: usize) -> usize {
        self.fields[field].components
    }

    /// `(start, count)` of the contiguous global DOF block owned by `field`
    /// under [`Layout::Block`] numbering.
    pub fn field_range(&self, field: usize) -> (usize, usize) {
        self.field_ranges[field]
    }

    /// All global DOFs owned by `node`.
    pub fn dofs_of(&self, node: usize) -> &[usize] {
        &self.node_dofs[node]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};

    fn two_node_line() -> Mesh {
        let mut b = MeshBuilder::new();
        let a = b.add_node(vec![0.0]);
        let c = b.add_node(vec![1.0]);
        b.add_element(CellType::Line, vec![a, c]);
        b.build()
    }

    #[test]
    fn two_fields_interleaved_matches_hand_computed() {
        let mesh = two_node_line();
        // field 0: 3 components, field 1: 1 component. Node-major layout.
        let map = MultiFieldDofMap::new(
            &mesh,
            &[FieldSpec::new("a", 3), FieldSpec::new("b", 1)],
            Layout::Interleaved,
        );
        assert_eq!(map.ndof, 2 * 4);
        // Node 0: a0 a1 a2 b0 -> 0 1 2 3
        assert_eq!(map.node_field_dof(0, 0, 0), 0);
        assert_eq!(map.node_field_dof(0, 0, 2), 2);
        assert_eq!(map.node_field_dof(0, 1, 0), 3);
        // Node 1: shift by total_components (4)
        assert_eq!(map.node_field_dof(1, 0, 0), 4);
        assert_eq!(map.node_field_dof(1, 1, 0), 7);
        assert_eq!(map.node_dofs[0], vec![0, 1, 2, 3]);
        assert_eq!(map.node_dofs[1], vec![4, 5, 6, 7]);
    }

    #[test]
    fn two_fields_block_matches_hand_computed() {
        let mesh = two_node_line();
        let map = MultiFieldDofMap::new(
            &mesh,
            &[FieldSpec::new("a", 3), FieldSpec::new("b", 1)],
            Layout::Block,
        );
        // Field 0 (3 comp) occupies dofs [0,3) for node0 and [3,6) for node1;
        // field 1 (1 comp) occupies [6,7) (node0) and [7,8) (node1).
        assert_eq!(map.node_field_dof(0, 0, 0), 0);
        assert_eq!(map.node_field_dof(1, 0, 0), 3);
        assert_eq!(map.node_field_dof(0, 1, 0), 6);
        assert_eq!(map.node_field_dof(1, 1, 0), 7);
        assert_eq!(map.field_range(0), (0, 6));
        assert_eq!(map.field_range(1), (6, 2));
    }

    #[test]
    fn single_field_degenerates_to_number_dofs() {
        let mesh = two_node_line();
        let single = MultiFieldDofMap::new(&mesh, &[FieldSpec::new("u", 2)], Layout::Block);
        let legacy = mesh.number_dofs(2);
        for n in 0..mesh.node_count() {
            assert_eq!(single.node_dofs[n], legacy.node_dofs[n]);
        }
        assert_eq!(single.ndof, legacy.ndof);
    }
}
