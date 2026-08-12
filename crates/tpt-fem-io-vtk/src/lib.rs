//! ParaView-compatible result export for `tpt-fem`.
//!
//! Wraps [`vtkio`] to write the linear element types of `tpt-fem-mesh`
//! (`Line2`, `Tri3`, `Quad4`, `Tet4`, `Hex8`) as an unstructured-grid `.vtk`
//! (legacy) or `.vtu` (XML) file, with optional per-node scalar fields such as
//! a computed temperature or displacement magnitude.
//!
//! # Example
//!
//! ```
//! use tpt_fem_io_vtk::PointData;
//! use tpt_fem_mesh::{CellType, MeshBuilder};
//!
//! let mut b = MeshBuilder::new();
//! let n0 = b.add_node(vec![0.0, 0.0]);
//! let n1 = b.add_node(vec![1.0, 0.0]);
//! let n2 = b.add_node(vec![0.0, 1.0]);
//! b.add_element(CellType::Tri, vec![n0, n1, n2]);
//! let mesh = b.build();
//!
//! let vtk = tpt_fem_io_vtk::mesh_to_vtk(&mesh, &[PointData::new("u", vec![0.0, 1.0, 1.0])]);
//! if let vtkio::model::DataSet::UnstructuredGrid { pieces, .. } = &vtk.data {
//!     if let vtkio::model::Piece::Inline(p) = &pieces[0] {
//!         assert_eq!(p.num_points(), 3);
//!     }
//! }
//! ```

use std::path::Path;

use tpt_fem_mesh::{CellType, Mesh};
use vtkio::model::{
    Attribute, Attributes, ByteOrder, Cells, DataArray, DataSet, ElementType, IOBuffer,
    UnstructuredGridPiece, Version, VertexNumbers, Vtk,
};

/// A named per-node scalar field to embed in the output file.
pub struct PointData {
    /// Field name (as shown in ParaView).
    pub name: String,
    /// One value per mesh node.
    pub values: Vec<f64>,
}

impl PointData {
    /// Create a new point-data field.
    pub fn new(name: impl Into<String>, values: Vec<f64>) -> Self {
        PointData {
            name: name.into(),
            values,
        }
    }
}

fn vtk_cell_type(cell: CellType) -> vtkio::model::CellType {
    match cell {
        CellType::Line => vtkio::model::CellType::Line,
        CellType::Tri => vtkio::model::CellType::Triangle,
        CellType::Quad => vtkio::model::CellType::Quad,
        CellType::Tet => vtkio::model::CellType::Tetra,
        CellType::Hex => vtkio::model::CellType::Hexahedron,
    }
}

/// Build a [`Vtk`] unstructured-grid model from a `tpt-fem` mesh and any
/// per-node scalar fields.
pub fn mesh_to_vtk(mesh: &Mesh, point_data: &[PointData]) -> Vtk {
    let mut points = Vec::with_capacity(mesh.node_count() * 3);
    for n in &mesh.nodes {
        for d in 0..3 {
            points.push(*n.coords.get(d).unwrap_or(&0.0));
        }
    }

    let mut vertices: Vec<u32> = Vec::new();
    let mut types = Vec::new();
    for e in &mesh.elements {
        vertices.push(e.nodes.len() as u32);
        for &nd in &e.nodes {
            vertices.push(nd as u32);
        }
        types.push(vtk_cell_type(e.cell_type));
    }

    let point_attrs: Vec<Attribute> = point_data
        .iter()
        .map(|pd| {
            Attribute::DataArray(DataArray {
                name: pd.name.clone(),
                elem: ElementType::Scalars {
                    num_comp: 1,
                    lookup_table: None,
                },
                data: IOBuffer::from(pd.values.clone()),
            })
        })
        .collect();

    Vtk {
        version: Version::new((4, 1)),
        byte_order: ByteOrder::BigEndian,
        title: String::from("tpt-fem mesh"),
        file_path: None,
        data: DataSet::inline(UnstructuredGridPiece {
            points: points.into(),
            cells: Cells {
                cell_verts: VertexNumbers::Legacy {
                    num_cells: mesh.elements.len() as u32,
                    vertices,
                },
                types,
            },
            data: Attributes {
                point: point_attrs,
                ..Default::default()
            },
        }),
    }
}

/// Write the mesh (no scalar fields) as a binary legacy `.vtk` file.
pub fn write_vtk(mesh: &Mesh, path: impl AsRef<Path>) -> Result<(), vtkio::Error> {
    mesh_to_vtk(mesh, &[]).export(path)
}

/// Write the mesh (no scalar fields) as an ASCII legacy `.vtk` file.
pub fn write_vtk_ascii(mesh: &Mesh, path: impl AsRef<Path>) -> Result<(), vtkio::Error> {
    mesh_to_vtk(mesh, &[]).export_ascii(path)
}

/// Write the mesh together with per-node scalar fields as a `.vtk` file.
pub fn write_vtk_with_data(
    mesh: &Mesh,
    point_data: &[PointData],
    path: impl AsRef<Path>,
) -> Result<(), vtkio::Error> {
    mesh_to_vtk(mesh, point_data).export(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fem_mesh::{CellType, MeshBuilder};
    use vtkio::model::{DataSet, Piece, Vtk};

    fn tri_mesh() -> Mesh {
        let mut b = MeshBuilder::new();
        let n0 = b.add_node(vec![0.0, 0.0]);
        let n1 = b.add_node(vec![1.0, 0.0]);
        let n2 = b.add_node(vec![0.0, 1.0]);
        b.add_element(CellType::Tri, vec![n0, n1, n2]);
        b.build()
    }

    #[test]
    fn builds_correct_grid() {
        let mesh = tri_mesh();
        let vtk = mesh_to_vtk(&mesh, &[PointData::new("u", vec![0.0, 0.5, 1.0])]);
        if let DataSet::UnstructuredGrid { pieces, .. } = &vtk.data {
            if let Piece::Inline(p) = &pieces[0] {
                assert_eq!(p.num_points(), 3);
                assert_eq!(p.cells.num_cells(), 1);
            }
        } else {
            panic!("expected unstructured grid");
        }
    }

    #[test]
    fn round_trips_through_file() {
        let mesh = tri_mesh();
        let dir = std::env::temp_dir();
        let path = dir.join("tpt_fem_io_vtk_test.vtk");
        write_vtk(&mesh, &path).expect("export");
        let imported = Vtk::import(&path).expect("import");
        if let DataSet::UnstructuredGrid { pieces, .. } = &imported.data {
            if let Piece::Inline(p) = &pieces[0] {
                assert_eq!(p.num_points(), 3);
            }
        } else {
            panic!("expected unstructured grid");
        }
        let _ = std::fs::remove_file(&path);
    }
}
