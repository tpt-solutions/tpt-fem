//! `tpt-fem` — umbrella crate for the `tpt-fem` core.
//!
//! This crate re-exports the `tpt-fem-*` crates behind Cargo features, so
//! downstream users can depend on a single crate and opt into only the pieces
//! they need:
//!
//! | Feature       | Re-exports                       | Backing crate          |
//! |---------------|---------------------------------|-----------------------|
//! | `quadrature`  | `tpt_fem_quadrature::*`         | `tpt-fem-quadrature`  |
//! | `element`     | `tpt_fem_element::*`            | `tpt-fem-element`     |
//! | `mesh`        | `tpt_fem_mesh::*`               | `tpt-fem-mesh`        |
//! | `sparse`      | `tpt_fem_sparse::*`             | `tpt-fem-sparse`      |
//! | `assembly`    | `tpt_fem_assembly::*`           | `tpt-fem-assembly`    |
//! | `thermal`     | `tpt_fem_thermal::*`            | `tpt-fem-thermal`     |
//! | `io-vtk`      | `tpt_fem_io_vtk::*`             | `tpt-fem-io-vtk`      |
//! | `elasticity`  | `tpt_fem_elasticity::*`         | `tpt-fem-elasticity`  |
//! | `solve`       | `tpt_fem_solve::*`              | `tpt-fem-solve`       |
//! | `eigen`       | `tpt_fem_eigen::*`              | `tpt-fem-eigen`       |
//! | `io-abaqus`   | `tpt_fem_io_abaqus::*`          | `tpt-fem-io-abaqus`   |
//! | `io-exodus`   | `tpt_fem_io_exodus::*`          | `tpt-fem-io-exodus`   |
//! | `mesh-gen`    | `tpt_fem_mesh_gen::*`           | `tpt-fem-mesh-gen`    |
//!
//! All features are enabled by default. Each constituent's public API is
//! documented in its own crate; see the module re-exports below.
//!
//! The end-to-end pipeline — reference-element shape functions and gradients,
//! quadrature, triplet assembly, and a sparse solve — is exercised by the
//! integration test `tests/patch_test.rs`.
#![cfg_attr(
    not(feature = "quadrature"),
    doc = "The `quadrature` feature is disabled; `tpt_fem_quadrature` is not re-exported."
)]
#![cfg_attr(
    not(feature = "element"),
    doc = "The `element` feature is disabled; `tpt_fem_element` is not re-exported."
)]
#![cfg_attr(
    not(feature = "mesh"),
    doc = "The `mesh` feature is disabled; `tpt_fem_mesh` is not re-exported."
)]
#![cfg_attr(
    not(feature = "sparse"),
    doc = "The `sparse` feature is disabled; `tpt_fem_sparse` is not re-exported."
)]
#![cfg_attr(
    not(feature = "assembly"),
    doc = "The `assembly` feature is disabled; `tpt_fem_assembly` is not re-exported."
)]
#![cfg_attr(
    not(feature = "thermal"),
    doc = "The `thermal` feature is disabled; `tpt_fem_thermal` is not re-exported."
)]
#![cfg_attr(
    not(feature = "io-vtk"),
    doc = "The `io-vtk` feature is disabled; `tpt_fem_io_vtk` is not re-exported."
)]
#![cfg_attr(
    not(feature = "elasticity"),
    doc = "The `elasticity` feature is disabled; `tpt_fem_elasticity` is not re-exported."
)]
#![cfg_attr(
    not(feature = "solve"),
    doc = "The `solve` feature is disabled; `tpt_fem_solve` is not re-exported."
)]
#![cfg_attr(
    not(feature = "eigen"),
    doc = "The `eigen` feature is disabled; `tpt_fem_eigen` is not re-exported."
)]
#![cfg_attr(
    not(feature = "io-abaqus"),
    doc = "The `io-abaqus` feature is disabled; `tpt_fem_io_abaqus` is not re-exported."
)]
#![cfg_attr(
    not(feature = "io-exodus"),
    doc = "The `io-exodus` feature is disabled; `tpt_fem_io_exodus` is not re-exported."
)]
#![cfg_attr(
    not(feature = "mesh-gen"),
    doc = "The `mesh-gen` feature is disabled; `tpt_fem_mesh_gen` is not re-exported."
)]

#[cfg(feature = "assembly")]
pub use tpt_fem_assembly::*;
#[cfg(feature = "eigen")]
pub use tpt_fem_eigen::*;
#[cfg(feature = "elasticity")]
pub use tpt_fem_elasticity::*;
#[cfg(feature = "element")]
pub use tpt_fem_element::*;
#[cfg(feature = "io-abaqus")]
pub use tpt_fem_io_abaqus::*;
#[cfg(feature = "io-exodus")]
pub use tpt_fem_io_exodus::*;
#[cfg(feature = "io-vtk")]
pub use tpt_fem_io_vtk::*;
#[cfg(feature = "mesh")]
pub use tpt_fem_mesh::*;
#[cfg(feature = "mesh-gen")]
pub use tpt_fem_mesh_gen::*;
#[cfg(feature = "quadrature")]
pub use tpt_fem_quadrature::*;
#[cfg(feature = "solve")]
pub use tpt_fem_solve::*;
#[cfg(feature = "sparse")]
pub use tpt_fem_sparse::*;
#[cfg(feature = "thermal")]
pub use tpt_fem_thermal::*;

/// A curated, feature-gated prelude.
///
/// `use tpt_fem::prelude::*;` pulls in the public API of every constituent
/// crate that is enabled by its Cargo feature, so a typical application can be
/// written against a single import. Like the crate-root re-exports, each
/// constituent is gated: only crates whose feature is on appear in the prelude.
pub mod prelude {
    #[cfg(feature = "assembly")]
    pub use tpt_fem_assembly::*;
    #[cfg(feature = "eigen")]
    pub use tpt_fem_eigen::*;
    #[cfg(feature = "elasticity")]
    pub use tpt_fem_elasticity::*;
    #[cfg(feature = "element")]
    pub use tpt_fem_element::*;
    #[cfg(feature = "io-abaqus")]
    pub use tpt_fem_io_abaqus::*;
    #[cfg(feature = "io-exodus")]
    pub use tpt_fem_io_exodus::*;
    #[cfg(feature = "io-vtk")]
    pub use tpt_fem_io_vtk::*;
    #[cfg(feature = "mesh")]
    pub use tpt_fem_mesh::*;
    #[cfg(feature = "mesh-gen")]
    pub use tpt_fem_mesh_gen::*;
    #[cfg(feature = "quadrature")]
    pub use tpt_fem_quadrature::*;
    #[cfg(feature = "solve")]
    pub use tpt_fem_solve::*;
    #[cfg(feature = "sparse")]
    pub use tpt_fem_sparse::*;
    #[cfg(feature = "thermal")]
    pub use tpt_fem_thermal::*;
}
