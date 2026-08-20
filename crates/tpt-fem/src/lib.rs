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
//! | `dofmap`     | `tpt_fem_dofmap::*`            | `tpt-fem-dofmap`     |
//! | `dynamic`    | `tpt_fem_dynamic::*`           | `tpt-fem-dynamic`     |
//! | `plasticity` | `tpt_fem_plasticity::*`        | `tpt-fem-plasticity`  |
//! | `hyperelastic` | `tpt_fem_hyperelastic::*`     | `tpt-fem-hyperelastic` |
//! | `composite`  | `tpt_fem_composite::*`         | `tpt-fem-composite`   |
//! | `porous`     | `tpt_fem_porous::*`            | `tpt-fem-porous`      |
//! | `contact`    | `tpt_fem_contact::*`           | `tpt-fem-contact`     |
//! | `fluid`      | `tpt_fem_fluid::*`             | `tpt-fem-fluid`       |
//! | `coupling`   | `tpt_fem_coupling::*`          | `tpt-fem-coupling`    |
//!
//! All Phase 1–4 features are enabled by default. The Phase 12 crates
//! (`dofmap`, `dynamic`, `plasticity`, `hyperelastic`, `composite`, `porous`,
//! `contact`, `fluid`, `coupling`) are experimental and **off by default** —
//! enable them explicitly with `--features`.
//!
//! The end-to-end pipeline — reference-element shape functions and gradients,
//! quadrature, triplet assembly, and a sparse solve — is exercised by the
//! integration test `tests/patch_test.rs`.
//!
//! # Error handling
//!
//! This crate also provides a unified [`Error`] type that aggregates the
//! per-crate error enums (plus `std::io::Error` and TOML deserialization
//! errors) behind `#[from]` conversions, so consumers such as the CLI can use a
//! single error type instead of threading `Box<dyn std::error::Error>`.

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

/// Unified error type for the `tpt-fem` umbrella and its consumers.
///
/// Each per-crate error type is forwarded via `#[from]` (gated behind the
/// corresponding Cargo feature), alongside `std::io::Error`, TOML config errors,
/// and a generic message variant for ad-hoc CLI diagnostics.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Held when a per-crate error type is unavailable (feature disabled).
    #[error("{0}")]
    Msg(String),
    /// I/O failure (file read/write, import/export).
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// TOML problem-config parse failure.
    #[error("config error: {0}")]
    Config(#[from] toml::de::Error),
    /// `tpt-fem-mesh` error (validation, import, selectors).
    #[cfg(feature = "mesh")]
    #[error(transparent)]
    Mesh(#[from] tpt_fem_mesh::MeshError),
    /// `tpt-fem-sparse` solver error.
    #[cfg(feature = "sparse")]
    #[error(transparent)]
    Sparse(#[from] tpt_fem_sparse::SparseError),
    /// `tpt-fem-io-abaqus` parser error.
    #[cfg(feature = "io-abaqus")]
    #[error(transparent)]
    Inp(#[from] tpt_fem_io_abaqus::InpError),
    /// `tpt-fem-io-exodus` codec error.
    #[cfg(feature = "io-exodus")]
    #[error(transparent)]
    Exodus(#[from] tpt_fem_io_exodus::ExodusError),
    /// `tpt-fem-io-vtk` reader/writer error.
    #[cfg(feature = "io-vtk")]
    #[error(transparent)]
    Vtk(#[from] tpt_fem_io_vtk::VtkError),
    /// `tpt-fem-solve` nonlinear/continuation error.
    #[cfg(feature = "solve")]
    #[error(transparent)]
    Newton(#[from] tpt_fem_solve::NewtonError),
}

#[cfg(feature = "assembly")]
pub use tpt_fem_assembly::*;
#[cfg(feature = "composite")]
pub use tpt_fem_composite::*;
#[cfg(feature = "contact")]
pub use tpt_fem_contact::*;
#[cfg(feature = "coupling")]
pub use tpt_fem_coupling::*;
#[cfg(feature = "dofmap")]
pub use tpt_fem_dofmap::*;
#[cfg(feature = "dynamic")]
pub use tpt_fem_dynamic::*;
#[cfg(feature = "eigen")]
pub use tpt_fem_eigen::*;
#[cfg(feature = "elasticity")]
pub use tpt_fem_elasticity::*;
#[cfg(feature = "element")]
pub use tpt_fem_element::*;
#[cfg(feature = "fluid")]
pub use tpt_fem_fluid::*;
#[cfg(feature = "hyperelastic")]
pub use tpt_fem_hyperelastic::*;
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
#[cfg(feature = "modal")]
pub use tpt_fem_modal::*;
#[cfg(feature = "plasticity")]
pub use tpt_fem_plasticity::*;
#[cfg(feature = "porous")]
pub use tpt_fem_porous::*;
#[cfg(feature = "quadrature")]
pub use tpt_fem_quadrature::*;
#[cfg(feature = "solve")]
pub use tpt_fem_solve::*;
#[cfg(feature = "sparse")]
pub use tpt_fem_sparse::*;
#[cfg(feature = "thermal")]
pub use tpt_fem_thermal::*;
#[cfg(feature = "topopt")]
pub use tpt_fem_topopt::*;

/// A curated, feature-gated prelude.
///
/// `use tpt_fem::prelude::*;` pulls in the public API of every constituent
/// crate that is enabled by its Cargo feature, so a typical application can be
/// written against a single import. Like the crate-root re-exports, each
/// constituent is gated: only crates whose feature is on appear in the prelude.
pub mod prelude {
    #[cfg(feature = "assembly")]
    pub use tpt_fem_assembly::*;
    #[cfg(feature = "composite")]
    pub use tpt_fem_composite::*;
    #[cfg(feature = "contact")]
    pub use tpt_fem_contact::*;
    #[cfg(feature = "coupling")]
    pub use tpt_fem_coupling::*;
    #[cfg(feature = "dofmap")]
    pub use tpt_fem_dofmap::*;
    #[cfg(feature = "dynamic")]
    pub use tpt_fem_dynamic::*;
    #[cfg(feature = "eigen")]
    pub use tpt_fem_eigen::*;
    #[cfg(feature = "elasticity")]
    pub use tpt_fem_elasticity::*;
    #[cfg(feature = "element")]
    pub use tpt_fem_element::*;
    #[cfg(feature = "fluid")]
    pub use tpt_fem_fluid::*;
    #[cfg(feature = "hyperelastic")]
    pub use tpt_fem_hyperelastic::*;
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
    #[cfg(feature = "modal")]
    pub use tpt_fem_modal::*;
    #[cfg(feature = "plasticity")]
    pub use tpt_fem_plasticity::*;
    #[cfg(feature = "porous")]
    pub use tpt_fem_porous::*;
    #[cfg(feature = "quadrature")]
    pub use tpt_fem_quadrature::*;
    #[cfg(feature = "solve")]
    pub use tpt_fem_solve::*;
    #[cfg(feature = "sparse")]
    pub use tpt_fem_sparse::*;
    #[cfg(feature = "thermal")]
    pub use tpt_fem_thermal::*;
    #[cfg(feature = "topopt")]
    pub use tpt_fem_topopt::*;
}
