//! `tpt-fem` — umbrella crate for the `tpt-fem` core.
//!
//! This crate re-exports the Phase 1 `tpt-fem-*` core crates behind Cargo
//! features, so downstream users can depend on a single crate and opt into only
//! the pieces they need:
//!
//! | Feature       | Re-exports                       | Backing crate          |
//! |---------------|---------------------------------|-----------------------|
//! | `quadrature`  | `tpt_fem_quadrature::*`         | `tpt-fem-quadrature`  |
//! | `element`     | `tpt_fem_element::*`            | `tpt-fem-element`     |
//! | `mesh`        | `tpt_fem_mesh::*`               | `tpt-fem-mesh`        |
//! | `sparse`      | `tpt_fem_sparse::*`             | `tpt-fem-sparse`      |
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

#[cfg(feature = "element")]
pub use tpt_fem_element::*;
#[cfg(feature = "mesh")]
pub use tpt_fem_mesh::*;
#[cfg(feature = "quadrature")]
pub use tpt_fem_quadrature::*;
#[cfg(feature = "sparse")]
pub use tpt_fem_sparse::*;
