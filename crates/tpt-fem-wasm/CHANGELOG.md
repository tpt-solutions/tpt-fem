# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This crate is `publish = false` (an in-browser demo, built with wasm-pack), so it
is not versioned on crates.io.

## [Unreleased]

### Added

- In-browser (WebAssembly) topology-optimisation demo: `tpt-fem-topopt`'s
  SIMP cantilever optimiser exposed through `wasm-bindgen`, with a
  `www/index.html` page that runs it and draws the density field.
- CI job `wasm-demo` building the crate for `wasm32-unknown-unknown`.
