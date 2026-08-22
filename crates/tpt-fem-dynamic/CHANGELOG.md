# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `modal_frequency_response` — harmonic frequency-response analysis by modal
  superposition: solves the `K x = ω² M x` eigenproblem with
  `tpt-fem-eigen`'s generalized shift-invert Lanczos and evaluates the
  damped modal sum over a frequency sweep, returning real/imaginary
  displacement amplitudes (`ModalFrequencyResponse`). The frequency-domain
  counterpart of stepping `newmark` to steady state.
- `DynamicError::Sparse` / `DynamicError::InvalidInput` variants for the new
  workflow's failure modes.

## [0.1.0] - 2026-08-20

### Added

- `coo_scale` — scale a `Coo` matrix by a scalar, returning a new `Coo`.
- `coo_add` — add two `Coo` matrices, returning a new `Coo`.
- `coo_matvec` — sparse matrix-vector product `y = A x` for a `Coo`.
- `NewmarkOptions` — `dt`, `β`, `γ` options for the implicit integrator (with defaults).
- `newmark` — implicit Newmark-beta integration of `M·ü + C·v + K·u = f(t)`.
- `CentralOptions` — `dt` option for the explicit integrator (with default).
- `central_difference` — explicit central-difference integration (mass lumped internally).

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-dynamic-0.1.0
