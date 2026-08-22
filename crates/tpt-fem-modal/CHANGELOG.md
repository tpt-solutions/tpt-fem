# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned

- Mode-acceleration and residual-vector corrections for truncated bases.
- Per-mode damping ratios (modal/Caughey damping) instead of a uniform ζ.

## [0.1.0] - 2026-08-13

### Added

- `modal_analysis`: generalized Lanczos eigensolve of `K φ = ω² M φ`
  returning `ModalData` (ascending natural frequencies `ω_i`, mode shapes,
  modal masses, and an assumed uniform modal damping ratio `ζ`).
- `ModalData::frequency_response`: steady-state harmonic response
  `U(Ω)` under sinusoidal loading via modal superposition, including
  viscous modal damping (`HarmonicResponse` complex amplitudes).
- `ModalData::modal_superposition`: time-domain response history by
  integrating each independent modal equation with `tpt-fem-dynamic`'s
  Newmark scheme and recombining into the physical DOFs.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-modal-0.1.0
