# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-13

### Added

- `newton` Newton–Raphson solver with `NewtonOptions` / `NewtonError`.
- `continuation` parameter-continuation driver (warm-started Newton).
- `arc_length_continuation` arc-length (Ramm) continuation with
  `ArcLengthOptions` / `ArcLengthError`.
- Physics-agnostic residual/jacobian callback interface.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-solve-0.1.0
