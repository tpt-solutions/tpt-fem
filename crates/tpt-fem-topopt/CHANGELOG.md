# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-23

### Added

- `SimpProblem` / `SimpOptions` / `SimpResult` and `simp_optimize` — SIMP
  compliance minimisation with a sensitivity density filter, optimality-
  criteria update (bisection on the Lagrange multiplier), move limits, and
  under-relaxation.
- `cantilever_problem` benchmark builder (clamped quad grid, tip load).

[Unreleased]: https://github.com/tpt-solutions/tpt-fem/compare/HEAD
