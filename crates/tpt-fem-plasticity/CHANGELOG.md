# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `PlasticityParams` — material parameters `young`, `poisson`, `yield_stress`, `iso_hardening`, `kin_hardening`, with the `steel()` preset (E = 200 GPa, ν = 0.3, σ_y = 250 MPa, perfect plasticity).
- `PlasticState` — per-integration-point state holding the plastic strain `plastic_strain` (Voigt `6`), the kinematic `back_stress` (Voigt `6`), and the equivalent plastic strain `eps_eq` (`ε̄ᵖ`).
- `elastic_stiffness` — isotropic elastic stiffness as a Voigt `6×6` matrix from Young's modulus and Poisson's ratio.
- `j2_return_mapping` — the 3-D associative J2 (von Mises) radial-return update as incremental flow theory: the elastic predictor `σ_tr = C:(ε − εᵖ)` is projected back onto the yield surface along the relative-stress normal, updating `plastic_strain`, `back_stress` and `eps_eq`. Supports combined isotropic + kinematic hardening (the kinematic modulus appears in the consistency denominator as well as translating the surface).
- `plastic_1d` — the scalar 1-D axial analogue of `j2_return_mapping`: `(total axial strain, prev_eps_eq) → (σ, eps_eq, plastic)`, consistent with the tensor law for the uniaxial-stress configuration.
- `uniaxial_tangent_1d` — the consistent tangent `dσ/dε` for the 1-D law (elastic `E` before yield, radial-return tangent once plastic).
- `solve_elastic_plastic_rod` — a force-controlled 1-D rod solver (built from `Line2` elements) driven through a sequence of end loads with [`tpt-fem-solve`]'s Newton loop, carrying a per-element equivalent-plastic-strain state.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-plasticity-0.1.0
