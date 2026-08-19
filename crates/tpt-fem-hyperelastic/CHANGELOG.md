# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- `Mat3` — a row-major `3×3` matrix type used throughout the crate.
- `mat_mul` / `mat_det` / `mat_inv` / `mat_transpose` — small `3×3` linear-algebra helpers (matrix product, determinant, inverse, transpose).
- `neo_hookean_piola` — the first Piola–Kirchhoff stress `P = μ F − p F^{-T}` for the incompressible Neo-Hookean model.
- `mooney_rivlin_piola` — the first Piola–Kirchhoff stress `P = 2 c₁ F + 2 c₂ (I₁ F − C F) − p F^{-T}` for the incompressible Mooney–Rivlin model (`C = FᵀF`, `I₁ = tr(C)`).
- `OgdenTerm` — a single Ogden term `(μ_i, α_i)`.
- `ogden_piola` — the first Piola–Kirchhoff stress for the incompressible (multi-term) Ogden model, evaluated from the principal stretches `λ` and principal directions (columns of `N`): `P = Σ μᵢ αᵢ λᵢ^{αᵢ−1} nᵢnᵢᵀ − p F^{-T}`.
- `neo_hookean_1d` — the scalar incompressible neo-Hookean nominal (PK1) stress `P = μ(F − F⁻²)` for a uniaxial stretch `F = λ`.
- `solve_hyperelastic_bar` — a 1-D bar solver (built from `Line2` elements) driven to a target end stretch with [`tpt-fem-solve`]'s Newton loop using the incompressible neo-Hookean response.

[0.1.0]: https://github.com/tpt-solutions/tpt-fem/releases/tag/tpt-fem-hyperelastic-0.1.0
