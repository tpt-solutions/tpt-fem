# Changelog

All notable changes to this crate are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This crate is `publish = false` (opt-in, excluded from the Cargo workspace).

## [Unreleased]

### Added

- `GpuContext` (wgpu adapter/device + compute pipelines), `GpuContext::spmv`
  (CSR matvec) and `solve_cg_gpu`: GPU-resident Jacobi-PCG for SPD systems with
  `f64` mixed-precision iterative refinement. CG scalars are computed on the
  GPU; the CPU syncs once per 16 iterations.
- `GpuError`, `GpuCgOptions`, `GpuCgSolution`.
- Tests (skip when no adapter) and a `bench_cg` example; measured 3.6x-7.7x
  faster than the CPU PCG from 250k to 1.4M unknowns on an RTX 3050.
