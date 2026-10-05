# tpt-fem-gpu

Optional GPU acceleration for `tpt-fem`, built on [wgpu](https://wgpu.rs)
compute shaders (Vulkan, DX12, Metal). Dev/opt-in crate: `publish = false`,
excluded from the Cargo workspace like `tpt-fem-py`, so the core stays free of
GPU dependencies.

| Item | What it does |
|------|--------------|
| `GpuContext::new()` | Opens the highest-performance adapter; `adapter_name()`, `backend()`. |
| `GpuContext::spmv(&Csr, &[f64])` | CSR matrix–vector product on the GPU. |
| `solve_cg_gpu(&ctx, &Coo, &[f64], &GpuCgOptions)` | Jacobi-preconditioned CG for SPD systems, fully GPU-resident, with `f64` mixed-precision iterative refinement. |

From Python: build `tpt-fem-py` with `--features gpu` to get `gpu_solve_cg`.

Shaders compute in `f32`; `solve_cg_gpu` recovers double precision by solving
`A d = r` on the GPU and updating `x` / recomputing `r = b - A x` in `f64`
until `‖r‖/‖b‖ ≤ tol` (default `1e-10`). The CG scalars (`alpha`, `beta`) are
computed on the GPU and the CPU synchronises once per 16 iterations, which is
what makes it fast. Very ill-conditioned systems (condition number ≳ 1e6) can
stall in `f32`; use the CPU solvers in `tpt-fem-sparse` there.

## Measured (RTX 3050 laptop GPU, Vulkan, 2-D Laplacian, release build)

| unknowns | GPU PCG | CPU PCG (f64) | CPU skyline Cholesky |
|---------:|--------:|--------------:|---------------------:|
| 40 000 | 0.14 s | 0.11 s | 0.25 s |
| 250 000 | 0.76 s | 2.7 s | 8.4 s |
| 640 000 | 2.3 s | 12.6 s | 53 s |
| 1 440 000 | 6.4 s | 49 s | — |

Break-even is around 40–50 k unknowns; beyond that the GPU wins by 3.6× → 7.7×.
Reproduce with `cargo run --release --example bench_cg -- <grid size>`.

```sh
cargo test            # GPU tests skip (pass) when no adapter is present
cargo run --release --example bench_cg -- 500
```
