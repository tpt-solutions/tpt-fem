//! Optional GPU acceleration for `tpt-fem` (wgpu compute: Vulkan / DX12 / Metal).
//!
//! The crate offloads the two kernels that dominate large iterative solves:
//!
//! * [`GpuContext::spmv`] — CSR sparse matrix–vector product `y = A x`.
//! * [`solve_cg_gpu`] — Jacobi-preconditioned conjugate gradients for symmetric
//!   positive-definite systems, run **entirely on the GPU** (matrix and Krylov
//!   vectors stay resident; only two scalars per iteration cross the bus).
//!
//! Shaders compute in `f32` (portable across every adapter). To still reach
//! double-precision answers, [`solve_cg_gpu`] wraps the GPU solve in
//! **mixed-precision iterative refinement**: the inner `f32` CG solves
//! `A d = r`, the outer loop updates `x += d` and recomputes the residual
//! `r = b - A x` in `f64` on the CPU until `‖r‖/‖b‖ ≤ tol`. Systems whose
//! condition number is beyond what `f32` can resolve (≳ 1e6) may stall; use the
//! CPU solvers in `tpt-fem-sparse` for those.
//!
//! ```no_run
//! use tpt_fem_gpu::{solve_cg_gpu, GpuContext, GpuCgOptions};
//! use tpt_fem_sparse::Coo;
//!
//! let mut a = Coo::new();
//! for i in 0..1000 {
//!     a.push(i, i, 2.0);
//!     if i + 1 < 1000 { a.push(i, i + 1, -1.0); a.push(i + 1, i, -1.0); }
//! }
//! let b = vec![1.0; 1000];
//! let gpu = GpuContext::new().expect("a GPU adapter");
//! println!("running on {}", gpu.adapter_name());
//! let sol = solve_cg_gpu(&gpu, &a, &b, &GpuCgOptions::default()).unwrap();
//! assert!(sol.relative_residual <= 1e-10);
//! ```

use std::sync::mpsc;

use bytemuck::{Pod, Zeroable};
use tpt_fem_sparse::{Coo, Csr};
use wgpu::util::DeviceExt;

/// Errors from the GPU layer.
#[derive(Debug)]
pub enum GpuError {
    /// No usable GPU adapter was found.
    NoAdapter,
    /// The adapter refused to create a device.
    Device(String),
    /// Mapping a GPU buffer back to the CPU failed.
    Readback(String),
    /// The input system is unusable (non-square, wrong rhs length, non-positive
    /// diagonal, too large for 32-bit GPU indices, ...).
    InvalidInput(String),
    /// The solver did not reach the tolerance within the iteration budget.
    NotConverged {
        /// Final relative residual `‖b - A x‖ / ‖b‖` (computed in `f64`).
        relative_residual: f64,
        /// Outer refinement steps taken.
        refinements: usize,
    },
}

impl std::fmt::Display for GpuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GpuError::NoAdapter => write!(f, "no usable GPU adapter found"),
            GpuError::Device(m) => write!(f, "GPU device creation failed: {m}"),
            GpuError::Readback(m) => write!(f, "GPU readback failed: {m}"),
            GpuError::InvalidInput(m) => write!(f, "invalid GPU solve input: {m}"),
            GpuError::NotConverged {
                relative_residual,
                refinements,
            } => write!(
                f,
                "GPU conjugate gradients did not converge after {refinements} refinement \
                 step(s) (relative residual {relative_residual:e})"
            ),
        }
    }
}

impl std::error::Error for GpuError {}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Params {
    n: u32,
    alpha: f32,
    beta: f32,
    _pad: u32,
}

const WG: u32 = 64;
const DOT_WG: u32 = 256;

const SPMV: &str = r#"
struct P { n: u32, alpha: f32, beta: f32, pad: u32 };
@group(0) @binding(0) var<storage, read> row_ptr: array<u32>;
@group(0) @binding(1) var<storage, read> col_ind: array<u32>;
@group(0) @binding(2) var<storage, read> vals: array<f32>;
@group(0) @binding(3) var<storage, read> x: array<f32>;
@group(0) @binding(4) var<storage, read_write> y: array<f32>;
@group(0) @binding(5) var<uniform> p: P;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    if (row >= p.n) { return; }
    var s = 0.0;
    let end = row_ptr[row + 1u];
    for (var k = row_ptr[row]; k < end; k = k + 1u) {
        s = s + vals[k] * x[col_ind[k]];
    }
    y[row] = s;
}
"#;

const DOT: &str = r#"
struct P { n: u32, alpha: f32, beta: f32, pad: u32 };
@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read_write> partial: array<f32>;
@group(0) @binding(3) var<uniform> p: P;
var<workgroup> sh: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>,
        @builtin(local_invocation_id) lid: vec3<u32>,
        @builtin(workgroup_id) wid: vec3<u32>) {
    var v = 0.0;
    if (gid.x < p.n) { v = a[gid.x] * b[gid.x]; }
    sh[lid.x] = v;
    workgroupBarrier();
    var s = 128u;
    loop {
        if (s == 0u) { break; }
        if (lid.x < s) { sh[lid.x] = sh[lid.x] + sh[lid.x + s]; }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (lid.x == 0u) { partial[wid.x] = sh[0]; }
}
"#;

// Scalars live in a storage buffer: [0]=rz, [1]=pap, [2]=alpha, [3]=beta.
//
// d += alpha p;  r -= alpha Ap;  z = minv .* r
const UPDATE: &str = r#"
struct P { n: u32, alpha: f32, beta: f32, pad: u32 };
@group(0) @binding(0) var<storage, read_write> d: array<f32>;
@group(0) @binding(1) var<storage, read_write> r: array<f32>;
@group(0) @binding(2) var<storage, read> p_vec: array<f32>;
@group(0) @binding(3) var<storage, read> ap: array<f32>;
@group(0) @binding(4) var<storage, read> minv: array<f32>;
@group(0) @binding(5) var<storage, read_write> z: array<f32>;
@group(0) @binding(6) var<storage, read> scal: array<f32>;
@group(0) @binding(7) var<uniform> p: P;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= p.n) { return; }
    let alpha = scal[2];
    d[i] = d[i] + alpha * p_vec[i];
    let ri = r[i] - alpha * ap[i];
    r[i] = ri;
    z[i] = minv[i] * ri;
}
"#;

// p = z + beta p
const PUPDATE: &str = r#"
struct P { n: u32, alpha: f32, beta: f32, pad: u32 };
@group(0) @binding(0) var<storage, read> z: array<f32>;
@group(0) @binding(1) var<storage, read_write> p_vec: array<f32>;
@group(0) @binding(2) var<storage, read> scal: array<f32>;
@group(0) @binding(3) var<uniform> p: P;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= p.n) { return; }
    p_vec[i] = z[i] + scal[3] * p_vec[i];
}
"#;

// Single-workgroup reductions of the per-workgroup dot partials into the
// scalar buffer. `p.n` carries the number of partials.
const REDUCE: &str = r#"
struct P { n: u32, alpha: f32, beta: f32, pad: u32 };
@group(0) @binding(0) var<storage, read> partial: array<f32>;
@group(0) @binding(1) var<storage, read_write> scal: array<f32>;
@group(0) @binding(2) var<uniform> p: P;
var<workgroup> sh: array<f32, 256>;

fn total(lid: u32) -> f32 {
    var v = 0.0;
    for (var k = lid; k < p.n; k = k + 256u) { v = v + partial[k]; }
    sh[lid] = v;
    workgroupBarrier();
    var s = 128u;
    loop {
        if (s == 0u) { break; }
        if (lid < s) { sh[lid] = sh[lid] + sh[lid + s]; }
        workgroupBarrier();
        s = s >> 1u;
    }
    return sh[0];
}

// pap = p . Ap ; alpha = rz / pap (0 if pap is not positive).
@compute @workgroup_size(256)
fn reduce_pap(@builtin(local_invocation_id) lid: vec3<u32>) {
    let pap = total(lid.x);
    if (lid.x == 0u) {
        scal[1] = pap;
        if (pap > 0.0) { scal[2] = scal[0] / pap; } else { scal[2] = 0.0; }
    }
}

// rz' = r . z ; beta = rz' / rz (0 if rz is not positive); rz <- rz'.
@compute @workgroup_size(256)
fn reduce_rz(@builtin(local_invocation_id) lid: vec3<u32>) {
    let rz_new = total(lid.x);
    if (lid.x == 0u) {
        let rz = scal[0];
        if (rz > 0.0) { scal[3] = rz_new / rz; } else { scal[3] = 0.0; }
        scal[0] = rz_new;
    }
}
"#;

/// A GPU device plus the compiled compute pipelines.
pub struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    info: wgpu::AdapterInfo,
    spmv: wgpu::ComputePipeline,
    dot: wgpu::ComputePipeline,
    update: wgpu::ComputePipeline,
    pupdate: wgpu::ComputePipeline,
    reduce_pap: wgpu::ComputePipeline,
    reduce_rz: wgpu::ComputePipeline,
}

fn pipeline(device: &wgpu::Device, label: &str, src: &str) -> wgpu::ComputePipeline {
    pipeline_entry(device, label, src, "main")
}

fn pipeline_entry(
    device: &wgpu::Device,
    label: &str,
    src: &str,
    entry: &str,
) -> wgpu::ComputePipeline {
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::Wgsl(src.into()),
    });
    device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some(label),
        layout: None,
        module: &module,
        entry_point: Some(entry),
        compilation_options: Default::default(),
        cache: None,
    })
}

impl GpuContext {
    /// Open the highest-performance adapter available.
    pub fn new() -> Result<Self, GpuError> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
        }))
        .ok_or(GpuError::NoAdapter)?;
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("tpt-fem-gpu"),
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(),
                memory_hints: wgpu::MemoryHints::Performance,
            },
            None,
        ))
        .map_err(|e| GpuError::Device(e.to_string()))?;
        Ok(GpuContext {
            spmv: pipeline(&device, "spmv", SPMV),
            dot: pipeline(&device, "dot", DOT),
            update: pipeline(&device, "update", UPDATE),
            pupdate: pipeline(&device, "pupdate", PUPDATE),
            reduce_pap: pipeline_entry(&device, "reduce_pap", REDUCE, "reduce_pap"),
            reduce_rz: pipeline_entry(&device, "reduce_rz", REDUCE, "reduce_rz"),
            device,
            queue,
            info,
        })
    }

    /// Human-readable adapter name (e.g. `"NVIDIA GeForce RTX 3050"`).
    pub fn adapter_name(&self) -> &str {
        &self.info.name
    }

    /// The graphics backend in use (`Vulkan`, `Dx12`, `Metal`, ...).
    pub fn backend(&self) -> String {
        format!("{:?}", self.info.backend)
    }

    fn storage(&self, data: &[u8], writable: bool) -> wgpu::Buffer {
        let mut usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
        if writable {
            usage |= wgpu::BufferUsages::COPY_SRC;
        }
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: data,
                usage,
            })
    }

    fn empty_storage(&self, floats: usize) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: (floats.max(1) * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        })
    }

    fn read_f32(&self, src: &wgpu::Buffer, floats: usize) -> Result<Vec<f32>, GpuError> {
        let bytes = (floats.max(1) * 4) as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: bytes,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, &staging, 0, bytes);
        self.queue.submit([enc.finish()]);
        self.map_staging(&staging, floats)
    }

    fn map_staging(&self, staging: &wgpu::Buffer, floats: usize) -> Result<Vec<f32>, GpuError> {
        let (tx, rx) = mpsc::channel();
        staging.slice(..).map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .map_err(|e| GpuError::Readback(e.to_string()))?
            .map_err(|e| GpuError::Readback(e.to_string()))?;
        let data = staging.slice(..).get_mapped_range();
        let out: Vec<f32> = bytemuck::cast_slice(&data)[..floats].to_vec();
        drop(data);
        staging.unmap();
        Ok(out)
    }

    fn uniform(&self, p: Params) -> wgpu::Buffer {
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: bytemuck::bytes_of(&p),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            })
    }

    fn bind(&self, pipe: &wgpu::ComputePipeline, buffers: &[&wgpu::Buffer]) -> wgpu::BindGroup {
        let entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .enumerate()
            .map(|(i, b)| wgpu::BindGroupEntry {
                binding: i as u32,
                resource: b.as_entire_binding(),
            })
            .collect();
        self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipe.get_bind_group_layout(0),
            entries: &entries,
        })
    }

    fn dispatch(
        enc: &mut wgpu::CommandEncoder,
        pipe: &wgpu::ComputePipeline,
        group: &wgpu::BindGroup,
        threads: usize,
        wg: u32,
    ) {
        let mut pass = enc.begin_compute_pass(&Default::default());
        pass.set_pipeline(pipe);
        pass.set_bind_group(0, group, &[]);
        pass.dispatch_workgroups((threads as u32).div_ceil(wg).max(1), 1, 1);
    }

    /// `y = A x` on the GPU (single precision internally; the result is widened
    /// back to `f64`).
    pub fn spmv(&self, a: &Csr, x: &[f64]) -> Result<Vec<f64>, GpuError> {
        let m = GpuMatrix::upload(self, a)?;
        if x.len() < a.ncols {
            return Err(GpuError::InvalidInput(format!(
                "x has {} entries, matrix has {} columns",
                x.len(),
                a.ncols
            )));
        }
        let xf: Vec<f32> = x.iter().map(|v| *v as f32).collect();
        let xb = self.storage(bytemuck::cast_slice(&xf), false);
        let yb = self.empty_storage(a.nrows);
        let params = self.uniform(Params {
            n: a.nrows as u32,
            alpha: 0.0,
            beta: 0.0,
            _pad: 0,
        });
        let bg = self.bind(
            &self.spmv,
            &[&m.row_ptr, &m.col_ind, &m.vals, &xb, &yb, &params],
        );
        let mut enc = self.device.create_command_encoder(&Default::default());
        Self::dispatch(&mut enc, &self.spmv, &bg, a.nrows, WG);
        self.queue.submit([enc.finish()]);
        Ok(self
            .read_f32(&yb, a.nrows)?
            .into_iter()
            .map(f64::from)
            .collect())
    }
}

/// A CSR matrix resident in GPU memory.
struct GpuMatrix {
    row_ptr: wgpu::Buffer,
    col_ind: wgpu::Buffer,
    vals: wgpu::Buffer,
}

impl GpuMatrix {
    fn upload(ctx: &GpuContext, a: &Csr) -> Result<Self, GpuError> {
        if a.nnz() > u32::MAX as usize || a.nrows > u32::MAX as usize {
            return Err(GpuError::InvalidInput(
                "matrix is too large for 32-bit GPU indices".into(),
            ));
        }
        let row_ptr: Vec<u32> = a.row_ptrs.iter().map(|&v| v as u32).collect();
        let col_ind: Vec<u32> = a.col_ind.iter().map(|&v| v as u32).collect();
        let vals: Vec<f32> = a.values.iter().map(|&v| v as f32).collect();
        Ok(GpuMatrix {
            row_ptr: ctx.storage(bytemuck::cast_slice(&row_ptr), false),
            col_ind: ctx.storage(bytemuck::cast_slice(&col_ind), false),
            vals: ctx.storage(bytemuck::cast_slice(&vals), false),
        })
    }
}

/// Options for [`solve_cg_gpu`].
#[derive(Clone, Copy, Debug)]
pub struct GpuCgOptions {
    /// Target relative residual `‖b - A x‖ / ‖b‖`, measured in `f64`.
    pub tol: f64,
    /// Maximum outer (`f64` residual) refinement steps.
    pub max_refinements: usize,
    /// Maximum inner GPU CG iterations per refinement (`0` means `10 · n`).
    pub max_inner: usize,
    /// Inner-solve reduction target (relative, in the preconditioned norm).
    pub inner_tol: f64,
}

impl Default for GpuCgOptions {
    fn default() -> Self {
        GpuCgOptions {
            tol: 1e-10,
            max_refinements: 40,
            max_inner: 0,
            inner_tol: 1e-4,
        }
    }
}

/// Result of [`solve_cg_gpu`].
#[derive(Clone, Debug)]
pub struct GpuCgSolution {
    /// The solution (double precision).
    pub x: Vec<f64>,
    /// Outer refinement steps taken.
    pub refinements: usize,
    /// Total inner GPU CG iterations.
    pub gpu_iterations: usize,
    /// Final relative residual in `f64`.
    pub relative_residual: f64,
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Solve the symmetric positive-definite system `A x = b` with
/// Jacobi-preconditioned conjugate gradients on the GPU, wrapped in `f64`
/// mixed-precision iterative refinement (see the crate docs).
pub fn solve_cg_gpu(
    ctx: &GpuContext,
    coo: &Coo,
    rhs: &[f64],
    opts: &GpuCgOptions,
) -> Result<GpuCgSolution, GpuError> {
    let a = coo.to_csr();
    let n = a.nrows;
    if a.ncols != n {
        return Err(GpuError::InvalidInput(format!(
            "matrix must be square, got {n} x {}",
            a.ncols
        )));
    }
    if rhs.len() != n {
        return Err(GpuError::InvalidInput(format!(
            "rhs length {} does not match matrix dimension {n}",
            rhs.len()
        )));
    }
    if n == 0 {
        return Ok(GpuCgSolution {
            x: Vec::new(),
            refinements: 0,
            gpu_iterations: 0,
            relative_residual: 0.0,
        });
    }
    let bnorm = norm(rhs);
    if bnorm == 0.0 {
        return Ok(GpuCgSolution {
            x: vec![0.0; n],
            refinements: 0,
            gpu_iterations: 0,
            relative_residual: 0.0,
        });
    }

    // Jacobi preconditioner.
    let mut minv = vec![0.0f32; n];
    for (r, slot) in minv.iter_mut().enumerate() {
        let mut d = 0.0;
        for c in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            if a.col_ind[c] == r {
                d += a.values[c];
            }
        }
        if !(d.is_finite() && d > 0.0) {
            return Err(GpuError::InvalidInput(format!(
                "solve_cg_gpu requires a positive diagonal (row {r} has {d})"
            )));
        }
        *slot = (1.0 / d) as f32;
    }

    let m = GpuMatrix::upload(ctx, &a)?;
    let minv_b = ctx.storage(bytemuck::cast_slice(&minv), false);
    let d_b = ctx.empty_storage(n);
    let r_b = ctx.empty_storage(n);
    let z_b = ctx.empty_storage(n);
    let p_b = ctx.empty_storage(n);
    let ap_b = ctx.empty_storage(n);
    let n_wg = n.div_ceil(DOT_WG as usize).max(1);
    let partial_b = ctx.empty_storage(n_wg);
    let params = ctx.uniform(Params {
        n: n as u32,
        alpha: 0.0,
        beta: 0.0,
        _pad: 0,
    });
    let reduce_params = ctx.uniform(Params {
        n: n_wg as u32,
        alpha: 0.0,
        beta: 0.0,
        _pad: 0,
    });
    // Scalars: [0]=rz, [1]=pap, [2]=alpha, [3]=beta.
    let scal_b = ctx.empty_storage(4);

    let bg_spmv = ctx.bind(
        &ctx.spmv,
        &[&m.row_ptr, &m.col_ind, &m.vals, &p_b, &ap_b, &params],
    );
    let bg_dot_pap = ctx.bind(&ctx.dot, &[&p_b, &ap_b, &partial_b, &params]);
    let bg_dot_rz = ctx.bind(&ctx.dot, &[&r_b, &z_b, &partial_b, &params]);
    let bg_update = ctx.bind(
        &ctx.update,
        &[&d_b, &r_b, &p_b, &ap_b, &minv_b, &z_b, &scal_b, &params],
    );
    let bg_pupdate = ctx.bind(&ctx.pupdate, &[&z_b, &p_b, &scal_b, &params]);
    let bg_reduce_pap = ctx.bind(&ctx.reduce_pap, &[&partial_b, &scal_b, &reduce_params]);
    let bg_reduce_rz = ctx.bind(&ctx.reduce_rz, &[&partial_b, &scal_b, &reduce_params]);

    let max_inner = if opts.max_inner == 0 {
        10 * n
    } else {
        opts.max_inner
    };
    // Iterations recorded per submission; the scalars (alpha, beta) are
    // computed on the GPU, so the CPU only syncs once per batch to test
    // convergence.
    const BATCH: usize = 16;

    let mut x = vec![0.0f64; n];
    let mut r64 = rhs.to_vec();
    let mut gpu_iterations = 0usize;
    let mut rel = 1.0f64;
    for refinement in 0..=opts.max_refinements {
        rel = norm(&r64) / bnorm;
        if rel <= opts.tol {
            return Ok(GpuCgSolution {
                x,
                refinements: refinement,
                gpu_iterations,
                relative_residual: rel,
            });
        }
        if refinement == opts.max_refinements {
            break;
        }

        // Inner solve A d = r on the GPU, scaled to unit norm to stay in f32 range.
        let scale = norm(&r64);
        let r32: Vec<f32> = r64.iter().map(|v| (v / scale) as f32).collect();
        let z32: Vec<f32> = r32.iter().zip(&minv).map(|(r, m)| r * m).collect();
        let rz0: f64 = r32
            .iter()
            .zip(&z32)
            .map(|(a, b)| f64::from(*a) * f64::from(*b))
            .sum();
        ctx.queue.write_buffer(&r_b, 0, bytemuck::cast_slice(&r32));
        ctx.queue.write_buffer(&z_b, 0, bytemuck::cast_slice(&z32));
        ctx.queue.write_buffer(&p_b, 0, bytemuck::cast_slice(&z32));
        ctx.queue
            .write_buffer(&d_b, 0, bytemuck::cast_slice(&vec![0.0f32; n]));
        ctx.queue.write_buffer(
            &scal_b,
            0,
            bytemuck::cast_slice(&[rz0 as f32, 0.0f32, 0.0, 0.0]),
        );
        let mut done = 0usize;
        // The first iteration must not apply a stale beta to p (= z already).
        let mut first = true;
        while done < max_inner {
            let k = BATCH.min(max_inner - done);
            let mut enc = ctx.device.create_command_encoder(&Default::default());
            for _ in 0..k {
                if !first {
                    GpuContext::dispatch(&mut enc, &ctx.pupdate, &bg_pupdate, n, WG);
                }
                first = false;
                GpuContext::dispatch(&mut enc, &ctx.spmv, &bg_spmv, n, WG);
                GpuContext::dispatch(&mut enc, &ctx.dot, &bg_dot_pap, n, DOT_WG);
                GpuContext::dispatch(&mut enc, &ctx.reduce_pap, &bg_reduce_pap, 256, 256);
                GpuContext::dispatch(&mut enc, &ctx.update, &bg_update, n, WG);
                GpuContext::dispatch(&mut enc, &ctx.dot, &bg_dot_rz, n, DOT_WG);
                GpuContext::dispatch(&mut enc, &ctx.reduce_rz, &bg_reduce_rz, 256, 256);
            }
            ctx.queue.submit([enc.finish()]);
            done += k;
            gpu_iterations += k;
            let sc = ctx.read_f32(&scal_b, 4)?;
            let (rz, pap) = (f64::from(sc[0]), f64::from(sc[1]));
            if !rz.is_finite()
                || pap.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
                || rz <= opts.inner_tol * opts.inner_tol * rz0
            {
                break;
            }
        }

        // x += scale * d ; r = b - A x  (double precision, on the CPU).
        let d32 = ctx.read_f32(&d_b, n)?;
        for i in 0..n {
            x[i] += scale * f64::from(d32[i]);
        }
        let ax = a.matvec(&x);
        for i in 0..n {
            r64[i] = rhs[i] - ax[i];
        }
    }
    Err(GpuError::NotConverged {
        relative_residual: rel,
        refinements: opts.max_refinements,
    })
}
