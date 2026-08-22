# tpt-fem-modal

Modal analysis and frequency-response for the `tpt-fem` finite-element core.

This crate sits at the intersection of [`tpt-fem-eigen`] (generalized
eigensolver) and [`tpt-fem-dynamic`] (time integration) and provides two
standard vibration workflows:

- **Frequency response** — the steady-state harmonic response
  `u(Ω)` of a structure under a sinusoidal load `f(t) = Re(F·e^{iΩt})`,
  resolved through a truncated set of eigenmodes. Used for vibration
  qualification and fatigue screening.
- **Modal superposition** — the time-domain complement: the eigenbasis
  reduces the second-order system to independent single-DOF modal
  equations, each integrated with Newmark, then recombined into `u(t)`.

## Example

```rust
use tpt_fem_modal::modal_analysis;
use tpt_fem_sparse::Coo;

// SDOF: M=1, K=4 (ω=2).
let m = Coo { rows: vec![0], cols: vec![0], vals: vec![1.0] };
let k = Coo { rows: vec![0], cols: vec![0], vals: vec![4.0] };
let data = modal_analysis(&k, &m, 0.0, 1, 8, 0.02).unwrap();
assert!((data.omega[0] - 2.0).abs() < 1e-6);
```

License: MIT OR Apache-2.0.
