---
title: Overview
---

# Overview

MOSE-GPU is the GPU version of [MOSE](https://github.com/open-hydra/MOSE) (**M**ultipurpose **O**pen **S**olver for
ideal-gas **E**quations), the compressible-flow solver of the Hydra suite. It solves the same equations with the same
numerical methods, reads the same input files, and runs every kernel of the time step on NVIDIA GPUs, one GPU per MPI rank.

The reference is MOSE commit `b29130b`, the head of `open-hydra/MOSE` when this version was released. On the test cases of
that commit the two codes give the same results to round-off; the [Verification & Validation](vv/index.md) section shows
the comparison case by case.

---

## Hydra CFD Suite

| Component | Role | Status |
|-----------|------|--------|
| [**ATLAS**](https://github.com/open-hydra/ATLAS) | Pre-processor: mesh, initial and boundary conditions, thermodynamic and chemical data | Separate package |
| [**MOSE**](https://github.com/open-hydra/MOSE) | Solver: compressible Euler and Navier–Stokes with finite-rate chemistry, on CPUs | Separate package |
| **MOSE-GPU** | The same solver on NVIDIA GPUs | This package |

!!! info "The same input files as MOSE"
    A case prepared for MOSE, by ATLAS or by hand, runs on MOSE-GPU as it is. The code translates the input keys at
    start-up and stops with a message, before the first step, on any option the GPU version does not have. See
    [Input File](user/input.md).

---

## What runs on the GPU

The [governing equations](https://open-hydra.github.io/MOSE/theory/governing-equations/) and the
[numerical methods](https://open-hydra.github.io/MOSE/theory/numerics/) are those of MOSE. On the GPU:

- **Convective fluxes.** MUSCL reconstruction with the minmod, van Leer, van Albada, MC and superbee limiters, or first
  order. Every Riemann solver of MOSE: HLLC, HLLE, HLLE++, HLLC+Tramel, HLLC+Chen, HLLC-PC, SLAU, SLAU2, AUSM+, AUSM+M,
  LLF/Rusanov, exact, LMRoe and MiczekRoe. The Tramel and Chen shock sensors.
- **Diffusive fluxes.** Constant Schmidt number or mixture-averaged multicomponent diffusion.
- **Time integration.** Explicit Euler, RK2 and RK3, global or local time step, implicit residual smoothing. Conservative,
  primitive or preconditioned integration variables.
- **Turbulence.** Spalart–Allmaras with the R, RC and compressibility corrections, Menter SST and SST-RC, Wilcox
  k-ω 2006, with the point-implicit treatment of the destruction terms.
- **Chemistry.** Finite-rate kinetics integrated cell by cell on the device with RADAU5 or RODAS, for the mechanisms
  compiled into the device code.
- **Boundary conditions and ghost cells.** Walls, symmetry, extrapolation, inlets and outlets, block connections,
  rotational periodicity and chimera.
- **Rotating frame.** Rigid-body rotation.

The host reads the input, builds the mesh and the boundary tables, and writes the output. Everything inside the time step
runs on the device: see [GPU Execution Model](theory/gpu-execution.md).

---

## What MOSE-GPU does not have

These options of MOSE stop the run with a message before the first step:

- grid sequencing (`levels > 1`);
- chemical equilibrium and soot;
- mechanisms other than the device ones, and ODE solvers other than RADAU5 and RODAS;
- ablating and gas–surface-interaction walls;
- inflow boundaries in a rotating frame.

The full list is in [Supported Options](user/options.md).

---

## How it is validated

Every change is checked in two ways:

1. **No unintended change.** A table of MOSE test cases, one iteration each, must keep the same field bit for bit.
2. **Same results as MOSE.** Both codes run the same input files, and the fields are compared value by value.

The method is described in [Testing](development/testing.md). The results are in
[Comparison with MOSE](vv/mose-comparison.md). The few places where the two codes still differ, and why, are in
[Known Differences](vv/known-differences.md).
