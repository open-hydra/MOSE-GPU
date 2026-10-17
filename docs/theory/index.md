# Theoretical Guide

MOSE-GPU solves the equations of MOSE with the numerical methods of MOSE. The theory is therefore the one of the
[MOSE Theoretical Guide](https://open-hydra.github.io/MOSE/theory/):

| Topic | MOSE page |
|---|---|
| Governing equations | [Governing Equations](https://open-hydra.github.io/MOSE/theory/governing-equations/) |
| Spatial discretization, reconstruction and limiters | [Spatial Discretization](https://open-hydra.github.io/MOSE/theory/numerics/) |
| Riemann solvers | [Riemann Solvers](https://open-hydra.github.io/MOSE/theory/riemann-solvers/) |
| Time integration, residual smoothing, preconditioning | [Time Integration](https://open-hydra.github.io/MOSE/theory/time-integration/) |
| Turbulence models | [Turbulence Modelling](https://open-hydra.github.io/MOSE/theory/turbulence/) |
| Rotating frame | [Rotating Frame](https://open-hydra.github.io/MOSE/theory/rotating-frame/) |
| Thermodynamic and transport properties | [Thermodynamic and Transport Properties](https://open-hydra.github.io/MOSE/theory/thermo/) |
| Finite-rate kinetics | [Finite-Rate Kinetics](https://open-hydra.github.io/MOSE/theory/kinetics/) |

Chemical equilibrium and gas–surface interaction are described there too, but are not available on the GPU.

## What is specific to the GPU version

- **Execution.** The time step runs on the device: see [GPU Execution Model](gpu-execution.md).
- **Chemistry integrators.** The stiff integration of the chemical source runs cell by cell on the device, with device
  versions of RADAU5 and RODAS. Their differences from the reference implementations are listed in
  [Build Flags](../development/build-flags.md#the-device-radau5).
- **Corrections of MOSE defects.** Where MOSE reads a value that is never assigned, the GPU version computes the value the
  method defines. These places are listed in [Known Differences](../vv/known-differences.md).
