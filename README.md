<p align="center">
  <h1 align="center">MOSE-GPU</h1>
  <p align="center"><b>MOSE on NVIDIA GPUs</b></p>
</p>

<p align="center">
  <a href="https://open-hydra.github.io/MOSE-GPU/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran%20%7C%20OpenACC%20%7C%20CUDA%20Fortran-734f96.svg" alt="Language: Fortran, OpenACC, CUDA Fortran">
  <a href="https://github.com/open-hydra/MOSE-GPU/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
</p>

---

MOSE-GPU is the GPU version of [MOSE](https://github.com/open-hydra/MOSE), the open-source compressible-flow solver of the
Hydra suite for the Euler and Navier–Stokes equations on multi-block structured grids. It solves the same equations with
the same numerical methods, reads the same input files, and runs every kernel of the time step on NVIDIA GPUs, one GPU per
MPI rank. On the test cases of MOSE it gives the same results as MOSE.

## Features

- **The physics and numerics of MOSE** — every Riemann solver, limiter, shock sensor, time scheme and set of integration
  variables of MOSE, including low-Mach preconditioning.
- **Turbulence modeling** — Spalart–Allmaras with the R, RC and compressibility corrections, Menter SST and SST-RC,
  Wilcox k-ω 2006.
- **Finite-rate chemistry** — stiff integration of every cell on the device, with RADAU5 or RODAS.
- **Multi-GPU** — one GPU per MPI rank, with the halo exchange on the device.
- **MOSE input files** — a MOSE case runs as it is; options the GPU version does not have stop the run with a message.
- **Automatic builds** — the `mose` launcher reads a case, builds the binary it needs once, and submits the run.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **NVIDIA HPC SDK** | `nvfortran` 24.11, CUDA 12.6 |
| **GPU** | NVIDIA, compute capability 7.0 or 8.0 |
| **CMake** | ≥ 3.23 |
| **MPI** | CUDA-aware |
| **SLURM** | for the `mose` launcher |

### Build and run

```bash
git clone --recurse-submodules https://github.com/open-hydra/MOSE-GPU.git

# from a case folder (input.ini + INPUT/): build the binary the case needs, then run it
cd my_case
/path/to/MOSE-GPU/mose check .
/path/to/MOSE-GPU/mose submit . -p <gpu partition>
```

To build by hand: `./install.sh build --production`; the executable is placed in `bin/MOSE`.

See the [Installation Guide](https://open-hydra.github.io/MOSE-GPU/getting-started/installation/) and the
[Quick Start](https://open-hydra.github.io/MOSE-GPU/getting-started/quick-start/).

## Dependencies

| Library | Role | In this repository |
|---|---|---|
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Thermodynamics, transport and chemistry, with the device integrators | `lib/FLINT`, included as files |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | ODE solver library used by FLINT | `lib/FLINT/lib/OSlo`, included as files |
| [ORION](https://github.com/MarcoGrossi92/ORION) | I/O: Tecplot (including TecIO), VTK, Plot3D | `lib/ORION`, included as files |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI parser | `lib/third_party/FiNeR`, Git submodule |

FLINT and ORION are included as files because MOSE-GPU uses commits of them that are not yet in their own repositories.

## Project Structure

```
MOSE-GPU/
├── src/
│   ├── app/           # Main application
│   ├── lib/           # Solver library, in the directories of MOSE
│   │   └── gpu/       # The flat device state and the shared kernel bodies
│   └── test/          # Unit tests
├── lib/               # FLINT, ORION, FiNeR
├── test/              # Regression table on the test cases of MOSE
├── scripts/           # Boundary-file checks
├── docs/              # MkDocs documentation source
├── cmake/             # CMake modules
├── mose               # The launcher
├── install.sh         # Build script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/MOSE-GPU](https://open-hydra.github.io/MOSE-GPU/)**, covering:

- Installation and quick start
- The launcher, the input file and the supported options
- The comparison with MOSE, test by test, and the known differences
- The GPU execution model and the code structure

## License

MOSE-GPU is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
