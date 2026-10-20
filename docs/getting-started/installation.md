# Installation

This page describes how to obtain MOSE-GPU and how its binaries are built.

## Prerequisites

| Requirement | Details |
|---|---|
| **NVIDIA HPC SDK** | `nvfortran` 24.11, with CUDA 12.6 |
| **GPU** | NVIDIA, compute capability 7.0 (V100) or 8.0 (A30, A100) |
| **CMake** | 3.23 or newer |
| **MPI** | a CUDA-aware MPI, such as the one shipped with the HPC SDK |
| **SLURM** | used by the `mose` launcher to build and run |
| **Python** | 3.6 or newer, standard library only, for the launcher |

## Source layout

```bash
git clone --recurse-submodules https://github.com/open-hydra/MOSE-GPU.git
cd MOSE-GPU
```

| Path | Content |
|---|---|
| `lib/FLINT` | thermodynamics, transport and chemistry, with the device integrators |
| `lib/FLINT/lib/OSlo` | ODE solvers used by FLINT |
| `lib/ORION` | I/O: Tecplot ASCII and binary (TecIO), VTK, Plot3D |
| `lib/third_party/FiNeR` | INI parser (Git submodule) |

`FLINT` and `ORION` are included as files: the GPU port needs commits of these libraries that are not yet published in
their own repositories. `FiNeR` is a Git submodule, as in MOSE.

## Building with the launcher

From a case folder, `mose solve` builds the binary the case needs and runs it:

```bash
cd my_case                       # input.ini + INPUT/ + MESH/
/path/to/MOSE-GPU/mose solve .
```

To type only `mose`, and to complete it with Tab, add to `~/.bashrc`:

```bash
export PATH=/path/to/MOSE-GPU:$PATH
source /path/to/MOSE-GPU/scripts/mose_completion.sh
```

The launcher reads three compile-time constants from the case:

| Constant | Read from |
|---|---|
| species count | the lines after the first of `INPUT/phase.txt` |
| turbulence equations | `equations` and `turbulence` in `[MOSE-Physics]`: 0 for Euler or laminar, 1 for Spalart–Allmaras, 2 for SST or Wilcox |
| chemistry integrator | `ode-solver` in `[MOSE-Chemistry]`: RADAU5 by default, RODAS on request |

Each combination is built once, by a short build job, and kept read-only in `$MOSE_BUILDS`. The run job waits for the
build job when the binary does not exist yet. See [The mose Launcher](../user/using.md).

## Building by hand

`install.sh` is the build script the launcher calls:

```bash
# the reference recipe: one binary for V100 and A30, 7 species, 2 turbulence equations
./install.sh build --production

# the same recipe for another case type
./install.sh build --production --nsc=4 --nrans=1

# incremental rebuild of the last configuration
./install.sh compile
```

The binary is written to `bin/MOSE`. The options are listed by `./install.sh --help` and described in
[Build Flags](../development/build-flags.md).

## Adapting the launcher to your cluster

The launcher was developed on one SLURM cluster, and its defaults are those of that machine. Change them with these
variables and options:

| Setting | Default | How to change it |
|---|---|---|
| binary cache | `/scratch/$USER/mose_builds` | `MOSE_BUILDS` |
| build partition | `free` | `MOSE_BUILD_PARTITION` |
| run partition | `fatgpu` | `mose solve -p <partition>` |
| HPC SDK location | `/opt/share/libs/nvidia/hpc_sdk/Linux_x86_64/24.11` | `NVHPC_ROOT` |
| modules loaded by the jobs | `nvidia/nvhpc/24.11`, `nvidia/cuda-12.6.2` | the job templates in `mose` |

## Troubleshooting

!!! warning "Failed to configure TecIO: Unable to find MPI installation"
    The first build of a clean tree also builds TecIO, inside `lib/ORION`, and its configuration looks for MPI on its own.
    If another MPI is loaded, for example by the default modules of a login shell, it finds the wrong one. The build jobs of
    `mose` start with `module purge` for this reason. When you run `install.sh` yourself, run it in a shell where only the
    NVIDIA HPC SDK provides MPI.
