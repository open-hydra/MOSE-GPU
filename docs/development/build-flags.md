# Build Flags

Every performance option of the port that proved itself is now the only code path: no compile-time switch is left to
choose a kernel variant. What remains is the mechanism of the binary, the chemistry integrator and the infrastructure.

## The reference recipe

`./install.sh build --production` builds `bin/MOSE` with the recipe of the binaries the regression checks are run with:

- `-DMASTER=None`, OpenACC, CUDA Fortran, MPI and TecIO, NVIDIA HPC SDK 24.11;
- `-gpu=cc70,cc80,cuda12.6,ptxinfo`: one binary for V100 and A30. `ptxinfo` prints the register use of every kernel in the
  build log;
- 7 species and 2 turbulence equations, RADAU5 on the device; other counts with `--nsc` and `--nrans`;
- `CMAKE_Fortran_FLAGS='-cpp -Mbackslash -tp=x86-64-v3'`: generic AVX2 host code, so the binary is not tied to the CPU
  of the build node.

!!! warning "Compare outputs, not binaries"
    The build is not bit-reproducible in one routine, the device RADAU5: register spills can differ between two builds of
    the same source. The arithmetic does not change, but the binary md5 can. Regression checks therefore compare the md5
    of the output field.

## Automatic builds

You normally do not call `install.sh`: `mose submit` reads the case and builds the matching binary with the recipe above,
once. The case sets the compile-time constants:

| Constant | Read from |
|---|---|
| `MOSE_NSC` | the lines after the first of `INPUT/phase.txt` |
| `MOSE_NRANS` | the turbulence model: 2 for SST and Wilcox, 1 for Spalart–Allmaras, 0 otherwise |
| `MOSE_ODE_RODAS` | `ode-solver = H-rodas` (on) or `H-radau5` (off, the default) |

## The options

| Option | Default | What it selects | `install.sh` |
|---|---|---|---|
| `MOSE_NSC`, `MOSE_NRANS` | 4, 2 | the mechanism of the binary: species and turbulence counts are compile-time constants of the flux kernels, the device chemistry and the thermodynamic and transport loops | `--nsc=N --nrans=R` |
| `MOSE_ODE_RODAS` | off | device chemistry integrator: off is RADAU5, on is RODAS (about half the per-thread matrix storage; RODAS and RADAU5 differ by about 3e-8) | `--ode-rodas` |
| `USE_MPI`, `USE_TECIO`, `USE_OPENMP`, `MASTER`, GPU architectures | MPI and TecIO on | infrastructure | `--no-mpi`, `--no-tecio`, `--use-openmp`, `--master`, `--gpu-arch` |

Run-time switches, read from the environment:

| Variable | Default | Meaning |
|---|---|---|
| `MOSE_BC_NQ` | 2 | number of asynchronous queues of the boundary kernels |
| `MOSE_GPU_AWARE_MPI` | 1 | 0 stages the halo exchange through host memory, for an MPI that is not CUDA-aware |
| `MOSE_BIND` | 1 | 0 disables the CPU-core placement of `mose_bind.sh` |

## The device RADAU5

The device integrator (`lib/FLINT/src/lib/Lib_Radau5_dev.f90`) differs from Hairer's RADAU5 in four deliberate ways;
everything else follows the reference:

1. **Sizing.** The size of the system is a compile-time constant.
2. **Tolerances.** The tolerance transform is computed once at device initialisation instead of per cell. This is
   bit-exact.
3. **Right-hand side.** It calls the right-hand side and the Jacobian directly, because procedure pointers cannot be
   called on the GPU.
4. **Cold start.** The first Newton iteration reuses the right-hand side already computed at the start of the step,
   because the three stage evaluations equal it there. This is bit-exact.
