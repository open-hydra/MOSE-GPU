# GPU Execution Model

MOSE-GPU keeps the structure of MOSE: the same modules, the same block arrays, the same order of operations in a time
step. What changes is where the work runs. Every kernel of the time step runs on the GPU, written in OpenACC or CUDA Fortran.
The host reads the input, builds the mesh and the boundary tables, and writes the output.

## The state on the device

MOSE keeps the primitive state of each block in an array with ghost layers, `blk(b)%P(nprim,i,j,k)`. Next to it are the
state at the start of the step, `PO`, and the residual, `R`. The device works on a flat copy of the same numbers:
`P_flat(cell, var)`, `PO_flat` and `R_flat`. The cells of all the local blocks are concatenated with their ghost layers,
variable last, so that neighbouring threads read neighbouring memory.

| Phase | Where the state lives |
|---|---|
| set-up, restart, output | the block arrays `P` on the host |
| inside the time step | the flat arrays on the device |

The per-cell thermodynamic quantities, such as density, temperature, gas constant and viscosity, are cached the same way.
The conversion routines are in `src/lib/gpu/MOSE_Field_m.f90`.

## One time step

1. **Strang splitting of the chemistry.** When chemistry is active, a step is half a step of transport, a full step of the
   stiff chemical source, and another half step of transport.
2. **Runge–Kutta stages.** Each stage runs on the device, in order:
   1. the ghost cells: physical boundaries, connections, periodicity and chimera on the device, then the halo exchange
      between ranks;
   2. the thermodynamic cache of every cell;
   3. the convective fluxes, one CUDA Fortran kernel per direction, with the reconstruction, the limiter and the
      Riemann solver inlined;
   4. the diffusive fluxes, with the transport properties evaluated at the faces;
   5. the boundary fluxes, one kernel per block over its boundary records;
   6. the turbulence source terms;
   7. the state update and the check of the state, which stops the run on a non-physical state as MOSE does.
3. **The time step.** The local or global time step is computed on the device from the CFL and viscous limits of MOSE.

The fluxes of a face are accumulated into the residuals of its two cells without atomic operations, in a fixed order. A
run on a given number of ranks is therefore reproducible bit for bit.

## Compile-time constants

The number of species and of turbulence equations are compile-time constants of the device kernels. Loops over species
then have fixed trip counts, and the per-thread arrays fit in registers. The price is one binary per combination, which the
[`mose` launcher](../user/using.md) builds and caches.

## Where to read the code

`src/lib/numerics/time/explicit/Mod_Explicit.f90` holds the time step. The [Code Structure](../development/structure.md)
page maps every directory.
