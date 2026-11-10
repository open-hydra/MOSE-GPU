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

You normally do not call `install.sh`: `mose solve` reads the case and builds the matching binary with the recipe above,
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
| `MOSE_BC_LIST` | 1 | explicit step: the boundary fluxes on the record lists of every local block, one launch per block face over all the blocks (`BC_Fluxes_list`, as a multirate list level); `0` keeps the per-block launches. Off by itself with a shock detector or 3D rotational periodic records |
| `MOSE_FLUX_PACK` | 1 | explicit step: the flux phase on block packs (`MOSE_Pack_m`): one launch per kernel type over each batch of blocks (convective, face thermo, diffusive, pass 2), geometry in flat device pools; `0` keeps the per-block launches. Off by itself with a shock detector or a host source (rotating frame, soot) |
| `MOSE_FLUX_PACK_MB` | 2048 | memory budget (MB) of the face buffers of one batch of the block pack: fewer, larger batches with more memory |
| `MOSE_MR_XCH` | 1 | multirate on more than one rank: each stage ships only the inter-rank halo entries its level reads (the receiver's level boxes grown by gc), as compact per-level sub-lists on the resident buffers; `0` keeps the full exchange at every stage. Needs the device-bound halo path (CUDA-aware MPI). |
| `MOSE_PG_EDGE` | 1 | more than one rank: the transverse stencil of the connection records (`bc_pg(:,3:6)`) is rebuilt by the receiver from its own ghost cells (the same numbers the donor would pack) and shipped only for the entries at the face edges; `0` ships it for every entry. Needs the device-bound halo path. |
| `MOSE_PART` | (unset) | more than one rank: the blocks go to the ranks by contiguity along the connection faces (`partition_blocks_connectivity`, a deterministic multi-start search that minimises the inter-rank connection cells within the load tolerance), replacing the size-only LPT map when it cuts fewer connection cells; `lpt` keeps the LPT map (A/B). An `INPUT/block_owner.txt` map overrides both. Ownership only: the field is bit-identical under any map. |
| `MOSE_PART_TOL` | 1.02 | load tolerance of the contiguous partition: the cells of every rank stay within this factor of the average. |
| `MOSE_MR_LIST_DEV` | 1 | multirate list levels: the face and cell lists of each level are built on the device at every full partition (`mr_build_lists_dev`: per block a count kernel, then flag / prefix-sum / gather kernels per direction and for the cells, in the lexicographic order of the host build; the lists are created on the device without a host copy); `0` keeps the host build and its upload. Same lists to the element. |
| `MOSE_MR_THERMO_EXACT` | 1 | multirate list levels: the thermo-cache lists of a level are the cells its stage reads in the cache (the ten cells of every face, the level's cells, the interior cell of every ghost record), built on the device, instead of the level's boxes grown by gc. Same result to the bit (every cached value is a function of the cell's primitives; Compute_dt refreshes every interior cell at each macro step). `0` = the boxes. |
| `MOSE_MR_GHOST_EXACT` | 1 | multirate list levels, one rank: the ghost records of a level are those whose ghost cells the stage reads (convective and diffusive stencils of its faces, its cells and their six neighbours, a 5 x 3 x 3 box around each boundary record), marked on the device, instead of the boxes grown by gc. With more ranks the boxes are kept. `0` = the boxes. |
| `MOSE_BC_SORT` | 1 | boundary record lists (explicit and multirate): the records of each launch grouped by boundary type, stable, so that a warp follows one branch of the record body. Within a launch every cell receives at most one record: the same result to the bit. `0` keeps the face order. |
| `MOSE_EXPL_ASYNC` | 1 | explicit stage: the per-block kernels (ghost fills, thermo cache, wall ghosts of the turbulence, turbulence sources, update) on the stage queue in program order, a wait before the fluxes (synchronous stream) and one at the end of the stage; the thermodynamic guard is read after the first. Same kernels, same order: the same result to the bit. `0` = synchronous. Off by itself with the SA and Wilcox2006 sources and with a host source (rotating frame, soot): their kernels are not on the queue. |
| `MOSE_MR_ASYNC_BLK` | 1 | multirate: the block-path levels (level 0) on the stage queue as well, the fluxes bracketed on the synchronous stream; the queue is reset to synchronous after the macro step. `0` = level 0 synchronous. |
| `MOSE_MR_ASYNC` | 1 | multirate list levels: every kernel of a stage is queued on one async queue in program order (the CUDA Fortran list launches on its stream) and the host waits once at the end of the stage, instead of before and after each launch; the thermodynamic guard is read after that wait. Same kernels, same order: the same result to the bit. `0` keeps every stage kernel synchronous. Off by itself with more than one OpenMP thread. |
| `MOSE_MR_ASYNC_FLUX` | 1 | multirate block-path stages (level 0): the per-block flux kernels (beta, convective and diffusive pass 1 with the face thermo, telescoping, pass 2) on the stage queue as well, in program order on its one stream, instead of synchronous launches after a join of the queue. Same kernels, same order: the same result to the bit. `0` = synchronous. Synchronous by itself with a host source. |
| `MOSE_MR_REUSE` | 1 | multirate full partition: when no cell of any rank changed level since the last build and the ladder is the same, the level structures are kept; they depend on the levels and the geometry only, so the stages are the same to the bit. `0` rebuilds them at every full partition. |
| `MOSE_MR_BUILD_FAST` | 1 | multirate build of the level structures: the in-chunk prefix sums of the list build by sub-chunks, the counts of the list build and the changed cells of the nesting rule summed on the device with one read back. Integers only: the same levels, lists and positions. `0` = the previous scan (A/B timing). |
| `MOSE_MR_FLUX_FANOUT` | 1 | multirate list levels with `MOSE_MR_ASYNC`: the independent work of a stage on six side queues forked from the stage queue and joined back to it -- the pass-1 flux chains of each direction (convective; face thermo + diffusive: each writes only its own direction's arrays) and the coarser levels' halo interpolations (each writes only the cells of its own level). Each kernel computes what it did on the one queue, from the same inputs: the same result to the bit. `0` = one queue. |
| `MOSE_PG_EARLY` | 1 | multirate on more than one rank, with `MOSE_MR_XCH` and `MOSE_MR_ASYNC`: the receives of the transverse stencil (`bc_pg(:,3:6)`) are posted with the halo receives, and a rank whose shipped transverse entries of the level read no ghost cell another rank fills packs and sends them right after its local ghost fill, in flight with the halo, instead of after the halo unpack. Decided per rank and per level at every full partition; the same copies of the same cells: the same result to the bit. `0` = after the unpack. |
| `MOSE_MR_BLOCK_CENSUS` | 0 | diagnostic: `1` appends, at every build of the level structures, the cells per level of each local block to `OUTPUT/mr_block_census_r<rank>.txt` (full partition number, iteration, block, rank, levels 0..8), for a block-to-rank map that balances each level (`INPUT/block_owner.txt`). Prints only. |
| `MOSE_CONN_FAST` | 2 | the block-connection records of the boundary fluxes: `0` the original routines, `1` the face state of the two viscous flux calls computed once, `2` as `1` and the cell values read where they are used (fewer per-thread arrays: less local memory and fewer spills). The same operations on the same values: the same result to the bit |
| `MOSE_MR_FORCE_REDO` | 0 | test of `mr-nonlocal = redo`: `1` runs every accepted macro step once more, with the same step and map, from its restored start state; the result must not change |
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
