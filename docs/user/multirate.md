# Multirate Time Integration

With an explicit scheme the whole domain advances at the step of its stiffest cell, usually a wall row of a RANS boundary
layer. The multirate integrator gives every cell a step close to its own limit. The cells are sorted into levels: over the
macro step `DT` a cell of level `l` takes `2^l` Heun steps of `DT / 2^l`.

The scheme is that of branch `multirate` of DIMA-THRUST/MOSE (88de097, `doc/archive/MULTIRATE_PLAN.md`):

- **Levels.** `lev = clamp(ceil(log2(s DT / dt_cell)), 0, L)`, where `dt_cell` is the cell's explicit limit from
  `Compute_dt` and `s` is the safety margin.
- **Level rules.**
  - Two face neighbours differ by at most one level.
  - Across a block connection the first two rows on each side share one level, so connection faces are never level
    interfaces.
- **Step cap.** `DT <= min over the cells of dt_cell 2^lev / s`. A cell the levels do not resolve shortens the macro
  step; it never makes the scheme unstable.
- **Recursion.**
  - The stage 1 of a level runs before the finer levels and its stage 2 after them.
  - The finer levels read the coarser cells next to them interpolated linearly in time.
- **Flux telescoping.** At the faces between two levels the coarse cell receives exactly the flux the fine cells took
  over its step: mass, momentum and energy are conserved.

This version keeps the scheme but none of the launch machinery of that branch:

- The stage kernels sweep boxes of the stage's level and update only its cells.
- The ghosts, boundary fluxes and sources are skipped on the blocks without cells of that level.
- The boundary flux kernels run over lists of the boundary records of the level, rebuilt with the levels.

## Keys

In `[MOSE-scheme]` (GPU dialect) or `[MOSE-Numerics]` (MOSE_open dialect):

| key | default | meaning |
|---|---|---|
| `multirate` | `f` | `t` integrates with the levels, `census` builds and prints them while the step stays explicit |
| `mr-ratio` | 16 | requested macro step in units of the smallest cell step, a power of 2 |
| `mr-safety` | 2 | margin: a cell's sub-step is at most its own limit divided by this factor |
| `mr-wall-safety` | 0 | margin of the cells within `mr-wall-rows` of a viscous wall (0: `mr-safety`) |
| `mr-wall-rows` | 10 | rows of the wall margin |
| `mr-lmax` | -1 | finest level (-1: from `mr-ratio` and the margins) |
| `mr-repart` | 20 | the levels are rebuilt at least every `mr-repart` macro steps |
| `mr-repart-hcap` | 0.8 | ... or when the cap of the frozen levels falls below this share of the step they gave at the full partition, scaled with the requested step since |
| `mr-print` | 200 | interval, in macro steps, of the `[MR]` lines |
| `mr-profile` | `f` | `t` times every call of the stages per level (a device wait after each call: a diagnostic, slower) |
| `mr-boxes` | 8 | at most this many cell boxes per block and level the stage kernels sweep (1 to 64) |
| `mr-list-level` | -1 | the levels from this one up run the whole stage (ghost fill, thermo cache, flux phase, boundary fluxes, wall ghosts, SST source, state update, state copy) on face, cell and record lists concatenated over the blocks, one launch per kernel and level (-1: all on the boxes; needs no shock detector, a viscous run and no 3D rotational periodic seam) |
| `mr-list-records-level` | 0 | with `mr-list-level` >= 0, the levels from the smaller of the two up run the ghost fill, the boundary fluxes and the wall ghosts on the record lists of their level; below `mr-list-level` their faces and cells stay on the boxes (-1: as `mr-list-level`) |
| `mr-ghost-near` | t | on a list level only the boundary and connection records that write a ghost the stage reads are filled, and only the cells of the coarser halos the stage reads are interpolated to its time (f: every record of the level's blocks and every coarser halo, as the boxes) |
| `mr-land` | t | the macro step is cut so that it ends exactly on the next output instant `k sol-dtime` or on `time-threshold`; the time is assigned, not summed, and the time-based output follows that landing |
| `mr-cons-halo` | t | the coarser halo is interpolated in time in the conservative variables, then converted back with the Newton iteration of the state update (f: in the primitive variables) |
| `mr-depth` | `ratio` | `auto`: at every full partition the depth comes from the spread of the admissible steps, `L = floor(log2(dt_cap,max / dt_min))` (at most `mr-lmax` and 8), and the macro step is `2^L dt_min / mr-safety`; `mr-ratio` is then not used |
| `mr-margin` | f | t: a cell's level and the step cap also take the smallest admissible step of its face neighbours in the block, divided by `kappa = max(cfl (1 + theta), 1.1)` (`theta` of the limiter: 1/2 minmod, 1 van Leer, MC and superbee, sqrt(2)/2 van Albada) |
| `mr-nonlocal` | `off` | `count`: the certificate of condition (4) of Gnedin, Semenov and Kravtsov (2018) is checked at every sub-step and counted (it reads the field and changes nothing); `redo`: a macro step with a violation, or with a state that is not physical, is run again from its start state with a deeper map or half the step (forces `mr-margin`) |
| `mr-octaves` | f | t: at the first full partition and at every 10th, the cells per octave of their admissible step and the ideal work of other macro steps (a diagnostic, prints only) |

Supported:

- `time-scheme = RK2`, `integration-variables = cons`, time-accurate runs;
- frozen chemistry, any number of ranks.

The run stops with an error on:

- IRS, preconditioning, `omega-source-relax`, a rotating frame or soot;
- a turbulence model other than SST (the sources of SA and Wilcox2006 are not restricted to the level);
- chimera faces, more than one multigrid level.

An iteration (`Iter =` in the log, `iter-threshold`) is one macro step.
With `state-check = t` the check runs once per macro step, after its last stage, on every cell. With
`mr-nonlocal = redo` a state that is not physical does not stop the run: the macro step is redone with half the step.

## Restart

A restart gives the result of the uninterrupted run only from a lossless backup:

- `bck-format = native binary` writes `OUTPUT/field.bin` and, next to it, `OUTPUT/field_state_r<rank>.bin`: the time,
  the iteration counters, the time-based dumps taken so far, the number of the output files (the dumps after the
  restart go on as `field<n+1>.tec`, not from 1) and the multirate partition of the step;
- the restart reads it with `rst-format = native binary` (default: `bck-format`). The state file is used only if its
  time is the time of the field, to the bit; without it the counters and the partition start afresh (`[restart]` line);
- a native file holds no grid: the grid comes from `mesh-file`. For the same result to the bit it must be the grid the
  run started from; a tecplot binary file (`.szplt`, e.g. the start field) is accepted.

The tecplot backups lose digits (ascii writes 15, binary single precision): a restart from them is close, not equal.

## Reading the `[MR]` lines

```
 [MR] iter 1  dt_min 1.144E-10  R_req 16.0  L 5  H_cap/H_req 1.000  R_eff 16.00  C_p 1.5  C_p_raw 1.4  passes 1  rule rounds 2
 [MR]   cells per level %:  L0 82.2  L1 9.7  L2 5.2  L3 2.3  L4 0.6  L5 0.0
```

- **`R_eff`:** the macro step in units of `dt_min`.
- **`C_p`:** the mean number of sub-steps per cell, `sum 2^lev / N`. The ideal gain against the explicit scheme is
  `R_eff / C_p`.
- **`light` line:** a macro step that kept the levels of the last full partition.
- **`wall s` line:** the wall since the start of the macro steps, of their stages and of the partition. With
  `mr-profile = t` a table follows: per level, the mean milliseconds of each call of a stage.

The certificate and the redo add lines to the `[MR]` census:

```
 [MR]   certificate (mr-nonlocal): violations 404  in macro steps 3
 [MR]   redo (mr-nonlocal): macro steps redone 3  rejected attempts 3  deeper maps 0
```

## Gates

- `multirate = t` with `mr-ratio = 1`, `mr-lmax = 0`, `mr-safety = 1` must give the explicit result bit for bit.
- With levels, compare a multirate run against the explicit run from the same field over the same time.
- `mr-nonlocal = count` must give the result of `off` to the bit; `redo` with `MOSE_MR_FORCE_REDO=1` (every accepted
  macro step run once more from the restored start state) the result of `redo`.
- A run split in two by a native restart must give the result of the run in one go to the bit.
