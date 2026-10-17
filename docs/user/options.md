# Supported Options

Every option of MOSE, with what MOSE-GPU does with it. The option names are those of the
[MOSE input reference](https://open-hydra.github.io/MOSE/user/registry/). An option marked *stops* ends the run before the
first step, with a message; `mose check` says it before the queue.

## Numerics

| Key | Values on the GPU |
|---|---|
| `time-scheme` | `euler`, `RK2`, `RK3` |
| `time-accurate` | `true` (global time step), `false` (local time step) |
| `cfl`, `vnn`, `cfl-rise-threshold` | as in MOSE |
| `irs`, `irs-beta` | implicit residual smoothing, as in MOSE |
| `integration-variables` | `cons`, `prim`; `prec` with `time-accurate = false` only, as in MOSE |
| `preconditioning-Uref`, `-Mach`, `-eps-min`, `-Ur-min`, `-Ur-factor`, `-Ur-smooth` | as in MOSE |
| `space-reconstruction` | `MUSCL`, `first-order` |
| `flux-limiter` | `minmod`, `vanleer`, `vanalbada`, `mc`, `superbee` |
| `shock-detector` | `Tramel`, `Chen`; by default the one the Riemann solver asks for |
| `riemann-solver` | `HLLC`, `HLLE`, `HLLE++`, `HLLC+Tramel`, `HLLC+Chen`, `HLLC-PC` (with `prec`), `SLAU`, `SLAU2`, `AUSM+`, `AUSM+M`, `LLF`, `Rusanov`, `exact`, `LMRoe`, `MiczekRoe` |
| `riemann-options-Mco` | as in MOSE; required by `AUSM+M`, `HLLC+Chen` and `LMRoe` |
| `levels` greater than 1 | stops: grid sequencing is not on the GPU |

!!! note "First order with a flux limiter"
    With `space-reconstruction = first-order`, MOSE-GPU runs first order whether or not a `flux-limiter` is given. MOSE
    `b29130b` prints the same intention, but keeps the limiter and runs second order when the key is present: see
    [Known Differences](../vv/known-differences.md).

## Physics

| Key | Values on the GPU |
|---|---|
| `equations` | `euler`, `navier-stokes` |
| `turbulence` | `none`, `SA`, `SA-R`, `SA-RC`, `SAcomp`, `SST`, `SST-RC`, `Wilcox2006`; `SSGLRR` stops |
| `point-implicit` | as in MOSE, default `true` |
| `omega-wall-bc` | `practical`, `asymptotic` |
| `sst-production` | `incompressible`, `compressible` |
| `k-energy-coupling` | as in MOSE |
| `Prt`, `Sct` | as in MOSE |
| `schmidt` | a positive value: constant Schmidt number; `0` or absent with several species: mixture-averaged multicomponent diffusion |
| `chemistry` | `frozen`, `finite-rate`; `equilibrium` stops |
| mechanisms with `finite-rate` | `ONERA-7` and `Frolov_nopressure`, the mechanisms compiled for the device; any other stops |
| `ode-solver` | `H-radau5`, `H-rodas`; any other stops |
| `soot-generation` | stops |
| `rotational-frame` | `none`, `rigid-body`; an inflow boundary in a rotating frame stops |

## Boundary conditions

The boundary types are listed in [Boundary Conditions](boundary-conditions.md). Ablating and gas–surface-interaction walls
(503–506) stop.

## Output

| Key | Values on the GPU |
|---|---|
| `sol-format` | as in MOSE: `tecplot ascii`, `tecplot binary` (szplt), `vtk ascii`, `vtk raw` |
| `sol-diter`, `shell-diter`, `ini-diter` | as in MOSE |
| probes | as in MOSE: `[MOSE-Probes]` and one section per probe with `index-position` or `position` |

!!! warning "Probe position key"
    The probe cell is given by `index-position` (block, i, j, k). The key `index-location`, found in some MOSE test inputs,
    does not exist in MOSE: it is ignored there, and every such probe lands on the cell nearest the origin.
