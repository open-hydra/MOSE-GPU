# Input File

MOSE-GPU reads two dialects of `input.ini`:

- **The MOSE dialect.** This is the file of a MOSE case, as written for MOSE.
- **The native GPU dialect.** This is the format of the GPU code before it was aligned with MOSE.

The code recognises the dialect at start-up. A file with a `[MOSE-Parameters]`, a `[MOSE-Numerics]` or a `[MOSE-Physics]`
section with `equations` is read as a MOSE file. Anything else is read as the GPU dialect, unchanged.

## The MOSE dialect

The meaning of every key is the one of the [MOSE input reference](https://open-hydra.github.io/MOSE/user/input/). At
start-up the code translates the keys of the file into the keys the GPU code reads. It adds them in memory and never
rewrites `input.ini`. The result is saved in `OUTPUT/input.translated.ini`, one line per key, with the MOSE key it comes
from.

The translation follows three rules:

1. **The defaults of MOSE.** A key absent from the file takes the default of MOSE, also where the GPU code had another
   default.
2. **The conventions of MOSE where they differ from the GPU code.** Without a `schmidt` key, a mixture of several species
   gets mixture-averaged multicomponent diffusion, which needs `INPUT/diffusion.dat`. The turbulence sources are treated
   point-implicitly unless `point-implicit = false`. The time step has no species-diffusion limit.
3. **A stop, never an approximation.** An option the GPU version does not have stops the run before the first step, with
   a message that names the key. `mose check` applies the same rules before the queue.

The options and what happens to each are listed in [Supported Options](options.md).

## Files read from `INPUT/`

The files are those of MOSE:

| File | Content |
|---|---|
| `bc.txt` | boundary conditions; see [Boundary Conditions](boundary-conditions.md) |
| `ic.tec` | the grid and the initial state, or the restart field |
| `phase.txt`, `thermo.dat` | species and their thermodynamics |
| `chemistry-info.txt` | the reaction mechanism, with `chemistry = finite-rate` |
| `diffusion.dat` | the binary diffusion table, with multicomponent diffusion |
| `block_owner.txt` | optional: the block-to-rank map; see [Parallel Execution](parallel.md) |

## The native GPU dialect

The GPU code also keeps its own dialect, with the sections `[MOSE-parameters]`, `[MOSE-scheme]`, `[MOSE-RANS]` and
`[MOSE-ODE]`. `mose example` prints an annotated example:

```ini
--8<-- "docs/user/input.gpu.example.ini"
```

New cases should use the MOSE dialect: it is the one validated against MOSE.
