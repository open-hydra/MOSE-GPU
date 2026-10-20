# The mose Launcher

`mose` is a Python script at the root of the repository. It reads a case folder, decides which binary the case needs,
builds it if it does not exist yet, and submits the run to SLURM. It needs Python 3.6 or newer and nothing else.

```text
mose check  [CASE]                  what the case needs and whether the GPU build runs it
mose build  [CASE]                  make sure the matching binary exists (submits its build job if not)
mose solve  [CASE] [-n N] [-p PART] [-t TIME] [-J NAME] [--mem MEM] [--nodelist NODE] [--dry-run]
mose list                           the binaries already built
mose example                        print an annotated input.ini of the native GPU dialect
```

`CASE` is the case folder, by default the current one. It holds `input.ini` and `INPUT/`, and `MESH/` when the case has one.

## What it reads from the case

| Compile-time constant | Read from |
|---|---|
| species count | the lines after the first of `INPUT/phase.txt` |
| turbulence equations | `[MOSE-Physics] equations` and `turbulence` |
| chemistry integrator | `[MOSE-Chemistry] ode-solver`: `H-radau5` (default) or `H-rodas` |

These are compile-time constants of the device kernels, so each combination has its own binary. Nothing is set by hand.

## check

`mose check` applies the same rules as the code, before the queue. It prints the binary the case needs and stops with a
message on any option the GPU version does not have, for example a chemistry mechanism that is not compiled for the
device or an inflow boundary in a rotating frame.

## build and solve

`mose solve` runs `check`, then looks for the binary in the cache. If it is missing, it submits a build job and makes the
run job wait for it. A second solve of the same combination reuses the binary, or waits for the build job already in the
queue.

| Option | Meaning |
|---|---|
| `-n N` | MPI ranks, one GPU each; default: the ranks of `INPUT/block_owner.txt`, else 1 |
| `--gpus G` | GPUs, by default one per rank |
| `-p PART` | SLURM partition of the run, the GPU nodes it goes to; default `fatgpu` |
| `-t TIME` | wall-time limit |
| `-J NAME` | job name |
| `--mem MEM` | memory per node |
| `--dry-run` | show what would be done, submit nothing |

## Suggestions and Tab completion

A wrong or incomplete command is answered with its usage, what was probably meant and a few examples: `mose solv`
suggests `solve`, `--partiton` suggests `--partition`, and `-p` with a partition that does not exist or has no GPUs
lists the partitions with GPUs. Run from a folder without `input.ini`, mose lists the case folders below it.
`mose submit` is the old name of `mose solve` and still works.

For Tab completion of the commands, the options, the case folders and the partitions, add to `~/.bashrc`:

```bash
source /path/to/MOSE-GPU/scripts/mose_completion.sh
```

## The binary cache

Binaries are kept read-only in `$MOSE_BUILDS`, one per combination and per state of the source tree. The source state is
the commit, the uncommitted changes and the untracked files, so a binary is never reused for different sources. The
build copies the source tree at submit time, so you can keep editing while a build waits in the queue.

## Provenance

The first lines of every job log `mose_<id>.out` name the binary, its md5 and the source commits it was built from; the
last line gives the exit code and the wall time. Keep the log with the results.

## Environment

| Variable | Default | Meaning |
|---|---|---|
| `MOSE_SRC` | the tree `mose` is in | source tree to build |
| `MOSE_BUILDS` | `/scratch/$USER/mose_builds` | binary cache |
| `MOSE_BUILD_PARTITION` | `free` | partition of the build jobs |
| `NVHPC_ROOT` | `/opt/share/libs/nvidia/hpc_sdk/Linux_x86_64/24.11` | NVIDIA HPC SDK |

The defaults are those of the cluster where MOSE-GPU was developed; see
[Adapting the launcher to your cluster](../getting-started/installation.md#adapting-the-launcher-to-your-cluster).
