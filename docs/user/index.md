# User Guide

MOSE-GPU runs the cases of MOSE. This guide covers what is specific to the GPU version: how a case is built and launched,
how it runs on several GPUs, which input keys are read, and which options are available. For the meaning of the physical
and numerical options, the [MOSE User Guide](https://open-hydra.github.io/MOSE/user/) applies as it is.

| Page | Content |
|---|---|
| [The mose Launcher](using.md) | check a case, build its binary, submit it |
| [Parallel Execution](parallel.md) | MPI ranks, GPUs and the block-to-rank map |
| [Input File](input.md) | the MOSE input file, its translation, and the native GPU dialect |
| [Supported Options](options.md) | every option of MOSE, with what the GPU version does with it |
| [Boundary Conditions](boundary-conditions.md) | the boundary types of MOSE on the GPU |
