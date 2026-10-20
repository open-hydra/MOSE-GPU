# Parallel Execution

MOSE-GPU runs one MPI rank per GPU. Each rank owns whole blocks of the multi-block grid and advances them on its GPU; the
ghost cells between blocks on different ranks are exchanged every stage.

## The block-to-rank map

`INPUT/block_owner.txt` gives the owner of every block: one integer per block, the rank, starting from 0. Without the file
the blocks are distributed by size: largest first, each to the least loaded rank. `mose solve` takes the number of ranks
from the file when `-n` is not given, and refuses a `-n` that does not match it.

Two constraints come from the device boundary conditions:

- **Chimera blocks.** The donor cells of a chimera boundary must be on the same rank as the receiving block.
- **Manifold connections.** A manifold connection (MOSE type 501) runs on a single rank.

The run stops at start-up with a message when either is violated.

## GPUs and CPU cores

Node-local rank *k* uses CUDA device *k* modulo the number of visible devices. The script `mose_bind.sh`, which the
launcher puts in front of the binary, also pins each rank to the CPU cores close to its GPU. Set `MOSE_BIND=0` to switch
the pinning off. Results do not depend on the placement.

## The halo exchange

The ghost cells are packed and unpacked on the device, and the messages go from GPU memory to GPU memory through a
CUDA-aware MPI. On a machine whose MPI is not CUDA-aware, set `MOSE_GPU_AWARE_MPI=0` to stage the messages through host
memory.

## Reproducibility

The same binary with the same input and the same number of ranks gives the same field bit for bit from run to run. Runs
with a different number of ranks are not bit-identical: compare them by value, not by checksum. Two builds of the same
source are not always bit-identical either, because the compiler may allocate the registers of the device chemistry
integrator differently. Regression checks therefore compare the output, not the binary: see
[Testing](../development/testing.md).
