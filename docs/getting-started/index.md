# Getting Started

MOSE-GPU is built with the NVIDIA HPC SDK and runs on NVIDIA GPUs through SLURM. Two pages get you from a clone to a
first result:

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Requirements, the library layout, and how a binary is built for a case.

    [:octicons-arrow-right-24: Installation](installation.md)

-   :material-rocket-launch:{ .lg .middle } __Quick Start__

    ---

    Run the Sod shock tube of MOSE on a GPU and compare it with MOSE.

    [:octicons-arrow-right-24: Quick Start](quick-start.md)

</div>

!!! tip "One binary per case type"
    The number of species and of turbulence equations are compile-time constants of the device kernels. You do not build
    by hand: the [`mose` launcher](../user/using.md) reads them from the case, builds the matching binary once, keeps it,
    and reuses it for every later run of the same kind.
