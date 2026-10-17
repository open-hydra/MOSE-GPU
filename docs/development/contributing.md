# Contributing

Thank you for your interest in contributing to MOSE-GPU. This page describes the workflow and what a change must show
before it is merged.

## Workflow

1. **Branch** from `main`.
2. **Implement** the change in small, focused commits.
3. **Check** it: the regression table, and the comparison with MOSE when results are meant to change. See
   [Testing](testing.md).
4. **Push** the branch and open a pull request against `main`.
5. **Review.** Address the feedback; maintainers merge when ready.

```bash
git clone --recurse-submodules https://github.com/open-hydra/MOSE-GPU.git
cd MOSE-GPU
git switch -c feature/my-change
```

## What a pull request shows

- **A change that must not alter results** keeps the one-iteration table identical to `test/reference_1it.md5`. Say so
  in the description, with the table.
- **A change that alters results** brings a comparison with MOSE on the same input: one step and to the end of a test
  that exercises the change, with the largest difference per variable.
- **A change to a device kernel** keeps the register and spill lines of the build log, from `ptxinfo`, no worse for the
  production kernels.

## Physics and numerics

MOSE-GPU follows MOSE. A new method, or a fix to an existing one, goes into MOSE as well, so that the two codes keep
giving the same results. When the GPU code finds a defect in MOSE, report it to MOSE with the file and the line, a run
that shows it, and the proposed fix.

## Coding conventions

- Fortran 2008 with OpenACC and CUDA Fortran, built with `nvfortran`.
- `implicit none` and argument intents in every new routine.
- No line longer than 132 characters.
- A device kernel and its host twin, when there is one, share one body through an included file, so that they cannot
  drift apart.

## Reporting bugs

Open an issue with:

- the commit;
- the GPU and the HPC SDK version;
- the job log, whose first lines name the binary and the source commits;
- the input file and the smallest case that shows the problem.

## License

MOSE-GPU is released under the GNU General Public License v3.0. Contributions are accepted under the same license.
