# Verification & Validation

MOSE-GPU is validated against MOSE itself. The physical validation of the methods, against analytical solutions and
experiments, is the one of the [MOSE V&V pages](https://open-hydra.github.io/MOSE/vv/). Here the question is narrower:
does the GPU code give the same result as MOSE on the same input?

## Method

- **Same input, same files.** Each test of MOSE commit `b29130b` runs on both codes from its own `input.ini` and the files
  its `MOSE.sh` links. MOSE is built with the Intel compilers on CPU; MOSE-GPU runs on one GPU.
- **Same choices where the test leaves one open.** For example the 400-cell grid of the 1-D tests, or an iteration cap
  where the test has no stop criterion. Both codes get the same choice.
- **Measure.** For each variable of the final field, the relative L2 difference ‖GPU − MOSE‖ / ‖MOSE‖ over all cells,
  and the largest difference relative to the largest value. Velocities are measured against the velocity magnitude.

Two fields that differ only by round-off of double precision agree to 1e-13 or better, after the same number of steps.

## Summary

| | Tests of MOSE `b29130b` |
|---|---|
| total | 31 |
| run on the GPU | 24 |
| stop before the first step, by design | 7 |
| agree with MOSE | 22 of the 24 |
| differ, for a reason found and verified | 2 |

The seven tests that stop use options the GPU version does not have: other chemical mechanisms and integrators,
chemical equilibrium, grid sequencing, ablating walls, an inflow in a rotating frame.

The two differences are explained on their own pages:

- **Toro shock tube with HLLE++.** MOSE fills the ghost cells of a one-cell direction wrongly, and the shock sensor reads
  them. With that fill corrected in a test copy of MOSE, the two codes agree to 6e-14. See
  [Known Differences](known-differences.md).
- **Hypersonic cylinder.** The flow never reaches a steady state with the second-order scheme of the test: a few cells of
  the bow shock flip between two states every iteration. Each code differs from itself by as much as from the other when
  the CFL number changes by one part in 10^12. At first order both converge, and they agree to 1e-13. See
  [Hypersonic Cylinder](hypersonic-cylinder.md).

The test-by-test numbers, and the variants that exercise the options no test of MOSE uses, are in
[Comparison with MOSE](mose-comparison.md).
