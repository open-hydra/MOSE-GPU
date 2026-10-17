# Testing

Every change to MOSE-GPU is checked in two ways before it is committed:

- **Regression.** Where the change must not alter results, the fields stay identical to the bit.
- **Comparison with MOSE.** Where the change is meant to alter results, MOSE-GPU must agree with MOSE on the same input.

## Regression: the one-iteration table

The test cases of MOSE that run on the GPU are run for one iteration each, and the md5 of the field is compared with a
reference table. One iteration is enough to touch every kernel a case uses: the ghost cells, the thermodynamic cache, the
fluxes, the boundary fluxes, the sources and the update. The md5 changes if a single bit of the field changes. A
refactoring or a performance change must leave the whole table unchanged.

`test/run_mose_tests.py` stages the cases from a MOSE checkout, submits them with `mose`, and compares the md5:

```bash
# stage the cases of MOSE b29130b, one iteration each
python3 test/run_mose_tests.py stage --mose /path/to/MOSE --work /scratch/$USER/mose_tests --one-iteration

# submit them: mose builds the binaries the cases need
python3 test/run_mose_tests.py submit --work /scratch/$USER/mose_tests -p <gpu partition>

# when the jobs have ended: the md5 table, and the comparison with the reference
python3 test/run_mose_tests.py md5 --work /scratch/$USER/mose_tests > new.md5
python3 test/run_mose_tests.py compare test/reference_1it.md5 new.md5
```

The cases are staged as the validation of [Comparison with MOSE](../vv/mose-comparison.md) staged them. The files are those
the `MOSE.sh` of each test links, the 1-D tests use the 400-cell grid, and the isentropic vortex and the manufactured
solution use 64 × 64. The two shock / boundary-layer cases start from a field of the coarse grid levels that MOSE computes,
since the GPU has no grid sequencing; they are staged only when `--swbli-seed` points to those fields.

`test/reference_1it.md5` holds the reference table. It was made with binaries built by `install.sh build --production`,
on V100 and A30, which give the same bits.

## Comparison with MOSE

When a change adds or alters a feature, both codes run the same input, and the fields are compared value by value:

- **one step**, which isolates the new term;
- **to the end of the test**, which shows whether a difference grows.

The case must exercise the new term: a shock for a Riemann solver, curved streamlines for a rotation correction, a wall
for a turbulence model. Two fields that differ only by round-off agree to 1e-13 or better, after the same number of steps.

When MOSE reads an undefined value, compare against a test copy of MOSE with that single line fixed, and report the defect
to MOSE: see [Known Differences](../vv/known-differences.md).

## Kernel resources

Building with `ptxinfo` writes the registers and the spills of every device kernel to the build log. A change to a kernel
should not add spills to the production kernels; compare the lines before and after.
