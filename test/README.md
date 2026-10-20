# Tests

The test cases of [MOSE](https://github.com/open-hydra/MOSE) that run on the GPU, ready to run, at their path in `test/` of
MOSE (`1D/Sod79`, `2D/euler/chang-nozzle`, ...), and the regression table.

| Path | Content |
|---|---|
| `1D/`, `2D/`, `3D/` | 21 test cases of MOSE b29130b: `input.ini`, `INPUT/`, `MESH/` |
| `run_mose_tests.py` | stages the cases, submits them with `mose`, collects and compares md5 tables; `import` takes the cases from a MOSE checkout |
| `reference_1it.md5` | the regression table: md5 of the field after one iteration of each case |

To run one case, `mose case` copies it to `/scratch/$USER`:

```bash
/path/to/MOSE-GPU/mose case Sod79
cd /scratch/$USER/1D_Sod79 && /path/to/MOSE-GPU/mose submit . -p <gpu partition>
```

The regression table:

```bash
python3 test/run_mose_tests.py stage --work /scratch/$USER/mose_tests --one-iteration
python3 test/run_mose_tests.py submit --work /scratch/$USER/mose_tests -p <gpu partition>
python3 test/run_mose_tests.py md5 --work /scratch/$USER/mose_tests > new.md5
python3 test/run_mose_tests.py compare test/reference_1it.md5 new.md5
```

The cases were written by `run_mose_tests.py import` as the validation against MOSE staged them: the 1-D tests on the
400-cell grid, the isentropic vortex and the manufactured solution on 64 × 64 (the fields of the other grids are left
out), the nozzles capped at 5000 iterations and pressure-centrifugal at 2000. The two shock / boundary-layer cases (swbli
SA and SST) are not here: they restart from the fields of the coarse grid levels computed by MOSE, and are staged from a
MOSE checkout with `--mose /path/to/MOSE --swbli-seed DIR`, DIR holding `SA/field.tec` and `SST/field.tec`; without these
options they are skipped.

See [Testing](https://open-hydra.github.io/MOSE-GPU/development/testing/) for the method and for the comparison with MOSE.
