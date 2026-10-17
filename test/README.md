# Tests

MOSE-GPU is tested on the test cases of MOSE: the cases are not copied here, they are taken from a checkout of
[MOSE](https://github.com/open-hydra/MOSE).

| File | Content |
|---|---|
| `run_mose_tests.py` | stages the cases of MOSE, submits them with `mose`, collects and compares md5 tables |
| `reference_1it.md5` | the regression table: md5 of the field after one iteration of each case |

```bash
python3 test/run_mose_tests.py stage --mose /path/to/MOSE --work /scratch/$USER/mose_tests --one-iteration
python3 test/run_mose_tests.py submit --work /scratch/$USER/mose_tests -p <gpu partition>
python3 test/run_mose_tests.py md5 --work /scratch/$USER/mose_tests > new.md5
python3 test/run_mose_tests.py compare test/reference_1it.md5 new.md5
```

The reference table holds MOSE commit `b29130b`. The two shock / boundary-layer cases need `--swbli-seed`, the fields of
the coarse grid levels computed by MOSE; without it they are skipped.

See [Testing](https://open-hydra.github.io/MOSE-GPU/development/testing/) for the method and for the comparison with MOSE.
