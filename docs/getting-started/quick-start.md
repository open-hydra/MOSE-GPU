# Quick Start

This walkthrough runs the Sod shock tube of MOSE on a GPU and compares the result with MOSE.

## 1. Take a MOSE test case

MOSE-GPU reads the input files of MOSE as they are. Take the case from a MOSE checkout:

```bash
git clone https://github.com/open-hydra/MOSE.git
mkdir sod && cd sod
cp ../MOSE/test/1D/Sod79/input.ini .
mkdir -p INPUT
cp ../MOSE/test/1D/Sod79/INPUT/ic_x4.tec INPUT/ic.tec
cp ../MOSE/test/1D/Sod79/INPUT/bc_x4.txt INPUT/bc.txt
cp ../MOSE/test/common/Air/phase.txt ../MOSE/test/common/Air/thermo.dat INPUT/
```

`ic_x4.tec` holds the 400-cell grid of the test together with the initial state.

## 2. Check the case

```bash
/path/to/MOSE-GPU/mose check .
```

`mose check` says which binary the case needs and stops, before the queue, on any option the GPU version does not have.

## 3. Run it

```bash
/path/to/MOSE-GPU/mose submit . -n 1 -p <gpu partition>
```

The first run of a case type also builds its binary; the run job starts when the build ends. The job log `mose_<id>.out`
starts with the binary, its md5 and the source commits, and ends with the exit code and the wall time.

## 4. Look at the result

The output is written to `OUTPUT/`, as in MOSE:

| File | Content |
|---|---|
| `field.tec` | the final field |
| `residual-history.dat` | the residual of every iteration |
| `input.translated.ini` | the keys the code actually read |

Run the same case with MOSE: the two `field.tec` agree to round-off, about 1e-14, after the same 195 steps. The
[Comparison with MOSE](../vv/mose-comparison.md) page has every test of MOSE.
