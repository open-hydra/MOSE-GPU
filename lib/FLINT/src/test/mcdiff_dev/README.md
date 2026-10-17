# co_DS_expr_dev unit test (branch gpu-port-mcdiff)

Checks that the device routine `co_DS_expr_dev` (Lib_ThermoTransport_dev.f90) reproduces the host routine `co_DS_expr`
(Lib_ThermoTransport.f90, from main 6e95b37). Both compute mixture-averaged diffusion coefficients:
D_k = (1 − X_k) / Σ_{j≠k} X_j / D_kj(T), rescaled by p_ref/p.

## What the test does
- Reads `phase.txt`, `thermo.dat` and `diffusion.dat` from an input directory with `read_idealgas_thermo` and `read_idealgas_diffusion`.
- Sends the tables to the device through `flint_acc_upload_thermo`, the same upload the solver uses.
- Evaluates both routines on the same states:
  - the host routine in a plain loop;
  - the device routine in an `acc parallel loop`.
- Each state gets the same ρ, (Tint, Tdiff) and p in both routines, so only the two implementations are compared.
- Writes both results to a binary file and prints, per species, the number of non-identical values and the largest relative difference.
- `check_py.py` recomputes D_k independently in numpy and compares it with both results.
- `split_check.py` separates real states from synthetic ones and expresses each difference in units of ε/(1 − X_k) (see "Threshold").

## Files
| file | role |
|---|---|
| `CMakeLists.txt` | stand-alone build with the flags of the MOSE_GPU production device build (MOSEL target, RELEASE, USE_OPENACC, MOSE_CT with MOSE_TT_NS = 7) |
| `test_co_DS_dev.F90` | driver: `test_co_DS_dev <input dir> <states.bin> <out.bin>` |
| `make_states.py` | writes `states.bin` from THOR pocket slices (`pocket_mr*.npz`) plus edge and random cases |
| `make_input.py` | cuts a `thermo.dat` at a given Tmax (the 7-species `diffusion.dat` stops at 5000 K) |
| `check_py.py` | independent numpy implementation, compared with host and device |
| `split_check.py` | real vs synthetic states; differences normalised by the conditioning of (1 − X_k) |

## Build
The configure step needs three external pieces:
- the ORION source tree the solver builds with;
- the OSlo modules and libraries from an existing MOSE_GPU build (OSlo is not rebuilt here);
- the variant options.
```
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DORION_SRC=<ORION tree> \
      -DOSLO_MODDIR=<dir with oslo.mod + interface_definitions.mod> -DOSLO_LIBS="<libOSlo.a;MKL;iomp5;libfatode.a;libiode;-lstdc++;-lnvcpumath>" \
      [-DTT_NS=7|OFF] [-DACC=ON|OFF] [-DGPU_OPTS=cc70,cc80,cuda12.6]
```
Variants:
- `TT_NS=7`: the compile-time build used in production;
- `TT_NS=OFF`: the universal build;
- `ACC=OFF`: host-only, where the device routine runs on the CPU. It allows the Python check without a GPU.

## Inputs used on 25/09/2026
| input | origin | md5 |
|---|---|---|
| `diffusion.dat` | 7-species ONERA table (HyShotII copy), T = 1…5000 K, Pref 101325 Pa | 4e21b371 |
| `thermo.dat` | THOR `thermo.dat` (T = 1…10000 K) cut at 5000 K with `make_input.py` | 64df8365 |
| `phase.txt` | THOR (O2, H2O, H2, H, O, OH, N2) | d7bd227f |
| `states.bin` | 1 234 016 states: 1 224 000 THORBL2 pocket cells (x 12 and 25 mm, 1800 × 68, fields +396/576/612/648/684 µs; T 142–3034 K, p 0.29–31 bar), plus 10 016 synthetic states | 4b3f182e |

The synthetic states cover:
- pure O2, H2, H and N2 (the near-pure fallback);
- air at T from 0.3 to 12000 K (clamp of the temperature bracket);
- air at p from 1e3 to 5e6 Pa;
- 10 000 random Dirichlet compositions (α 0.3) over 200–4500 K and 1e4–3e6 Pa.

**Prerequisite for the solver side:**
- `read_idealgas_diffusion` requires the diffusion table to end at the same Tmax as the thermo table (otherwise ios = 3).
- The THOR thermo table ends at 10000 K, so the solver needs a diffusion table to 10000 K.
- ATLAS GPB writes it with `write_diffusion_properties(name, T_low, T_max, …)`.

## Threshold
- The threshold was set before the run: host against device ≤ 1e-14 relative on every value.
- That threshold did not account for the conditioning of the numerator. D_k is proportional to (1 − X_k), so one ulp of difference in X_k gives a relative difference of about ε/(1 − X_k) in D_k. For a species close to pure, that is well above 1e-14.
- **Real states:** the 1e-14 threshold holds.
- **Synthetic near-pure states:** the difference is judged in units of ε/(1 − X_k) and must stay at a few units (round-off amplified by the conditioning, not an implementation error).
- Host and device values need not be bit-identical in the compile-time build. There the device routine inlines and unrolls its own copies of the helper routines (`f_molecularWeight_m`), so FMA contraction and summation order differ from the generic host routine.

## Results (25/09/2026)
Host-only CT 7 build (job 302752):
- 0 NaN.
- Python against host and against the device routine on the CPU: largest relative difference 3.9e-13, worst species H.
- **Real states:**
  - host against device ≤ **3.0e-15**;
  - 0.11 % of the values not bit-identical;
  - ≤ 3.6 ε/(1 − X_k).
- **Synthetic states:** up to 7.3e-14 at X_H = 0.997, ≤ 5.3 ε/(1 − X_k).

Device runs on one A30 (job 302753; binaries: CT 7 90daac7b, universal c0dc2fe5):

| | CT 7 (production) | universal |
|---|---|---|
| kernel launches of the test loop (NV_ACC_NOTIFY) | 1 | 1 |
| NaN host / device | 0 / 0 | 0 / 0 |
| real states: largest relative difference host vs device | **3.3e-15** | **3.3e-15** |
| real states: values not bit-identical | 12.8 % | 12.8 % |
| real states: largest difference in units of ε/(1 − X_k) | 4.1 | 4.1 |
| synthetic states: largest relative difference | 2.6e-13 (X_H = 0.99916) | 2.6e-13 |
| synthetic states: largest difference in units of ε/(1 − X_k) | 4.4 | 4.4 |
| numpy vs device | 3.9e-13 (H) | 3.9e-13 (H) |
| md5 of the output file | ccd66f29 | ccd66f29 |

- The two device builds give bit-identical results.
- On the GPU, 12.8 % of the real-state values differ from the host in the last bits, because the compiler contracts to FMA.
- **Verdict:**
  - real states within 1e-14;
  - synthetic near-pure states within a few ε/(1 − X_k);
  - the device routine reproduces the host routine to round-off.
