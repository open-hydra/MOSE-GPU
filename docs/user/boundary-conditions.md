# Boundary Conditions

MOSE-GPU reads `INPUT/bc.txt` in the format of MOSE, the three-digit boundary types written by ATLAS. It also reads the
older format of the GPU code; the format is recognised from the first typed record. Each MOSE type is mapped at read time
onto a closure of the GPU code.

## Types on the GPU

| MOSE type | Boundary | Notes |
|---|---|---|
| 0, 100 | no boundary condition | 100 marks the blanked faces of a chimera overlap |
| 101, 201 | block connection | |
| 102 | chimera | the donors must be on the rank of the receiving block |
| 200 | rotational periodicity | |
| 300 | symmetry | |
| 301 | wall, prescribed heat flux | roughness and emissivity must be zero |
| 302 | wall, prescribed temperature | roughness and emissivity must be zero |
| 400 | extrapolation | |
| 401, 402 | total-pressure inlet | `alpha`, `beta` as angles or `normal` |
| 403 | mass-flux inlet with total temperature | |
| 404 | mass-flux inlet with static temperature | `rel_fac = 1` only |
| 405 | supersonic inlet | |
| 406 | back-pressure outlet | |
| 407 | inflow or outflow by the sign of the normal velocity | |
| 420 | nozzle inlet | |
| 501 | manifold connection | one rank only |

## Types that stop the run

| MOSE type | Boundary |
|---|---|
| 103 | coupled multi-solver wall |
| 303, 304, 305–309 | walls with radiation, pyrolysis and the other wall models |
| 408 | prescribed-velocity inlet |
| 409, 410 | other inlets and mapped states |
| 502 | solid-rocket-motor grain |
| 503–506 | ablating and gas–surface-interaction walls |
| any inlet | inside a rotating frame |

## Differences inherited from the GPU closures

These differences are known and do not appear in any test of MOSE:

- **Over-pressurised stagnation inlet.** With `rel_fac = 1`, MOSE-GPU closes the inlet for that iteration. MOSE iterates
  the relaxation up to 1000 times, then sets the discriminant to zero.
- **Subsonic outflow.** The guard on the exit pressure uses the total Mach number in MOSE-GPU and the face-normal Mach number
  in MOSE.
- **Mass-flux inlets with outward flow.** MOSE-GPU keeps injecting; MOSE closes the face.

The inflow direction `normal` uses the same formula as MOSE, including its defect on faces without an x component: see
[Known Differences](../vv/known-differences.md).

## Checking a boundary file

Three scripts in `scripts/` check a boundary file without running the solver:

| Script | Purpose |
|---|---|
| `bc_open_lint.py` | walk a MOSE `bc.txt` with the record structure of the reader, and report its types and columns |
| `bc_chimera_check.py` | check that the chimera weights sum to one and that every donor exists |
| `bc_legacy_to_open.py` | convert a file of the old GPU format into the MOSE format |
