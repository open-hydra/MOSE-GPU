# Known Differences

The places where MOSE-GPU and MOSE commit `b29130b` do not give the same result, or would not on some input. Each item was
checked in two ways:

- **Reading the code.** The defect is located in MOSE `b29130b`.
- **A run that exposes it.** A test copy of MOSE sets a chosen value in one variable, or applies the proposed fix, next to a
  control run where the effect must not appear.

Without the switch, the test copy gives the same bits as the reference binary on the Wilcox flat plate. Line numbers refer
to MOSE `b29130b`.

| | Item | Effect | Priority |
|---|---|---|---|
| 1 | ghost cells of extrapolation faces in a one-cell direction | thread-dependent results; NaN with HLLC+Chen | high |
| 2 | inflow direction `normal` on faces without an x component | tangential inflow velocity up to 95 % of the normal one | high |
| 3 | `u2Roe` never assigned in LMRoe | results depend on undefined memory | high |
| 4 | six entries of Ŝ never assigned in Wilcox 2006 | undefined values read in 3-D or in a rotating frame | medium |
| 5 | reference velocity of the preconditioning left at zero in ghost cells | HLLC-PC uses zero at block connections | medium |
| 6 | `first-order` with a `flux-limiter` key stays second order | the requested scheme is silently not applied | medium |
| 7 | `level1-iter` replaces `iter-threshold` with grid sequencing | `iter-threshold` ignored without warning | low |
| 8 | eddy viscosity of the preconditioning floor without the velocity gradient | no effect measured | low |
| 9 | no species-diffusion limit on the time step | not tested | to assess |

## 1. Ghost cells of extrapolation faces in a one-cell direction

**Where.** Type 400 falls in the default branch of `Fill_Ghost_Cell`, which is `Ghost_Extrapolate` in
`numerics/space/Lib_Ghost.f90`. That routine uses P0 = 3P1 − 3P2 + P3. The state, ghost cells included, is zeroed in
`driver/Mod_Allocate_Data.f90:363`.

**What happens.** In a direction with one cell, P2 and P3 are the ghost cells of the opposite face. With one thread the
far face keeps its initial values G2, G3, which are zero, and the near face gets 3c − 3G2 + G3 = 3c, then 6c, at every fill.
With several threads the order of the two faces changes from column to column, under `schedule(dynamic, 64)`.

**Evidence.** Sod tube with HLLC+Chen and the lateral faces of the test:

- **Ghost pressures.** A print in the Chen sensor shows 0.3 on one side and −2.8e-17 on the other, next to p = 0.1, in all
  200 cells of the right state. The sensor h = min p_C/p_k is −2.8e-16, and g = 1 − h^M (`Lib_Riemann_HLL.f90:206`) is
  not defined for h < 0.
- **The NaN.** The state check stops at the first bad cell of each thread: cell 201 with one thread, 201 and 301 with
  four, 201, 251, 301 and 351 with eight.
- **Thread dependence.** Two identical runs of the Toro test with eight threads differ by 1.5e-5 in pressure.

**The fix, tested.** With zero-gradient ghosts in directions with fewer than three cells, MOSE becomes deterministic: one
and eight threads give the same bits. It also agrees with MOSE-GPU, whose fill gives the same ghost values:

| Case | Steps, MOSE / GPU | Difference from the GPU |
|---|---|---|
| Sod with HLLE++ | 194 / 194 | 1.1e-14 |
| Sod with HLLC+Tramel | 194 / 194 | 2.7e-14 |
| Sod with HLLC+Chen | 194 / 194, no NaN | 3.9e-14 |
| Toro, as in the test | 244 / 244 | 6.0e-14 |

**Proposed fix in MOSE.** In `Ghost_Extrapolate`, use `Ghost_ZG_Extrapolate` when the face direction has fewer than three
cells.

## 2. Inflow direction `normal` on faces without an x component

**Where.** `Compute_Inflow_Direction` in `numerics/fluxes/bc/Lib_BC_Fluxes_Inflow.f90:559-576` sets
alpha = atan(n2/(n1 + 1e-20)) and beta = atan(n3/(n1 + 1e-20)). The tangential components E1 and E2 then follow from
tan(alpha) and tan(beta), at lines 121-128, 205-210 and 284-289.

**What happens.** With n1 = 0 both angles saturate at π/2, and tan is 1.63e16 for both. The ratio of n2 to n3 is lost,
and E1, E2, which must be zero for an inflow along the normal, come out of order one.

**Evidence.** The same formulas in double precision, on cylindrical faces around the x axis, normal (0, cos θ, sin θ):

- **Axial face edge.** E2 reaches 0.95.
- **Circumferential face edge.** E1 reaches 0.95.
- **Control.** A face with n1 = 0.6 gives 1e-17.

**Consequence.** An inflow `normal` on faces parallel to the x axis enters with a tangential velocity up to 95 % of the
normal one, varying from face to face. MOSE-GPU has the same formula, so the two codes agree, and both are wrong.

**Proposed fix.** Build the direction from the unit normal: for `normal`, the tangential components are zero, without
going through the angles. To be fixed in both codes.

## 3. `u2Roe` never assigned in LMRoe

**Where.** In `numerics/fluxes/riemann/Lib_Riemann_Roe.f90`:

- `u2Roe` is local to `Riemann_LMRoe`, declared at line 47;
- it is read at line 120, in the tangential speed, and at line 185, in the entropy eigenvector;
- it is never assigned;
- the comment at line 164 says "u2Roe computed above";
- `Riemann_MiczekRoe`, in the same file, assigns it as |u_Roe|² at line 371.

**Evidence.** Isentropic vortex, one step, LMRoe with M_co = 0.005:

- **The value matters.** Setting `u2Roe` = 1e3 instead of the undefined value changes the field by 23 % on T.
- **A large value.** Setting it to 1e6 stops the run at step 1.
- **Control.** The same line in MiczekRoe, where `u2Roe` is assigned before use, leaves the field identical to the bit.
- **Two builds.** The reference binary and the test copy without the switch are two builds of the same code. They differ
  by 0.8 % on T after one step.

**Proposed fix.** `u2Roe = uRoe*uRoe + vRoe*vRoe + wRoe*wRoe` after `unRoe`. With that line in a test copy of MOSE, the two
codes agree to 1e-14 on the vortex and 4e-10 at the Sod shock. MOSE-GPU assigns it so.

## 4. Six entries of Ŝ never assigned in Wilcox 2006

**Where.** In `physics/turbulence/Lib_Wilcox2006.f90`:

- `Sij_hat(3,3)` is local, at line 73;
- only its diagonal is assigned, at lines 141-143;
- the loop of the invariant χ_ω reads all nine entries, at line 154.

**What happens.** The off-diagonal entries weigh ω_i ω_k with i ≠ k, so they matter as soon as the vorticity has two
components: in 3-D, or in a rotating frame whose axis is not normal to the plane.

**Evidence.** Wilcox flat plate, one step:

- **Plain 2-D.** Values 1 and 1e3 in the unassigned entries give fields identical to the bit: they weigh zero there. NaN
  stops the run: they are read anyway.
- **Rotating frame with axis (1, 1, 0)/√2.** Values 1 and 1e3 give different fields, 2e-4 on ρk after one step. MOSE as it
  is differs from the fixed version by 5e-7 on ρk.

**Proposed fix.** `Sij_hat = Sij` before the three diagonal lines. MOSE-GPU uses all nine entries.

## 5. Reference velocity of the preconditioning left at zero in ghost cells

**Where.**

- **Zeroing.** `blk%Ur` is zeroed with its ghost cells (`Mod_Allocate_Data.f90:347`).
- **Computation.** It is computed in the interior cells only (`Lib_Preconditioning.f90:210`), and no ghost fill touches it.
- **Use.** At a block connection the Riemann solver receives the Ur of the ghost cell (`Lib_BC_Fluxes_Connection.f90:82`).

**Evidence.** HLLC-PC with preconditioning, one step:

- **Periodic isentropic vortex, 256 connection faces.** Ghost values 0, 1 and 10 give three different fields.
- **Control: Prandtl–Meyer, no connection faces.** The field is identical to the bit for any value.

**Consequence.** HLLC-PC uses V_r = min(Ur, a) = 0 on the ghost side, so connection and periodic faces do not behave like
interior faces, and the result depends on how the grid is cut into blocks. MOSE-GPU reproduces this, to give the same
results as MOSE.

**Proposed fix.** Fill the ghost cells of Ur with the exchange of the primitive variables, or recompute it there from the
filled ghost state.

## 6. `first-order` with a `flux-limiter` key stays second order

**Where.** `Setup_Space_Scheme` in `numerics/space/Mod_Space.f90:19-33`. The first-order branch warns "Flux limiter
specified but MUSCL reconstruction is not enabled. Ignoring flux limiter.", then calls `Assign_Limiter` with the limiter of
the input.

**Evidence.** Hypersonic cylinder with `first-order` and `flux-limiter = vanalbada`: the log says "First-order", but the
residual history equals, digit for digit for 5000 iterations, the one of MUSCL with van Albada. Without the limiter key
the run is first order and converges.

**Consequence.** No test of MOSE is affected: Noh and the Gresho vortex use first order without a limiter key. MOSE-GPU
runs first order, as the warning says. With this combination of keys the two codes differ until MOSE is fixed.

**Proposed fix.** Call `Assign_Limiter('none')` in the first-order branch.

## 7. `level1-iter` replaces `iter-threshold` with grid sequencing

**Where.**

- `driver/Wrap_Setup.f90:57` sets `domain(1)%itermax = iter_threshold`.
- `Setup_Multigrid`, called later, sets `domain(m)%itermax` from the level keys for every level. That includes level 1,
  the fine grid (`numerics/multigrid/Mod_Multigrid.f90:52-55`).

**Consequence.** With more than one level, `iter-threshold` is ignored without warning. If this is intended, the input
reference should say so. Grid sequencing is not on the GPU.

## 8. Eddy viscosity of the preconditioning floor without the velocity gradient

**Where.**

- The viscous floor of Ur computes the eddy viscosity with `blk%vel_gradient` (`Lib_Preconditioning.f90:176-182`).
- Only `Compute_Velocity_Gradient` writes that array (`physics/turbulence/Lib_SpalartShur.f90:34`), and only the SA and
  SST sources with the RC corrections call it.
- Otherwise the gradient stays zero: for SST and Wilcox the floor uses ρk/ω, without the strain limiter.

**Evidence.** SST flat plate with preconditioning: 300 steps as is and with the gradient computed before Ur give the same
bits.

**Consequence.** An inconsistency without a measured effect on this case. Low priority.

## 9. No species-diffusion limit on the time step

`Compute_Dt` (`numerics/time/Mod_dt.f90`) has the CFL limit and the viscous limit, with the laminar plus turbulent
viscosity. With mixture-averaged diffusion, a light species such as H₂ diffuses much faster than momentum: its effective
Schmidt number is about 0.2. In that case the viscous limit alone may not keep the explicit species diffusion stable. Not
tested.

The native dialect of MOSE-GPU has such a limit (`vnn-diff`); with a MOSE input file it is switched off, so that the time
step is the one of MOSE.

## Not differences

- **Ablating wall test.** Its `INPUT` and `MESH` are not in the repository because `build_case.py` builds them, as its
  README says.
- **One-dimensional tests.** The `solve` block of their `MOSE.sh` links `phase.txt` and `thermo.dat` only; grid, initial
  state and boundary conditions are linked by the `test` block, one resolution at a time. `MOSE.sh solve` does not start
  on a clean checkout. This is an inconvenience, not a computational error.
