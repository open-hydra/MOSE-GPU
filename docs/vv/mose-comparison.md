# Comparison with MOSE

Every test of MOSE commit `b29130b`, run on both codes from the same input. *Difference* is the largest relative L2
difference over the variables of the final field, with that variable in brackets; *steps* are those of MOSE-GPU and MOSE.

## Tests of MOSE

| Test | Difference | Steps | Verdict |
|---|---|---|---|
| 1D Sod shock tube | 4.9e-15 (u) | 195 / 195 | agree |
| 1D Einfeldt | 3.0e-14 (p) | 184 / 184 | agree |
| 1D Noh | 5.2e-10 (T) | 821 / 821 | agree: 4 cells of 400 at the shock fronts, at 5e-4 K below the thermodynamic tables |
| 1D Shu–Osher | 6.0e-15 (p) | 336 / 336 | agree |
| 1D Toro, HLLE++ | 3.7e-2 (T) | 244 / 241 | differ: a defect of MOSE, see [Known Differences](known-differences.md) |
| 1D unity-Lewis diffusion | 4.2e-9 (u) | 452 873 / 452 873 | agree |
| 1D multicomponent diffusion | 1.2e-8 (u) | 438 639 / 438 639 | agree |
| 2D Chang nozzle | 3.7e-15 (T) | 4698 / 4698 | agree |
| 2D Gresho vortex | 1.1e-8 (u) | 16 959 / 16 959 | agree |
| 2D hypersonic cylinder | 2.4e-2 (T) | 5000 / 5000 | not steady, see [Hypersonic Cylinder](hypersonic-cylinder.md) |
| 2D isentropic vortex | 1.6e-15 (u) | 369 / 369 | agree |
| 2D oblique shock | 1.9e-14 (ρ) | 3000 / 3000 | agree |
| 2D Prandtl–Meyer expansion | 3.2e-15 (p) | 385 / 385 | agree |
| 2D ramp channel | 7.9e-14 (T) | 10 000 / 10 000 | agree; the flow does not reach a steady state in either code |
| 2D rocket nozzle, frozen | 2.6e-15 (ρ) | 3103 / 3103 | agree |
| 2D supersonic forward step | 4.4e-14 (p) | 1658 / 1658 | agree |
| 2D laminar flat plate | 7.8e-8 (T) | 100 000 / 100 000 | agree: on 3 of the 64 cells of the ambient faces the normal velocity is zero to round-off, and the two codes pick inflow or outflow differently |
| 2D turbulent flat plate, Spalart–Allmaras | 2.8e-9 (v) | 200 000 / 200 000 | agree |
| 2D turbulent flat plate, SST | 1.9e-13 (v) | 100 000 / 100 000 | agree |
| 2D turbulent flat plate, Wilcox 2006 | 5.3e-14 (ρk) | 100 000 / 100 000 | agree |
| 2D manufactured solution, without the source | 1.0e-15 (T) | 9400 / 9400 | agree |
| 2D shock / boundary-layer interaction, SA | 3.6e-4 (p) | 40 000 from the same field | agree within the sensitivity of the flow |
| 2D shock / boundary-layer interaction, SST | 1.4e-4 (p) | 40 000 from the same field | agree within the sensitivity of the flow |
| 3D pressure in a centrifugal field | 9.4e-11 (v) | 2000 / 2000 | agree |

The seven tests that stop on the GPU: Fer14 and the premixed flame (mechanism, integrator, grid sequencing), the rocket
nozzle with equilibrium and with finite-rate chemistry (equilibrium, mechanism), the ablating wall, the Coriolis channel
(inflow in a rotating frame) and HyShot II (mechanism, grid sequencing).

!!! note "Choices shared by both codes"
    - Einfeldt, Noh and Toro run on the 400-cell grid, as Sod.
    - The isentropic vortex and the manufactured solution run on the 64 × 64 grid.
    - The manufactured solution runs without its source term, which needs a separate driver of MOSE.
    - The nozzles stop at 5000 iterations, the centrifugal case at 2000.
    - The shock / boundary-layer interaction runs the coarse grid levels on MOSE, which has grid sequencing, then both codes
      for 40 000 iterations on the fine grid from the same field.
    - The SST variant uses `point-implicit = true`: as shipped, MOSE itself stops on it at iteration 300 of the coarsest level.

## The shock / boundary-layer interaction

The flow does not converge in either code: the interaction oscillates and the residual stays at a plateau. The difference
between the codes is the same size as the difference of each code from itself when the CFL number changes by one part in
10^12:

| | GPU against MOSE | GPU against itself | MOSE against itself |
|---|---|---|---|
| SA, pressure | 3.6e-4 | 3.6e-4 | 4.5e-4 |
| SST, pressure | 1.4e-4 | 2.6e-4 | 1.6e-4 |

## Riemann solvers on the Sod shock tube

Each solver of MOSE on the Sod tube, changing only the solver. The lateral faces of the test are one cell wide; they are
run as in the test, with extrapolation, and with symmetry.

| Solver | Lateral symmetry faces | Lateral extrapolation faces |
|---|---|---|
| HLLE++ | 2.4e-15 | differs: [Known Differences](known-differences.md), item 1 |
| HLLC+Tramel | 8.6e-15 | differs: item 1 |
| HLLC+Chen | 1.1e-14 | MOSE gives NaN at step 1: item 1 |
| LLF, Rusanov | 1.6e-14 | 1.6e-14 |
| SLAU2 | 4.0e-15 | 4.0e-15 |
| exact | 1.7e-14 | 9.4e-15 |
| MiczekRoe, M_co = 0.5 | 1.2e-13 | 2.2e-13 |
| LMRoe | differs: a defect of MOSE, item 3 | differs: item 3 |
| AUSM+, AUSM+M | both codes stop at step 1, in the same cell with the same state | the same |

With M_co = 0.005, MiczekRoe stops at step 1 in both codes, in the same cell with the same state: with the fluid at rest
the preconditioning is at its maximum on a 10:1 pressure jump.

With the one-cell ghost fill of MOSE corrected in a test copy, HLLE++, HLLC+Tramel and HLLC+Chen with extrapolation faces
agree with the GPU to 1e-14 or better, in the same number of steps. With `u2Roe` assigned in a test copy of MOSE, LMRoe
agrees to 1e-14 on the isentropic vortex and 4e-10 at the Sod shock.

## Integration variables and turbulence models

The options no test of MOSE uses, each on a test of the matrix with only the keys changed:

| Variant | What it exercises | Difference | Steps |
|---|---|---|---|
| `prim`, isentropic vortex, 1 step | the primitive residual | 5.9e-16 | 1 |
| `prim`, rocket nozzle, 13 species, 100 steps | the species rows | 1.7e-15 | 100 |
| `prim`, SST flat plate, 1 step | the turbulence rows | 1.4e-15 | 1 |
| `prim`, centrifugal case with residual smoothing, 1 step | the smoothed residual | 8.2e-13 | 1 |
| `prec` + HLLC-PC, Gresho vortex, 2000 steps | the reference velocity, the time step, the preconditioned flux | 4.2e-12 | 2000 |
| `prec` + HLLC-PC, isentropic vortex, 2000 steps | the periodic faces | 2.1e-15 | 2000 |
| `prec`, laminar flat plate, 1 step | the viscous floor of the reference velocity | 3.5e-16 | 1 |
| `prec`, SST flat plate, 1 step | the eddy viscosity in that floor | 2.0e-15 | 1 |
| SA-RC, shock / boundary-layer interaction, 200 steps | rotation and curvature correction | 6.8e-11 | 200 |
| SST-RC, shock / boundary-layer interaction, 200 steps | rotation and curvature correction | 4.9e-14 | 200 |
| SAcomp, shock / boundary-layer interaction, 200 steps | compressibility correction | 6.0e-15 | 200 |
| Wilcox 2006 with `k-energy-coupling`, flat plate, 1 step | the k term in the energy | 5.9e-16 | 1 |

The largest SA-RC difference sits in the free stream near the inlet, where the turbulence variable is about 1e-6. There,
strain and vorticity are zero to round-off, and the correction is a ratio of two round-off numbers. Where the variable is
above a thousandth of its maximum, the two codes agree to 5e-12.

## Hypersonic cylinder at first order

With `space-reconstruction = first-order` and no flux limiter, both codes converge, stop at the same iteration (6748)
on the residual threshold, and agree to 1.1e-13.
