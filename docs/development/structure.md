# Code Structure

MOSE-GPU is laid out as MOSE (commit `b29130b`), with the device port folded into the modules of MOSE. One physics, one execution path: every kernel runs on the device
(OpenACC or CUDA Fortran); the host reads the input, builds the mesh and the boundary tables, owns the block arrays
at the I/O boundaries and writes the output. Build options: [Build Flags](build-flags.md); validation against MOSE:
[Verification & Validation](../vv/index.md).

## The state: block arrays and their flat twins

MOSE keeps the primitive state in per-block arrays, `blk(b)%P(nprim,i,j,k)` (with `PO` the state at the start of
the step and `R` the residual). The device works on a flat, var-last copy of the same numbers: `P_flat(cell, var)`,
`PO_flat`, `R_flat` (`src/lib/gpu/MOSE_Field_m.f90`), the cells of every local block concatenated with their ghost
layers, so that neighbouring threads read neighbouring doubles. `P` is the truth at setup and at the I/O boundary of
a step (restart, probes, wall files, residual); `P_flat` is the truth inside the step. `pack_P_to_flat` /
`unpack_flat_to_P` convert; `field_update_P_flat_device` / `_host` move the flat field. The per-cell thermodynamic
caches (`rho_c`, `T_c`, `Rg_c`, `mil_c`, ...) have the same two forms (`blk(b)%rho_c` and `rho_c_f`).

The kernels that receive the flat arrays as arguments keep the dummy names `w`, `wo`, `wr` of the port; the call
sites pass `P_flat`, `PO_flat`, `R_flat`.

## src/lib — the solver (the directories of MOSE)

**driver/** — `Wrap_Setup.f90` (reads everything, `Setup_Device_Residency` creates the device copies, the setup
ghost fill `Fill_Ghost_Cell_Setup`), `Wrap_Solve.f90`, `Wrap_Postprocess.f90`, `Mod_Allocate_Data.f90`,
`Procedures_m.f90`.

**base/** — `Global_m.f90` (the compile-time dimensions `mose_ct_dims.inc`: species, RANS variables),
`MOSE_Types_m.f90` (block and domain types, the configuration), `General_Types_m.f90`, `Parameters_m.f90`.

**numerics/time/** — `Mod_dt.f90` (`Compute_dt` and its device kernel `Compute_dt_blk`, `Set_Global_dt`),
`explicit/Mod_Explicit.f90` (the step: Strang splitting, RK substages, the halo, the thermo cache, the fluxes, the
boundary fluxes, the state update, chemistry once per step; read this to understand an iteration),
`explicit/Lib_Newstate.f90` (`Newstate_Conservative` / `Newstate_Primitive` and their per-block kernels, the state
check), `explicit/Lib_IRS.f90` (residual smoothing: `IRS_apply` on the device residual), `explicit/Lib_Preconditioning.f90`
(integration-variables = prec: the reference velocity Ur of the cells, `Update_Derived_Variables`), `explicit/Lib_RK.f90`,
`explicit/Lib_Strang.f90`.

**numerics/fluxes/** — `Mod_Fluxes.f90` (the flux driver: `Fluxes` → `Fluxes_blk` per block, the fused pass 2 that
accumulates the face fluxes into `R_flat`), `Mod_Convective.f90` (the CUDA-Fortran convective kernels, one per
direction, with the reconstruction, the limiter and the Riemann solver inlined from `gpu/conv_physics.inc`;
`Convective_Flux` launches them), `Mod_Diffusive.f90` (the same for the diffusive flux, `gpu/diff_physics.inc`;
`Diffusive_Face_Thermo_Fill` evaluates the transport properties at the faces; `Compute_Diffusive_Flux` is
MOSE's per-face routine, used by the viscous connection BC), `Lib_Riemann.f90` (`Assign_RS`: the Riemann
solvers the device has — HLLC Batten, HLLE, SLAU, AUSM+, AUSM+M, HLLE++, HLLC+Tramel, LLF (Rusanov), SLAU2, HLLC+Chen,
exact, LMRoe, MiczekRoe, HLLC-PC (with integration-variables = prec): every solver of MOSE — and the device dispatch
`Riemann`), `Lib_Riemann_HLL.f90`, `riemann_open_body.inc` (eleven of MOSE's solvers, ids 4..14), `Lib_Shock_Detector.f90`.

**numerics/fluxes/bc/** — `Mod_BC.f90` (`BC_Fluxes`: one kernel per block over its boundary records,
`bc_record_body.inc` / `bc_record_apply.inc`; the manifold update), `Lib_BC.f90` (the boundary conditions
themselves, one routine per type, `!$acc routine`), `Lib_Extrapolation.f90` (`Fill_Ghost_Cell`, MOSE's entry:
the ghost cells of every local block — walls, symmetry, extrapolation, rotational periodicity, connections, chimera
— on one rank by `Fill_Ghost_Cell_local`, on more by `Fill_Ghost_Cell_mr` around the inter-rank exchange;
`Extrapolate_Wall`), `Mod_BC_Q2D.f90`, `Mod_BC_Time.f90`.

**numerics/space/** — `Lib_Metrics.f90`, `Lib_Limiters.f90` (`Assign_Limiter`, the device dispatch `rlimiter`),
`Lib_Interpolation.f90`, `Lib_KDTree.f90`; **numerics/multigrid/** — `Lib_Multigrid.f90`.

**physics/** — `Lib_Chemistry.f90` (`Chemistry_Newstate_blk`: the stiff integration of every cell on the device,
RADAU5 or RODAS of `lib/FLINT`), `Mod_Thermo_Cache.f90` (the per-cell caches and their kernels), `Lib_RotatingFrame.f90`,
`Lib_Soot.f90`, `Mod_GSI.f90`; **physics/turbulence/** — `Mod_RANS.f90` (`Set_RANS_Model`), `Lib_RANS.f90`,
`Mod_SST.f90` (`SST_Blk` + `gpu/sst_body.inc`), `Mod_Spalart.f90` (`SA_Source_Blk` + `gpu/sa_body.inc`),
`Mod_Wilcox2006.f90` (`Wilcox_Blk` + `gpu/wilcox_body.inc`), `Lib_SpalartShur.f90` (the velocity gradient and the rotation/curvature terms of SA-RC and SST-RC; `SST_RC_Blk`, `SA_RC_Source_Blk` are the -RC twins), `Mod_SSGLRR.f90`, `Mod_Fluid.f90`, `Mod_QCR.f90`, `Lib_Omega_Relax.f90` (the point-implicit
omega source, `omega-source-relax`).

**parallel/** — `Mod_MPI.f90` (ranks, the block → rank map, the device binding), `Mod_BlockSplit.f90`,
`Mod_GhostExchange.f90` (the inter-rank halo on the device: pack / unpack kernels, persistent requests, the
transverse ghost layers `Pg`).

**io/** — `IO_Ini.f90` (the GPU input dialect; `Stop_On_Removed_Keys`: multirate, j-implicit and the other
features this version does not have stop the run), `IO_Ini_Open.f90` (the MOSE input dialect, translated
key by key), `IO_BC.f90` / `IO_BC_Open.f90` (the two bc.txt dialects), `IO_Chemistry.f90`, `IO_RANS.f90`,
`IO_Solution.f90` (ascii and szplt fields), `IO_ThermoTransport.f90`, `IO_Wall.f90`, `Mod_Probes.f90`.

**diagnostic/** — `Mod_Diagnostic.f90` (the residual norms, on the device).

## src/lib/gpu — what MOSE does not have

`MOSE_Field_m.f90` (the flat field, above) and the included kernel bodies: `conv_decls.inc` / `conv_physics.inc`
and `diff_decls.inc` / `diff_physics.inc` / `diff_cuf_kernels.inc` (the CUDA-Fortran flux kernels, instantiated
once per direction and, for the diffusive flux, once more for the multicomponent-diffusion closure),
`newstate_cons_body.inc`, `dt_blk.inc`, `pass2_body.inc`, `sst_body.inc`, `sa_body.inc`.

## Outside src

| Path | Content |
|---|---|
| `lib/FLINT` | thermodynamics, transport and chemistry, with the device code: RADAU5 and RODAS in `Lib_Radau5_dev.f90`, the thermodynamic and transport routines in `Lib_ThermoTransport_dev.f90`; its ODE library `lib/FLINT/lib/OSlo` |
| `lib/ORION` | I/O, including TecIO |
| `lib/third_party/FiNeR` | INI files (Git submodule) |
| `mose` | the launcher: reads the case, builds the binary the case needs, submits it |
| `mose_bind.sh` | per-rank placement of the CPU cores next to each GPU |
| `install.sh`, `CMakeLists.txt`, `cmake/` | the build |
| `test/` | the regression table on the test cases of MOSE: see [Testing](testing.md) |
| `scripts/` | boundary-file checks, Tab completion of `mose` (`mose_completion.sh`) |
| `docs/` | this documentation |
