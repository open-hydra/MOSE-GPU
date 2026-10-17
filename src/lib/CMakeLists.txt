# Source layout (MOSE_open taxonomy; see doc/SOURCE_MAP.md for the full per-file map):
#   base/ driver/ io/ diagnostic/ parallel/ physics[/turbulence]/
#   numerics/{fluxes[/bc],space,time[/explicit,/implicit],multigrid}/
#   gpu/  -- the var-last device port (GPU production path): W twins, CUF flux
#            kernels and their physics .inc bodies.
# All of it is compiled into the same library/TU set -- purely organizational.
file(GLOB_RECURSE MOSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
add_library(${lib} ${MOSE_SRC} ${FLINT_SRC})
target_link_libraries(${lib} PRIVATE FiNeR::FiNeR ORION)
link_FLINT_dependencies(${lib})

target_include_directories(${lib} INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
# Fortran INCLUDE search paths: the generated mose_ct_dims.inc lands in the build dir
# (never in the source tree); gpu/ holds the CUF physics .inc bodies included from
# gpu/*_cuf.f90 (same-directory resolution also works, the -I is belt and braces).
target_include_directories(${lib} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/gpu)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${lib} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(${lib} PUBLIC USE_MPI)
    target_link_libraries(${lib} PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)

# GPU Phase 2 / Milestones 1+2: run the CONVECTIVE (M1, Mod_Convective_W) and DIFFUSIVE
# The var-last flat SoA w pipeline (fhat two-pass gather-difference kernels, full
# device residency) IS the GPU path -- HARDWIRED to _OPENACC since 2026-07-11 (the
# old MOSE_FIELD_W flag is gone; source guards use _OPENACC directly). CPU builds (no
# OpenACC) compile the original var-first AoS P path (Mod_Convective/Mod_Diffusive,
# host Fluxes_blk) -- generality across turbulent/laminar/Euler and CPU/GPU is
# runtime + USE_OPENACC, not a tuning flag.

# GPU Phase 2 / M2: diffusive face thermo comes from the RESIDENT CELL CACHE
# (mil_c/kl_c/cp_c + rho_c/T_c, face-averaged) -- UNCONDITIONAL since 2026-07-11.
# Validated THOR solve -9.7% (2026-06-24). The per-face Wilke RECOMPUTE path
# (Diffusive_Flux_w_blk) was DELETED: it had been broken since Res var-last Stage 4
# (silently dropped diffusive fluxes; doc/facex_closure_REPORT_2026-07-07.md 6.1)
# and required a 32 KB CUDA stack. Git history keeps it (pre-2026-07-11).

# Face-exact diffusive closure is the SOLE closure (facex-hardwire pass, 2026-07-17).
# The CUF diffusive kernels always read the tface buffer (Diffusive_Face_Thermo_Fill),
# which recomputes rho/T/mil/kl/cp from the FACE-AVERAGED primitive state -- the CPU
# original ("Option B", mu(1/2(s1+s2))). The old cell-cache-closure path
# (1/2(mu(s1)+mu(s2))) and the MOSE_DIFF_FACE_EXACT compile flag were removed; no source
# references the macro any more, so there is no option here.

# CUF flux kernels: convective + diffusive pass-1 as CUDA-Fortran attributes(global)
# kernels (launch_bounds) -- the PRODUCTION GPU flux path (bit-exact vs the old OpenACC
# scatter, ~21.5% faster). Now UNCONDITIONAL on the var-last GPU path (no user flag): the
# Mod_Convective_W_CUF / Mod_Diffusive_W_CUF modules are gated by _OPENACC.
# They contain CUDA-Fortran, so -cuda is added library-wide (coexists with -acc, proven).
# The diffusive CUF reads the resident thermo cache (unconditional).
IF(USE_OPENACC)
    target_compile_options(${lib} PUBLIC -cuda)
    target_link_options(${lib} PUBLIC -cuda)
ENDIF()

# One binary per mechanism: the species and RANS-variable counts are compile-time parameters (the generated
# mose_ct_dims.inc); the mose launcher sets MOSE_NSC / MOSE_NRANS from the case.
SET(MOSE_NSC 4 CACHE STRING "Chemical species count (compile-time)")
SET(MOSE_NRANS 2 CACHE STRING "RANS variable count (compile-time)")
MATH(EXPR MOSE_NU "${MOSE_NSC}+1")
MATH(EXPR MOSE_NV "${MOSE_NSC}+2")
MATH(EXPR MOSE_NW "${MOSE_NSC}+3")
MATH(EXPR MOSE_NP "${MOSE_NSC}+4")
MATH(EXPR MOSE_NT "${MOSE_NSC}+5")
MATH(EXPR MOSE_NPRIM "${MOSE_NSC}+4+${MOSE_NRANS}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mose_ct_dims.inc.in ${CMAKE_CURRENT_BINARY_DIR}/mose_ct_dims.inc @ONLY)

# The chemistry chain of FLINT specialized to the mechanism: the radau5_dev ODE size + tolerance-transform hoist
# (MOSE_CHEM_NZ), the same-module inlined chemistry chain (MOSE_CHEM_INLINE), its literal species trip counts
# (MOSE_CHEM_NS) and the thermo/transport trip counts (MOSE_TT_NS). Guards against a mismatch: the launcher stops on
# nsc+1 /= MOSE_CHEM_NZ, FLINT checks the table upload.
IF(USE_OPENACC)
    MATH(EXPR MOSE_NZ "${MOSE_NSC}+1")
    target_compile_definitions(${lib} PUBLIC MOSE_CHEM_NZ=${MOSE_NZ})
    target_compile_definitions(${lib} PUBLIC MOSE_CHEM_INLINE)
    target_compile_definitions(${lib} PUBLIC MOSE_CHEM_NS)
    target_compile_definitions(${lib} PUBLIC MOSE_TT_NS=${MOSE_NSC})
    # nvfortran 24.11: the HOST inliner (-Minline) on the merged chemistry module ICEs
    # the NVVM device backend (bisected 2026-07-02, job 293908). Compile that one file
    # without it; the same-TU NVVM inliner still inlines the _m chain on device.
    set_source_files_properties(${CMAKE_SOURCE_DIR}/lib/FLINT/src/lib/Lib_Radau5_dev.f90
                                PROPERTIES COMPILE_OPTIONS "-Mnoinline")
ENDIF()


# GPU chemistry integrator selection. Default = radau5_dev (RADAU5). MOSE_ODE_RODAS=ON
# swaps in rodas_dev (device RODAS4 Rosenbrock) via the ODE_DEV macro in Lib_Chemistry
# -- compile-time (like the mechanism binaries) so only the chosen integrator is
# compiled into the cell kernel (its register/stack footprint alone -> the occupancy
# lever; RODAS carries ~half RADAU5's per-thread matrix storage). V0-validated
# faithful (host port vs Hairer RODAS ~1e-10; RODAS vs RADAU5 ~3e-8).
OPTION(MOSE_ODE_RODAS "GPU chemistry uses rodas_dev (device RODAS4) instead of radau5_dev" OFF)
IF(MOSE_ODE_RODAS AND USE_OPENACC)
    target_compile_definitions(${lib} PUBLIC MOSE_ODE_RODAS)
ENDIF()


# BC boundary-flux: opposite face pairs {1,2}{3,4}{5,6} fused into one launch each
# (6->3 per block, ~2x SM fill) -- UNCONDITIONAL since 2026-07-11 (bit-exact,
# +1.47% solve; per-face fallback for degenerate n_d=1 stays as a runtime guard).

# Per-block BC kernels on concurrent async queues -- UNCONDITIONAL since 2026-07-11
# (bit-exact across schedules by construction; NQ default 2 = measured A30 plateau,
# runtime env MOSE_BC_NQ overrides).


# Device-side residual on the resident var-last w/wo -- UNCONDITIONAL since
# 2026-07-11 for the nres=7 RANS layout (field bit-exact; probes refresh from w).
# Other layouts (laminar/Euler, nres/=7) take the HOST path via the RUNTIME guard
# in Mod_Explicit (res_needs_host_P) -- generality preserved, not a compile switch.

