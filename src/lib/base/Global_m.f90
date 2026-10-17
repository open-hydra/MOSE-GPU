module MOSE_Global_m
  use MOSE_Parameters_m

  implicit none

  character(clen) :: MOSE_phase_prefix = ''

  include 'mose_ct_dims.inc'   ! nsc/nu/nv/nw/np/nt/nprim/nrans = compile-time parameters
  integer :: nc, nsoot   ! Index of first soot variable, Number of soot variables
  integer :: npass=0,nps ! Number of Passive scalar variables
  !> Device-build upper bound for nprim-sized flux work arrays:
  !> nprim = nsc+4+nrans; device-supported ONERA-7 (nsc=7) + SST (nrans=2) -> 13.
  !> Fixed compile-time size avoids per-call device-heap malloc; bump if a
  !> bigger device mechanism/turbulence model is added.
  integer, parameter :: NPRIM_DEV = max(16, nprim)   ! one binary per species count: grows with nsc+4+nrans (16 up to nsc = 12)
  integer :: nres        ! Number of variables whose residuals are printed
  integer :: ndir        ! Number of dimensions of computational frame
  integer :: gc=2        ! ghost cells
  integer :: nrc         ! Number of reactions

  logical :: RANS = .false., VISC = .true.
  !> k-omega (SST, Wilcox2006) omega at a smooth wall, rho*omega_w = omega_wall_c(1) * mu / d^2 (MOSE_open Mod_RANS omega_wall_coef): 800 =
  !> 60/beta_1, Menter's practical form (default), or 80 = 6/beta_1, the asymptotic y -> 0 limit ([MOSE-RANS] omega-wall-bc).
  !> An array of one: nvfortran creates declare-create module arrays at start-up (not scalars); IO_RANS updates the device.
  real(8) :: omega_wall_c(1) = 8d2
  !$acc declare create(omega_wall_c)
  !> SST production ([MOSE-RANS] sst-production, MOSE_open Lib_SST): incompressible P = mu_t S^2 (SST-2003m, default) or
  !> compressible P = tau_ij du_i/dx_j (SST-2003e, the dilatation terms of OpenFOAM kOmegaSST). Same device pattern.
  logical :: sst_comp_c(1) = .false.
  !$acc declare create(sst_comp_c)
  logical :: CHIMERA = .false., WALL = .false.
  logical :: soot = .false.
  logical :: ROTATING_FRAME = .false.
  !> RANS/VISC are read inside the device flux kernels (`if (RANS)`/`if (VISC)`);
  !> make them device-resident. Synced in Set_RANS_Model alongside the devirt ids.
  !$acc declare create(RANS, VISC)
  !> The prim-vector indices are referenced directly inside the `acc routine seq`
  !> flux callees (Convective_Flux/Diffusive_Flux/SST). Device-resident; synced
  !> once in Allocate_Data after they are computed from the mechanism.
  !> gc (ghost-cell count) is referenced in the device BC twins' array bounds
  !> (1-gc:n+gc) — an acc routine using a module var requires it declare-create'd.
  !$acc declare create(gc)

  contains

    function Is_Wall(type) result (ans)
      implicit none
      integer, intent(in) :: type
      logical :: ans

      select case (type)
        case (5, 6, 8, 9, 10, 12, 13, 1000)
          ans = .true.
        case default
          ans = .false.
      end select
    end function Is_Wall


end module MOSE_Global_m
