module MOSE_Lib_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
  !> Concrete procedures pointing to one of the subroutine realizations.
  procedure(Eddy_Viscosity_if), pointer, public :: Eddy_Viscosity
  procedure(RANS_diffusive_flux_if), pointer, public :: RANS_Diffusive_Flux
  procedure(RANS_source_terms_if), pointer, public :: RANS_Source_Terms
  procedure(RANS_Set_Wall_Values_if), pointer, public :: RANS_Set_Wall_Values
  procedure(RANS_Set_Blowing_Wall_if), pointer, public :: RANS_Set_Blowing_Wall
  procedure(RANS_Extrapolate_Wall_if), pointer, public :: RANS_Extrapolate_Wall
  procedure(RANS_Enforce_Realizability_if), pointer, public :: RANS_Enforce_Realizability

  !> Devirtualized turbulence-model ids for the GPU device build (procedure
  !> pointers can't run on the GPU). Set in Set_RANS_Model.
  !>   eddy_dev / ransdiff_dev : 1 = SST ; 2 = SA / SA-R ; 3 = Wilcox2006 ; 0 = unsupported (device hard-stop in Set_RANS_Model)
  !>   stress_dev              : 1 = Std ; 2 = 2-eq (k-energy-coupling with SST or Wilcox2006)
  integer, public :: eddy_dev = 0, ransdiff_dev = 0, stress_dev = 1
  !$acc declare create(eddy_dev, ransdiff_dev, stress_dev)

  !> Device-resident copies of the diffusion non-dimensional numbers. obj_rans
  !> is a derived type, so its scalar members can't be read inside a GPU acc
  !> routine; these plain scalars are captured + pushed in Set_RANS_Model.
  real(R8), public :: sc_dev = 0.7d0, sct_dev = 0.9d0, prt_dev = 0.85d0
  real(R8), public :: prl_dev = 0.0d0   !> [MOSE-Physics] prandtl on the device (0 = Wilke conductivity from the table)
  !$acc declare create(prl_dev)
  !$acc declare create(sc_dev, sct_dev, prt_dev)

  !> DIFF_MC: .true. with [MOSE-Physics] schmidt <= 0 -> mixture-averaged multicomponent D_k (co_DS_expr) + mu_t/(rho Sct)
  !> in place of (mu/Sc + mu_t/Sct)/rho. Device copy pushed in Set_RANS_Model (and copied in on the frozen path).
  logical, public :: mcdiff_dev = .false.
  !$acc declare create(mcdiff_dev)

  !> Abstract interface relative to the eddy viscosity computation procedure.
  abstract interface
    subroutine Eddy_Viscosity_if ( mut, rans_variables, mul, rho, vel_gradient, walldist )
      use MOSE_Global_m
      implicit none
      real(kind=8), intent(in), dimension(nRANS) :: rans_variables  ! : vector of RANS eqs variables in form (rho*R)
      real(kind=8), intent(in)                   :: mul             ! : Laminar viscosity
      real(kind=8), intent(in)                   :: rho             ! : Density
      real(kind=8), intent(in), dimension(3,3)   :: vel_gradient    ! : Velocity gradient
      real(kind=8), intent(in)                   :: walldist        ! : Distance nearest wall
      real(kind=8), intent(out)                  :: mut             ! : Eddy viscosity
    end subroutine Eddy_Viscosity_if
  end interface

  !> Abstract interface relative to the computation of conservative diffusive fluxes for RANS equations.
  abstract interface
    subroutine RANS_diffusive_flux_if(flux, rans_variables, mul, vel_gradient, rans_gradient, rho, area, normal, dist)
      use MOSE_Global_m
      implicit none
      real(kind=8), intent(in), dimension(nrans)   :: rans_variables         ! : vector of RANS eqs variables in form (rho*R)
      real(kind=8), intent(in)                     :: mul                    ! : laminar viscosity
      real(kind=8), intent(in), dimension(3,3)     :: vel_gradient           ! : Velocity gradient
      real(kind=8), intent(in), dimension(nrans,3) :: rans_gradient          ! : Gradient in form dRi/dxj
      real(kind=8), intent(in)                     :: rho                    ! : Density
      real(kind=8), intent(in)                     :: area, normal(3), dist  ! : Metrics (face area, face normal, distance from wall)
      real(kind=8), intent(out), dimension(nrans)  :: flux                   ! : Diffusive flux
    end subroutine RANS_diffusive_flux_if
  end interface

  !> Abstract interface relative to the computation of RHS terms for RANS equations.
  abstract interface
    subroutine RANS_source_terms_if(domain)
      use MOSE_Types_m
      implicit none
      type(MOSE_domain_type), intent(inout) :: domain
    end subroutine RANS_source_terms_if
  end interface

  !> Abstract interface relative to the computation of wall values for RANS variables.
  abstract interface
    subroutine RANS_Set_Wall_Values_if(mil, rans_variables, dist)
      use MOSE_Global_m
      implicit none
      real(kind=8), intent(out), dimension(nRANS)  :: rans_variables   ! : ... at wall
      real(kind=8), intent(in)                     :: mil              ! : laminar viscosity at wall
      real(kind=8), intent(in)                     :: dist             ! : cell center-wall distance
    end subroutine RANS_Set_Wall_Values_if
  end interface

  !> Abstract interface relative to the computation of blowing wall values for RANS variables.
  abstract interface
    subroutine RANS_Set_Blowing_Wall_if(rho, mil, tau, rans_variables, mdot, dist)
      use MOSE_Global_m
      implicit none
      real(kind=8), intent(out), dimension(nRANS)  :: rans_variables
      real(kind=8), intent(in)                     :: rho, mil, mdot, dist
      real(kind=8), intent(in), dimension(3)       :: tau
    end subroutine RANS_Set_Blowing_Wall_if
  end interface

  abstract interface
    subroutine RANS_Extrapolate_Wall_if ( prim, wall_rans, rho, wall_rho, ghost_rans )
      use MOSE_Global_m
      implicit none
      real(kind=8), intent(in),  dimension(nRANS) :: prim, wall_rans  ! : boundary cell and wall variables
      real(kind=8), intent(in) :: rho, wall_rho ! : boundary cell density and wall density
      real(kind=8), intent(out), dimension(nRANS) :: ghost_rans ! : ghost cell RANS variables
    end subroutine RANS_Extrapolate_Wall_if
  end interface

  abstract interface
    subroutine RANS_Enforce_Realizability_if ( var )
      implicit none
      real(kind=8), intent(inout),  dimension(:) :: var
    end subroutine RANS_Enforce_Realizability_if
  end interface

contains


  ! Setting the value of turbulent stresses at symmetry interfaces and in ghost cells
  ! NOTE: symmetry interface MUST be either an x-plane, y-plane or z-plane

end module MOSE_Lib_RANS