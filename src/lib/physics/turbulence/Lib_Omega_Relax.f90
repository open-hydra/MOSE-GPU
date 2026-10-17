!> omega-source-relax: the two per-cell SST relaxations of the production RDE runs, with the THOR constants.
!>   - The point-implicit SST destruction on the residual: the k and omega rows of the resident residual R_flat are
!>     divided by (1 + dt beta* omega) and (1 + dt 2 beta omega) where either factor exceeds SRC_EPS, only outside
!>     the reacting band (T < tign, or Y_fuel*Y_ox < yeps, or a heat-release measure hr_c below hreps): the cold
!>     fresh fill and the burned gas relax, an igniting cell keeps its explicit destruction.
!>   - A physical omega floor on the updated state, after Newstate: omega >= rho * omfloor. A convective undershoot
!>     on a sub-micron SST wall cell clips omega to ~0, nu_t = k/omega explodes and the viscous time step collapses;
!>     the floor at the ambient scale bounds nu_t <= k/omfloor and never touches a cell whose omega is meaningful.
!> SST only: [MOSE-scheme] omega-source-relax = t | f (default f). The species of the reacting-band test are the
!> fuel and the oxidizer of ONERA-7 (H2 = 3, O2 = 1).
module MOSE_Lib_Omega_Relax
  use iso_fortran_env, only: R8 => real64

  implicit none
  private
  public :: obj_omega_relax

  type :: omega_relax_type
    logical  :: active  = .false.
    real(R8) :: omfloor = 1.0e2_R8     ! omega floor [1/s]
    real(R8) :: tign    = 1100.0_R8    ! reacting band: T >= tign ...
    real(R8) :: yeps    = 1.0e-4_R8    ! ... and Y_fuel*Y_ox >= yeps ...
    real(R8) :: hreps   = 0.05_R8      ! ... and |dT|/T of the last chemistry step >= hreps
    integer  :: ifu     = 3            ! fuel species (ONERA-7: H2)
    integer  :: iox     = 1            ! oxidizer species (ONERA-7: O2)
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type omega_relax_type

  type(omega_relax_type) :: obj_omega_relax

  public :: Omega_Relax_Source, Omega_Relax_Floor
  ! cells where both implicit factors are below SRC_EPS are left untouched
  real(R8), parameter :: SRC_EPS = 0.05_R8

contains


  subroutine read_ini ( self, fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    use MOSE_Mod_MPI, only: mpi_is_root
    use MOSE_IO_RANS, only: obj_rans
    implicit none
    class(omega_relax_type) :: self
    type(file_ini) :: fini
    integer :: error

    ! the jimp-* keys of the omega relaxation became omega-source-relax = t (THOR constants) | f
    call jimp_obsolete_keys ( fini, trim(codename)//'-scheme', [character(len=20) :: 'jimp-src', 'jimp-omega-floor', &
         'jimp-chem-tign', 'jimp-chem-yeps', 'jimp-r-on', 'jimp-mask-power', 'jimp-chem-hr', 'jimp-chem-fuel', 'jimp-chem-ox'] )
    call fini%get(section_name=trim(codename)//'-scheme', option_name='omega-source-relax', val=self%active, error=error)
    if (error/=0) self%active = .false.
    ! k-omega destruction and omega floor: SST physics, refused with any other model
    if ( self%active .and. index(obj_rans%rans_name, 'SST') == 0 ) &
      error stop 'omega-source-relax needs the SST model: not with Wilcox2006, SA, laminar or Euler'
    if ( self%active .and. mpi_is_root ) write(*,'(A,F7.1,A,ES9.2,A,F5.2,A,I0,A,I0,A,ES9.2)') &
      ' [omega-source-relax] on: tign ', self%tign, '  yeps ', self%yeps, '  chem-hr ', self%hreps, &
      '  fuel ', self%ifu, '  ox ', self%iox, '  omega floor ', self%omfloor
  end subroutine read_ini


  !> The old spelling of the relaxation: stop if one of its keys is still in the input (it would be ignored and the
  !> run would silently lose or change the omega relaxation).
  subroutine jimp_obsolete_keys ( fini, sec, keys )
    use Finer, only: file_ini
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    type(file_ini)               :: fini
    character(len=*), intent(in) :: sec
    character(len=*), intent(in) :: keys(:)
    character(len=256) :: sval
    integer :: i, error
    do i = 1, size(keys)
      call fini%get(section_name=sec, option_name=trim(keys(i)), val=sval, error=error)
      if (error /= 0) cycle
      if (mpi_is_root) write(*,'(5A)') ' ERROR: input key ', trim(keys(i)), ' in [', sec, &
        '] no longer exists -- use omega-source-relax = t (THOR constants) | f (same section)'
      error stop 'obsolete input key of the omega relaxation (see the ERROR line)'
    end do
  end subroutine jimp_obsolete_keys


  !> The point-implicit SST destruction on the resident residual R_flat: one kernel per local block.
  !> Called inside `!$omp single`.
  subroutine Omega_Relax_Source ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: nrans
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    if (nrans < 2) return
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Omega_Relax_Source_blk ( b, domain%blk(b)%dtlocal, domain%blk(b)%rho_c, &
             domain%blk(b)%Rg_c, domain%blk(b)%dim, obj_omega_relax%tign, obj_omega_relax%yeps, obj_omega_relax%ifu, &
             obj_omega_relax%iox, domain%blk(b)%hr_c, obj_omega_relax%hreps )
    end do
  end subroutine Omega_Relax_Source


  subroutine Omega_Relax_Source_blk ( b, dtlocal, rho_c, Rg_c, n, tign, yeps, ifu, iox, hr_c, hreps )
    use MOSE_Global_m, only: gc, np, nt
    use MOSE_Field_m,  only: P_flat, R_flat, field_block_index_params, rho_c_f, Rg_c_f
    use MOSE_Mod_SST,  only: beta_star, beta_2
    implicit none
    integer,  intent(in) :: b, n(3), ifu, iox
    real(R8), intent(in) :: tign, yeps, hreps
    real(R8), dimension(1:n(1), 1:n(2), 1:n(3)), intent(in) :: dtlocal
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: rho_c, Rg_c, hr_c
    ! Local
    integer  :: i, j, k, c0, boff, pi, pj, nb1_, nb2_, nb3_
    real(R8) :: rho, om, dtc, fk, fw, T
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    call field_block_index_params(b, boff, pi, pj)

    !$acc data copyin(dtlocal, rho_c, Rg_c, hr_c)
    !$acc parallel loop collapse(3) gang vector present(P_flat, R_flat, dtlocal, rho_c, Rg_c, hr_c, rho_c_f, Rg_c_f) &
    !$acc&  private(c0, rho, om, dtc, fk, fw, T)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0  = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      ! omega = (rho omega)/rho from the state and the start-of-stage density cache
      rho = rho_c_f(c0)
      om  = P_flat(c0, nt+1) / rho
      dtc = dtlocal(i,j,k)
      fk  = dtc * beta_star * om
      fw  = dtc * 2d0 * beta_2 * om
      if ( fk > SRC_EPS .or. fw > SRC_EPS ) then
        ! only outside the reacting band
        T = P_flat(c0, np) / ( rho * Rg_c_f(c0) )
        if ( ( T < tign ) .or. ( P_flat(c0, ifu) * P_flat(c0, iox) < yeps * rho * rho ) .or. &
             ( hreps > 0d0 .and. hr_c(i,j,k) < hreps ) ) then
          R_flat(c0, nt  ) = R_flat(c0, nt  ) / ( 1d0 + max(fk, 0d0) )
          R_flat(c0, nt+1) = R_flat(c0, nt+1) / ( 1d0 + max(fw, 0d0) )
        end if
      end if
    end do
    end do
    end do
    !$acc end data
  end subroutine Omega_Relax_Source_blk


  !> omega >= rho * omfloor on the updated state (rho from the updated species): one kernel per local block.
  !> Called inside `!$omp single`.
  subroutine Omega_Relax_Floor ( domain )
    use MOSE_Types_m
    use MOSE_Global_m, only: nrans
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    integer :: b
    if (nrans < 2) return
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      call Omega_Relax_Floor_blk ( b, domain%blk(b)%dim, obj_omega_relax%omfloor )
    end do
  end subroutine Omega_Relax_Floor


  subroutine Omega_Relax_Floor_blk ( b, n, omfloor )
    use MOSE_Global_m, only: gc, nsc, nt
    use MOSE_Field_m,  only: P_flat, field_block_index_params
    implicit none
    integer,  intent(in) :: b, n(3)
    real(R8), intent(in) :: omfloor
    ! Local
    integer  :: i, j, k, s, c0, boff, pi, pj, nb1_, nb2_, nb3_
    real(R8) :: rho
    nb1_ = n(1); nb2_ = n(2); nb3_ = n(3)
    call field_block_index_params(b, boff, pi, pj)

    !$acc parallel loop collapse(3) gang vector present(P_flat) private(c0, rho, s)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      c0  = boff + (i+gc-1) + (j+gc-1)*pi + (k+gc-1)*pi*pj + 1
      rho = 0d0
      do s = 1, nsc
        rho = rho + P_flat(c0, s)
      end do
      P_flat(c0, nt+1) = max( P_flat(c0, nt+1), rho * omfloor )
    end do
    end do
    end do
  end subroutine Omega_Relax_Floor_blk

end module MOSE_Lib_Omega_Relax
