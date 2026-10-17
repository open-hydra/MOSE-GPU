module MOSE_IO_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use MOSE_Parameters_m

  implicit none

  type :: rans_type
    real(R8) :: Sc  = 0.7d0  !> Schmidt laminare
    real(R8) :: Sct = 0.90d0 !> Schmidt turbolento
    real(R8) :: Prt = 0.85d0 !> Prandtl turbolento
    character(llen) :: rans_name
    logical :: QCR2000 = .false., SpalartShur = .false., SAR = .false., k_energy_coupling = .false., &
               blowing_corr = .false., RSM = .false., SAcomp = .false.   ! the flags of Setup_RANS_Model, false unless set
    logical :: point_implicit = .false.   !> [MOSE-RANS] point-implicit: Patankar treatment of the turbulence destruction (MOSE_open default t)
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type

  type(rans_type), public :: obj_rans

contains

  subroutine read_ini(self, fini)
    use MOSE_Parameters_m
    use MOSE_Global_m, only: omega_wall_c, sst_comp_c
    use Finer, only: file_ini

    implicit none
    class(rans_type) :: self
    type(file_ini) :: fini
    ! Local
    integer :: error, i
    character(len=:), allocatable :: list(:) 
    character(len=32) :: wkey


    self%rans_name = 'inviscid'

    ! Look for RANS section
    call fini%get_sections_list(list=list)
    do i = 1, size(list,1)
      if (index(list(i),'RANS')>0) then

        ! Model name
        call fini%get(section_name=trim(codename)//'-RANS', option_name='turbulence-model', val=self%rans_name, error=error)
        if (error/=0) self%rans_name = 'laminar'

        ! Nondimensional parameters : Prt, Sct, Sc
        call fini%get(section_name=trim(codename)//'-RANS', option_name='Prt', val=self%Prt, error=error)
        call fini%get(section_name=trim(codename)//'-RANS', option_name='Sct', val=self%Sct, error=error)
        call fini%get(section_name=trim(codename)//'-RANS', option_name='Sc', val=self%Sc, error=error)

        ! Submodels
        if ( index ( trim(self%rans_name), '-R' ) > 0 ) then
          if ( index ( trim(self%rans_name), '-RC' ) > 0 ) then
            self%SpalartShur = .true.
            self%SAR = .false.
          else
            self%SAR = .true.
            self%SpalartShur = .false.
          end if
        end if

        if ( index ( trim(self%rans_name), '-QCR2000' ) > 0 ) then
          self%QCR2000 = .true.
        else
          self%QCR2000 = .false.
        end if

        if ( index ( trim(self%rans_name), '-blowcorr' ) > 0 ) then
          self%blowing_corr = .true.
        else
          self%blowing_corr = .false.
        end if

        call fini%get(section_name=trim(codename)//'-RANS', option_name='k-coupling', val=self%k_energy_coupling, error=error)
        if (error/=0) call fini%get(section_name=trim(codename)//'-RANS', option_name='k-energy-coupling', &
                                    val=self%k_energy_coupling, error=error)   ! MOSE_open's name of the same switch
        if (error/=0) self%k_energy_coupling = .false.
        call fini%get(section_name=trim(codename)//'-RANS', option_name='point-implicit', val=self%point_implicit, error=error)
        if (error/=0) self%point_implicit = .false.

        ! k-omega (SST, Wilcox2006) omega at a smooth wall (MOSE_open omega-wall-bc, Mod_RANS): practical, 60/beta_1 = 800 (default), or
        ! asymptotic, 6/beta_1 = 80, the y -> 0 limit for wall-resolved meshes
        wkey = 'practical'
        call fini%get(section_name=trim(codename)//'-RANS', option_name='omega-wall-bc', val=wkey, error=error)
        if (error/=0) wkey = 'practical'
        select case ( trim(adjustl(wkey)) )
        case ('practical');  omega_wall_c(1) = 8d2
        case ('asymptotic'); omega_wall_c(1) = 8d1
        case default
          write(*,'(3A)') ' [ERROR] [MOSE-RANS] omega-wall-bc = ', trim(wkey), ': practical or asymptotic'
          error stop 'omega-wall-bc must be practical or asymptotic (see the [ERROR] line)'
        end select

        ! SST production (MOSE_open sst-production): incompressible (default) or compressible
        wkey = 'incompressible'
        call fini%get(section_name=trim(codename)//'-RANS', option_name='sst-production', val=wkey, error=error)
        if (error/=0) wkey = 'incompressible'
        select case ( trim(adjustl(wkey)) )
        case ('incompressible'); sst_comp_c(1) = .false.
        case ('compressible');   sst_comp_c(1) = .true.
        case default
          write(*,'(3A)') ' [ERROR] [MOSE-RANS] sst-production = ', trim(wkey), ': incompressible or compressible'
          error stop 'sst-production must be incompressible or compressible (see the [ERROR] line)'
        end select

      endif
    enddo
    !$acc update device(omega_wall_c, sst_comp_c)

    ! SAcomp: the model name exactly (MOSE_open Setup_RANS_Model)
    self%SAcomp = ( trim(self%rans_name) == 'SAcomp' )

  end subroutine read_ini


end module MOSE_IO_RANS