module MOSE_Lib_Shock_Detector
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

  real(R8), parameter, private :: DELTA=20.d0

  logical, public :: SD = .false.

  !> The shock sensor of MOSE_open on the device (Mod_Fluxes Fluxes_blk + Lib_Shock_Detector): 0 none (beta = 1),
  !> 1 Tramel, 2 Chen. Set by Assign_SD from [MOSE-scheme] shock-detector or, when absent, from the Riemann solver
  !> as MOSE_open Assign_Setup does (HLLE++ / HLLC+Tramel -> Tramel, AUSM+M -> Chen). With a sensor on, beta of a
  !> cell scales the reconstruction slopes and goes to the Riemann solver of its faces (MOSE_open: the face's left
  !> cell; connections: the record's own cell; symmetry: 1).
  integer, public :: sd_dev = 0
  !$acc declare create(sd_dev)

contains

  !> MOSE_open Assign_Setup (shock detector): the key, else the Riemann solver's own sensor. rs_name is the GPU
  !> riemann key (AUSM+M, HLLE++, HLLC+Tramel ...). Sets sd_dev and, with a sensor, SD_space_order (the slopes are
  !> scaled by beta as in MOSE_open Reconstruction).
  subroutine Assign_SD ( sd_key, rs_name )
    use MOSE_Lib_Limiters, only: SD_space_order
    implicit none
    character(len=*), intent(in) :: sd_key, rs_name
    character(len=len(sd_key)) :: k
    integer :: i, c
    k = adjustl(sd_key)
    do i = 1, len_trim(k)
      c = iachar(k(i:i)); if (c >= iachar('A') .and. c <= iachar('Z')) k(i:i) = achar(c + 32)
    end do
    select case ( trim(k) )
      case ('tramel'); sd_dev = 1
      case ('chen');   sd_dev = 2
      case ('none');   sd_dev = 0
      case ('')
        sd_dev = 0
        if ( trim(rs_name) == 'HLLE++' .or. trim(rs_name) == 'HLLC+Tramel' ) sd_dev = 1   ! descriptions with 'Tramel'
        if ( trim(rs_name) == 'AUSM+M' .or. trim(rs_name) == 'HLLC+Chen' ) sd_dev = 2        ! descriptions with 'Chen'
      case default
        write(*,'(3A)') ' [ERROR] [MOSE-scheme] shock-detector = ', trim(sd_key), ': tramel, chen or none'
        error stop 'Assign_SD: unknown shock-detector'
    end select
    if ( sd_dev > 0 ) SD_space_order = .true.
    !$acc update device(sd_dev, SD_space_order)   ! created on the device by IO_Chemistry, which runs first
  end subroutine Assign_SD


  !> The sensor of one cell from its pressure p and the pressures of its six face neighbours (i-1, i+1, j-1, j+1,
  !> k-1, k+1): id 1 = MOSE_open SD_Tramel, 2 = SD_Chen, otherwise 1.
  pure function SD_beta ( id, p, p1, p2, p3, p4, p5, p6 ) result ( beta )
    !$acc routine seq
    implicit none
    integer,  intent(in) :: id
    real(R8), intent(in) :: p, p1, p2, p3, p4, p5, p6
    real(R8) :: beta
    real(R8) :: kp1, kp2, kp3, sw, phi
    select case (id)
    case (1)
      kp1 = abs( (p2 - 2d0*p + p1) / (p2 + 2d0*p + p1) )
      kp2 = abs( (p4 - 2d0*p + p3) / (p4 + 2d0*p + p3) )
      kp3 = abs( (p6 - 2d0*p + p5) / (p6 + 2d0*p + p5) )
      sw = max(kp1,kp2,kp3)
      if (sw<1d0/DELTA) then
        phi = max(sw/DELTA,0d0)
        beta = 1d0-tanh(10d0*phi*phi*phi)
      else
        beta = 0d0
      end if
    case (2)
      beta = 1.d0
      beta = min(beta, min(p, p1) / max(p, p1))
      beta = min(beta, min(p, p2) / max(p, p2))
      beta = min(beta, min(p, p3) / max(p, p3))
      beta = min(beta, min(p, p4) / max(p, p4))
      beta = min(beta, min(p, p5) / max(p, p5))
      beta = min(beta, min(p, p6) / max(p, p6))
    case default
      beta = 1.d0
    end select
  end function SD_beta


  pure function SD_Tramel(pmatrix) result(beta)
    !$acc routine seq
    use MOSE_Global_m
    implicit none
    real(R8), intent(in)  :: pmatrix(0:2,0:2,0:2)
    real(R8)              :: beta
    ! Local
    real(R8) :: kp1, kp2, kp3, sw, phi
    real(R8) :: p1, p2, p3, p4, p5, p6, p

    p  = pmatrix(1,1,1)
    p1 = pmatrix(0,1,1)
    p2 = pmatrix(2,1,1)
    p3 = pmatrix(1,0,1)
    p4 = pmatrix(1,2,1)
    p5 = pmatrix(1,1,0)
    p6 = pmatrix(1,1,2)

    ! Jameson-type sensor
    kp1 = abs( (p2 - 2d0*p + p1) / (p2 + 2d0*p + p1) )
    kp2 = abs( (p4 - 2d0*p + p3) / (p4 + 2d0*p + p3) )
    kp3 = abs( (p6 - 2d0*p + p5) / (p6 + 2d0*p + p5) )
    sw = max(kp1,kp2,kp3)

    if (sw<1d0/DELTA) then
      phi = max(sw/DELTA,0d0)
      beta = 1d0-tanh(10d0*phi*phi*phi)
    else
      beta = 0d0
    end if

  end function SD_Tramel

end module MOSE_Lib_Shock_Detector
