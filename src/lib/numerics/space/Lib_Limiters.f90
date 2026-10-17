module MOSE_Lib_Limiters
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: assign_limiter
  public :: rlimiter_MINMOD   ! exposed for the devirtualized device convective flux
  public :: rlimiter      ! the device limiter dispatch (lim_dev = 1 minmod, 2 vanleer, 3 vanalbada, 4 mc, 5 superbee)

  logical, public :: SD_space_order
  !> Read inside the device convective flux (`if (SD_space_order)`); make it
  !> device-resident. Synced at the end of Assign_Limiter.
  !$acc declare create(SD_space_order)

  !> The limiter id of the device kernels (rlimiter): 1 Min-Mod, 2 Van Leer, 3 Van Albada, 4 MC, 5 Super-Bee,
  !> 6 first order; 0 = a limiter the device kernels do not have (Assign_RS stops).
  integer, public :: lim_dev = 0
  !$acc declare create(lim_dev)

contains

  subroutine Assign_Limiter(limiter_word)
    implicit none
    character(len=*), intent(inout) :: limiter_word

    lim_dev = 0
    if ( limiter_word(len_trim(limiter_word)-1:len_trim(limiter_word)) == 'SD' ) then
      limiter_word = trim(limiter_word(1:len_trim(limiter_word)-3))
      SD_space_order = .true.
    else
      SD_space_order = .false.
    endif

    select case (limiter_word)
      case ('IORD')
        lim_dev = 6   ! first order on the device: zero slope (rlimiter case 6)
      case ('MINMOD','MinMod','Min-Mod')
        limiter_word = 'Min-Mod'
        lim_dev = 1   ! device-supported limiter
      case ('VANLEER','VanLeer','Van Leer')
        limiter_word = 'Van Leer'
        lim_dev = 2
      case ('VANALBADA','VanAlbada','Van Albada')
        limiter_word = 'Van Albada'
        lim_dev = 3
      case ('MC')
        limiter_word = 'MC'
        lim_dev = 4
      case ('LIMO3')
        limiter_word = 'LIMO3'
      case ('SUPERBEE','Super-Bee','Super Bee')
        limiter_word = 'Super-Bee'
        lim_dev = 5
      case default
        limiter_word = 'IORD'
    end select

    if (SD_space_order) limiter_word = trim(limiter_word)//' with shock detection'

    !$acc update device(SD_space_order, lim_dev)

  end subroutine Assign_Limiter

  !> The device limiter: one function dispatched by lim_dev (acc routine for the OpenACC BC kernels and the host side of
  !> the device build; the CUF convective kernels carry their own copy rlimiter_d). Case 1 is rlimiter_MINMOD verbatim.
  pure function rlimiter(lim, x, y) result(result)
    !$acc routine seq
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    integer,  intent(in) :: lim
    real(R8), intent(in) :: x, y
    real(R8) :: ru, phi, result
    select case (lim)
    case (1)   ! Min-Mod: the statements of rlimiter_MINMOD (bit-identical)
      result = 0.d0
      if (x*y<0.d0) return
      result = x
      if (abs(x)>abs(y)) result = y
    case (2)   ! Van Leer (rlimiter_VANLEER)
      ru = y/(x + 1d-30)
      phi = (ru+abs(ru))/(1d0+abs(ru))
      result = phi*x
    case (3)   ! Van Albada (rlimiter_VANALBADA)
      ru = y/(x+1d-30)
      phi = (ru*ru+ru)/(1d0+ru*ru)
      result = phi*x
    case (4)   ! MC (rlimiter_MC: maxval/minval of the same values)
      ru = y/(x+1d-30)
      phi = max(0d0, min(2d0*ru, 0.5d0*(1d0+ru), 2d0))
      result = phi*x
    case (5)   ! Super-Bee (rlimiter_SB)
      ru = y/(x+1d-30)
      phi = max(0d0, min(2d0*ru,1d0), min(ru,2d0))
      result = phi*x
    case default   ! 6 = IORD (first order): zero slope
      result = 0d0
    end select
  end function rlimiter

  !> Min-Mod flux limiter
  pure function rlimiter_MINMOD(x,y) result(result)
    !$acc routine seq
    use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
    implicit none
    real(R8), intent(in) :: x, y
    real(R8) :: ru
    real(R8) :: phi
    real(R8) :: result
    
    result = 0.d0
    if (x*y<0.d0) return
    result = x
    if (abs(x)>abs(y)) result = y

  end function rlimiter_MINMOD

end module MOSE_Lib_Limiters