module MOSE_Lib_RK
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  real(R8), parameter, private :: RKcoeff(3,3) = reshape([  1._R8,   0._R8,      0._R8,  &
                                                            1._R8, 0.50_R8,      0._R8,  &
                                                            1._R8, 0.25_R8,2._R8/3._R8 ],&
                                                            shape(RKcoeff), order=[2,1])

  type :: rk_type
    integer :: n_RK
    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type rk_type

  type(rk_type), public :: obj_rk
  private :: read_ini

contains


  subroutine read_ini ( self, fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename

    implicit none
    class(rk_type) :: self
    type(file_ini) :: fini
    integer :: error

    call fini%get(section_name=trim(codename)//'-scheme', option_name='RK', val=self%n_RK, error=error)
    if (error/=0) self%n_RK = 1
    
  end subroutine



  !> The scalar RK-stage coefficient RKcoeff(n_RK, irk). Exposed so the
  !> devirtualized device state-update (Newstate_Primitive_acc) can inline the
  !> RK_Stage combination without a procedure call / nprim automatic.
  pure function rk_coeff ( irk ) result ( c )
    implicit none
    integer, intent(in) :: irk
    real(R8) :: c
    c = rkcoeff(obj_rk%n_RK, irk)
  end function rk_coeff
      

end module MOSE_Lib_RK
