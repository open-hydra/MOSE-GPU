module MOSE_Lib_Strang
  implicit none

  integer :: N_strang     !> 1 w/o strang;   2 w strang  
  real(8) :: strangcoeff  !> 1 w/o strang; 1/2 w strang  

contains
  

  subroutine set_StrangSplitting(switch_on)
    implicit none
    logical, intent(in) :: switch_on
    
    strangcoeff = 1.d0
    N_strang = 1
    
    if (switch_on) then
      strangcoeff = 0.5d0
      N_strang = 2
    endif

  end subroutine set_StrangSplitting


end module MOSE_Lib_Strang