module MOSE_Wrap_Solve

  implicit none

contains


  subroutine MOSE_solve ( simulation, External_Function )
    use MOSE_types_m
    use MOSE_Mod_Explicit
    use MOSE_Lib_Interpolation
    
    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    external :: External_Function
    ! Local
    logical :: new_domain

    ! 1. Calculate solution on local grid
    call Explicit_Step( simulation%domain, simulation%config, simulation%config%MGL, simulation%config%MG_level, External_Function )

    ! 2. Interpolate solution on finer grid
    if (simulation%config%change_MG) then
      call Prolongation ( Fine=simulation%domain(simulation%config%MG_level-1), Coarse=simulation%domain(simulation%config%MG_level) )
    endif
    
  end subroutine MOSE_solve


end module MOSE_Wrap_Solve