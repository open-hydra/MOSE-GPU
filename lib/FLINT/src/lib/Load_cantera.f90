module FLINT_cantera_load
# if defined (CANTERA)
  use cantera
  implicit none

contains

  subroutine load_phase(phase, filename, phase_name)
    implicit none
    type(phase_t), intent(out) :: phase
    character(len=*), intent(in) :: filename
    character(len=*), intent(inout), optional :: phase_name
    ! Import the phase from the specified file and name
    if (.not. present(phase_name)) then
      phase_name = ''
    end if
    phase = importPhase(filename, phase_name)
  end subroutine load_phase

# endif
end module FLINT_cantera_load
