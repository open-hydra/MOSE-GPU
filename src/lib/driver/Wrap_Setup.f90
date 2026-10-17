module MOSE_Wrap_Setup

  implicit none

contains


  !> The setup fill: P (the host truth at setup) packed into the resident P_flat, MOSE_open's Fill_Ghost_Cell on the device,
  !> w back into P, ghost cells included.
  subroutine Fill_Ghost_Cell_Setup ( domain )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Lib_Extrapolation, only: Fill_Ghost_Cell
    use MOSE_Field_m, only: pack_P_to_flat, unpack_flat_to_P, field_update_P_flat_device, field_update_P_flat_host
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    call pack_P_to_flat ( domain )
    call field_update_P_flat_device ()
    call Fill_Ghost_Cell ( domain )
    !$acc wait
    call field_update_P_flat_host ()
    call unpack_flat_to_P ( domain )
  end subroutine Fill_Ghost_Cell_Setup

  subroutine MOSE_setup ( simulation )
    use MOSE_Types_m
    use MOSE_Parameters_m
    use MOSE_Global_m
    use MOSE_IO_Ini
    use MOSE_IO_Table
    use MOSE_IO_Chemistry
    use MOSE_IO_RANS
    use MOSE_Mod_RANS
    use MOSE_IO_Solution
    use MOSE_Mod_Allocate_Data
    use MOSE_Lib_Multigrid
    use MOSE_IO_BC
    use MOSE_IO_Wall
    use MOSE_Lib_Metrics
    use MOSE_Lib_RK
    use MOSE_Mod_Probes
    use MOSE_Mod_Diagnostic, only : unitRES, unitRES_format
    use MOSE_Lib_Extrapolation
    use MOSE_Lib_Interpolation
    use MOSE_Lib_Soot
    use MOSE_Mod_BlockSplit, only: split_domain_for_mpi, rebuild_IOfield_after_split
    use MOSE_Mod_MPI, only: mpi_is_root, is_local_block
    use MOSE_Mod_GhostExchange, only: build_ghost_schedule, build_local_bc_index
    use MOSE_Mod_BC, only: Setup_BC_Device, Check_Device_Boundary_Types
    use MOSE_Mod_Explicit, only: Setup_Device_Residency
    use MOSE_IO_BC_Open, only: Resolve_Open_BC_Direction

    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    ! Local
    integer :: m, ios


    ! Print Header
    if (mpi_is_root) then
    call Print_Header()
    endif
    
    ! Read input.ini
    call Read_Inifile ( simulation%config )

    ! Allocate container objects Domain and IOfield.
    allocate ( simulation%domain  ( simulation%config%MGL ) )
    allocate ( simulation%IOfield ( simulation%config%MGL ) )
    simulation%domain(1)%itermax = simulation%config%itermax

    ! Read file for initial solution.
    call Read_IC ( simulation%IOfield(1), simulation%config%newrun )

    ! Allocation of data structures and copy solution from IOfield.
    call Setup_Data_Structure ( simulation%domain(1), simulation%IOfield(1), simulation%config%MGL )
    allocate(simulation%config%residuotot(nres))
    
    ! With multigrid, allocate Grid and IOfield 2,...,MGL
    if ( simulation%config%MGL > 1 ) then
      call Setup_Multigrid ( simulation )
    end if
    
    ! Read Boundary Conditions file
    call Setup_BC ( simulation%domain, simulation%config%MGL )
    do m = 1, simulation%config%MGL
      call Check_Device_Boundary_Types ( simulation%domain(m) )
    end do

    ! MPI: split blocks if needed to ensure enough blocks for all ranks
    do m = 1, simulation%config%MGL
      !call split_domain_for_mpi(simulation%domain(m))
      !call rebuild_IOfield_after_split(simulation%IOfield(m), simulation%domain(m))
      call build_ghost_schedule(simulation%domain(m))
      call build_local_bc_index(simulation%domain(m))
    end do

   ! Read grid file and setup metrics in each grid level.
    do m = 1, simulation%config%MGL
      if (mpi_is_root) write(*,'(A28,I1,A3)') ' Computing metrics for Grid ', m, '...'
      call Setup_Metrics ( simulation%domain(m) )
      if (mpi_is_root) write(*,'(A8)') ' ...done'
      ! MOSE_open-dialect inflow records written 'normal': the angles from the face
      ! normal, now that the metrics exist (no-op for a legacy bc file).
      call Resolve_Open_BC_Direction ( simulation%domain(m) )
    end do

    ! Flatten the boundary-face list into plain device-mappable
    ! arrays (after Setup_BC has filled bc(:) and n_bf). CPU-only no-op otherwise.
    do m = 1, simulation%config%MGL
      call Setup_BC_Device ( simulation%domain(m) )
    end do

    ! Print simulation info onto the logfile/shell.
    if (mpi_is_root) then
    call Print_Shell_Info ( simulation%config )
    endif

    ! Setup probes location, output files, and variables to be printed.
    call Setup_Probes ( simulation%domain(1) )

    ! Solution setup
    call Output_Solution_Setup ( simulation%IOfield, simulation%config%MGL )

    ! Initialize Wall file
    if (WALL .and. mpi_is_root) call Initialize_Wall_File ( simulation%domain(1), simulation%IOfield(1)%tec%extension)
  
    ! Residual history file (only root writes)
    if (mpi_is_root) then
      open(newunit=unitRES,file='OUTPUT/'//trim(MOSE_phase_prefix)//'residual-history.dat',status='unknown',form='formatted')
      if ( .not. simulation%config%newrun ) then
        ios = 0
        do while ( ios == 0 )
          read (unitRES, *, iostat=ios)
        enddo
        backspace (unitRES)
      endif
    end if

    ! Setting the format for writing residual_history file
    write(unitRES_format,'(A11,I0,A7)') '(I8,E20.10,', nres, 'E20.10)'

    ! Initialize run variables depending on scheme
    select case (trim(simulation%config%solver_name))

      case ('explicit')

        simulation%config%MG_level = simulation%config%MGL

        do m = 1, simulation%config%MGL
          simulation%domain(m)%iter = 0

          if (simulation%config%newrun) then
            if (simulation%config%time_accurate) then
              simulation%domain(m)%time = 0.0d0
            else
              simulation%domain(m)%time = -1.0d0
            endif
            simulation%config%iter_general = 0
          else
            if (simulation%config%time_accurate) then
              simulation%domain(m)%time = simulation%IOfield(m)%solutiontime
              simulation%config%iter_general = 0
            else
              simulation%config%iter_general = int(simulation%IOfield(m)%solutiontime)
              simulation%domain(m)%time = -1.d0
            endif
          endif
          simulation%config%time_from_call = simulation%domain(m)%time
          simulation%config%iter_from_call = 0
          simulation%config%dump_dt = obj_io%sol_dtime   ! dump-clip period (a runtime re-read of sol-dtime is NOT followed)
        enddo
        call Setup_Device_Residency ( simulation%domain(1) )
        call Fill_Ghost_Cell_Setup ( simulation%domain(1) )

        ! Interpolate initial solution on coarser domains
        if ( simulation%config%MGL > 1 ) then
          do m = 2, simulation%config%MGL
            call Restriction ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
            call Fill_Ghost_Cell_Setup ( simulation%domain(m) )
          end do
        endif

    end select

    ! Deallocate heavy arrays from non-local blocks to save memory.
    ! Must run AFTER Setup_Metrics (needs node/M/dl on all blocks),
    ! Setup_Probes (needs node on all blocks), and solver init
    ! (Restriction needs vol/P on all blocks).
    do m = 1, simulation%config%MGL
      call deallocate_remote_computation_data(simulation%domain(m))
    end do

    ! Calculate time at beginning of simulation
    call Cpu_Time ( simulation%config%cputime(1) )


    contains


      subroutine Print_Header()

        write(*,*)
        write(*,'(A89)')'_/\\\\____________/\\\\________/\\\\\___________/\\\\\\\\\\\_____/\\\\\\\\\\\\\\\__      '
        write(*,'(A89)')'_\/\\\\\\________/\\\\\\______/\\\///\\\_______/\\\/////////\\\__\/\\\///////////___     '
        write(*,'(A89)')' _\/\\\//\\\____/\\\//\\\____/\\\/__\///\\\____\//\\\______\///___\/\\\______________    '
        write(*,'(A89)')'  _\/\\\\///\\\/\\\/_\/\\\___/\\\______\//\\\____\////\\\__________\/\\\\\\\\\\\______   '
        write(*,'(A89)')'   _\/\\\__\///\\\/___\/\\\__\/\\\_______\/\\\_______\////\\\_______\/\\\///////_______  '
        write(*,'(A89)')'    _\/\\\____\///_____\/\\\__\//\\\______/\\\___________\////\\\____\/\\\______________ '
        write(*,'(A89)')'     _\/\\\_____________\/\\\___\///\\\__/\\\______/\\\______\//\\\___\/\\\______________'
        write(*,'(A89)')'      _\/\\\_____________\/\\\_____\///\\\\\/______\///\\\\\\\\\\\/____\/\\\\\\\\\\\\\\\_'
        write(*,'(A89)')'       _\///______________\///________\/////__________\///////////______\///////////////_'
        write(*,*)
        write(*,*)
        write(*,'(A89)') "  MO•dern  /'mɒd.ən/  n. 1. Relating to the present or recent times as opposed to the    "
        write(*,'(A89)') '  remote past.                                                                           '
        write(*,'(A89)') '  2. Characterized by or using the most up-to-date techniques, ideas, or equipment.      '
        write(*,'(A89)') '  3. Denoting the form of language that is currently used as opposed to any earlier form.'
        write(*,'(A89)') "  SE•quel  /'si:kw(ə)l/  n. 1. A published, broadcast, or recorded work that continues   "
        write(*,'(A89)') '  the story or develops the theme of an earlier one.                                     '
        write(*,'(A89)') '  2. Something that takes place after or as a result of an earlier event.                '
        write(*,*)
        write(*,'(A89)') '                                                               Oxford English Dictionary '
        write(*,*)
                                           
      end subroutine Print_Header


      subroutine Print_Shell_Info( config )
        use MOSE_Global_m
        implicit none
        type(config_type), intent(in) :: config
        ! Local
        character(llen) :: eosword

        write(*,*)
        write(*,'(A)')   ' IO specifications:'
        write(*,'(2A)')  ' - Initial Conditions --> ',trim(obj_io%nameinit)
        write(*,'(4A)')  ' - Backup format      --> ',trim(obj_io%bck_format(1)),' ',trim(obj_io%bck_format(2))
        write(*,'(4A)')  ' - Solution format    --> ',trim(obj_io%sol_format(1)),' ',trim(obj_io%sol_format(2))

        eosword = 'Ideal'

        write(*,*)
        write(*,'(A)')   ' Gas model:'
        write(*,'(2A)')  ' - Equation of state  --> ', trim(eosword)
        write(*,'(A,I4)')' - Number of species  -->', nsc
        write(*,'(A,I4)')' - Momentum equations -->', nw - nsc
        write(*,*)

        write(*,'(A)')   ' RANS specifications:'
        write(*,'(2A)')  ' - Turbulence model   --> ', trim(obj_rans%rans_name)
        write(*,*)

        if (obj_chem%reactive) then
          write(*,'(A)') ' Chemistry specifications:'
          write(*,*) ' - Mechanism          --> ', trim(obj_chem%mechanism_name)
          write(*,*) ' - Strang operator splitting'
          write(*,*) ' - ODE solver         --> ', trim(obj_chem%ode_name)
           if ( size(obj_chem%no_chem_list) == 0 ) then
            write(*,*) ' - Chemistry active in all blocks'
          else
            write(*,'(A, I0, *( " ", I0))') ' - Blocks without chemistry --> ', obj_chem%no_chem_list
          endif
          write(*,*)
        endif

        if (obj_soot%is_on) then
          write(*,'(A)') ' Soot model specifications:'
          if (trim(obj_soot%model_name)=='LL91') then
            write(*,*) ' - Mechanism          --> Leung (1991)'
          elseif (trim(obj_soot%model_name)=='LIN') then
            write(*,*) ' - Mechanism          --> Lindsted (1994)'
          end if
          if (nsoot==1) then
            write(*,*) ' - C(gr) => Heavy Gas'
          endif
          write(*,*)
        end if

        if (npass>0) then
            write(*,'(A,I5)')' - Number of Passive scalars  --> ', npass
        endif

        write(*,'(A)')' Numerical scheme:'
        if (index(trim(config%limiter_name),'IORD')>0) then
          write(*,'(A)')' - Space  --> I order'
        else
          write(*,'(2A)')' - Space  --> II order with ',trim(config%limiter_name)
        endif
        if (simulation%config%solver_name=='explicit') then
          if (obj_rk%n_RK==1) then
            write(*,'(A)',advance='no') ' - Time   --> explicit Euler'
          else
            write(*,'(2A,I2)',advance='no') ' - Time   --> ', trim(config%solver_name), obj_rk%n_RK
          endif
          write(*,'(3A)')' over the ',trim(config%intvar_name),' variables'
        else
          write(*,'(2A,I2)') ' - Time   --> implicit Euler Newton-Krylov'
        endif
        write(*,'(2A)')' - Fluxes --> ',trim(config%riemann_name)
        write(*,*)

        if ( simulation%config%MGL > 1 ) then
          write(*,'(A21,I1,A8)') ' Full Multigrid with ', config%MGL, ' levels.'
          write(*,*)
        endif

    end subroutine Print_Shell_Info

  
  end subroutine MOSE_setup


end module MOSE_Wrap_Setup
