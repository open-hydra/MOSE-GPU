module MOSE_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none

  integer:: id_stampa=0
  character(81) :: RANS_shell_format  =    "('MOSE | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(69) :: URANS_shell_format =    "('MOSE | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"
  character(96) :: RANS_shell_format_MG  = "('MOSE Grid Level', i2, ' | Iter =', i9, ' | Global iter =', i9, ' | Density residual =', E13.6)"
  character(84) :: URANS_shell_format_MG = "('MOSE Grid Level', i2, ' | Iter =', i9, ' | Time =', E13.6,  ' | Delta t =', E13.6)"

contains


  subroutine MOSE_postprocess ( simulation )
    use IR_precision
    use MOSE_Types_m
    use MOSE_IO_Ini
    use MOSE_Global_m
    use MOSE_Parameters_m
    use MOSE_IO_Solution
    use MOSE_IO_Wall
    use MOSE_Mod_Diagnostic
    use MOSE_Mod_Probes
    use MOSE_Lib_Interpolation
    use MOSE_Lib_RotatingFrame, only: RF_Torque_Walls
    use MOSE_Mod_MPI, only: mpi_is_root

    implicit none
    type(MOSE_simulation_type), intent(inout) :: simulation
    ! Local
    real(R8)        :: sim_time
    character(llen) :: solfile, bckfile, dgsfile, wallfile, mgfile
    integer         :: m, level
    logical         :: skip_bck_
    
    level = simulation%config%MG_level
    skip_bck_ = .false.

    ! Update IOfield solutiontime
    if (simulation%config%time_accurate) then
      simulation%IOfield(1)%solutiontime = simulation%domain(1)%time
      IOwall%solutiontime = simulation%domain(1)%time
    else
      simulation%IOfield(1)%solutiontime = -real(simulation%config%iter_general,8)
      IOwall%solutiontime = -real(simulation%config%iter_general,8)
    endif

    if (level == 1) then
      solfile  = '/'//trim(MOSE_phase_prefix)//'field'
      bckfile  = '/'//trim(MOSE_phase_prefix)//'field'
      wallfile = '/'//trim(MOSE_phase_prefix)//'wall'
      dgsfile  = '/'//trim(MOSE_phase_prefix)//'diagnostic'
    else
      solfile  = '/'//trim(MOSE_phase_prefix)//'field-level'//trim(str(.true.,level))
      bckfile  = '/'//trim(MOSE_phase_prefix)//'field-level'//trim(str(.true.,level))
      mgfile   = '/'//trim(MOSE_phase_prefix)//'field-prolongated'
    endif

    ! END OF SIMULATION
    if ( simulation%config%TODO == 3 ) then

      if (mpi_is_root) then
        ! Residual history
        if ( mod (simulation%domain(1) % iter, obj_io%res_diter) == 0d0 ) then
          write (unitRES,trim(unitRES_format)) simulation%config%iter_general, simulation%domain(1)%time, simulation%config%residuotot
          if (ROTATING_FRAME .and. level == 1) &
            call RF_Torque_Walls( simulation%domain(1), simulation%config%iter_general, simulation%domain(1)%time )
        endif

        ! Shell
        if (simulation%config%time_accurate) then
          if (simulation%config%MGL > 1 ) then
            write (*,URANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          else
            write (*,URANS_shell_format) simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
          endif
        else
          if (simulation%config%MGL > 1 ) then
            write (*,RANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
          else
            write (*,RANS_shell_format) simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
          endif
        endif

        ! Calculate time at end of simulation
        call Cpu_Time ( simulation%config%cputime(2) )

        sim_time = ( simulation%config%cputime(2) - simulation%config%cputime(1) ) / simulation%config%nthreads
        write(*,*)
        write(*,*) '  Time of operation was', sim_time/60, 'min'
      end if

      ! MOSE_STATE_HASH: the device state hash before anything is written (2 = no final backup: gates without a field file;
      ! 3 = also at every iteration, see hash_print_)
      block
        character(len=8) :: hev_
        integer :: hst_, hmode_
        hmode_ = 0
        call get_environment_variable ( 'MOSE_STATE_HASH', hev_, status=hst_ )
        if ( hst_ == 0 ) read ( hev_, *, iostat=hst_ ) hmode_
        if ( hst_ /= 0 ) hmode_ = 0
        if ( hmode_ > 0 ) call hash_print_ ( simulation%domain(1)%nb, 'end', simulation%config%iter_from_call )
        skip_bck_ = ( hmode_ == 2 )
      end block
      ! Write output solution
      if ( .not. skip_bck_ ) call Write_Bck ( simulation%domain(1), simulation%IOfield(1), bckfile )
      call Write_Solution ( simulation%domain(1), simulation%IOfield(1), solfile )
      if (WALL) call Write_Wall_Solution( simulation%domain(1), wallfile, obj_io%sol_format )
      if (.not. simulation%config%time_accurate) call Write_Diagnostic ( simulation%domain(1), simulation%IOfield(1), dgsfile, obj_io%sol_format )
        

    
    ! INTERMEDIATE SOLUTION EVALUATION
    elseif ( simulation%config%TODO == 2 ) then

      if ( level == 1 ) then   ! MOSE_STATE_HASH=3: state hash at every iteration (A6 restart diagnosis)
        block
          character(len=8) :: hev_
          integer :: hst_, hmode_
          hmode_ = 0
          call get_environment_variable ( 'MOSE_STATE_HASH', hev_, status=hst_ )
          if ( hst_ == 0 ) read ( hev_, *, iostat=hst_ ) hmode_
          if ( hst_ == 0 .and. hmode_ >= 3 ) call hash_print_ ( simulation%domain(1)%nb, 'iter', simulation%config%iter_from_call )
        end block
      end if

      if (level == 1) then
        id_stampa = id_stampa + 1
        if (.not.obj_io%sol_overwrite) solfile = trim(solfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) dgsfile = trim(dgsfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%bck_overwrite) bckfile = trim(bckfile)//trim(str(.true.,id_stampa))
        if (.not.obj_io%sol_overwrite) wallfile = trim(wallfile)//trim(str(.true.,id_stampa))
      endif

      ! Probes (each rank writes its own local probes)
      if (level == 1) call Write_Probes_Data( simulation%config%iter_general, simulation%domain(1)%time)

      if (mpi_is_root) then
        ! Residuals history
        if ( mod(simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (unitRES,trim(unitRES_format)) simulation%config%iter_general, simulation%domain(level)%time, simulation%config%residuotot
          if (ROTATING_FRAME .and. level == 1) &
            call RF_Torque_Walls( simulation%domain(1), simulation%config%iter_general, simulation%domain(1)%time )
        endif

        ! Shell
        if ( mod (simulation%domain(level) % iter, simulation%config%diter_shell) == 0d0 ) then
          if (simulation%config%time_accurate) then
            if (simulation%config%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (simulation%config%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
            else
              write (*,RANS_shell_format) simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
            endif
          endif
        endif
      end if

      ! Iter-based solution
      if ( mod(simulation%domain(level) % iter, obj_io%sol_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )

        if (simulation%config%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = simulation%config%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgfile )
        endif
        if (level == 1) then
          if (WALL) call Write_Wall_Solution( simulation%domain(level), wallfile, obj_io%sol_format )
          if (.not. simulation%config%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile, obj_io%sol_format )
        endif

      ! Time-based solution
      elseif ( simulation%domain(level) % time >= simulation%config % time_from_call + obj_io%sol_dtime ) then
        if ( obj_io%sol_dtime_grid ) then   ! sol-dtime-grid: next threshold on the fixed grid, the offset does not accumulate
          simulation%config % time_from_call = simulation%config % time_from_call + obj_io%sol_dtime * &
            real( floor( ( simulation%domain(level) % time - simulation%config % time_from_call ) / obj_io%sol_dtime ), 8 )
        else
          simulation%config % time_from_call = simulation%domain(level) % time
        end if
        obj_io%n_tdump = obj_io%n_tdump + 1
        if ( mod(obj_io%n_tdump, obj_io%sol_field_every) == 0 ) then   ! sol-field-every: 3D field every N-th dump
        if (mpi_is_root) write(*,*) ' ... writing time-based solution'
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        else
        if (mpi_is_root) write(*,*) ' ... writing time-based wall file only (sol-field-every)'
        end if
        if (simulation%config%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = simulation%config%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgfile )
        endif
        if (level == 1) then
          if (WALL) call Write_Wall_Solution ( simulation%domain(level), wallfile, obj_io%sol_format )
          if (.not. simulation%config%time_accurate) call Write_Diagnostic ( simulation%domain(level), simulation%IOfield(level), dgsfile, obj_io%sol_format )
        endif
      endif

      ! Save backup solution (only iter-based)
      if ( mod(simulation%domain(level) % iter, obj_io%bck_diter) == 0d0 ) then
        if (mpi_is_root) write(*,*) ' ... writing iter-based backup'
        call Write_Bck ( simulation%domain(level), simulation%IOfield(level), bckfile )
        if (simulation%config%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = simulation%config%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Bck (  simulation%domain(1), simulation%IOfield(1), mgfile )
        endif
      endif

      if (simulation%config%change_MG) then
        if (mpi_is_root) write(*,*) ' ... writing solution of grid-level ', level
        call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
        call Write_Bck ( simulation%domain(level), simulation%IOfield(level), bckfile )
        if (simulation%config%MGL > 1 .and. level > 1) then
          if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
          do m = simulation%config%MG_level, 2, -1
            call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
          enddo
          call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgfile )
          call Write_Bck (  simulation%domain(1), simulation%IOfield(1), mgfile )
        endif
        simulation%config%MG_level = simulation%config%MG_level - 1
        simulation%config%iter_from_call = 0
        simulation%config%time_from_call = 0.0
      endif


    ! AUXILIARY OPERATIONS AT EACH ITERATION
    elseif ( simulation%config%TODO == 1 ) then

      ! Probes (each rank writes its own local probes)
      if (level == 1) call Write_Probes_Data( simulation%config%iter_general, simulation%domain(1)%time )

      if (mpi_is_root) then
        if ( mod (simulation%domain(level) % iter, obj_io%res_diter) == 0d0 ) then
          write (unitRES,trim(unitRES_format)) simulation%config%iter_general, simulation%domain(level)%time, simulation%config%residuotot
          if (ROTATING_FRAME .and. level == 1) &
            call RF_Torque_Walls( simulation%domain(1), simulation%config%iter_general, simulation%domain(1)%time )
        end if

        ! Shell
        if ( mod (simulation%domain(level) % iter, simulation%config%diter_shell) == 0d0 ) then
          if (simulation%config%time_accurate) then
            if (simulation%config%MGL > 1 ) then
              write (*,URANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            else
              write (*,URANS_shell_format) simulation%config%iter_from_call, simulation%domain(1)%time, simulation%domain(1)%dtglobal
            endif
          else
            if (simulation%config%MGL > 1 ) then
              write (*,RANS_shell_format_MG) level, simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
            else
              write (*,RANS_shell_format) simulation%config%iter_from_call, simulation%config%iter_general, simulation%config%residuotot(1)
            endif
          endif
        endif
      end if

    endif

    ! Update input data
    if ( mod (simulation%domain(level) % iter, simulation%config%diter_ini) == 0d0 ) then
      call Read_Inifile_Runtime ( simulation%config )
    end if
    

    ! PRINT FOR HYDRA-COUPLING PROBLEMS
    if ( simulation%config%HYDRA_postprocess ) then
      if (mpi_is_root) write(*,*) ' ... writing HYDRA-coupling solution'
      call Write_Solution (  simulation%domain(level), simulation%IOfield(level), solfile )
      if (simulation%config%MGL > 1 .and. level > 1) then
        if (mpi_is_root) write(*,*) ' ... writing prolongated solution to grid-level 1'
        do m = simulation%config%MG_level, 2, -1
          call Prolongation ( Fine=simulation%domain(m-1), Coarse=simulation%domain(m) )
        enddo
        call Write_Solution (  simulation%domain(1), simulation%IOfield(1), mgfile )
      endif
      if (WALL .and. level == 1) call Write_Wall_Solution( simulation%domain(1), wallfile, obj_io%sol_format )
    endif

   
  end subroutine MOSE_postprocess


  !> MOSE_STATE_HASH: global (MPI-summed) and per-rank hash of the interior device state w, tagged (A6, 16/09).
  subroutine hash_print_ ( nb, tag, iter )
    use MOSE_Field_m, only: field_state_hash, field_is_allocated
    use MOSE_Mod_MPI, only: mpi_is_root, mpi_rank_
#ifdef USE_MPI
    use mpi
#endif
    integer, intent(in) :: nb, iter
    character(len=*), intent(in) :: tag
    integer :: ierr_
    integer(8) :: h1_, h2_, n_, g1_, g2_, gn_
    if ( .not. field_is_allocated() ) return
    call field_state_hash ( nb, h1_, h2_, n_ )
    g1_ = h1_; g2_ = h2_; gn_ = n_
#ifdef USE_MPI
    call MPI_ALLREDUCE ( h1_, g1_, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr_ )
    call MPI_ALLREDUCE ( h2_, g2_, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr_ )
    call MPI_ALLREDUCE ( n_,  gn_, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr_ )
#endif
    if ( mpi_is_root ) write(*,'(A,A,A,I0,A,I0,A,I0,A,I0)') ' [HASH] state w interior ', trim(tag), ' ', iter, ': n ', gn_, ' lo ', g1_, ' hi ', g2_
    if ( tag == 'end' ) write(*,'(A,I0,A,I0,A,I0,A,I0)') ' [HASH] rank ', mpi_rank_, ' n ', n_, ' lo ', h1_, ' hi ', h2_
  end subroutine hash_print_


end module MOSE_Wrap_Postprocess