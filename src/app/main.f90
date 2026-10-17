program MOSE_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use MOSE_Procedures_m
  use MOSE_Types_m
  use MOSE_Mod_MPI
#ifdef USE_MPI
  use MOSE_Mod_GhostExchange, only: cleanup_ghost_schedule
#endif

  implicit none
  type(MOSE_type)            :: MOSE
  type(MOSE_simulation_type) :: simulation

  ! --- startup/IO instrumentation (root rank wall clock) ---
  integer(8) :: tk0, tk1, tkrate
  real(8)    :: t_setup, t_solve, t_post

  ! Initialize MPI environment (no-op if USE_MPI is not defined)
  call mpi_init_env()
  call system_clock(count_rate=tkrate)

#if defined (_OPENMP)
  !$omp parallel
  simulation%config%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  if (mpi_is_root) then
    write(*,'(A)')   ' OpenMP:'
    write(*,'(A,I4)')' -  Number of threads --> ', simulation%config%nthreads
  end if
#else
  if (mpi_is_root) write(*,'(A)')   ' Serial execution'
  simulation%config%nthreads = 1
#endif

#ifdef USE_MPI
  if (mpi_is_root) then
    write(*,'(A)')      ' MPI:'
    write(*,'(A,I4)')   ' -  Number of ranks   --> ', mpi_size_
  end if
#endif

  ! Solving with MOSE
  call system_clock(tk0)
  call MOSE%setup( simulation )
  call system_clock(tk1)
  t_setup = real(tk1-tk0,8)/real(tkrate,8)
  if (mpi_is_root) write(*,'(A,F12.3,A)') ' [TIMING] setup (mesh+ic+bc read) = ', t_setup, ' s'

  t_solve = 0d0; t_post = 0d0
  simulation%config%TODO = 1
  do while ( simulation%config%TODO <= 2 )
    call system_clock(tk0)
    call MOSE%solve( simulation, Dummy_Function )
    call system_clock(tk1); t_solve = t_solve + real(tk1-tk0,8)/real(tkrate,8)
    if ( simulation%config%TODO <= 2 ) then
      call system_clock(tk0)
      call MOSE%postprocess( simulation )
      call system_clock(tk1); t_post = t_post + real(tk1-tk0,8)/real(tkrate,8)
    end if
  enddo

  call system_clock(tk0)
  call MOSE%postprocess( simulation )
  call system_clock(tk1); t_post = t_post + real(tk1-tk0,8)/real(tkrate,8)
  if (mpi_is_root) then
    write(*,'(A,F12.3,A)') ' [TIMING] solve  (compute loop)   = ', t_solve, ' s'
    write(*,'(A,F12.3,A)') ' [TIMING] post   (output writes)  = ', t_post,  ' s'
  end if
#ifdef USE_MPI
  ! Per-rank loss-decomposition probes (multigpu_plan M0): the solve spread across
  ! ranks is the imbalance share, halo-wait is the true comm-wait share.
  if (mpi_size_ > 1) write(*,'(A,I4,A,F12.3,A,F12.3,A)') &
    ' [TIMING] rank ', mpi_rank_, ' solve = ', t_solve, ' s  halo-wait = ', halo_wait_time, ' s'
#endif

  ! Free persistent MPI requests before finalizing
#ifdef USE_MPI
  call cleanup_ghost_schedule()
#endif

  ! Finalize MPI environment
  call mpi_finalize_env()


contains

    subroutine Dummy_Function

    end subroutine Dummy_Function


end program MOSE_program
