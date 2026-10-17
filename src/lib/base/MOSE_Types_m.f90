module MOSE_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Lib_ORION_data
  use MOSE_General_Types_m
  use MOSE_Parameters_m
  use MOSE_Mod_BC_Time

  implicit none
  
  private :: read_ini

  type :: block_type
    integer                                :: dim(3)           ! Number of cells in i-j-k (ghost not included)
    real(R8), allocatable                  :: vol(:,:,:)       ! Cell volume
    type(MOSE_vector_3D_type), allocatable :: node(:,:,:)      ! Mesh grid points (including ghost)
    type(MOSE_tensor_3D_type), allocatable :: M(:,:,:)         ! Metric transformation tensor
    ! Flat SoA mirror of M%c for device-mappable, coalesced access (GPU Phase 2,
    ! P2.0). m_f(:,:,i,j,k) == M(i,j,k)%c ; populated in Setup_Metrics.
    real(R8), allocatable                  :: m_f(:,:,:,:,:)   ! (3,3,i,j,k) flat metric tensor
    type(MOSE_vector_3D_type), allocatable :: dl(:,:,:)        ! Average cell length (in i/j/k direction). eg: dl%c(1) is sqrt(dx**2+dy**2+dz**2) of the cell in the i direction
    real(R8), allocatable                  :: dl_f(:,:,:,:)    ! (3,i,j,k) flat mirror of dl%c
    type(d_metrics_type)                   :: dir(3)           ! Direction object. Contains: i-faces, j-faces, k-faces; eg: dir(1)%face(i,j,k)%n
    ! Flat SoA mirror of dir(d)%f: fm<d>(1:3,..)=face normal,
    ! fm<d>(4,..)=face area. Block-level (not dir% component) so they map cleanly
    ! to device. Asymmetric per-direction bounds match dir(d)%f.
    real(R8), allocatable                  :: fm1(:,:,:,:)     ! (4, 0:ni, nj, nk)
    real(R8), allocatable                  :: fm2(:,:,:,:)     ! (4, ni, 0:nj, nk)
    real(R8), allocatable                  :: fm3(:,:,:,:)     ! (4, ni, nj, 0:nk)
    real(R8), allocatable                  :: yn(:,:,:)        ! Nearest wall distance
  end type block_type


  type :: bc_type
    integer                    :: i, j, k, b, f                             ! ijk coordinates, block and face in which the boundary element is located
    integer                    :: type                                      ! BC type (1,9)
    integer                    :: bs, is, js, ks, fs, d11, d12, d21, d22    ! BC 1 (connection) specifications
    type(MOSE_tensor_3D_type)  :: Mg(2)                                     ! Ghost cell metric tensor
    type(MOSE_vector_3D_type)  :: dlg(2)                                    ! Ghost cell average cell length
    real(R8)                   :: volg(2)                                   ! Ghost cell volume
    real(R8), allocatable      :: Pg(:,:)                                   ! Ghost cell primitive stencil
    integer                    :: ni(2)                                     ! BC chimera
    integer, allocatable       :: donorID(:,:)                              ! BC chimera
    real(R8), allocatable      :: volume_fraction(:)                        ! BC chimera
    real(R8), allocatable      :: ext_flux(:)                               ! Multi-Solver Coupling
  end type bc_type


  ! Perfect gas extension for a block.
  type, extends(block_type) :: MOSE_block_type
    real(R8), dimension(:,:,:,:), allocatable :: P, PO                 ! Primitive variables at time n and n-1: { rho(s) vel p rho*r }, r generic RANS variable
    real(R8), dimension(:,:,:,:), allocatable :: R                     ! Residuals
    real(R8), dimension(:,:,:,:), allocatable :: TE                    ! Truncation error
    real(R8), dimension(:,:,:,:), allocatable :: RC                    ! Conservation form residuals (temporary storage)
    real(R8), dimension(:,:,:,:), allocatable :: RS1, RS2              ! Implicit smoothing residuals (temporary storage)
    real(R8), dimension(:,:,:), allocatable   :: dtlocal               ! Local time step
    real(R8), dimension(:,:,:), allocatable   :: beta                  ! Shock detector flag
    type(MOSE_tensor_3D_type), allocatable    :: vel_gradient(:,:,:)   ! Gradient of velocity
    real(R8), allocatable                     :: rc_term1(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
    real(R8), allocatable                     :: rc_term2(:,:,:)       ! Spalart-Shur rotation/curvature correction terms
    ! Per-cell thermo scratchpad (refreshed once per RK stage by Mod_Thermo_Cache).
    ! Avoids redundant co_rotot_Rtot / T = p/(rho R) / f_laminarViscosity calls
    ! across Mod_Diffusive, Mod_dt, Mod_SST, Lib_Newstate.
    real(R8), dimension(:,:,:), allocatable   :: rho_c, T_c, Rg_c      ! density, temperature, gas constant per cell
    real(R8), dimension(:,:,:), allocatable   :: hr_c                ! |dT|/T of the last chemistry step per cell (heat-release measure, j-implicit guard)
    real(R8), dimension(:,:,:), allocatable   :: mil_c, kl_c, cp_c     ! laminar viscosity, conductivity, mixture cp per cell (VISC)
    logical                                   :: no_chem = .false.     ! Per disattivare la chimica su un unico blocco
  end type MOSE_block_type


  type, extends(bc_type) :: MOSE_bc_type
    real(R8)              :: qw, Tw, Taw, hg, qrad             ! BC viscous wall specifications   
    real(R8)              :: T0, p0, alpha, beta, mach, pamb   ! BC 4 (inflow/outflow) specifications
    real(R8)              :: mdot                              ! BC 4 (inflow/outflow) specifications
    real(R8)              :: un                                ! BC 25 (prescribed-velocity inlet) normal velocity
    real(R8)              :: rel_fac                           ! BC 4 (inflow/outflow) specifications
    real(R8), allocatable :: ci(:)                             ! BC 4 (inflow/outflow) specifications
    real(R8)              :: SF                                ! BC ? specifications
    real(R8)              :: aCoeff, n, pRef, Taf, haf         ! BC 14 specification (SRM grain)
    type(time_bc_type)    :: p0time
    type(time_bc_type), allocatable :: q2d_map(:)              ! BC 667 (Q2D mapped) time-varying primitive data (one per primitive)
    character(len=256) :: q2d_file = ' '                       ! BC 667 (Q2D mapped) Tecplot filename from bc.txt
    logical            :: q2d_periodic = .false.                ! BC 667 (Q2D mapped) periodic time signal flag
    ! MOSE_open (ATLAS) dialect provenance, filled by MOSE_IO_BC_Open when the bc
    ! file uses 3-digit ids. Zero/false for a legacy file. Kept for diagnostics and
    ! for the post-metrics checks that need the original record semantics.
    integer            :: open_type = 0                        ! original MOSE_open id (0 = legacy dialect)
    logical            :: alpha_from_normal = .false.          ! the 'normal' token, per component (MOSE_open parse_dir_tok):
    logical            :: beta_from_normal  = .false.          ! resolved from the face normal after the metrics (Resolve_Open_BC_Direction)
  end type MOSE_bc_type


  type :: MOSE_domain_type
    real(R8)                                         :: time             ! Solution time
    real(R8)                                         :: dtglobal         ! Global dt for time accurate simulation
    integer                                          :: iter, itermax    ! Iteration number
    integer                                          :: nb, nbound       ! Number of blocks, number of boundary faces
    integer, dimension(:,:), allocatable             :: n_bf             ! Number of bc elements per faces per block
    type(MOSE_block_type), dimension(:), allocatable :: blk              ! Allocatable block type
    type(MOSE_bc_type), dimension(:), allocatable    :: bc               ! Allocatable bc object
    ! MPI: pre-filtered BC index lists for local rank (avoid scanning all nbound)
    integer                                          :: n_local_bc = 0   ! Number of BC entries with local Bm
    integer, dimension(:), allocatable               :: local_bc_idx     ! Indices into bc(:) where Bm is local
    integer                                          :: n_local_bs = 0   ! Number of type-1 BC entries with local Bs
    integer, dimension(:), allocatable               :: local_bs_idx     ! Indices into bc(:) where Bs is local (for Pg computation)
    ! Flat SoA mirror of the boundary-face list. bc(:) is a
    ! derived type with allocatable components (the device-mapping blocker, same
    ! AoS problem as the geometry); the device BC kernel reads these plain int
    ! arrays instead. Populated once at setup by Setup_BC_Device, in the same
    ! (block,face) order the host BC_Fluxes driver walks. bc_blk_lo/hi give each
    ! block's contiguous index range so the per-block device kernel can launch a
    ! face loop without scanning all nbound.
    integer, dimension(:), allocatable               :: bc_bi, bc_bj, bc_bk  ! per-face cell ijk (size nbound)
    integer, dimension(:), allocatable               :: bc_bf, bc_typ        ! per-face local face id (1..6) and BC type
    integer, dimension(:), allocatable               :: bc_blk_lo, bc_blk_hi ! per-block [lo,hi] range into the above (size nb)
    integer, dimension(:), allocatable               :: bc_perm             ! sorted position i -> source bc(:) index (block-major,face-major; see Setup_BC_Device)
    ! Flat per-face inflow/outflow (type 4) scalars + the
    ! tangential face vector (precomputed from node geometry). Filled for type-4
    ! faces by Setup_BC_Device; bc_p0 is refreshed per step for time-varying inlets.
    real(R8), dimension(:),   allocatable            :: bc_mach, bc_t0, bc_p0, bc_relfac
    real(R8), dimension(:),   allocatable            :: bc_alpha, bc_beta, bc_pamb
    real(R8), dimension(:,:), allocatable            :: bc_ci    ! (nsc,   nbound) inflow species mass fractions
    real(R8), dimension(:,:), allocatable            :: bc_rans  ! (nrans, nbound) inflow RANS variables
    real(R8), dimension(:,:), allocatable            :: bc_tvec  ! (3,     nbound) unit tangential face vector
    real(R8), dimension(:),   allocatable            :: bc_tw    ! prescribed wall temperature [K] (type-6 faces)
    ! Flat mirror of bc(:)%Pg (donor ghost stencil
    ! for type-1 connection faces). Dynamic — refreshed each step in BC_Fluxes
    ! from the host ghost-fill, since it carries the donor block's solution.
    real(R8), dimension(:,:,:), allocatable          :: bc_pg    ! Flat per-entry donor/receiver
    ! index tables for type-1 connections with BOTH blocks local (single-rank), grouped
    ! by ordered (Bm,Bs) block pair so ONE kernel per pair fills cross-block (reads
    ! blk(Bs)%P, writes blk(Bm)%P ghost + bc_pg) with both block P arrays present-matched.
    integer                                          :: n_conn = 0, n_conn_pair = 0
    integer                                          :: pair_xhits = -1      ! transverse reads of a ghost written by ANOTHER connection pair (-1 = census not run)
    integer, dimension(:), allocatable               :: conn_im, conn_jm, conn_km, conn_fm   ! receiver face cell + face (size n_conn)
    integer, dimension(:), allocatable               :: conn_is, conn_js, conn_ks, conn_fs   ! donor face cell + face
    integer, dimension(:), allocatable               :: conn_d11, conn_d12, conn_d21, conn_d22 ! transverse Pg transforms
    integer, dimension(:), allocatable               :: conn_pg                              ! sorted bc_pg slot for this entry
    integer, dimension(:), allocatable               :: conn_pair_bm, conn_pair_bs           ! the two blocks of pair p (size n_conn_pair)
    integer, dimension(:), allocatable               :: conn_pair_lo, conn_pair_hi           ! [lo,hi] range into the conn_* entry arrays
    integer, dimension(:), allocatable               :: conn_all, conn_all_bs                ! every entry 1..n_conn + its donor block (pass 2)
    ! Chimera (666) device ghost-fill tables (single-rank, donors local). All cell
    ! references are GLOBAL 1-based flat indices into the var-last w arena (cidx), so the
    ! gather kernel needs no per-block offsets and donors may span any number of blocks.
    ! Built by Setup_BC_Device_Flat (needs the field arena allocated).
    integer                                          :: n_chim = 0
    integer, dimension(:), allocatable               :: chim_g1, chim_g2      ! ghost row 1 / row 2 cell (size n_chim)
    integer, dimension(:), allocatable               :: chim_rc               ! receptor boundary cell (T guess)
    integer, dimension(:), allocatable               :: chim_dlo, chim_dmid, chim_dhi  ! donor ranges: row1=[dlo,dmid], row2=[dmid+1,dhi]
    integer, dimension(:), allocatable               :: chim_slot             ! sorted bc_pg slot of this entry
    integer, dimension(:), allocatable               :: chim_df               ! donor cells, concatenated (P_flat idx)
    real(R8), dimension(:), allocatable              :: chim_vf               ! donor volume fractions (matches chim_df)
    integer, dimension(:,:), allocatable             :: chim_tn               ! (4,n_chim) receptor's transverse neighbours -> bc_pg(:,3:6)
    integer, dimension(:,:), allocatable             :: chim_cg               ! (4,n_chim) corner ghost-of-ghost cells at block edges, copied from ghost row 1 (0 = none)
    ! Manifold (1001) mean-connection device tables: one entry per contiguous run of
    ! 1001 slots sharing (block,face,bs,fs). The per-substage source-face reduction
    ! (BC_Manifold_Update) refreshes bc_t0/bc_p0/bc_pamb for the slot range
    ! [mf_lo,mf_hi]; static inflow params are preset in Setup_BC_Device. Source-face
    ! cells are stored as GLOBAL P_flat indices (mf_cf) + their face areas (mf_af,
    ! from host dir%f%A -- the exact values host BC_Manifold sums), in the host loop
    ! order (ks,js,is; is fastest) so the device serial reduction is bit-ordered
    ! like the host.
    integer                                          :: n_mfold = 0
    integer, dimension(:), allocatable               :: mf_clo, mf_chi   ! per-face range into mf_cf/mf_af
    integer, dimension(:), allocatable               :: mf_lo, mf_hi     ! sorted slot range of the receptor entries
    real(R8), dimension(:), allocatable              :: mf_aM            ! receptor-face total area (static host dir%f%A sum)
    integer, dimension(:), allocatable               :: mf_cf            ! source-face cells (P_flat idx), concatenated
    real(R8), dimension(:), allocatable              :: mf_af            ! matching face areas
    real(R8), dimension(:,:), allocatable            :: mf_tmp           ! device scratch (max face cells, np+2): cons*A, T*A, A
  end type MOSE_domain_type


  type :: config_type
    logical         :: newrun           ! Flag di restart
    real(R8)        :: restop           ! Residuo min. per arresto esecuzione
    integer         :: itermax          ! Numero max iterate per arresto esecuzione
    real(R8)        :: timemax          ! Tempo max per arresto esecuzione 
    integer         :: diter_shell      ! Shell update
    integer         :: diter_ini        ! input.ini update
    integer         :: iter_general     ! Number of iteration - including all MG levels
    integer         :: iter_from_call
    real(R8)        :: time_from_call
    real(R8), allocatable :: residuotot(:)
    real(R8)        :: cfl              ! Parametro di stabilita' convettiva
    real(R8)        :: vnn              ! Parametro di stabilita' diffusiva
    integer         :: rampa_iter       ! Numero di iterazioni per rampa di vnn
    logical         :: time_accurate    ! Flag per integrazione time-accurate
    logical         :: time_clip = .false. ! Clip the last dt so the run lands exactly on time-threshold
    logical         :: dump_clip = .false. ! dump-clip: clip dt so the run lands exactly on every time-based output instant (A5, 16/09)
    real(R8)        :: dump_dt   = huge(1.0_R8) ! sol-dtime copied at setup (the output period used by dump-clip)
    integer         :: MGL              ! Number of multigrid levels
    integer         :: MG_level         ! Current MG-level being solved
    logical         :: change_MG        ! Logical variable to change MG-level
    character(llen) :: limiter_name     ! Solver description
    character(llen) :: riemann_name     ! Solver description
    character(llen) :: solver_name      ! Solver description
    character(llen) :: intvar_name      ! Solver description
    real(R8)        :: Minf             ! Solver description: free-stream Mach number (for AUSM+-up)
    real(R8)        :: Mco = 0.0_R8     ! riemann-options-Mco: the low-Mach cutoff of AUSM+M (MOSE_open obj_riemann%Mco)

    real(R8)        :: cputime(2)       ! Simulation time duration
    integer         :: nthreads         ! Number of threads for simulation
    integer         :: TODO             ! Decide solve and/or postprocess
    logical         :: HYDRA_time_accurate = .false.
    logical         :: HYDRA_postprocess   = .false.
    logical         :: HYDRA_MG = .false.

    contains

      procedure, pass(self) :: ini_setup => read_ini

  end type config_type


  type :: MOSE_simulation_type

    type(MOSE_domain_type), dimension(:), allocatable  :: domain
    type(ORION_data), dimension(:), allocatable        :: IOfield
    type(config_type)                                  :: config

  end type MOSE_simulation_type


  contains

    subroutine read_ini(self, fini)
      use Finer, only: file_ini
      use MOSE_Parameters_m, only: codename
      use MOSE_Global_m

      implicit none
      type(file_ini)     :: fini
      class(config_type) :: self
      ! Local
      integer :: error

      ! Newrun/restart
      call fini%get(section_name=trim(codename)//'-parameters', option_name='newrun', val=self%newrun, error=error)
      if (error/=0) then
        call fini%get(section_name=trim(codename)//'-parameters', option_name='restart', val=self%newrun, error=error)
        if (error==0) self%newrun = .not. self%newrun
      endif
      if (error/=0) self%newrun = .true.

      ! Stop criteria
      self%restop = 1d-10
      self%itermax = huge(1)
      self%timemax = huge(1d0)
      call fini%get(section_name=trim(codename)//'-parameters', option_name='res-threshold', val=self%restop, error=error)
      call fini%get(section_name=trim(codename)//'-parameters', option_name='time-threshold', val=self%timemax, error=error)
      call fini%get(section_name=trim(codename)//'-parameters', option_name='iter-threshold', val=self%itermax, error=error)

      ! Shell options
      self%diter_shell = 1
      call fini%get(section_name=trim(codename)//'-IO', option_name='shell-diter', val=self%diter_shell, error=error)
      self%diter_ini = 10000
      call fini%get(section_name=trim(codename)//'-IO', option_name='ini-diter', val=self%diter_ini, error=error)

      ! Stability coefficients and related options
      call fini%get(section_name=trim(codename)//'-parameters', option_name='cfl', val=self%cfl, error=error)
      if (error/=0) then
        call fini%get(section_name=trim(codename)//'-parameters', option_name='CFL', val=self%cfl, error=error)
        if (error/=0) then
          error stop ( '  CFL number not provided in ['//trim(codename)//'-parameters].' )
        end if
      end if
      self%vnn = 0.3d0
      call fini%get(section_name=trim(codename)//'-parameters', option_name='vnn', val=self%VNN, error=error)
      if (error/=0) then
        call fini%get(section_name=trim(codename)//'-parameters', option_name='VNN', val=self%VNN, error=error)
      end if
      self%rampa_iter = 0
      call fini%get(section_name=trim(codename)//'-parameters', option_name='cfl-rise-threshold', val=self%rampa_iter, error=error)

      ! Time-accurate switch
      self%time_accurate = .false.
      call fini%get(section_name=trim(codename)//'-parameters', option_name='time-accurate', val=self%time_accurate, error=error)

      ! time-clip: land exactly on time-threshold (last dt clipped). Default off so that
      ! existing references/oracles (never clipped) stay comparable step for step.
      self%time_clip = .false.
      call fini%get(section_name=trim(codename)//'-parameters', option_name='time-clip', val=self%time_clip, error=error)
      ! dump-clip (A5): land exactly on the time-based output instants too (time_from_call + sol-dtime). Default off.
      self%dump_clip = .false.
      call fini%get(section_name=trim(codename)//'-parameters', option_name='dump-clip', val=self%dump_clip, error=error)

      ! Solver-type
      self%solver_name = 'explicit'
      call fini%get(section_name=trim(codename)//'-scheme', option_name='solver-type', val=self%solver_name, error=error)
      if ( trim(self%solver_name) /= 'explicit' ) &
        error stop '[MOSE-scheme] solver-type: only explicit in this version (MOSE_open has no implicit solver)'

      ! Space discretization
      self%limiter_name = 'IORD'
      call fini%get(section_name=trim(codename)//'-scheme', option_name='limiter', val=self%limiter_name, error=error)

      ! Riemann solver
      self%riemann_name = 'HLLC Batten'
      call fini%get(section_name=trim(codename)//'-scheme', option_name='riemann', val=self%riemann_name, error=error)
      call fini%get(section_name=trim(codename)//'-scheme', option_name='riemann-options-Minf', val=self%Minf, error=error)
      self%Mco = 0.0_R8
      call fini%get(section_name=trim(codename)//'-scheme', option_name='riemann-options-Mco', val=self%Mco, error=error)
      if (error/=0) self%Mco = 0.0_R8
      if (error/=0) then
        if ( trim(self%riemann_name) == 'AUSM+-up' .or. trim(self%riemann_name) == 'AUSM+-up2' ) then
          write(*,*) '[ERROR] AUSM+-up flux selected. Minf must be defined in input.ini.'
          stop
        end if
      endif

      ! Integration variables
      self%intvar_name = 'cons'
      call fini%get(section_name=trim(codename)//'-scheme', option_name='integration-variables', val=self%intvar_name, error=error)

      ! Multigrid Parameters
      if (.not. self%HYDRA_MG) then
        call fini%get(section_name=trim(codename)//'-multigrid', option_name='levels', val=self%MGL, error=error)
        if (error/=0) self%MGL = 1
      endif
      
    end subroutine read_ini


end module MOSE_Types_m