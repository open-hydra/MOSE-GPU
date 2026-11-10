!> Multirate time integration: the state shared by the partition (Mod_Multirate_Partition), the driver (Mod_Multirate) and
!> the stage kernels. A cell of level l takes 2**l Heun sub-steps of h_l = DT / 2**l per macro step DT; the levels come from
!> each cell's own explicit limit (Compute_dt). The scheme is the one of branch multirate of DIMA-THRUST/MOSE (88de097,
!> doc/archive/MULTIRATE_PLAN.md sec. 1-3): the recursion of the levels, the coarse halo interpolated in time, the flux
!> telescoping at the level interfaces inside a block, levels equal across the block connections. This version keeps the
!> scheme and drops the launch machinery of that branch: the stage kernels run on boxes of the stage's level.
!>
!>   [MOSE-scheme] multirate      = f | t | census   (census: the levels are built and printed, the step stays explicit)
!>                 mr-ratio       = 16    macro step / smallest cell step, a power of 2
!>                 mr-safety      = 2     a cell's sub-step is at most its own limit / mr-safety
!>                 mr-wall-safety = 0     margin of the cells within mr-wall-rows of a viscous wall (0: mr-safety)
!>                 mr-wall-rows   = 10
!>                 mr-lmax        = -1    finest level (-1: from mr-ratio and the margins)
!>                 mr-repart      = 20    the levels are rebuilt at least every mr-repart macro steps,
!>                 mr-repart-hcap = 0.8   or when the cap of the frozen levels falls below this share of the requested step
!>                 mr-print       = 200   interval (macro steps) of the [MR] census line
!>                 mr-profile     = f     t: the wall of every call of the stages, per level (a device wait after each call)
!>                 mr-boxes       = 8     at most this many cell boxes per block and level (up to MR_NBOX_MAX)
!>                 mr-list-level  = -1    the levels >= this one run their flux phase on lists concatenated over the blocks
!>                                        (one launch per kernel and level; -1: every level on the boxes)
!>                 mr-list-records-level = 0  the levels >= this one (and >= mr-list-level) fill the ghosts, the boundary fluxes
!>                                        and the wall ghosts on the record lists (the coarser list-less levels keep their
!>                                        faces and cells on the boxes; -1: as mr-list-level)
!>                 mr-ghost-near  = t     a list level fills only the ghosts, and interpolates only the coarser halo cells, its
!>                                        stage reads (f: every ghost of its blocks, every coarser halo)
!>                 mr-land        = t     the macro step ends exactly on the output instants k sol-dtime and on time-threshold
!>                 mr-cons-halo   = t     the coarse halo interpolated in time in the conservative variables (f: primitives)
!>                 mr-depth       = ratio (auto: the depth from the octaves of the admissible steps, H = 2^L dt_min / mr-safety)
!>                 mr-margin      = f     t: a cell's level and the cap also see the admissible steps of its face neighbours in
!>                                        the block, over kappa = max(cfl (1 + theta of the limiter), 1.1)
!>                 mr-nonlocal    = off   count: the certificate of condition (4) is counted; redo: a violated macro step is
!>                                        redone with a deeper map or a shorter step (forces mr-margin)
!>                 mr-octaves     = f     t: the census of the admissible steps in octaves at the full partitions (prints only)
module MOSE_Lib_Multirate_Base
  use iso_fortran_env, only: R8 => real64
  implicit none
  private

  integer, parameter, public :: MR_LMAX_HARD = 8   !< the finest level the arrays hold
  integer, parameter, public :: MR_NBOX_MAX  = 64  !< boxes per (block, level) the arrays hold (mr-boxes: the ones used)
  integer, parameter, public :: MR_NCALL     = 8   !< the calls of a stage timed by mr-profile (Mod_Multirate MR_Stage)

  type, public :: mr_type
    ! keys
    logical  :: active = .false.          !< multirate = t or census: the partition runs in Compute_dt
    logical  :: drive  = .false.          !< multirate = t: the levels integrate
    real(R8) :: ratio  = 16.0_R8
    real(R8) :: sigma  = 2.0_R8
    real(R8) :: sigma_wall = 0.0_R8
    integer  :: wall_rows = 10
    integer  :: lmax = -1
    integer  :: repart = 20
    real(R8) :: repart_hcap = 0.8_R8
    integer  :: nprint = 200
    logical  :: profile = .false.
    integer  :: nboxmax = 8
    integer  :: listlev = -1
    integer  :: reclev = 0
    logical  :: ghost_near = .true.
    logical  :: land = .true.             !< mr-land (= open-hydra feature/multirate-integrator): outputs on k sol-dtime exactly
    logical  :: cons_halo = .true.        !< mr-cons-halo (= open-hydra feature/multirate-integrator): coarse halo in conservative variables
    logical  :: depth_auto = .false.      !< mr-depth = auto (= open-hydra feature/multirate-integrator): depth from the octaves
    logical  :: margin = .false.          !< mr-margin (= open-hydra interface-margin): the neighbours' step in the levels and the cap
    integer  :: nonlocal = 0              !< mr-nonlocal (= open-hydra nonlocal-cfl): 0 off, 1 count (the certificate), 2 redo
    logical  :: octaves = .false.         !< mr-octaves: the census of the admissible steps in octaves (prints only)
    integer  :: viol_step = 0             !< certificate: violations of the current macro step (this rank)
    integer(8) :: ncert_viol = 0, ncert_steps = 0   !< certificate: violations and macro steps with any, over the run
    logical  :: nonphys = .false.         !< redo: a state that is not physical at the end of the attempt (Lib_Newstate)
    integer(8) :: nredo_steps = 0, nredo_att = 0, nredo_deeper = 0   !< redo: macro steps redone, rejected attempts, deeper maps
    real(R8) :: kappa = 1.1_R8            !< its factor, max(cfl (1 + theta), 1.1), theta of the limiter (MR_Setup)
    real(R8) :: ratio_eff = 16.0_R8       !< H / dt_min requested: mr-ratio, or 2**L^ / mr-safety with mr-depth = auto
    ! partition
    integer  :: L = 0                     !< finest level of the ladder
    real(R8) :: dtmin = 0.0_R8, Hreq = 0.0_R8, Hcap = 0.0_R8, H = 0.0_R8, Hbuilt = 0.0_R8, Cp = 1.0_R8, Cp_raw = 1.0_R8
    real(R8) :: Hreqb = 1.0_R8           !< the requested step at the last full partition (Hbuilt is the step it gave)
    real(R8) :: hist(0:MR_LMAX_HARD) = 0.0_R8
    logical  :: have_part = .false.
    integer  :: nsince = 0, nlight = 0, nfull = 0, npass = 0
    integer  :: nrounds = 0, nseam = 0, nviol = 0, nwallc = 0
    integer  :: nifc = 0                  !< interface faces of this rank (telescoping slots)
    ! mr-reuse (Mod_Multirate_Partition mr_build_from_levels): the ladder the structures were built for, the cells whose level
    ! differs from that build at the last full partition (all ranks; -1: not compared), the full partitions that kept them
    integer  :: Lbuilt = -1, nchg = -1, nreuse = 0
    ! stage
    logical  :: in_stage = .false.        !< inside a stage of the ladder: the kernels take the boxes of level lcur
    integer  :: lcur = 0, irkcur = 1
    real(R8) :: tau_lev(0:MR_LMAX_HARD) = 0.0_R8   !< start (from the macro step start) of the current step of each level
    real(R8) :: h_lev(0:MR_LMAX_HARD) = 0.0_R8     !< its length
    ! census of the stages
    real(R8) :: tlev(0:MR_LMAX_HARD) = 0.0_R8
    integer  :: nlev(0:MR_LMAX_HARD) = 0
    real(R8) :: tcall(MR_NCALL, 0:MR_LMAX_HARD) = 0.0_R8   !< mr-profile: the wall of each call of the stages, per level
    real(R8) :: tmacro = 0.0_R8, tpart = 0.0_R8, tfull = 0.0_R8   !< the wall of the macro steps, of the partition (all, full)
    real(R8) :: tgh(8, 0:MR_LMAX_HARD) = 0.0_R8   !< mr-profile: the ghost fill of the stages by part (Fill_Ghost_Cell_local; 8 = stages with 410 faces)
    real(R8) :: tit(3) = 0.0_R8                    !< mr-profile: Compute_dt (partition included), the rest of the step, between steps
    real(R8) :: tpb(9) = 0.0_R8                    !< mr-profile: the full partitions by part (Mod_Multirate_Partition MR_Partition)
  contains
    procedure, pass(self) :: ini_setup => mr_read_ini
  end type mr_type

  type(mr_type), public :: obj_mr
  !> mr-land (open-hydra MOSE feature/multirate-integrator, Mod_Multirate.f90 619-638 and 785-790 there): Compute_dt cuts the
  !> macro step so that it ends exactly on the next output instant t_k = k sol-dtime (k integer, the grid anchored at t = 0)
  !> or on time-threshold, whichever comes first; a step that would end within 1e-12 H of it ends there. The time is then
  !> assigned (mr_t_land), not summed, and the time-based output is triggered by mr_output_landed, not by a comparison of
  !> times (k sol-dtime and (k-1) sol-dtime + sol-dtime can differ in the last place). Set by MR_Land (Mod_dt), reset at
  !> the start of every Explicit_Step.
  logical,  public :: mr_landed = .false., mr_output_landed = .false.
  real(R8), public :: mr_t_land = 0.0_R8

  !> mr-async (MOSE_MR_ASYNC, default on): the stage kernels of a list level are queued on the async queue mr_q in
  !> program order (one stream, the same order as the synchronous launches: the same result to the bit) and the host
  !> waits once, at the end of the stage, instead of before and after every launch (the launch census of 05/11: ~24
  !> launches and ~45 host round-trips per stage, 1.1 ms for 1 132 cells). Outside a list stage, and with the switch
  !> off, mr_q holds acc_async_sync (-1): the async clause then means a synchronous launch, the behaviour so far. The
  !> queue id is above every block index (the block path of the boundary fluxes queues on the block index) and is not
  !> the overlap queue of the inter-rank fill. One host thread only (a queue belongs to the thread that launched on it).
  integer, public :: mr_q = -1
  integer, public :: mr_qid = -1
  logical, public :: mr_async = .false.
  public :: mr_async_setup, mr_set_queue, mr_q_async, mr_set_queue_expl, mr_set_queue_sync, mr_flux_async, mr_q_flux
  public :: mr_fanout_on, mr_fanout_q



  !> The interface faces of one direction of a block: face index in that direction's fhat/dfhat, coarse level, slot in gacc.
  type, public :: mr_ifl_type
    integer :: n = 0
    integer, allocatable :: f(:), lc(:), s(:)
  end type mr_ifl_type

  !> Per local block: the level of every cell, the scratch of the rules, the margin class (1 = wall row), the cell limit
  !> of the last Compute_dt, the boxes of each level, the interface faces, and the boundary records of each level by
  !> launch of BC_Fluxes_blk (Mod_BC bc_launch_ranges): bcl(bco(g,l)+1 : bco(g,l)+bcn(g,l)) on the device.
  type, public :: mr_blk_type
    integer,  allocatable :: lev(:,:,:), lev2(:,:,:), cls(:,:,:)
    integer,  allocatable :: levb(:,:,:)   !< mr-reuse: host copy of the levels the structures were built from
    real(R8), allocatable :: dtcap(:,:,:)
    real(R8), allocatable :: dtc(:,:,:), dmin(:,:,:)                      !< mr-nonlocal: step of the stage, smallest seen
    real(R8), allocatable :: crho(:,:,:), crg(:,:,:), ct(:,:,:), cmil(:,:,:)   !< mr-nonlocal: its thermo scratch (with ghosts)
    integer :: nbox(0:MR_LMAX_HARD) = 0
    integer :: box(6, MR_NBOX_MAX, 0:MR_LMAX_HARD) = 0      !< (ilo, ihi, jlo, jhi, klo, khi)
    type(mr_ifl_type) :: ifl(3)
    integer, allocatable :: bcl(:)
    integer :: bcn(6, 0:MR_LMAX_HARD) = 0, bco(6, 0:MR_LMAX_HARD) = 0
  end type mr_blk_type
  type(mr_blk_type), allocatable, public :: mrb(:)

  !> The halo of level m: the level-m cells within two cells (any direction) of a finer cell of their block, as flat rows of
  !> P_flat, and the stash of their stage-1 predictor while the finer levels run.
  type, public :: mr_halo_type
    integer :: n = 0
    integer,  allocatable :: c(:)
    real(R8), allocatable :: ws(:,:)
  end type mr_halo_type
  type(mr_halo_type), allocatable, public :: mrh(:)   ! (0:L)
  !> ... and, per list level l and coarser level m < l (mr-ghost-near), the cells of the halo of m in the rows the stage
  !> of l reads (Mod_Multirate_Partition mr_build_halo_near), with their position in the halo list of m (its ws row).
  type, public :: mr_hn_type
    integer :: n = 0
    integer, allocatable :: c(:), pos(:)
  end type mr_hn_type
  type(mr_hn_type), allocatable, public :: mrhn(:,:)   ! (0:L, 0:L) = (m, l)

  !> List path (mr-list-level): the interior faces of direction d of the cells of one level, concatenated over the local
  !> blocks (each face once), with what the face kernels read from the block arrays gathered by list position: the
  !> flat index c of the face's left cell (the c0 of the convective, the c1 of the diffusive kernel), the flat strides
  !> along the face normal (sa) and the two other axes in the order of the diffusive stencil (sb, sc), the widths
  !> dl(d, -1..2) of the convective stencil, the face metric fm(1:4), the cell metrics M of the two cells and their wall
  !> distances yn; the face buffers fh, dfh, tf (the twins of fhat, dfhat, tface by list position); the interface faces of
  !> the level in the list (position, coarse level, slot), the records MR_Iface_Update works on.
  type, public :: mr_fl_type
    integer :: n = 0, ni = 0
    integer,  allocatable :: c(:), sa(:), sb(:), sc(:)
    real(R8), allocatable :: dl(:,:), fm(:,:), mm(:,:,:,:), yn(:,:)   ! (4,n), (4,n), (3,3,2,n), (2,n)
    real(R8), allocatable :: fh(:,:), dfh(:,:), tf(:,:)               ! (n,nprim), (n,nprim), (n,5[+nsc])
    integer,  allocatable :: if_(:), ilc(:), isl(:)
  end type mr_fl_type
  !> ... and the cells of one level, concatenated over the local blocks: flat index, and the positions of their six faces
  !> in the face lists in the order of pass 2 (i+, i-, j+, j-, k+, k-; 0 = a boundary face, no interior flux); for the
  !> cell kernels (Newstate, SST) the block-local indices ijk, the block's flat offset bo and row lengths pp (pi, pj),
  !> and the block arrays they read at the cell, gathered: volume, metric, wall distance.
  type, public :: mr_cl_type
    integer :: n = 0
    integer,  allocatable :: c(:), pf(:,:), ijk(:,:), bo(:), pp(:,:)
    real(R8), allocatable :: vol(:), met(:,:,:), wd(:)
  end type mr_cl_type
  type(mr_fl_type), allocatable, public :: mrfl(:,:)   ! (3, 0:MR_LMAX_HARD)
  !> ... and the ghost fill of a list level (Lib_Extrapolation Fill_Ghost_Cell_local): the physical boundary records of the
  !> blocks with cells of the level (record, block, the face normal the symmetry / mirror ghosts read), the connection
  !> records of the pairs into those blocks from another block (record, receiving block, donor block), and the pairs of a
  !> block with itself, which keep their own launch. mr-ghost-near: only the records that write a ghost the stage reads (in
  !> the boxes of the level grown by gc, or a donor ghost the transverse gather of a level record reads), and the gather of
  !> the connection records whose cell is of the level (pbr/pbs: record, donor block); otherwise pbr/pbs = cr/cs.
  type, public :: mr_gl_type
    integer :: np = 0, nc = 0, nself = 0, npb = 0, nq = -1
    integer,  allocatable :: pr(:), pb(:), cr(:), cm(:), cs(:), sp(:), pbr(:), pbs(:)
    integer,  allocatable :: qk(:)   ! mr-ghost-near: the 410 slots (Mod_BC_Q2D) whose ghosts the stage reads (nq = -1: all)
    real(R8), allocatable :: pn(:,:), pyn(:)          ! + the wall distance of the record's cell (Extrapolate_Wall)
  end type mr_gl_type
  type(mr_gl_type), allocatable, public :: mrgl(:)     ! (0:MR_LMAX_HARD)
  !> ... and the boundary fluxes of a list level (Mod_BC BC_Fluxes_list): per block face f = 1..6 the records whose interior
  !> cell is of the level, concatenated over the local blocks (record, block), with the values the per-block record body
  !> (bc_record_body.inc) takes from the block arrays gathered by position: normal and area (gn), the half widths and the
  !> four reconstruction widths (gdl: Bound_Dl, Int_Dl, dlc4), the cell metrics of the record and ghost cells (gm), their
  !> wall distances and beta (gy: ynrec, yngho, beta).
  type, public :: mr_bl_type
    integer :: n = 0
    integer,  allocatable :: rec(:), rb(:)
    real(R8), allocatable :: gn(:,:), gdl(:,:), gm(:,:,:,:), gy(:,:)   ! (4,n), (6,n), (3,3,2,n), (3,n)
  end type mr_bl_type
  type(mr_bl_type), allocatable, public :: mrbl(:,:)   ! (6, 0:MR_LMAX_HARD)
  !> ... and the thermo cache of a list level (Mod_Thermo_Cache Refresh_Thermo_Cache): the cells the box refresh sweeps,
  !> each once, as flat rows -- ca: the boxes of the level grown by gc, ghosts included (rho/Rg/T); cb: their interior part
  !> (rho/Rg/T and the viscosity).
  type, public :: mr_tl_type
    integer :: na = 0, nb = 0
    integer, allocatable :: ca(:), cb(:)
  end type mr_tl_type
  type(mr_tl_type), allocatable, public :: mrtl(:)   ! (0:MR_LMAX_HARD)
  integer, allocatable, public :: mr_dn(:,:)            ! (3, nb): the block dimensions on the device (list kernels)
  type(mr_cl_type), allocatable, public :: mrcl(:)     ! (0:MR_LMAX_HARD)

  real(R8), allocatable, public :: gacc(:,:), gsave(:,:)   ! (nifc, nprim): the telescoping accumulators

  public :: mr_stage_on, mr_block_active, mr_nbox, mr_get_box, mr_dt_sel, mr_check_now, mr_list_on, mr_ptoc, mr_rec_on, &
            mr_rec_level, mr_setup_dn, mr_gather_bc_list
  public :: MR_Iface_Update, MR_Iface_Update_List

contains


  subroutine mr_read_ini ( self, fini )
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    implicit none
    class(mr_type) :: self
    type(file_ini) :: fini
    character(len=32) :: v
    integer :: error, i
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='multirate', val=v, error=error)
    do i = 1, len_trim(v)
      if ( v(i:i) >= 'A' .and. v(i:i) <= 'Z' ) v(i:i) = achar( iachar(v(i:i)) + 32 )
    end do
    select case ( trim(adjustl(v)) )
      case ('f', 'false', '.false.', 'no'); self%active = .false.; self%drive = .false.
      case ('t', 'true', '.true.', 'yes');  self%active = .true.;  self%drive = .true.
      case ('census');                      self%active = .true.;  self%drive = .false.
      case default
        write(*,'(3A)') ' [MR] ERROR: [MOSE-scheme] multirate = ', trim(v), ': t, f or census'
        error stop 'invalid multirate value (see the [MR] ERROR line)'
    end select
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-ratio', val=self%ratio, error=error)
    if (error /= 0) self%ratio = 16.0_R8
    self%ratio_eff = self%ratio
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-safety', val=self%sigma, error=error)
    if (error /= 0) self%sigma = 2.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-wall-safety', val=self%sigma_wall, error=error)
    if (error /= 0) self%sigma_wall = 0.0_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-wall-rows', val=self%wall_rows, error=error)
    if (error /= 0) self%wall_rows = 10
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-lmax', val=self%lmax, error=error)
    if (error /= 0) self%lmax = -1
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-repart', val=self%repart, error=error)
    if (error /= 0) self%repart = 20
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-repart-hcap', val=self%repart_hcap, error=error)
    if (error /= 0) self%repart_hcap = 0.8_R8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-print', val=self%nprint, error=error)
    if (error /= 0) self%nprint = 200
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-profile', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('t', 'T', 'true', '.true.', 'yes'); self%profile = .true.
      case default;                             self%profile = .false.
    end select
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-boxes', val=self%nboxmax, error=error)
    if (error /= 0) self%nboxmax = 8
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-list-level', val=self%listlev, error=error)
    if (error /= 0) self%listlev = -1
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-list-records-level', val=self%reclev, error=error)
    if (error /= 0) self%reclev = 0
    v = 't'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-ghost-near', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('f', 'F', 'false', '.false.', 'no'); self%ghost_near = .false.
      case default;                              self%ghost_near = .true.
    end select
    v = 't'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-land', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('f', 'F', 'false', '.false.', 'no'); self%land = .false.
      case default;                              self%land = .true.
    end select
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-margin', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('t', 'T', 'true', '.true.', 'yes'); self%margin = .true.
      case default;                             self%margin = .false.
    end select
    v = 'off'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-nonlocal', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('count');  self%nonlocal = 1
      case ('redo');   self%nonlocal = 2; self%margin = .true.   ! as the branch: redo forces the interface margin
      case default;    self%nonlocal = 0
    end select
    v = 'ratio'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-depth', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('auto', 'octaves'); self%depth_auto = .true.
      case default;             self%depth_auto = .false.
    end select
    v = 'f'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-octaves', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('t', 'T', 'true', '.true.', 'yes'); self%octaves = .true.
      case default;                             self%octaves = .false.
    end select
    v = 't'
    call fini%get(section_name=trim(codename)//'-scheme', option_name='mr-cons-halo', val=v, error=error)
    select case ( trim(adjustl(v)) )
      case ('f', 'F', 'false', '.false.', 'no'); self%cons_halo = .false.
      case default;                              self%cons_halo = .true.
    end select
    if ( .not. self%active ) return
    if ( self%ratio < 1.0_R8 .or. abs(log(self%ratio)/log(2.0_R8) - nint(log(self%ratio)/log(2.0_R8))) > 1.0e-9_R8 ) &
      error stop '[MOSE-scheme] mr-ratio must be a power of 2 (1, 2, 4, ...)'
    if ( self%sigma < 1.0_R8 ) error stop '[MOSE-scheme] mr-safety must be >= 1'
    if ( self%sigma_wall /= 0.0_R8 .and. self%sigma_wall < self%sigma ) &
      error stop '[MOSE-scheme] mr-wall-safety must be 0 (off) or >= mr-safety'
    if ( self%wall_rows < 1 ) self%wall_rows = 1
    if ( self%repart < 1 ) self%repart = 1
    if ( self%lmax > MR_LMAX_HARD ) error stop '[MOSE-scheme] mr-lmax: at most 8'
    if ( self%nboxmax < 1 .or. self%nboxmax > MR_NBOX_MAX ) error stop '[MOSE-scheme] mr-boxes: 1 to 64'
  end subroutine mr_read_ini


  !> mr-profile: t += the wall since the clock ct_ (after the device work queued so far), ct_ := now.
  subroutine mr_ptoc ( t, ct_ )
    real(R8),   intent(inout) :: t
    integer(8), intent(inout) :: ct_
    integer(8) :: c1_, cr_
    !$acc wait
    call system_clock ( c1_, cr_ )
    t = t + real(c1_ - ct_, R8) / real(cr_, R8)
    ct_ = c1_
  end subroutine mr_ptoc

  !> True inside a stage of the ladder: the stage kernels then work on the boxes of level obj_mr%lcur.
  logical function mr_stage_on ()
    mr_stage_on = obj_mr%drive .and. obj_mr%in_stage
  end function mr_stage_on

  !> Block b holds cells of the stage's level (always true outside a stage).
  logical function mr_block_active ( b )
    integer, intent(in) :: b
    mr_block_active = .true.
    if ( mr_stage_on() ) mr_block_active = mrb(b)%nbox(obj_mr%lcur) > 0
  end function mr_block_active

  !> The number of cell boxes of block b at this stage (1, the whole block, outside a stage).
  integer function mr_nbox ( b )
    integer, intent(in) :: b
    mr_nbox = 1
    if ( mr_stage_on() ) mr_nbox = mrb(b)%nbox(obj_mr%lcur)
  end function mr_nbox

  !> Box ib of block b at this stage: cells lo(d) .. hi(d); the whole block outside a stage.
  subroutine mr_get_box ( b, ib, n, lo, hi )
    integer, intent(in)  :: b, ib, n(3)
    integer, intent(out) :: lo(3), hi(3)
    if ( mr_stage_on() ) then
      lo(1) = mrb(b)%box(1,ib,obj_mr%lcur); hi(1) = mrb(b)%box(2,ib,obj_mr%lcur)
      lo(2) = mrb(b)%box(3,ib,obj_mr%lcur); hi(2) = mrb(b)%box(4,ib,obj_mr%lcur)
      lo(3) = mrb(b)%box(5,ib,obj_mr%lcur); hi(3) = mrb(b)%box(6,ib,obj_mr%lcur)
    else
      lo = 1; hi = n
    end if
  end subroutine mr_get_box

  !> The level predicate of the cell kernels: inside a stage of level l the sub-step of that level, DT / 2**l (exactly the
  !> dtlocal MR_Set_Local_dt gives its cells: the same expression), else -1 (every cell).
  real(R8) function mr_dt_sel ()
    mr_dt_sel = -1.0_R8
    if ( mr_stage_on() ) mr_dt_sel = obj_mr%h_lev(obj_mr%lcur)
  end function mr_dt_sel

  !> The flux phase of this stage runs on the lists of its level (mr-list-level, Mod_Fluxes Fluxes_List).
  logical function mr_list_on ()
    mr_list_on = .false.
    if ( mr_stage_on() .and. obj_mr%listlev >= 0 ) mr_list_on = ( obj_mr%lcur >= obj_mr%listlev )
  end function mr_list_on

  !> The first level whose ghost fill, boundary fluxes and wall ghosts run on the record lists (mr-list-records-level, at
  !> most mr-list-level; -1 without lists).
  integer function mr_rec_level ()
    mr_rec_level = -1
    if ( obj_mr%listlev < 0 ) return
    mr_rec_level = obj_mr%listlev
    if ( obj_mr%reclev >= 0 ) mr_rec_level = min( obj_mr%reclev, obj_mr%listlev )
  end function mr_rec_level

  !> The ghost fill, the boundary fluxes and the wall ghosts of this stage run on the record lists of its level.
  logical function mr_rec_on ()
    mr_rec_on = .false.
    if ( mr_stage_on() .and. obj_mr%listlev >= 0 ) mr_rec_on = ( obj_mr%lcur >= mr_rec_level() )
  end function mr_rec_on

  !> The state check follows this update: always outside a stage; inside, only after stage 2 of level 0, the last stage of
  !> the macro step (every level at its end).
  logical function mr_check_now ()
    mr_check_now = .true.
    if ( mr_stage_on() ) mr_check_now = ( obj_mr%lcur == 0 .and. obj_mr%irkcur == 2 )
  end function mr_check_now


  !> The flux telescoping of block b (plan sec. 3.2), between the face kernels and pass 2 of a stage of level lcur, on the
  !> interface faces of the block (|level difference| = 1). With g = fhat - dfhat, the net flux pass 2 gives the two cells:
  !>   coarse stage 1 : gsave = g                    (the coarse cell uses its own flux at the start of its step)
  !>   fine stages    : gacc += g / 4                (4 = 2 sub-steps x 2 stages of the fine side per coarse step)
  !>   coarse stage 2 : fhat = 2 gacc - gsave, dfhat = 0, then gacc = gsave = 0
  !> so over its step the coarse cell receives exactly the impulse the fine cells took through the face.
  subroutine MR_Iface_Update ( b, fh1, fh2, fh3, dfh1, dfh2, dfh3, nfmax, visc )
    use MOSE_Global_m, only: nprim
    implicit none
    integer,  intent(in)    :: b, nfmax
    logical,  intent(in)    :: visc
    real(R8), intent(inout) :: fh1(nfmax, nprim), fh2(nfmax, nprim), fh3(nfmax, nprim)
    real(R8), intent(inout) :: dfh1(nfmax, nprim), dfh2(nfmax, nprim), dfh3(nfmax, nprim)
    if ( .not. mr_stage_on() ) return
    if ( obj_mr%L == 0 ) return
    if ( mrb(b)%ifl(1)%n > 0 ) call mr_iface_dir ( fh1, dfh1, nfmax, mrb(b)%ifl(1)%f, mrb(b)%ifl(1)%lc, mrb(b)%ifl(1)%s, &
                                                   mrb(b)%ifl(1)%n, obj_mr%lcur, obj_mr%irkcur, visc )
    if ( mrb(b)%ifl(2)%n > 0 ) call mr_iface_dir ( fh2, dfh2, nfmax, mrb(b)%ifl(2)%f, mrb(b)%ifl(2)%lc, mrb(b)%ifl(2)%s, &
                                                   mrb(b)%ifl(2)%n, obj_mr%lcur, obj_mr%irkcur, visc )
    if ( mrb(b)%ifl(3)%n > 0 ) call mr_iface_dir ( fh3, dfh3, nfmax, mrb(b)%ifl(3)%f, mrb(b)%ifl(3)%lc, mrb(b)%ifl(3)%s, &
                                                   mrb(b)%ifl(3)%n, obj_mr%lcur, obj_mr%irkcur, visc )
  end subroutine MR_Iface_Update

  !> The flux telescoping of a stage of level lcur on its face lists (list path): the same records as MR_Iface_Update,
  !> addressed by their positions in the level's lists.
  subroutine MR_Iface_Update_List ( fl, visc )
    implicit none
    type(mr_fl_type), intent(inout) :: fl(3)
    logical, intent(in) :: visc
    integer :: d
    if ( .not. mr_stage_on() ) return
    if ( obj_mr%L == 0 ) return
    do d = 1, 3
      if ( fl(d)%ni > 0 ) call mr_iface_dir ( fl(d)%fh, fl(d)%dfh, fl(d)%n, fl(d)%if_, fl(d)%ilc, fl(d)%isl, fl(d)%ni, &
                                              obj_mr%lcur, obj_mr%irkcur, visc )
    end do
  end subroutine MR_Iface_Update_List

  subroutine mr_iface_dir ( fh, dfh, nfmax, f, lc, sl, n, lmin, irk, visc )
    use MOSE_Global_m, only: nprim
    implicit none
    integer,  intent(in)    :: nfmax, n, lmin, irk, f(n), lc(n), sl(n)
    logical,  intent(in)    :: visc
    real(R8), intent(inout) :: fh(nfmax, nprim), dfh(nfmax, nprim)
    integer  :: e, v, ff, ss
    real(R8) :: g
    !$acc parallel loop collapse(2) gang vector present(fh, dfh, f, lc, sl, gacc, gsave) private(g, ff, ss) async(mr_q)
    do v = 1, nprim
      do e = 1, n
        ff = f(e); ss = sl(e)
        if ( visc ) then
          g = fh(ff, v) - dfh(ff, v)
        else
          g = fh(ff, v)
        end if
        if ( lmin == lc(e) ) then
          if ( irk == 1 ) then
            gsave(ss, v) = g
          else
            fh(ff, v) = 2.0_R8*gacc(ss, v) - gsave(ss, v)
            if ( visc ) dfh(ff, v) = 0.0_R8
            gacc(ss, v) = 0.0_R8; gsave(ss, v) = 0.0_R8
          end if
        else if ( lmin == lc(e) + 1 ) then
          gacc(ss, v) = gacc(ss, v) + 0.25_R8*g
        end if
      end do
    end do
  end subroutine mr_iface_dir

  !> The block dimensions on the device (mr_dn), once: the list kernels of the stages and the boundary fluxes of the
  !> explicit step on its record lists (Mod_BC) read them.
  subroutine mr_setup_dn ( domain )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer :: b
    if ( allocated(mr_dn) ) return
    allocate ( mr_dn(3, domain%nb) )
    mr_dn = 1
    do b = 1, domain % nb
      if (is_local_block(b)) mr_dn(:,b) = domain%blk(b)%dim
    end do
    !$acc enter data copyin(mr_dn)
  end subroutine mr_setup_dn

  !> The boundary-flux lists bl(1:6) (mr_bl_type) of level l_ -- or, with l_ < 0, of every record (the explicit step,
  !> Mod_BC bc_list_setup): per face f the records of every local block whose interior cell is of the level, in block
  !> order, with the values of bc_record_body.inc gathered from the host block arrays by the same formulas (0.5 dl is
  !> exact, so Bound_Dl and Int_Dl are the device's doubles). A 3D rotational periodic record (type 2, nk > 1) is not on
  !> the list path: with l_ < 0 nothing is built and ok = .false. (the caller keeps the per-block launches), with l_ >= 0
  !> it stops the run.
  subroutine mr_gather_bc_list ( domain, l_, bl, ok )
    use MOSE_Types_m, only: MOSE_domain_type
    use MOSE_Global_m, only: gc
    use MOSE_Mod_MPI,  only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer, intent(in) :: l_
    type(mr_bl_type), intent(inout) :: bl(:)
    logical, intent(out) :: ok
    integer :: b, f, idx, e, lo, nrec(6), Im, Jm, Km, Fm, dir, modfm, gstep, Fi, Fj, Fk, Ii, Ij, Ik
    integer :: Ig_, Jg_, Kg_, o1_, o2_, o3_, o4_, t
    integer :: ef(6), fc, sl
    ok = .true.
    if ( l_ < 0 ) then
      do b = 1, domain % nb
        if ( .not. is_local_block(b) ) cycle
        if ( domain%blk(b)%dim(3) == 1 ) cycle
        do idx = domain % bc_blk_lo(b), domain % bc_blk_hi(b)
          if ( domain%bc_typ(idx) == 2 ) then
            ok = .false.; return
          end if
        end do
      end do
    end if
    ! The six lists are the launches, not the faces: the records of the two opposite faces of a direction d go to
    ! list 2d-1 when the block is not degenerate in d (n_d >= 2: the two faces write disjoint planes, so one launch; the
    ! per-block kernel fuses them the same way, bc_launch_ranges), and the even face of a degenerate block (n_d = 1: both
    ! faces write the same cells) goes to list 2d, launched after. A corner cell on faces 1 and 3 still accumulates
    ! 1 then 3 (list 1 before list 3), on faces 2 and 3 likewise (list 1 holds face 2), so the order of the six
    ! per-face launches is kept to the last bit; three launches per stage instead of six.
    nrec = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( domain % bc_blk_hi(b) < domain % bc_blk_lo(b) ) cycle
      lo = domain % bc_blk_lo(b)
      do f = 1, 6
        sl = slot_of(b, f)
        do idx = lo, lo + domain % n_bf(b,f) - 1
          if ( of_level(b, idx) ) nrec(sl) = nrec(sl) + 1
        end do
        lo = lo + domain % n_bf(b,f)
      end do
    end do
    do f = 1, 6
      bl(f)%n = nrec(f)
      e = max(nrec(f), 1)
      allocate ( bl(f)%rec(e), bl(f)%rb(e), bl(f)%gn(4,e), bl(f)%gdl(6,e), bl(f)%gm(3,3,2,e), &
                 bl(f)%gy(3,e) )
      bl(f)%rec = 1; bl(f)%rb = 1; bl(f)%gn = 0.0_R8; bl(f)%gdl = 0.0_R8; bl(f)%gm = 0.0_R8
      bl(f)%gy = 0.0_R8
    end do
    ef = 0
    do b = 1, domain % nb
      if (.not. is_local_block(b)) cycle
      if ( domain % bc_blk_hi(b) < domain % bc_blk_lo(b) ) cycle
      lo = domain % bc_blk_lo(b)
      associate ( bk => domain%blk(b) )
      do fc = 1, 6
        sl = slot_of(b, fc)
        do idx = lo, lo + domain % n_bf(b,fc) - 1
          Im = domain%bc_bi(idx); Jm = domain%bc_bj(idx); Km = domain%bc_bk(idx); Fm = domain%bc_bf(idx)
          if ( .not. of_level(b, idx) ) cycle
          f = sl
          t = domain%bc_typ(idx)
          ! rotational periodicity: on an nk = 1 block (2D) the record reads only its own cell; the 3D seam reads dl of the
          ! cells k-1, k, 1, 2 of the block, not gathered
          if ( t == 2 .and. bk%dim(3) /= 1 ) &
            error stop 'mr-list-level: 3D rotational periodic records (type 2) are not on the list path'
          ef(f) = ef(f) + 1; e = ef(f)
          bl(f)%rec(e) = idx; bl(f)%rb(e) = b
          dir = (Fm + 1)/2; modfm = mod(Fm, 2); gstep = 2*modfm - 1
          Fi = Im; Fj = Jm; Fk = Km; Ii = Im; Ij = Jm; Ik = Km
          if ( dir == 1 ) then
            Fi = Im - modfm; Ii = Im + gstep; bl(f)%gn(1:4,e) = bk%fm1(1:4,Fi,Fj,Fk)
          else if ( dir == 2 ) then
            Fj = Jm - modfm; Ij = Jm + gstep; bl(f)%gn(1:4,e) = bk%fm2(1:4,Fi,Fj,Fk)
          else
            Fk = Km - modfm; Ik = Km + gstep; bl(f)%gn(1:4,e) = bk%fm3(1:4,Fi,Fj,Fk)
          end if
          bl(f)%gdl(1,e) = 0.5d0 * bk%dl_f(dir, Im, Jm, Km)
          bl(f)%gdl(2,e) = 0.5d0 * bk%dl_f(dir, Ii, Ij, Ik)
          Ig_ = Im; Jg_ = Jm; Kg_ = Km
          if ( dir == 1 ) Ig_ = Im - gstep
          if ( dir == 2 ) Jg_ = Jm - gstep
          if ( dir == 3 ) Kg_ = Km - gstep
          o1_ = modfm*(-gstep*gc) + (1-modfm)*( gstep)
          o2_ = modfm*(-gstep)
          o3_ = (1-modfm)*(-gstep)
          o4_ = (1-modfm)*(-gstep*gc) + modfm*( gstep)
          if ( dir == 1 ) then
            bl(f)%gdl(3:6,e) = (/ bk%dl_f(1,Im+o1_,Jm,Km), bk%dl_f(1,Im+o2_,Jm,Km), bk%dl_f(1,Im+o3_,Jm,Km), &
                                       bk%dl_f(1,Im+o4_,Jm,Km) /)
          else if ( dir == 2 ) then
            bl(f)%gdl(3:6,e) = (/ bk%dl_f(2,Im,Jm+o1_,Km), bk%dl_f(2,Im,Jm+o2_,Km), bk%dl_f(2,Im,Jm+o3_,Km), &
                                       bk%dl_f(2,Im,Jm+o4_,Km) /)
          else
            bl(f)%gdl(3:6,e) = (/ bk%dl_f(3,Im,Jm,Km+o1_), bk%dl_f(3,Im,Jm,Km+o2_), bk%dl_f(3,Im,Jm,Km+o3_), &
                                       bk%dl_f(3,Im,Jm,Km+o4_) /)
          end if
          bl(f)%gy(3,e) = 1.0d0                        ! beta_ (1: no shock detector on the list path)
          if ( t == 1 .or. t == 666 ) then
            bl(f)%gy(3,e) = bk%beta(Im,Jm,Km)
            bl(f)%gm(:,:,1,e) = bk%m_f(:,:,Im,Jm,Km); bl(f)%gy(1,e) = bk%yn(Im,Jm,Km)
            bl(f)%gm(:,:,2,e) = bk%m_f(:,:,Ig_,Jg_,Kg_); bl(f)%gy(2,e) = bk%yn(Ig_,Jg_,Kg_)
          else if ( t == 667 ) then
            bl(f)%gy(3,e) = bk%beta(Im,Jm,Km)
            bl(f)%gm(:,:,1,e) = bk%m_f(:,:,Im,Jm,Km); bl(f)%gy(1,e) = bk%yn(Im,Jm,Km)
            bl(f)%gm(:,:,2,e) = bk%m_f(:,:,Ig_,Jg_,Kg_); bl(f)%gy(2,e) = bl(f)%gy(1,e)
          else if ( t == 3 .or. t == 5 .or. t == 6 ) then
            bl(f)%gm(:,:,1,e) = bk%m_f(:,:,Im,Jm,Km); bl(f)%gy(1,e) = bk%yn(Im,Jm,Km)
          end if
        end do
        lo = lo + domain % n_bf(b,fc)
      end do
      end associate
    end do
    if ( any(ef /= nrec) ) error stop 'mr_gather_bc_list: count and fill differ'
    ! bc-sort (MOSE_BC_SORT, default on): the records of each launch grouped by boundary type (stable: the original order
    ! within a type), so that the 32 threads of a warp follow one branch of the record body. Within one launch every cell
    ! receives at most one record (opposite faces are fused only when the block has two or more cells across), so the
    ! residual of every cell gets the same additions in the same order: the same result to the bit.
    if ( bc_sort_on() ) then
      do f = 1, 6
        if ( bl(f)%n > 1 ) call sort_by_type ( bl(f) )
      end do
    end if
    call print_types ( )
    do f = 1, 6
      !$acc enter data copyin(bl(f)%rec, bl(f)%rb, bl(f)%gn, bl(f)%gdl, bl(f)%gm, bl(f)%gy)
    end do
  contains
    logical function bc_sort_on ()
      logical, save :: checked = .false., on = .true.
      character(len=8) :: env
      integer :: st, ln
      if ( .not. checked ) then
        checked = .true.
        call get_environment_variable ( 'MOSE_BC_SORT', env, ln, st )
        if ( st == 0 .and. ln > 0 ) then
          if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
        end if
      end if
      bc_sort_on = on
    end function bc_sort_on

    subroutine sort_by_type ( lst )
      type(mr_bl_type), intent(inout) :: lst
      integer :: n, e, k, t, tmin, tnext, nt
      integer,  allocatable :: perm(:), irec(:), irb(:)
      real(R8), allocatable :: tgn(:,:), tgdl(:,:), tgm(:,:,:,:), tgy(:,:)
      n = lst%n
      allocate ( perm(n) )
      k = 0; nt = 0
      tmin = -huge(0)
      do   ! the types in ascending order, each group in the original order
        tnext = huge(0)
        do e = 1, n
          t = domain%bc_typ(lst%rec(e))
          if ( t > tmin .and. t < tnext ) tnext = t
        end do
        if ( tnext == huge(0) ) exit
        nt = nt + 1
        do e = 1, n
          if ( domain%bc_typ(lst%rec(e)) == tnext ) then
            k = k + 1; perm(k) = e
          end if
        end do
        tmin = tnext
      end do
      if ( k /= n ) error stop 'mr_gather_bc_list: the type sort lost records'
      if ( nt > 1 ) then
        allocate ( irec(n), irb(n), tgn(4,n), tgdl(6,n), tgm(3,3,2,n), tgy(3,n) )
        do e = 1, n
          irec(e) = lst%rec(perm(e)); irb(e) = lst%rb(perm(e))
          tgn(:,e) = lst%gn(:,perm(e)); tgdl(:,e) = lst%gdl(:,perm(e))
          tgm(:,:,:,e) = lst%gm(:,:,:,perm(e)); tgy(:,e) = lst%gy(:,perm(e))
        end do
        lst%rec(1:n) = irec; lst%rb(1:n) = irb; lst%gn(:,1:n) = tgn; lst%gdl(:,1:n) = tgdl
        lst%gm(:,:,:,1:n) = tgm; lst%gy(:,1:n) = tgy
        deallocate ( irec, irb, tgn, tgdl, tgm, tgy )
      end if
      deallocate ( perm )
    end subroutine sort_by_type

    !> The records per launch and per boundary type, once: the explicit lists, and each level at the first full partition.
    subroutine print_types ( )
      use MOSE_Mod_MPI, only: mpi_is_root
      logical, save :: done_expl = .false.
      integer :: ff, e, tt, ntyp, typ(16), cnt(16), j
      if ( .not. mpi_is_root ) return
      if ( l_ < 0 ) then
        if ( done_expl ) return
        done_expl = .true.
      else
        if ( obj_mr%nfull /= 1 ) return
      end if
      do ff = 1, 6
        if ( bl(ff)%n == 0 ) cycle
        ntyp = 0; typ = 0; cnt = 0
        do e = 1, bl(ff)%n
          tt = domain%bc_typ(bl(ff)%rec(e))
          j = findloc ( typ(1:max(ntyp,1)), tt, dim=1 )
          if ( ntyp == 0 ) j = 0
          if ( j == 0 ) then
            if ( ntyp < 16 ) then
              ntyp = ntyp + 1; typ(ntyp) = tt; j = ntyp
            else
              cycle
            end if
          end if
          cnt(j) = cnt(j) + 1
        end do
        write(*,'(A,I0,A,I0,A,I0,A)',advance='no') ' [BC]   records of level ', l_, ' launch ', ff, ': ', bl(ff)%n, ' by type'
        do j = 1, ntyp
          write(*,'(A,I0,A,I0)',advance='no') '  ', typ(j), ':', cnt(j)
        end do
        write(*,'(A)') ''
      end do
    end subroutine print_types

    !> The list (launch) of face fc of block b: 2d-1 for the odd face and for the even face of a block with n_d >= 2,
    !> 2d for the even face of a block degenerate in d
    integer function slot_of ( b, fc )
      integer, intent(in) :: b, fc
      integer :: d
      d = (fc + 1) / 2
      slot_of = 2*d - 1
      if ( mod(fc, 2) == 0 .and. domain%blk(b)%dim(d) < 2 ) slot_of = 2*d
    end function slot_of
    !> The record idx of block b is on the list: of level l_, or any with l_ < 0 (mrb is not touched then)
    logical function of_level ( b, idx )
      integer, intent(in) :: b, idx
      of_level = .true.
      if ( l_ >= 0 ) of_level = ( mrb(b)%lev(domain%bc_bi(idx), domain%bc_bj(idx), domain%bc_bk(idx)) == l_ )
    end function of_level
  end subroutine mr_gather_bc_list

  !> mr-async: the switch once (MOSE_MR_ASYNC, default on; off with more than one OpenMP thread, where a queue of the
  !> launching thread is undefined for another), the queue id above the block count.
  subroutine mr_async_setup ( nb )
    use MOSE_Mod_MPI, only: mpi_is_root
#if defined (_OPENMP)
    use omp_lib, only: omp_get_max_threads
#endif
    integer, intent(in) :: nb
    character(len=8) :: env
    integer :: st, ln
    logical, save :: done = .false.
    if ( done ) return
    done = .true.
    mr_async = .true.
    call get_environment_variable ( 'MOSE_MR_ASYNC', env, ln, st )
    if ( st == 0 .and. ln > 0 ) then
      if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) &
        mr_async = .false.
    end if
#if defined (_OPENMP)
    if ( mr_async .and. omp_get_max_threads() > 1 ) then
      mr_async = .false.
      if ( mpi_is_root ) write(*,'(A)') ' [MR]   mr-async off: more than one OpenMP thread'
    end if
#endif
    mr_qid = nb + 64
    if ( mpi_is_root ) then
      if ( mr_async ) then
        write(*,'(A,I0,A)') ' [MR]   async stage chain on the list levels (queue ', mr_qid, '; MOSE_MR_ASYNC=0 turns it off)'
      else
        write(*,'(A)') ' [MR]   async stage chain off: every stage kernel synchronous'
      end if
    end if
  end subroutine mr_async_setup

  !> The queue of the stage of level l: the async queue on a list level with mr-async, else the synchronous one.
  subroutine mr_set_queue ( l )
    use openacc, only: acc_async_sync
    integer, intent(in) :: l
    logical, save :: checked = .false., blk_on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then   ! mr-async-blk (MOSE_MR_ASYNC_BLK, default on): the block-path levels on the queue too
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_ASYNC_BLK', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) blk_on = .false.
      end if
    end if
    mr_q = acc_async_sync
    if ( mr_async .and. obj_mr%listlev >= 0 ) then
      if ( l >= obj_mr%listlev ) mr_q = mr_qid
    end if
    ! a level on the block path (the coarse level): its ghost fills, thermo cache, boundary lists, wall ghosts, turbulence
    ! and update on the queue as well; its fluxes (the synchronous stream) bracketed by MR_Stage
    if ( mr_async .and. blk_on ) mr_q = mr_qid
  end subroutine mr_set_queue

  !> expl-async (MOSE_EXPL_ASYNC, default on): the queue of the explicit stage's per-block kernels is the async queue.
  subroutine mr_set_queue_expl ()
    use openacc, only: acc_async_sync
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_EXPL_ASYNC', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_q = acc_async_sync
    if ( mr_async .and. on ) mr_q = mr_qid
  end subroutine mr_set_queue_expl

  !> Back to the synchronous queue (outside the explicit stages).
  subroutine mr_set_queue_sync ()
    use openacc, only: acc_async_sync
    mr_q = acc_async_sync
  end subroutine mr_set_queue_sync

  !> The current stage queues its kernels asynchronously.
  logical function mr_q_async ()
    mr_q_async = mr_async .and. ( mr_q == mr_qid )
  end function mr_q_async

  !> mr-async-flux (MOSE_MR_ASYNC_FLUX, default on): the per-block flux kernels of a block-path stage (Mod_Fluxes
  !> Fluxes_blk: beta, the convective and diffusive pass 1 with the face thermo, the telescoping, pass 2) on the stage
  !> queue as well, in program order on its one stream (the same result to the bit), instead of synchronous launches after
  !> a join of the queue (nsys 315348: ~200 000 host waits in 30 macro steps, 102 637 of them in the face thermo fill).
  logical function mr_flux_async ()
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_ASYNC_FLUX', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_flux_async = on
  end function mr_flux_async

  !> The queue of the per-block flux kernels: the stage queue inside a multirate stage (the synchronous one there with
  !> MOSE_MR_ASYNC_FLUX=0, MR_Stage sets it), the synchronous one outside (the explicit step keeps its per-block fluxes).
  integer function mr_q_flux ()
    use openacc, only: acc_async_sync
    mr_q_flux = acc_async_sync
    if ( mr_stage_on() ) mr_q_flux = mr_q
  end function mr_q_flux

  !> mr-flux-fanout (MOSE_MR_FLUX_FANOUT, default on): on a list level with the async stage chain, the independent work of
  !> a stage on side queues forked from the stage queue and joined back to it -- the six pass-1 flux chains (convective of
  !> each direction; face thermo + diffusive of each direction: each writes only its own direction's arrays) and the
  !> coarser levels' halo interpolations (each writes only the cells of its own level). Each kernel computes what it did on
  !> the one queue, from the same inputs: the same result to the bit. On a small level the kernels fill a few
  !> multiprocessors each, so side by side they overlap instead of queueing.
  logical function mr_fanout_on ()
    logical, save :: checked = .false., on = .true.
    character(len=8) :: env
    integer :: st, ln
    if ( .not. checked ) then
      checked = .true.
      call get_environment_variable ( 'MOSE_MR_FLUX_FANOUT', env, ln, st )
      if ( st == 0 .and. ln > 0 ) then
        if ( env(1:1) == '0' .or. env(1:1) == 'n' .or. env(1:1) == 'N' .or. env(1:1) == 'f' .or. env(1:1) == 'F' ) on = .false.
      end if
    end if
    mr_fanout_on = on .and. mr_q_async()
  end function mr_fanout_on

  !> The k-th side queue of mr-flux-fanout (k = 1..6), above the stage queue mr_qid and the list-build queue mr_qid + 1.
  integer function mr_fanout_q ( k )
    integer, intent(in) :: k
    mr_fanout_q = mr_qid + 1 + k
  end function mr_fanout_q

end module MOSE_Lib_Multirate_Base
