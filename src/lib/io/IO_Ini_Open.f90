!> @brief The MOSE_open input dialect: a MOSE_open input.ini runs on the GPU code as it is.
!>
!> One entry point, MOSE_Load_Ini, replaces every `fini%load(filename='input.ini')` of the code. It loads the file
!> and, when the file is written the MOSE_open way ([MOSE-Parameters], [MOSE-Numerics] or [MOSE-Physics] equations
!> present), adds IN MEMORY the keys the GPU code reads ([MOSE-parameters], [MOSE-scheme], [MOSE-RANS], [MOSE-ODE],
!> [MOSE-nochem], [MOSE-probes], [MOSE-rotating-frame]) with the MOSE_open defaults (src/lib/config/Register_*.f90
!> of MOSE_open) and MOSE_open's own conventions where the GPU default differs (schmidt 0.0 = multicomponent
!> diffusion, point-implicit turbulence sources, no species-diffusion time-step limit). It stops with an
!> explicit message on every choice the GPU build does not have (each case below: limiters other than minmod,
!> vanleer, vanalbada, mc, superbee; the turbulence names that are not MOSE_open's (SSGLRR, QCR2000); multigrid, equilibrium
!> chemistry, soot,
!> stationary faces of the rotating frame). Nothing is written back to input.ini; the
!> translated keys are saved in OUTPUT/input.translated.ini. A GPU-dialect input is returned untouched, so those
!> runs stay bit-identical.
!>
!> The translation is idempotent (FiNeR `add` updates an existing key), so the periodic re-read of input.ini
!> (ini-diter) and the other loaders see the same keys; the report is printed once, by the root rank.
module MOSE_IO_Ini_Open
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Finer, only: file_ini
  use MOSE_Parameters_m, only: codename
  implicit none
  private
  public :: MOSE_Load_Ini

  !> .true. when the loaded input.ini is written in the MOSE_open dialect (set by MOSE_Load_Ini).
  logical, public, protected :: open_dialect = .false.
  !> [MOSE-Physics] chemistry = frozen (the MOSE_open default) while INPUT/chemistry-info.txt exists: IO_Chemistry
  !> must not read the mechanism (the GPU code otherwise goes reactive on the mere presence of the file).
  logical, public, protected :: open_force_frozen = .false.

  logical, save :: first = .true.
  integer, parameter :: vlen = 256
  integer, parameter :: nlog = 80
  integer :: nput = 0, nwarn = 0
  character(len=160) :: putlog(nlog), warnlog(nlog)

contains

  !> Load input.ini and translate the MOSE_open dialect into the keys the GPU code reads.
  subroutine MOSE_Load_Ini ( fini )
    use MOSE_Mod_MPI, only: mpi_is_root
    implicit none
    type(file_ini), intent(inout) :: fini
    ! MOSE_open sections
    character(len=vlen) :: sP, sN, sPH, sCH, sT, sIO, sPR, sMG, sRF
    ! GPU sections
    character(len=vlen) :: gp, gs, gr, go, gn, gpr, grf
    character(len=vlen) :: v, model, key
    logical :: navier_stokes
    integer :: i, nsc_file, ios, lev
    real(R8) :: x

    call fini%load(filename='input.ini')
    open_force_frozen = .false.
    nput = 0; nwarn = 0

    sP  = trim(codename)//'-Parameters';  sN  = trim(codename)//'-Numerics';  sPH = trim(codename)//'-Physics'
    sCH = trim(codename)//'-Chemistry';   sT  = trim(codename)//'-Turbulence'; sIO = trim(codename)//'-IO'
    sPR = trim(codename)//'-Probes';      sMG = trim(codename)//'-Multigrid';  sRF = trim(codename)//'-Rotating-Frame'
    gp  = trim(codename)//'-parameters';  gs  = trim(codename)//'-scheme';     gr  = trim(codename)//'-RANS'
    go  = trim(codename)//'-ODE';         gn  = trim(codename)//'-nochem';     gpr = trim(codename)//'-probes'
    grf = trim(codename)//'-rotating-frame'

    open_dialect = fini%has_section(section_name=trim(sP)) .or. fini%has_section(section_name=trim(sN)) .or. &
                   has(fini, sPH, 'equations')
    if (.not. open_dialect) return

    ! ---------------------------------------------------------------- [MOSE-Parameters] -> [MOSE-parameters]
    call copy(fini, sP, 'newrun',         gp, 'newrun')
    call copy(fini, sP, 'iter-threshold', gp, 'iter-threshold')
    call copy(fini, sP, 'time-threshold', gp, 'time-threshold')
    call copy(fini, sP, 'res-threshold',  gp, 'res-threshold')

    ! ---------------------------------------------------------------- [MOSE-Numerics]
    if (.not. has(fini, sN, 'cfl')) call stop_open('[' // trim(sN) // '] cfl is required')
    call copy(fini, sN, 'cfl',                gp, 'cfl')
    call copy(fini, sN, 'vnn',                gp, 'vnn')
    ! MOSE_open's time step has the CFL and the viscous VNN limits only (Mod_dt Compute_dt): the species-diffusion limit of
    ! the GPU (DIFF_MC, vnn-diff, default 0.8 with schmidt <= 0) is off unless the input sets it
    call copy(fini, sP, 'vnn-diff', gp, 'vnn-diff')
    call copy(fini, sN, 'vnn-diff', gp, 'vnn-diff')
    if (.not. has(fini, gp, 'vnn-diff')) &
      call put(fini, gp, 'vnn-diff', '0', 'MOSE_open has no species-diffusion time-step limit')
    call copy(fini, sN, 'cfl-rise-threshold', gp, 'cfl-rise-threshold')
    if (.not. has(fini, sN, 'time-accurate')) call stop_open('[' // trim(sN) // '] time-accurate is required')
    call copy(fini, sN, 'time-accurate',      gp, 'time-accurate')

    v = gets(fini, sN, 'time-scheme', '')
    select case (trim(v))
    case ('euler'); call put(fini, gs, 'RK', '1', 'time-scheme = euler')
    case ('RK2');   call put(fini, gs, 'RK', '2', 'time-scheme = RK2')
    case ('RK3');   call put(fini, gs, 'RK', '3', 'time-scheme = RK3')
    case ('');      call stop_open('[' // trim(sN) // '] time-scheme is required (euler, RK2 or RK3)')
    case default;   call stop_open('[' // trim(sN) // '] time-scheme = ' // trim(v) // ': euler, RK2 or RK3')
    end select

    v = gets(fini, sN, 'integration-variables', 'cons')
    select case (trim(v))
    case ('cons'); call put(fini, gs, 'integration-variables', 'cons', 'integration-variables')
    case ('prim'); call put(fini, gs, 'integration-variables', 'prim', 'integration-variables')
    case ('prec')   ! the Weiss-Smith preconditioning (Lib_Preconditioning) and its keys, the MOSE_open names
      if (as_logical(gets(fini, sN, 'time-accurate', 'false'))) &   ! MOSE_open Assign_Setup
        call stop_open('integration-variables=prec is currently supported only for steady/pseudo-time runs (as MOSE_open)')
      call put(fini, gs, 'integration-variables', 'prec', 'integration-variables')
      call copy(fini, sN, 'preconditioning-Uref', gs, 'preconditioning-Uref')
      call copy(fini, sN, 'preconditioning-Mach', gs, 'preconditioning-Mach')
      call copy(fini, sN, 'preconditioning-eps-min', gs, 'preconditioning-eps-min')
      call copy(fini, sN, 'preconditioning-Ur-min', gs, 'preconditioning-Ur-min')
      call copy(fini, sN, 'preconditioning-Ur-factor', gs, 'preconditioning-Ur-factor')
      call copy(fini, sN, 'preconditioning-Ur-smooth', gs, 'preconditioning-Ur-smooth')
    case default;  call stop_open('[' // trim(sN) // '] integration-variables = ' // trim(v) // ': cons, prim or prec')
    end select

    ! MOSE_open switches the implicit residual smoothing on with irs = true or irs-beta > 0 (Assign_Setup:67):
    ! [MOSE-scheme] IRS = t, beta = irs-beta; njacobi stays the GPU default 2 = MOSE_open Lib_IRS. Device: Mod_IRS_W.
    v = gets(fini, sN, 'irs-beta', '0.0')
    read(v, *, iostat=ios) x
    if (ios /= 0) x = 0.0d0
    if (x > 0.0d0 .or. as_logical(gets(fini, sN, 'irs', 'false'))) then
      call put(fini, gs, 'IRS',  't',     'irs / irs-beta = ' // trim(v) // ' (MOSE_open: smoothing on)')
      call put(fini, gs, 'beta', trim(v), 'irs-beta')
    end if

    v = gets(fini, sN, 'space-reconstruction', '')
    select case (trim(v))
    case ('MUSCL')
    case ('first-order'); call put(fini, gs, 'limiter', 'IORD', 'space-reconstruction = first-order')
    case ('');            call stop_open('[' // trim(sN) // '] space-reconstruction is required (MUSCL)')
    case default;         call stop_open('[' // trim(sN) // '] space-reconstruction = ' // trim(v) // ': MUSCL or first-order')
    end select

    v = gets(fini, sN, 'flux-limiter', '')
    if (trim(gets(fini, sN, 'space-reconstruction', '')) == 'first-order') v = 'first-order'   ! no limiter needed
    select case (trim(lower(v)))
    case ('first-order')
    case ('minmod');    call put(fini, gs, 'limiter', 'MINMOD', 'flux-limiter = minmod')
    case ('vanleer');   call put(fini, gs, 'limiter', 'VANLEER', 'flux-limiter = vanleer')
    case ('vanalbada'); call put(fini, gs, 'limiter', 'VANALBADA', 'flux-limiter = vanalbada')
    case ('mc');        call put(fini, gs, 'limiter', 'MC', 'flux-limiter = mc')
    case ('superbee');  call put(fini, gs, 'limiter', 'SUPERBEE', 'flux-limiter = superbee')
    case ('');          call stop_open('[' // trim(sN) // '] flux-limiter is required with MUSCL (minmod, vanleer, vanalbada, mc, superbee)')
    case default;       call stop_open('flux-limiter = ' // trim(v) // ': the GPU has minmod, vanleer, vanalbada, mc and superbee')
    end select

    ! shock sensor of MOSE_open (Tramel / Chen): the key, else the Riemann solver's own (Assign_SD, as Assign_Setup)
    v = gets(fini, sN, 'shock-detector', '')
    select case (trim(lower(v)))
    case ('');               ! the solver decides (HLLE++, HLLC+Tramel -> Tramel; AUSM+M -> Chen)
    case ('tramel', 'chen'); call put(fini, gs, 'shock-detector', trim(lower(v)), 'shock-detector')
    case default;            call stop_open('shock-detector = ' // trim(v) // ': Tramel or Chen (as MOSE_open)')
    end select

    v = gets(fini, sN, 'riemann-solver', 'HLLC')
    select case (trim(v))
    case ('HLLC');        call put(fini, gs, 'riemann', 'HLLC', 'riemann-solver = HLLC')
    case ('HLLE');        call put(fini, gs, 'riemann', 'HLLE', 'riemann-solver = HLLE')
    case ('SLAU');        call put(fini, gs, 'riemann', 'SLAU', 'riemann-solver = SLAU')
    case ('AUSM+');       call put(fini, gs, 'riemann', 'AUSM+', 'riemann-solver = AUSM+')
    case ('HLLE++');      call put(fini, gs, 'riemann', 'HLLE++', 'riemann-solver = HLLE++ (sensor Tramel unless given)')
    case ('HLLC+Tramel'); call put(fini, gs, 'riemann', 'HLLC+Tramel', 'riemann-solver = HLLC+Tramel (sensor Tramel unless given)')
    case ('AUSM+M');      call put(fini, gs, 'riemann', 'AUSM+M', 'riemann-solver = AUSM+M (sensor Chen unless given)')
    case ('LLF', 'Rusanov'); call put(fini, gs, 'riemann', 'LLF', 'riemann-solver = ' // trim(v))
    case ('SLAU2');       call put(fini, gs, 'riemann', 'SLAU2', 'riemann-solver = SLAU2')
    case ('HLLC+Chen');   call put(fini, gs, 'riemann', 'HLLC+Chen', 'riemann-solver = HLLC+Chen (sensor Chen unless given)')
    case ('exact');       call put(fini, gs, 'riemann', 'exact', 'riemann-solver = exact (Godunov)')
    case ('LMRoe');       call put(fini, gs, 'riemann', 'LMRoe', 'riemann-solver = LMRoe')
    case ('MiczekRoe');   call put(fini, gs, 'riemann', 'MiczekRoe', 'riemann-solver = MiczekRoe')
    case ('HLLC-PC')
      ! MOSE_open Assign_Setup: the preconditioned HLLC needs integration-variables = prec
      if (trim(gets(fini, sN, 'integration-variables', 'cons')) /= 'prec') &
        call stop_open('riemann-solver = HLLC-PC: integration-variables must be set to prec (as MOSE_open)')
      call put(fini, gs, 'riemann', 'HLLC-PC', 'riemann-solver = HLLC-PC')
    case default
      ! the message of the MOSE_open input check (Register_Numerics: the allowed values of riemann-solver)
      call stop_open('riemann-solver = ' // trim(v) // ': riemann-solver must be one of: HLLC, HLLC+Tramel, HLLC+Chen, &
                     &HLLC-PC, HLLE, HLLE++, SLAU, SLAU2, LMRoe, MiczekRoe, AUSM+, AUSM+M, LLF, Rusanov, exact')
    end select
    ! riemann-options-Mco (low-Mach cut-off): MOSE_open Assign_Setup requires it > 0 for AUSM+M, HLLC+Chen and LMRoe;
    ! MiczekRoe reads it too (floor 1e-6)
    key = gets(fini, sN, 'riemann-options-Mco', '')
    x = 0.0d0
    if (len_trim(key) > 0) then
      read(key, *, iostat=ios) x
      if (ios /= 0) call stop_open('riemann-options-Mco = ' // trim(key) // ' is not a number')
      call put(fini, gs, 'riemann-options-Mco', trim(key), 'riemann-options-Mco')
    end if
    select case (trim(v))
    case ('AUSM+M', 'HLLC+Chen', 'LMRoe')
      if (x <= 0.0d0) call stop_open('riemann-solver = ' // trim(v) // ' needs riemann-options-Mco > 0 (as MOSE_open)')
    end select

    ! MOSE_open Lib_Ghost fills the ghost cells of the extrapolation faces (400) by cubic extrapolation (400 is not in its
    ! zero-gradient list `case(0,401:407,420)` and falls in `case default`). The four keys below are also the defaults of
    ! the GPU code (IO_Ini): written here so that OUTPUT/input.translated.ini states them.
    if (.not. has(fini, gs, 'extrapolation-ghost')) &
      call put(fini, gs, 'extrapolation-ghost', 'cubic', 'MOSE_open Lib_Ghost: cubic ghost cells on the 400 faces')
    ! the same case default fills the ghosts of the viscous walls (301-304) by cubic extrapolation, and before the RANS
    ! sources Ghost_Wall_Extrapolation copies the boundary cell's turbulence into the wall ghost
    if (.not. has(fini, gs, 'wall-ghost')) &
      call put(fini, gs, 'wall-ghost', 'cubic', 'MOSE_open Lib_Ghost: cubic ghost cells on the viscous walls')
    if (.not. has(fini, gs, 'rans-wall-ghost')) &
      call put(fini, gs, 'rans-wall-ghost', 'copy', 'MOSE_open Ghost_Wall_Extrapolation: turbulence ghost = boundary cell')
    ! MOSE_open Newstate calls Check_And_Fix_State on every cell after every update: stop on p < 0 or a NaN
    if (.not. has(fini, gs, 'state-check')) &
      call put(fini, gs, 'state-check', 't', 'MOSE_open Check_And_Fix_State: stop on NaN or p < 0 after every update')

    ! ---------------------------------------------------------------- [MOSE-Multigrid]
    v = gets(fini, sMG, 'levels', '1')
    read(v, *, iostat=ios) lev
    if (ios /= 0) lev = 1
    if (lev > 1) call stop_open('[' // trim(sMG) // '] levels = ' // trim(v) // ': multigrid is not on the GPU')

    ! ---------------------------------------------------------------- [MOSE-Physics]
    v = gets(fini, sPH, 'equations', '')
    select case (trim(v))
    case ('euler');         navier_stokes = .false.
    case ('navier-stokes'); navier_stokes = .true.
    case ('');              call stop_open('[' // trim(sPH) // '] equations is required (euler or navier-stokes)')
    case default;           call stop_open('[' // trim(sPH) // '] equations = ' // trim(v) // ': euler or navier-stokes')
    end select

    if (navier_stokes) then
      model = gets(fini, sPH, 'turbulence', 'none')
      if (index(model, 'QCR') > 0) call stop_open('turbulence = ' // trim(model) // ': the QCR stress is not on the GPU')
      select case (trim(model))
      case ('none', ''); model = 'laminar'
      case ('SA', 'SA-R', 'SA-RC', 'SAcomp', 'SST', 'SST-RC', 'Wilcox2006')
      case ('SSGLRR')
        call stop_open('turbulence = ' // trim(model) // ' is not on the GPU: use SST, SST-RC, Wilcox2006, SA, SA-R, SA-RC, SAcomp or none')
      case default;    call stop_open('[' // trim(sPH) // '] turbulence = ' // trim(model) // ' is not a model name')
      end select
      ! [MOSE-RANS] exists = viscous for the GPU code; turbulence-model = laminar is Navier-Stokes without a model
      call put(fini, gr, 'turbulence-model', model, 'equations = navier-stokes, turbulence')
      if (trim(model) /= 'laminar') then
        call copy(fini, sT, 'Prt', gr, 'Prt')
        call copy(fini, sT, 'Sct', gr, 'Sct')
        if (as_logical(gets(fini, sT, 'point-implicit', 'true'))) then
          call put(fini, gr, 'point-implicit', 't', '[' // trim(sT) // '] point-implicit (MOSE_open default true)')
        else
          call put(fini, gr, 'point-implicit', 'f', '[' // trim(sT) // '] point-implicit')
        end if
        v = gets(fini, sT, 'omega-wall-bc', 'practical')
        select case (trim(v))
        case ('practical', 'asymptotic')
          call put(fini, gr, 'omega-wall-bc', trim(v), '[' // trim(sT) // '] omega-wall-bc')
        case default
          call stop_open('[' // trim(sT) // '] omega-wall-bc = ' // trim(v) // ': practical or asymptotic')
        end select
        v = gets(fini, sT, 'sst-production', 'incompressible')
        select case (trim(v))
        case ('incompressible', 'compressible')
          call put(fini, gr, 'sst-production', trim(v), '[' // trim(sT) // '] sst-production')
        case default
          call stop_open('[' // trim(sT) // '] sst-production = ' // trim(v) // ': incompressible or compressible')
        end select
        ! k in the total energy = [MOSE-RANS] k-coupling: the -2/3 rho k delta_ij of the stress (stress_dev = 2 in
        ! diff_physics.inc, Compute_Diffusive_Flux, Lib_BC) and the k source out of the mean-flow energy (sst_body.inc,
        ! wilcox_body.inc): the terms of MOSE_open Lib_Fluid Stress_Vector_2eq, Lib_SST and Lib_Wilcox2006 (A10 of 2026-10-16)
        if (as_logical(gets(fini, sT, 'k-energy-coupling', 'false'))) then
          call put(fini, gr, 'k-coupling', 't', '[' // trim(sT) // '] k-energy-coupling')
        else
          call put(fini, gr, 'k-coupling', 'f', '[' // trim(sT) // '] k-energy-coupling')
        end if
      end if
      ! laminar species diffusion: MOSE_open default schmidt = 0.0 = mixture-averaged multicomponent (INPUT/diffusion.dat)
      if (.not. has(fini, sPH, 'schmidt')) then
        nsc_file = species_count()
        if (nsc_file > 1) then
          call put(fini, sPH, 'schmidt', '0.0', 'MOSE_open default: multicomponent diffusion, needs INPUT/diffusion.dat')
        else
          call warn_open('schmidt absent with one species: no species diffusion, the GPU default Sc is irrelevant')
        end if
      end if
      ! [MOSE-Physics] prandtl: the same key in both dialects, read by IO_ThermoTransport (kl = mu cp / Pr when > 0)
    else
      if (fini%has_section(section_name=trim(sT))) &
        call warn_open('[' // trim(sT) // '] is ignored with equations = euler')
      if (has(fini, sPH, 'turbulence')) then
        v = gets(fini, sPH, 'turbulence', 'none')
        if (trim(v) /= 'none' .and. len_trim(v) > 0) &
          call warn_open('turbulence = ' // trim(v) // ' is ignored with equations = euler (as MOSE_open does)')
      end if
    end if

    v = gets(fini, sPH, 'chemistry', 'frozen')
    select case (trim(v))
    case ('frozen')
      if (file_exists('INPUT/chemistry-info.txt')) then
        open_force_frozen = .true.
        call note_open('chemistry = frozen: INPUT/chemistry-info.txt is present but not read (frozen run)')
      end if
    case ('finite-rate')
      if (.not. file_exists('INPUT/chemistry-info.txt')) &
        call stop_open('chemistry = finite-rate needs INPUT/chemistry-info.txt (the mechanism)')
    case ('equilibrium'); call stop_open('chemistry = equilibrium is not on the GPU (frozen or finite-rate)')
    case default;         call stop_open('[' // trim(sPH) // '] chemistry = ' // trim(v) // ': frozen, finite-rate or equilibrium')
    end select

    v = gets(fini, sPH, 'soot-generation', 'none')
    if (trim(v) /= 'none' .and. len_trim(v) > 0) call stop_open('soot-generation = ' // trim(v) // ' is not on the GPU')

    v = gets(fini, sPH, 'rotational-frame', 'none')
    select case (trim(v))
    case ('none', '')
    case ('rigid-body')
      call put(fini, grf, 'active', 't', 'rotational-frame = rigid-body')
      call copy(fini, sRF, 'omega',  grf, 'omega')
      call copy(fini, sRF, 'axis',   grf, 'axis')
      call copy(fini, sRF, 'origin', grf, 'origin')
      if (has(fini, sRF, 'stationary-face-1')) &
        call stop_open('[' // trim(sRF) // '] stationary-face-N is not on the GPU rotating frame')
    case default; call stop_open('[' // trim(sPH) // '] rotational-frame = ' // trim(v) // ': rigid-body or none')
    end select

    ! ---------------------------------------------------------------- [MOSE-Chemistry] -> [MOSE-ODE], [MOSE-nochem]
    v = gets(fini, sCH, 'ode-solver', '')
    select case (trim(v))
    case ('')
    case ('H-radau5', 'H-rodas'); call put(fini, go, 'solver', v, 'ode-solver')
    case default; call stop_open('ode-solver = ' // trim(v) // ': the GPU integrates the chemistry with H-radau5 or H-rodas')
    end select
    call copy(fini, sCH, 'ode-max-steps',                go, 'max-steps-ode')
    call copy(fini, sCH, 'ode-analytical-jacobian',      go, 'analytical-jacobian')
    call copy(fini, sCH, 'ode-relative-tol-species',     go, 'relative-tol-species')
    call copy(fini, sCH, 'ode-absolute-tol-species',     go, 'absolute-tol-species')
    call copy(fini, sCH, 'ode-relative-tol-temperature', go, 'relative-tol-pressure')
    call copy(fini, sCH, 'ode-absolute-tol-temperature', go, 'absolute-tol-pressure')
    v = gets(fini, sCH, 'exclude-blocks', 'none')
    if (trim(v) /= 'none' .and. len_trim(v) > 0) call put(fini, gn, 'blocks', v, 'exclude-blocks')

    ! ---------------------------------------------------------------- [MOSE-IO] (same section, two renamed keys)
    call copy(fini, sIO, 'ic-format', sIO, 'ini-format')
    ! a restart (newrun = false) reads OUTPUT/field in the solution format (MOSE_open Setup_Input_Solution: sol-format,
    ! default tecplot ascii); the GPU code reads it in bck-format, whose own default is native
    if (.not. has(fini, sIO, 'bck-format')) &
      call put(fini, sIO, 'bck-format', gets(fini, sIO, 'sol-format', 'tecplot ascii'), 'MOSE_open restarts from sol-format')
    call copy(fini, sIO, 'res-diter', sIO, 'diter-res')
    if (has(fini, sIO, 'timer-diter')) call warn_open('[' // trim(sIO) // '] timer-diter is not on the GPU: ignored')

    ! ---------------------------------------------------------------- [MOSE-Probes] -> [MOSE-probes]
    if (fini%has_section(section_name=trim(sPR))) then
      i = 0
      do
        i = i + 1
        write(key, '(A,I0)') 'probe', i
        if (.not. has(fini, sPR, key)) exit
        call copy(fini, sPR, key, gpr, key)
      end do
    end if

    ! ---------------------------------------------------------------- the report, once, root rank
    if (first .and. mpi_is_root) then
      ! what the GPU code actually reads: the translated ini, saved next to the outputs (the input file is never touched)
      call fini%save(filename='OUTPUT/input.translated.ini', iostat=ios)
      if (ios == 0) write(*,'(A)') ' [MOSE_open input] the translated keys are saved in OUTPUT/input.translated.ini'
      write(*,'(A)') ' [MOSE_open input] input.ini is written in the MOSE_open dialect: keys read by the GPU code'
      do i = 1, nput
        write(*,'(A)') '   ' // trim(putlog(i))
      end do
      do i = 1, nwarn
        write(*,'(A)') ' [MOSE_open input] WARNING: ' // trim(warnlog(i))
      end do
    end if
    first = .false.

  end subroutine MOSE_Load_Ini


  ! ------------------------------------------------------------------ helpers
  function has ( fini, sec, key ) result ( r )
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: sec, key
    logical :: r
    r = fini%has_option(section_name=trim(sec), option_name=trim(key))
  end function has


  !> The value of a key as a string, or the default when absent.
  function gets ( fini, sec, key, dflt ) result ( v )
    implicit none
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: sec, key, dflt
    character(len=vlen) :: v
    integer :: err, i
    v = dflt
    if ( fini%has_option(section_name=trim(sec), option_name=trim(key)) ) then
      call fini%get(section_name=trim(sec), option_name=trim(key), val=v, error=err)
      if (err /= 0) v = dflt
    end if
    do i = 1, len(v)   ! CRLF files and tabs: FiNeR keeps the CR and the tab in the value
      if (v(i:i) == achar(13) .or. v(i:i) == achar(9)) v(i:i) = ' '
    end do
    v = adjustl(v)
  end function gets


  !> Add (or update) a key the GPU code reads, and log it for the report.
  subroutine put ( fini, sec, key, val, from )
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=*), intent(in) :: sec, key, val, from
    integer :: err
    character(len=:), allocatable :: tv
    tv = trim(val)   ! an exact-length variable: a trim() temporary into the class(*) dummy of FiNeR add segfaults (nvfortran 24.11)
    call fini%add(section_name=trim(sec), option_name=trim(key), val=tv, error=err)
    if (err /= 0) call stop_open('cannot add [' // trim(sec) // '] ' // trim(key))
    if (nput < nlog) then
      nput = nput + 1
      putlog(nput) = '[' // trim(sec) // '] ' // trim(key) // ' = ' // trim(val) // '   <- ' // trim(from)
    end if
  end subroutine put


  subroutine copy ( fini, sfrom, kfrom, sto, kto )
    implicit none
    type(file_ini), intent(inout) :: fini
    character(len=*), intent(in) :: sfrom, kfrom, sto, kto
    if ( has(fini, sfrom, kfrom) ) &
      call put( fini, sto, kto, gets(fini, sfrom, kfrom, ''), '[' // trim(sfrom) // '] ' // trim(kfrom) )
  end subroutine copy


  subroutine stop_open ( msg )
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(A)') ' [MOSE_open input] ERROR: ' // trim(msg)
    error stop 'MOSE_open input: the GPU code cannot run this input.ini as it is (see the message above)'
  end subroutine stop_open


  subroutine warn_open ( msg )
    implicit none
    character(len=*), intent(in) :: msg
    if (nwarn < nlog) then
      nwarn = nwarn + 1
      warnlog(nwarn) = msg
    end if
  end subroutine warn_open


  subroutine note_open ( msg )
    implicit none
    character(len=*), intent(in) :: msg
    if (nput < nlog) then
      nput = nput + 1
      putlog(nput) = msg
    end if
  end subroutine note_open


  !> Fortran list-directed logical: the first letter after an optional dot.
  function as_logical ( s ) result ( r )
    implicit none
    character(len=*), intent(in) :: s
    logical :: r
    character(len=vlen) :: t
    t = adjustl(s)
    if (t(1:1) == '.') t = t(2:)
    r = ( t(1:1) == 't' .or. t(1:1) == 'T' )
  end function as_logical


  function lower ( s ) result ( r )
    implicit none
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i, c
    r = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) r(i:i) = achar(c + 32)
    end do
  end function lower


  function file_exists ( path ) result ( r )
    implicit none
    character(len=*), intent(in) :: path
    logical :: r
    inquire(file=path, exist=r)
  end function file_exists


  !> The species count the code reads: the lines after the first of INPUT/phase.txt (-1 when the file is missing).
  function species_count () result ( n )
    implicit none
    integer :: n, u, ios
    character(len=vlen) :: line
    logical :: ex
    n = -1
    inquire(file='INPUT/phase.txt', exist=ex)
    if (.not. ex) return
    open(newunit=u, file='INPUT/phase.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    n = -1
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n = n + 1
    end do
    close(u)
  end function species_count

end module MOSE_IO_Ini_Open
