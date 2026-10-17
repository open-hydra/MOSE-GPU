!> MOSE_open (ATLAS) bc.txt dialect support.
!>
!> MOSE_GPU keeps its legacy integer BC ids (1 = connection, 4 = inflow/outflow,
!> 11 = extrapolation, ...) all the way down to the flat device arrays
!> (domain%bc_typ) and the device kernel dispatch. MOSE_open/ATLAS writes the
!> same boundary-face list with 3-digit ids and a different property-line layout.
!>
!> This module is the translation layer: it detects the dialect, skips the right
!> number of property lines during the consistency prescan, and reads an open
!> record straight into the legacy fields of MOSE_bc_type, rewriting bc%type to
!> the legacy id. Nothing downstream of Read_BCfile changes -- the solver, the
!> ghost fill and the device path never see a 3-digit id.
!>
!> Mappings that are NOT exactly representable in the legacy field set are
!> refused with an explicit message (never silently approximated).
module MOSE_IO_BC_Open
  use MOSE_Types_m
  use MOSE_Global_m
  use MOSE_IO_Wall, only: viscous_flag, coupling_flag

  implicit none
  private

  public :: BC_DIALECT_LEGACY, BC_DIALECT_OPEN
  public :: Detect_BC_Dialect
  public :: Skip_BC_Record
  public :: Read_BC_Record_Open
  public :: Resolve_Open_BC_Direction

  integer, parameter :: BC_DIALECT_LEGACY = 0
  integer, parameter :: BC_DIALECT_OPEN   = 1

  !> Sentinel written into bc%pamb by the inlet translations. The legacy
  !> BC_Inflow_Outflow reads it as "pure inlet": when the local flow turns
  !> outward the face falls back to a reflective (symmetry) boundary, which is
  !> exactly what every MOSE_open BC_Inlet_* routine does (error = 1).
  real(R8), parameter :: PAMB_FORCE_INFLOW = -3.14d0

  !> Legacy mach-column encodings of the inflow modes (see BC_Subsonic_Inflow):
  !>   >= 0  : stagnation (p0, T0) inlet   -- open 401/402/407
  !>   -10   : mass flux g with TOTAL T    -- open 403
  !>   < 0, /= -10, /= -15 : mass flux g with STATIC T -- open 404
  real(R8), parameter :: MACH_STAGNATION = 0.0d0
  real(R8), parameter :: MACH_MASSFLUX_T0 = -10.0d0
  real(R8), parameter :: MACH_MASSFLUX_T  = -1.0d0

contains

  !> Decide which dialect a bc file is written in by walking its records.
  !>
  !> Both dialects agree that a type-0 (degenerate/none) record carries no
  !> property line, so records can be walked safely until the first typed record
  !> is met -- that record decides, since the two id sets are disjoint.
  function Detect_BC_Dialect ( unitfile ) result ( dialect )
    implicit none
    integer, intent(in) :: unitfile
    integer             :: dialect
    ! Local
    integer :: di(5), ti, ios

    dialect = BC_DIALECT_LEGACY          ! empty file: harmless, the prescan will complain

    rewind(unitfile)
    ios = 0
    do while (ios == 0)
      read( unitfile,*,iostat=ios ) di(1), di(2), di(3), di(4), di(5), ti
      if (ios /= 0) exit
      select case (ti)
        case (0)
          cycle                          ! no property line in either dialect: keep looking
        case (100,101,102,103,200,201,300,301:309,400:410,420,501,502)
          dialect = BC_DIALECT_OPEN
          exit
        case (1,2,3,4,5,6,8,9,10,11,12,13,14,22,23,666,667,999,1000,1001)
          dialect = BC_DIALECT_LEGACY
          exit
        case default
          write(*,*) '  [ERROR] Unknown boundary condition type in bc file: ', ti
          error stop 'IO_BC: unknown BC type'
      end select
    end do

    rewind(unitfile)

  end function Detect_BC_Dialect


  !> Consistency-prescan helper: consume the property line(s) of one record.
  !> Mirrors the read loop below (and MOSE_open's Check_BC) exactly -- a wrong
  !> line count here desynchronises the boundary-face count against the IC.
  subroutine Skip_BC_Record ( unitfile, ti, dialect, ios )
    implicit none
    integer, intent(in)    :: unitfile, ti, dialect
    integer, intent(inout) :: ios
    ! Local
    integer :: nd(2), s

    if (dialect == BC_DIALECT_LEGACY) then
      select case (ti)
        case (1,4,5,6,8,9,10,12,13,14,22,1000,1001,999,667)
          read( unitfile,*,iostat=ios )
        case (666)
          ! Chimera: variable-length record -- ni(1)+ni(2) donor lines follow the
          ! ni line (one per donor, matching the reader and MOSE_open's prescan).
          read( unitfile,*,iostat=ios ) nd(1), nd(2)
          do s = 1, nd(1)+nd(2)
            read( unitfile,*,iostat=ios )
          end do
      end select
    else
      select case (ti)
        case (101,103,201,301:309,401:408,410,420,501,502)
          read( unitfile,*,iostat=ios )
        case (102)
          read( unitfile,*,iostat=ios ) nd(1), nd(2)
          do s = 1, nd(1)+nd(2)
            read( unitfile,*,iostat=ios )
          end do
        ! 0, 200, 300, 400, 409: no property line
      end select
    endif

  end subroutine Skip_BC_Record


  !> Read the property line(s) of one MOSE_open record into the legacy fields of
  !> bc, and rewrite bc%type to the legacy id. The record header (b,i,j,k,f,type)
  !> has already been read by the caller.
  subroutine Read_BC_Record_Open ( unitfile, bc, ios )
    implicit none
    integer, intent(in)                 :: unitfile
    type(MOSE_bc_type), intent(inout)   :: bc
    integer, intent(inout)              :: ios
    ! Local
    integer            :: cc, s, ip
    integer            :: open_type
    real(R8)           :: k_rough, eps_wall, rho_grain, sf_geo, psub, psup
    character(len=32)  :: p0file, alpha_tok, beta_tok
    character(len=256) :: q2d_line

    open_type = bc % type
    bc % open_type = open_type

    select case (open_type)

      ! ── 0 / 100: face with no boundary condition. No property line. ─────────
      ! ATLAS writes 100 on faces it deliberately leaves unassigned (typically
      ! the blanked faces of a chimera overlap, interleaved with the 102
      ! records); MOSE_open has no case for it either, so no flux is ever
      ! computed there. Legacy type 0 is the same "no flux" face, with a
      ! zero-gradient ghost instead of MOSE_open's extrapolated one.
      case (0, 100)
        bc % type = 0

      ! ── 101/201: block connection (201 = translational periodic connection) ──
      ! MOSE_open dispatches both exactly like a connection (case(101,102,201)).
      case (101, 201)
        read( unitfile,*,iostat=ios ) &
          bc%bs, bc%is, bc%js, bc%ks, bc%fs, bc%d11, bc%d12, bc%d21, bc%d22
        allocate ( bc % Pg (nprim, 6) )
        bc % type = 1

      ! ── 102: chimera overlap ────────────────────────────────────────────────
      case (102)
        read( unitfile,*,iostat=ios ) (bc%ni(cc),cc=1,2)
        allocate ( bc%donorID(1:sum(bc%ni),1:4) )
        allocate ( bc%volume_fraction(1:sum(bc%ni)) )
        do s = 1, sum(bc%ni)
          read( unitfile,*,iostat=ios ) bc%donorID(s,1:4), bc%volume_fraction(s)
        end do
        allocate ( bc % Pg (nprim, 6) )
        bc % type = 666

      ! ── 103: coupled multi-solver wall ──────────────────────────────────────
      case (103)
        coupling_flag( bc%b , bc%f ) = .true.
        read( unitfile,*,iostat=ios ) &
          bc%bs, bc%is, bc%js, bc%ks, bc%fs, bc%d11, bc%d12, bc%d21, bc%d22
        allocate ( bc%ext_flux(nprim) )
        bc%ext_flux = 0.0d0
        allocate ( bc % Pg (1, 6) )
        bc % type = 1000

      ! ── 200: rotational (axisymmetric) periodicity. No property line. ───────
      case (200)
        bc % type = 2

      ! ── 300: Euler symmetry. No property line. ──────────────────────────────
      case (300)
        bc % type = 3

      ! ── 400: extrapolation. No property line. ───────────────────────────────
      case (400)
        bc % type = 11

      ! ── 301: viscous wall, prescribed heat flux (q, ks, eps) ────────────────
      case (301)
        viscous_flag( bc%b , bc%f ) = .true.
        read( unitfile,*,iostat=ios ) bc%qw, k_rough, eps_wall
        call Refuse_Unsupported_Wall_Extras ( open_type, k_rough, eps_wall )
        bc % type = 5

      ! ── 302: viscous wall, prescribed temperature (T, ks, eps) ──────────────
      case (302)
        viscous_flag( bc%b , bc%f ) = .true.
        read( unitfile,*,iostat=ios ) bc%Tw, k_rough, eps_wall
        call Refuse_Unsupported_Wall_Extras ( open_type, k_rough, eps_wall )
        bc % type = 6

      ! ── 303: viscous wall, temperature + radiative flux (melting) ───────────
      case (303)
        viscous_flag( bc%b , bc%f ) = .true.
        read( unitfile,*,iostat=ios ) bc%Tw, bc%qrad, k_rough
        call Refuse_Unsupported_Wall_Extras ( open_type, k_rough, 0.0d0 )
        bc % type = 8

      ! ── 304: viscous wall, radiative flux only (pyrolysis/ablation) ─────────
      ! MOSE_open passes 304 to BC_Wall_Pyrolysis, whose material select-case
      ! groups it with the legacy id 9 (HDPE / generic ablation).
      case (304)
        viscous_flag( bc%b , bc%f ) = .true.
        read( unitfile,*,iostat=ios ) bc%qrad, k_rough
        call Refuse_Unsupported_Wall_Extras ( open_type, k_rough, 0.0d0 )
        bc % type = 9

      ! ── 401: subsonic inlet, stagnation T0 + p0 ─────────────────────────────
      case (401)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%p0, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        bc%mach = MACH_STAGNATION
        bc%pamb = PAMB_FORCE_INFLOW
        bc % type = 4

      ! ── 402: subsonic inlet, stagnation T0 + time-varying p0 ────────────────
      case (402)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, p0file, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call bc%p0time%initialize(file=p0file, bar=.true.)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        bc%mach = MACH_STAGNATION
        bc%pamb = PAMB_FORCE_INFLOW
        bc % type = 4

      ! ── 403: inlet, mass flux g + TOTAL temperature ─────────────────────────
      ! Legacy carries g in the p0 column and selects the total-T closure with
      ! mach == -10 (identical static-T recovery formula).
      case (403)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%p0, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        bc%mach = MACH_MASSFLUX_T0
        bc%pamb = PAMB_FORCE_INFLOW
        bc % type = 4

      ! ── 404: inlet, mass flux g + STATIC temperature ────────────────────────
      ! Same closure as the legacy mach<0 (/= -10, /= -15) branch, EXCEPT that
      ! MOSE_open relaxes the imposed g towards the interior mass flux with
      ! rel_fac. The legacy branch has no such relaxation, so only rel_fac == 1
      ! (no relaxation) is an exact translation.
      case (404)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%p0, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        if (bc%rel_fac /= 1.0d0) then
          write(*,*) '  [ERROR] BC 404 (mass flux + static T) with rel_fac = ', bc%rel_fac
          write(*,*) '          MOSE_open relaxes g towards the interior mass flux; the MOSE_GPU'
          write(*,*) '          inflow closure has no such relaxation. Only rel_fac = 1 is supported.'
          error stop 'IO_BC (MOSE_open dialect): BC 404 with rel_fac /= 1 not supported'
        endif
        bc%mach = MACH_MASSFLUX_T
        bc%pamb = PAMB_FORCE_INFLOW
        bc % type = 4

      ! ── 405: supersonic inlet, Mach + static T + static p ───────────────────
      ! The legacy supersonic closure (mach > 1) takes its p0/T0 arguments as
      ! static values, so the open fields map straight across; the direction
      ! follows MOSE_open (Set_Inflow_Direction / Resolve_Open_BC_Direction).
      case (405)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%mach, bc%T0, bc%p0, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        if (bc%mach <= 1.0d0) then
          write(*,*) '  [ERROR] BC 405 (supersonic inlet) with Mach = ', bc%mach
          write(*,*) '          The legacy inflow dispatch only reaches the supersonic closure'
          write(*,*) '          for Mach > 1; a subsonic Mach would silently switch to the'
          write(*,*) '          stagnation closure and reinterpret p as a total pressure.'
          error stop 'IO_BC (MOSE_open dialect): BC 405 requires Mach > 1'
        endif
        bc%pamb = PAMB_FORCE_INFLOW
        bc % type = 4

      ! ── 407: ambient / far field, T0 + p0 (inflow) or back pressure (outflow)
      ! This is exactly the legacy type-4 dual behaviour: stagnation inlet when
      ! the flow enters, back-pressure outlet when it leaves.
      case (407)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%p0, bc%pamb, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        bc%mach = MACH_STAGNATION
        bc % type = 4

      ! ── 410: mapped (Q2D) state from a time-varying file ────────────────────
      ! Legacy id 667. The open record separates file name and options with a
      ! comma, the legacy one with a blank.
      case (410)
        allocate ( bc % Pg (nprim, 6) )
        read( unitfile,'(A)',iostat=ios ) q2d_line
        q2d_line = adjustl(q2d_line)
        ip = index(trim(q2d_line), ',')
        if (ip > 0) then
          bc%q2d_file = q2d_line(1:ip-1)
          if (index(q2d_line(ip:), 'periodic') > 0) bc%q2d_periodic = .true.
        else
          bc%q2d_file = trim(q2d_line)
        endif
        bc % type = 667

      ! ── 420: choked-nozzle inlet ────────────────────────────────────────────
      ! The legacy BC_Inflow_Nozzle carries psub/psup in the alpha/beta columns
      ! *in bar* (it multiplies both by 1e5) and the mass flux g in the pamb
      ! column; mach == 0 selects the g-driven supersonic branch.
      case (420)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%p0, psub, psup, bc%pamb, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        bc%alpha = psub * 1.0d-5
        bc%beta  = psup * 1.0d-5
        bc%mach  = 0.0d0
        bc % type = 22

      ! ── 501: manifold mean connection ───────────────────────────────────────
      case (501)
        read( unitfile,*,iostat=ios ) bc%bs, bc%fs
        bc % type = 1001

      ! ── 502: SRM grain boundary ─────────────────────────────────────────────
      ! Open order: Taf, a, n, pRef, rhoGrain, SF, massf(1:nsc). haf is no longer
      ! written by ATLAS (zero); rhoGrain/SF_geo have no legacy counterpart.
      case (502)
        viscous_flag( bc%b , bc%f ) = .true.
        allocate( bc % ci(1 : nsc) )
        bc%haf = 0.0d0
        read( unitfile,*,iostat=ios ) bc%Taf, bc%aCoeff, bc%n, bc%pRef, &
          rho_grain, sf_geo, (bc%ci(s), s = 1, nsc)
        bc%SF = sf_geo
        bc % type = 14

      ! ── 406: back-pressure outlet ───────────────────────────────────────────
      ! No legacy counterpart: the legacy outflow branch neither relaxes the exit
      ! pressure nor closes the face on reversed flow. Internal type 24
      ! (BC_Outlet) implements MOSE_open's contract; host path only, the device
      ! capability guard stops on it.
      case (406)
        read( unitfile,*,iostat=ios ) bc%pamb, bc%rel_fac
        bc % type = 24

      ! ── 408: inlet, prescribed normal velocity + static T ───────────────────
      ! No legacy counterpart (the legacy closures prescribe a mass flux or
      ! stagnation conditions, never a velocity). Internal type 25
      ! (BC_Inlet_Velocity); host path only, the device guard stops on it.
      case (408)
        allocate( bc % ci(1 : nsc+nrans) )
        read( unitfile,*,iostat=ios ) &
          bc%T0, bc%un, alpha_tok, beta_tok, bc%rel_fac, &
          (bc%ci(s), s = nrans+1, nrans+nsc), (bc%ci(s), s = 1, nrans)
        call Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
        bc % type = 25

      ! ── Types that the legacy field set cannot express exactly ──────────────
      case (409)
        write(*,*) '  [ERROR] BC 409 (forced outlet) is not available in MOSE_GPU.'
        error stop 'IO_BC (MOSE_open dialect): BC 409 not supported'

      case (305:309)
        write(*,*) '  [ERROR] Wall BC subtype ', open_type, ' has no MOSE_GPU counterpart.'
        error stop 'IO_BC (MOSE_open dialect): unsupported wall subtype'

      case default
        write(*,*) '  [ERROR] Unknown MOSE_open boundary condition type: ', open_type
        error stop 'IO_BC (MOSE_open dialect): unknown BC type'

    end select

    if (ios /= 0) write(*,*) '  Error in BC file (MOSE_open record, type ', open_type, ')'

  end subroutine Read_BC_Record_Open


  !> Post-metrics resolution of the inflow directions of a MOSE_open bc file (MOSE_open parse_dir_tok +
  !> Compute_Inflow_Direction: each component is an explicit angle, used as it is, or the token 'normal',
  !> alpha = atan(n2/(n1+1e-20)), beta = atan(n3/(n1+1e-20)) of the face normal; explicit 0, 0 is +x).
  !>
  !> The subsonic and mass-flux closures (401-404, 407 -> legacy 4, 408 -> 25) evaluate that same expression at
  !> every call when both columns are zero (the legacy rule). A record with both components 'normal' therefore
  !> keeps 0, 0 and the closure resolves it on the device at every call, as MOSE_open does (and as the production
  !> cases, which read ATLAS bc files, always ran). A 'normal' component next to an explicit one is resolved here;
  !> a pair that ends at 0, 0 on a face whose normal is not along x is refused: the legacy rule would turn it into
  !> the face normal, MOSE_open means +x. The supersonic closure (405) has no face-normal rule and takes the
  !> angles as absolute: its 'normal' components are always resolved here, explicit zeros are +x as in MOSE_open.
  subroutine Resolve_Open_BC_Direction ( domain )
    use MOSE_Lib_BC, only: Face_Index
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer  :: i, Bm, Fm, dir, fi, fj, fk, ot
    logical  :: an, bn, rule
    real(R8) :: normal(3), alpha_n, beta_n
    real(R8), parameter :: ANGLE_TOL = 1.0d-10

    do i = 1, domain % nbound
      ot = domain % bc(i) % open_type
      select case (ot)
      case (401:405, 407, 408)
      case default
        cycle
      end select
      an   = domain % bc(i) % alpha_from_normal
      bn   = domain % bc(i) % beta_from_normal
      rule = ( ot /= 405 )                       ! the closure applies the face-normal rule to 0, 0
      if ( rule .and. an .and. bn ) cycle        ! both 'normal': the closure's own per-call rule
      if ( .not. rule .and. .not. (an .or. bn) ) cycle   ! 405 with explicit angles: absolute, nothing to do
      Bm = domain % bc(i) % b
      if (.not. is_local_block(Bm)) cycle
      Fm = domain % bc(i) % f
      call Face_Index ( Fm, dir, domain % bc(i) % i, domain % bc(i) % j, domain % bc(i) % k, fi, fj, fk )
      normal  = domain % blk(Bm) % dir(dir) % f(fi,fj,fk) % n
      alpha_n = atan( normal(2) / ( normal(1) + 1d-20 ) )
      beta_n  = atan( normal(3) / ( normal(1) + 1d-20 ) )
      if (an) domain % bc(i) % alpha = alpha_n
      if (bn) domain % bc(i) % beta  = beta_n
      if ( rule .and. domain % bc(i) % alpha == 0.0d0 .and. domain % bc(i) % beta == 0.0d0 ) then
        if ( abs(alpha_n) > ANGLE_TOL .or. abs(beta_n) > ANGLE_TOL ) then
          write(*,*) '  [ERROR] BC ', ot, ' with flow direction 0, 0 (explicit) on a face whose normal is not along x'
          write(*,*) '          (block ', Bm, ' face ', Fm, ', angles of the normal ', alpha_n, beta_n, ').'
          write(*,*) '          MOSE_open reads 0, 0 as the +x direction; the GPU subsonic closure reads it as'
          write(*,*) '          "along the face normal". Write the direction as normal, normal or as explicit angles.'
          error stop 'IO_BC (MOSE_open dialect): explicit 0, 0 inflow direction on a face not normal to x'
        end if
      end if
    end do

  end subroutine Resolve_Open_BC_Direction


  !> Translate the open direction tokens into the legacy alpha/beta columns.
  !>
  !> MOSE_open writes alpha/beta as character tokens: a number (radians) or the
  !> word 'normal', resolved per face as atan(n2/n1), atan(n3/n1). The legacy
  !> subsonic/mass-flux closure applies that very same formula whenever both
  !> columns are zero, so 'normal' translates to a plain 0 -- exactly, not
  !> approximately.
  subroutine Set_Inflow_Direction ( bc, alpha_tok, beta_tok, open_type )
    implicit none
    type(MOSE_bc_type), intent(inout) :: bc
    character(len=*),   intent(in)    :: alpha_tok, beta_tok
    integer,            intent(in)    :: open_type

    ! MOSE_open parse_dir_tok + Compute_Inflow_Direction: each component is either an explicit angle
    ! (radians, used as it is) or the token 'normal', the angle of the face normal; the two are independent
    ! (mixing is allowed). 'normal' needs the metrics, so it is flagged here and Resolve_Open_BC_Direction
    ! decides per record: both 'normal' on a subsonic/mass-flux closure stay 0, 0 (the closure's own rule),
    ! the other cases are resolved there.
    bc%alpha_from_normal = Is_Normal_Token ( alpha_tok )
    bc%beta_from_normal  = Is_Normal_Token ( beta_tok )
    bc%alpha = 0.0d0
    bc%beta  = 0.0d0
    if (.not. bc%alpha_from_normal) bc%alpha = Parse_Angle_Token ( alpha_tok, open_type )
    if (.not. bc%beta_from_normal)  bc%beta  = Parse_Angle_Token ( beta_tok,  open_type )

  end subroutine Set_Inflow_Direction


  pure function Is_Normal_Token ( tok ) result ( is_normal )
    implicit none
    character(len=*), intent(in) :: tok
    logical                      :: is_normal
    character(len=len(tok))      :: tok_

    tok_ = adjustl(tok)
    is_normal = index(trim(tok_), 'normal') > 0

  end function Is_Normal_Token


  function Parse_Angle_Token ( tok, open_type ) result ( val )
    implicit none
    character(len=*), intent(in) :: tok
    integer,          intent(in) :: open_type
    real(R8)                     :: val
    ! Local
    integer                 :: ios_loc
    character(len=len(tok)) :: tok_

    tok_ = adjustl(tok)
    read(tok_, *, iostat=ios_loc) val
    if (ios_loc /= 0) then
      write(*,*) '  [ERROR] BC ', open_type, ': unreadable flow-direction token "'//trim(tok_)//'"'
      error stop 'IO_BC (MOSE_open dialect): bad flow-direction token'
    endif

  end function Parse_Angle_Token


  !> Wall records carry roughness / emissivity columns that the MOSE_GPU wall
  !> closures do not model. Zero values are the ATLAS default and are dropped
  !> silently; anything else would change the physics, so it is refused.
  subroutine Refuse_Unsupported_Wall_Extras ( open_type, k_rough, eps_wall )
    implicit none
    integer,  intent(in) :: open_type
    real(R8), intent(in) :: k_rough, eps_wall

    if (k_rough /= 0.0d0 .or. eps_wall /= 0.0d0) then
      write(*,*) '  [ERROR] BC ', open_type, ': wall roughness (', k_rough, ') / emissivity (', &
                 eps_wall, ') are not modelled by the MOSE_GPU wall closures.'
      error stop 'IO_BC (MOSE_open dialect): wall roughness/emissivity not supported'
    endif

  end subroutine Refuse_Unsupported_Wall_Extras

end module MOSE_IO_BC_Open
