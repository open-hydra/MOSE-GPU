!> @brief Single Rotating Frame (SRF) source terms for compressible Navier-Stokes.
!>
!> When activated, the full computational domain is assumed to rotate rigidly at
!> constant angular velocity omega [rad/s] about a fixed axis.  The velocity
!> stored in P(nu:nw) is the RELATIVE (rotating-frame) velocity
!>
!>   w = u_abs - Omega x r
!>
!> so that no changes are required to the convective or diffusive flux modules.
!>
!> --------------------------------------------------------------------------
!> GOVERNING EQUATIONS (source terms added to the RHS)
!> --------------------------------------------------------------------------
!>
!>   Continuity / species:  unchanged
!>
!>   Momentum (per unit volume):
!>     S_mom = -2 rho (Omega x w)         [Coriolis]
!>           + rho omega^2 r_perp          [centrifugal = -rho Omega x (Omega x r)]
!>     where  r_perp = r - (r.e) e  is the position vector perpendicular to the
!>     axis,  r = x_cell - x_origin.
!>
!>   Energy (per unit volume, E = cv T + 0.5 |w|^2):
!>     S_E = w . (rho omega^2 r_perp)
!>     [Coriolis does no work:  w . (-2 rho Omega x w) = 0 identically]
!>
!>   RANS equations: unchanged.
!>     For Spalart-Allmaras, activate the Spalart-Shur correction
!>     (turbulence-model = 'SA-RC') to remove spurious production due to
!>     solid-body rotation in the vorticity magnitude.
!>     SST / Wilcox k-omega use the strain-rate for production, which is
!>     already frame-invariant; no correction is needed.
!>
!> --------------------------------------------------------------------------
!> VISCOUS FLUXES
!> --------------------------------------------------------------------------
!>   tau = mu ( grad(w) + grad(w)^T - 2/3 div(w) I )
!>   Viscous energy flux = tau . w  ,  heat conduction q = -k grad(T)
!>   All quantities use the relative velocity w stored in P(nu:nw).
!>   Mod_Diffusive.f90 is ALREADY CORRECT; no changes required.
!>
!> --------------------------------------------------------------------------
!> BOUNDARY CONDITIONS
!> --------------------------------------------------------------------------
!>
!>   Inflow / outflow (type 4):
!>     T0, p0, alpha, beta and mach must be given as RELATIVE total conditions
!>     (rotating frame).  For inflow from a stationary region upstream:
!>
!>       w_inlet  = u_abs_inlet - Omega x r_inlet        (relative velocity)
!>       T0_rel   = T_static + |w_inlet|^2 / (2 cp)
!>       p0_rel   = p_static (T0_rel/T_static)^(gamma/(gamma-1))
!>       mach_rel = |w_inlet| / sqrt(gamma R T_static)
!>
!>     Outflow static-pressure (pamb) is frame-independent; no change needed.
!>
!>   Walls (type 5, 6, 8-14):
!>     No-slip condition w = 0 applies to the relative velocity.  The wall is
!>     stationary in the rotating frame, so the existing BC is correct.
!>
!>   Rotational periodicity (type 1 connection with pitch angle, or type 2):
!>     The relative velocity transforms between the two periodic faces by the
!>     same rotation matrix R(Dtheta) as the position vector.  For a sector with
!>     pitch angle Dtheta about the z-axis:
!>
!>       u_ghost = u_source cos(Dtheta) - v_source sin(Dtheta)
!>       v_ghost = u_source sin(Dtheta) + v_source cos(Dtheta)
!>       w_ghost = w_source   (axial)
!>
!>     This transformation must be embedded in the block connectivity (Mg tensor)
!>     when the mesh is generated.  The connection BC (type 1) then handles it
!>     automatically.  No code changes are needed provided the grid file encodes
!>     the correct rotational mapping between the two periodic faces.
!>
!> --------------------------------------------------------------------------
!> INITIALISATION (called from IO_Ini.f90 after reading input.ini)
!> --------------------------------------------------------------------------
!>   The input.ini section [MOSE-rotating-frame] is read entirely in
!>   IO_Ini.f90, which already owns the Finer dependency.  The raw fields
!>   of obj_rot are filled there, then RotatingFrame_Setup() is called to
!>   normalise the axis, compute Omega_vec, set the global flag, and print.
!>
!>   [MOSE-rotating-frame]
!>     active = .true.
!>     omega  = 1000.0        ! angular speed [rad/s]; sign follows right-hand rule
!>     axis   = 0.0 0.0 1.0   ! rotation-axis direction (need not be a unit vector)
!>     origin = 0.0 0.0 0.0   ! coordinates of one point on the axis [m]


!>   n-stationary-face  = 4
!>   stationary-face-1 = 1 4
!>   stationary-face-2 = 2 4
!>   stationary-face-3 = 3 4
!>   stationary-face-4 = 4 4
!>
module MOSE_Lib_RotatingFrame
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none


  type :: rotating_frame_type
    logical  :: active    = .false.
    real(R8) :: omega     = 0.0_R8                       !> Angular speed [rad/s]
    real(R8) :: axis(3)   = [ 0.0_R8, 0.0_R8, 1.0_R8 ] !> Unit rotation axis 
    real(R8) :: origin(3) = [ 0.0_R8, 0.0_R8, 0.0_R8 ] !> A point on the axis [m]
    real(R8) :: Omega_vec(3) = 0.0_R8                   !> omega * axis  [rad/s]
    !> Stationary wall patches in the lab frame, stored as (2, N) pairs:
    !>   stationary_face(1, k) = block index
    !>   stationary_face(2, k) = face direction (1-6)
    !> When ROTATING_FRAME is active, patches in this list receive
    !>   w_wall = -(Omega x r)  (casing/hub); all other type-5/6 walls use
    !>   w_wall = 0             (blades, default).
    !> If unallocated, all type-5/6 walls are treated as rotating (default).
    integer, allocatable :: stationary_face(:,:)

    contains
      procedure, pass(self) :: ini_setup => read_ini
  end type rotating_frame_type

  type(rotating_frame_type), public :: obj_rot
  private :: read_ini


contains


  ! ---------------------------------------------------------------------------
  !> @brief Finalise the rotating-frame configuration after IO_Ini has filled
  !>        the raw fields of obj_rot.
  !>
  !> Actions:
  !>   1. Normalise obj_rot%axis to a unit vector.
  !>   2. Compute obj_rot%Omega_vec = omega * axis.
  !>   3. Set the global flag ROTATING_FRAME.
  !>   4. Print a summary to stdout.

  subroutine read_ini(self, fini)
    use Finer, only: file_ini
    use MOSE_Parameters_m, only: codename
    use MOSE_Global_m, only: ROTATING_FRAME

    implicit none
    type(file_ini), intent(inout) :: fini
    class(rotating_frame_type), intent(inout) :: self
    ! Local
    integer            :: error, k, n_stat
    character(len=48)  :: key
    character(len=128) :: str
    real(R8)           :: norm

    obj_rot%active = .false.
    call fini%get( section_name = trim(codename)//'-rotating-frame', &
                   option_name  = 'active', val = obj_rot%active, error = error )
    if (error /= 0) obj_rot%active = .false.

    if ( obj_rot%active ) then
      call fini%get( section_name = trim(codename)//'-rotating-frame', &
                     option_name  = 'omega',  val = obj_rot%omega, error = error )
      str = ''
      call fini%get( section_name = trim(codename)//'-rotating-frame', &
                     option_name  = 'axis',   val = str, error = error )
      if (error == 0) read(str, *) obj_rot%axis
      str = ''
      call fini%get( section_name = trim(codename)//'-rotating-frame', &
                     option_name  = 'origin', val = str, error = error )
      if (error == 0) read(str, *) obj_rot%origin
    end if

    if ( .not. obj_rot%active ) then
      ROTATING_FRAME = .false.
      return
    end if

    norm = sqrt( sum(obj_rot%axis**2) )
    if ( norm < 1.0e-14_R8 ) error stop '[RotatingFrame] Rotation axis has zero length'
    obj_rot%axis      = obj_rot%axis / norm
    obj_rot%Omega_vec = obj_rot%omega * obj_rot%axis

    ROTATING_FRAME = .true.

    ! Optional list of stationary wall patches: (block, face-direction) pairs.
    ! Keys: n-stationary-face = N
    !       stationary-face-K = block_index face_direction   (both on the same line)
    n_stat = 0
    call fini%get( section_name = trim(codename)//'-rotating-frame', &
                   option_name  = 'n-stationary-face', val = n_stat, error = error )
    if ( error == 0 .and. n_stat > 0 ) then
      allocate( obj_rot%stationary_face(2, n_stat) )
      do k = 1, n_stat
        write(key, '(A,I0)') 'stationary-face-', k
        str = ''
        call fini%get( section_name = trim(codename)//'-rotating-frame', &
                       option_name  = trim(key), val = str, error = error )
        if ( error == 0 ) read(str, *) obj_rot%stationary_face(1,k), obj_rot%stationary_face(2,k)
      end do
    end if

    write(*,*)
    write(*,'(A)')         ' [Rotating Frame] SRF active'
    write(*,'(A,F14.4,A)') '   omega  = ', obj_rot%omega,   ' rad/s'
    write(*,'(A,3F11.5)')  '   axis   = ', obj_rot%axis
    write(*,'(A,3F11.5)')  '   origin = ', obj_rot%origin
    if ( allocated(obj_rot%stationary_face) ) then
      write(*,'(A,I0,A)') '   stationary patches = ', size(obj_rot%stationary_face,2), ':'
      do k = 1, size(obj_rot%stationary_face, 2)
        write(*,'(A,I0,A,I0)') '     block ', obj_rot%stationary_face(1,k), &
                                '  face-dir ', obj_rot%stationary_face(2,k)
      end do
    else
      write(*,'(A)') '   all type-5/6 walls rotate with the frame'
    end if

  end subroutine read_ini


  ! ---------------------------------------------------------------------------
  !> @brief Add Coriolis and centrifugal source terms to the momentum and
  !>        energy residuals for all blocks in the domain.
  subroutine RotatingFrame_Source_Terms ( domain )
    use MOSE_Types_m
    use MOSE_Mod_MPI, only: is_local_block

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: b

    do b = 1, domain%nb
      if (.not. is_local_block(b)) cycle
      call RF_Source_Blk( domain%blk(b)%P,    &
                          domain%blk(b)%R,    &
                          domain%blk(b)%node, &
                          domain%blk(b)%vol,  &
                          domain%blk(b)%dim   )
    end do

  end subroutine RotatingFrame_Source_Terms


  ! ---------------------------------------------------------------------------
  !> @brief Block-level kernel: loops over interior cells and accumulates the
  !>        rotating-frame source terms into the residual array R.
  !>
  !> Residual sign convention (same as throughout MOSE):
  !>   V * dU/dt = -R(U)   with  R = CONV - DIFF - SOURCE
  !> => to add a physical source S, we do  R -= S * vol
  subroutine RF_Source_Blk ( Prim, Res, node, vol, n )
    use MOSE_General_Types_m
    use MOSE_Global_m

    implicit none
    integer,  intent(in) :: n(3)
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: Prim
    real(R8), dimension(nprim, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(inout) :: Res
    type(MOSE_vector_3D_type), dimension(0:n(1), 0:n(2), 0:n(3)),  intent(in) :: node
    real(R8), dimension(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc), intent(in) :: vol
    ! Local
    integer  :: i, j, k
    real(R8) :: rho, w(3)
    real(R8) :: xc, yc, zc, r_vec(3), r_proj, r_perp(3)
    real(R8) :: Cor(3), Cen(3)
    ! Cache module-level data in local scalars to avoid repeated indirect lookups
    real(R8) :: Omega(3), ax(3), ox(3), omega2

    Omega  = obj_rot%Omega_vec
    ax     = obj_rot%axis
    ox     = obj_rot%origin
    omega2 = obj_rot%omega**2

    !$omp do collapse(3) &
    !$omp private( i, j, k, rho, w, xc, yc, zc, r_vec, r_proj, r_perp, Cor, Cen )
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)

      rho = sum( Prim(1:nsc, i, j, k) )
      w   = Prim(nu:nw, i, j, k)            ! relative velocity in rotating frame

      ! ---- Cell centre (average of 8 hexahedral corner nodes) ----
      ! Node (i,j,k) in [0..ni,0..nj,0..nk]:
      ! cell (i,j,k) occupies the hexahedron with corners at i-1..i, j-1..j, k-1..k.
      xc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(1) + node(i,j-1,k-1)%c(1) &
                       + node(i-1,j,  k-1)%c(1) + node(i,j,  k-1)%c(1) &
                       + node(i-1,j-1,k  )%c(1) + node(i,j-1,k  )%c(1) &
                       + node(i-1,j,  k  )%c(1) + node(i,j,  k  )%c(1) )
      yc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(2) + node(i,j-1,k-1)%c(2) &
                       + node(i-1,j,  k-1)%c(2) + node(i,j,  k-1)%c(2) &
                       + node(i-1,j-1,k  )%c(2) + node(i,j-1,k  )%c(2) &
                       + node(i-1,j,  k  )%c(2) + node(i,j,  k  )%c(2) )
      zc = 0.125_R8 * ( node(i-1,j-1,k-1)%c(3) + node(i,j-1,k-1)%c(3) &
                       + node(i-1,j,  k-1)%c(3) + node(i,j,  k-1)%c(3) &
                       + node(i-1,j-1,k  )%c(3) + node(i,j-1,k  )%c(3) &
                       + node(i-1,j,  k  )%c(3) + node(i,j,  k  )%c(3) )

      ! ---- Position vector from axis origin to cell centre ----
      r_vec = [ xc - ox(1), yc - ox(2), zc - ox(3) ]

      ! ---- Component of r perpendicular to the rotation axis ----
      !   r_perp = r_vec - (r_vec . axis) * axis
      r_proj = dot_product( r_vec, ax )
      r_perp = r_vec - r_proj * ax

      ! ---- Coriolis force per unit volume:  -2 rho (Omega x w) ----
      Cor(1) = -2.0_R8 * rho * ( Omega(2)*w(3) - Omega(3)*w(2) )
      Cor(2) = -2.0_R8 * rho * ( Omega(3)*w(1) - Omega(1)*w(3) )
      Cor(3) = -2.0_R8 * rho * ( Omega(1)*w(2) - Omega(2)*w(1) )

      ! ---- Centrifugal force per unit volume:  rho omega^2 r_perp ----
      !   (= -rho Omega x (Omega x r)  via the vector triple-product identity)
      Cen = rho * omega2 * r_perp

      ! ---- Accumulate into residuals:  R -= source * vol ----
      !   Momentum equations
      Res(nu, i, j, k) = Res(nu, i, j, k) - ( Cor(1) + Cen(1) ) * vol(i,j,k)
      Res(nv, i, j, k) = Res(nv, i, j, k) - ( Cor(2) + Cen(2) ) * vol(i,j,k)
      Res(nw, i, j, k) = Res(nw, i, j, k) - ( Cor(3) + Cen(3) ) * vol(i,j,k)

      !   Energy equation:  only centrifugal does work on w
      !   w . Coriolis = w . (-2 rho Omega x w) = 0  identically
      
      Res(np, i, j, k) = Res(np, i, j, k) - dot_product( w, Cen ) * vol(i,j,k)

    end do ; end do ; end do
    !$omp end do

  end subroutine RF_Source_Blk


  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> @brief Compute the centre of a block face as the average of its 4 corner nodes.
  !>

  subroutine RF_Face_Center(node, Im, Jm, Km, Fm, r_face)
    use MOSE_General_Types_m, only: MOSE_vector_3D_type

    implicit none
    type(MOSE_vector_3D_type), intent(in) :: node(0:,0:,0:)
    integer,  intent(in)  :: Im, Jm, Km, Fm
    real(R8), intent(out) :: r_face(3)
    ! Local
    integer :: Fi, Fj, Fk, d

    ! Face node index (same logic as Face_Index in Lib_BC.f90)
    select case (Fm)
    case(1,2)
      Fi = Im - mod(Fm,2) ; Fj = Jm ; Fk = Km
    case(3,4)
      Fi = Im ; Fj = Jm - mod(Fm,2) ; Fk = Km
    case(5,6)
      Fi = Im ; Fj = Jm ; Fk = Km - mod(Fm,2)
    end select

    ! 4-node average for each Cartesian component
    do d = 1, 3
      select case (Fm)
      case(1,2)   ! i-face: nodes span (Fi, Fj-1:Fj, Fk-1:Fk)
        r_face(d) = 0.25_R8 * ( node(Fi,Fj-1,Fk-1)%c(d) + node(Fi,Fj,Fk-1)%c(d) &
                               + node(Fi,Fj-1,Fk  )%c(d) + node(Fi,Fj,Fk  )%c(d) )
      case(3,4)   ! j-face: nodes span (Fi-1:Fi, Fj, Fk-1:Fk)
        r_face(d) = 0.25_R8 * ( node(Fi-1,Fj,Fk-1)%c(d) + node(Fi,Fj,Fk-1)%c(d) &
                               + node(Fi-1,Fj,Fk  )%c(d) + node(Fi,Fj,Fk  )%c(d) )
      case(5,6)   ! k-face: nodes span (Fi-1:Fi, Fj-1:Fj, Fk)
        r_face(d) = 0.25_R8 * ( node(Fi-1,Fj-1,Fk)%c(d) + node(Fi,Fj-1,Fk)%c(d) &
                               + node(Fi-1,Fj  ,Fk)%c(d) + node(Fi,Fj  ,Fk)%c(d) )
      end select
    end do

  end subroutine RF_Face_Center


  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> @brief Return .true. if the wall patch (block b, face direction f) is
  !>        stationary in the lab frame (casing/hub).
  !>
  !> Called from Mod_BC.f90 for type-5/6 wall BCs when ROTATING_FRAME is active.
  !> Patches listed in obj_rot%stationary_face receive w_wall = -(Omega x r);
  !> all other patches (blades) use the standard w_wall = 0.
  logical function RF_Is_Stationary_Face(b, f)
    implicit none
    integer, intent(in) :: b, f
    integer :: k

    RF_Is_Stationary_Face = .false.
    if ( .not. allocated(obj_rot%stationary_face) ) return
    do k = 1, size(obj_rot%stationary_face, 2)
      if ( obj_rot%stationary_face(1,k) == b .and. &
           obj_rot%stationary_face(2,k) == f ) then
        RF_Is_Stationary_Face = .true.
        return
      end if
    end do

  end function RF_Is_Stationary_Face


  ! ---------------------------------------------------------------------------
  !> @brief Estimate the pressure torque on the rotating wall surfaces.
  !>
  !> Loops over every BC entry that is a wall face (types 5,6,8-10,12,13) and
  !> is NOT listed in obj_rot%stationary_face (i.e., the rotating blade/helix
  !> surface).  The pressure force on each face is
  !>
  !>   F_wall = p * (outward-from-fluid unit normal) * area
  !>
  !> and the axial torque contribution is
  !>
  !>   dT = axis . ( (r_face - origin) x F_wall )
  !>
  !> Only pressure contributions are included; viscous shear is neglected
  !> (adequate for a first estimate on a turbine/propeller profile).
  !>
  !> Results are written to OUTPUT/torque.dat at each call, at the same
  !> frequency as the residual history (controlled by the caller).
  !>
  !> @param[in]  domain   MOSE domain (read-only)
  !> @param[in]  iter     Current global iteration (for the output column)
  !> @param[in]  time     Current solution time [s or non-dim]
  subroutine RF_Torque_Walls(domain, iter, time)
    use MOSE_Types_m
    use MOSE_Global_m, only: np, Is_Wall

    implicit none
    type(MOSE_domain_type), intent(in) :: domain
    integer,  intent(in) :: iter
    real(R8), intent(in) :: time
    ! Local
    integer,  save :: unit_torque = -1
    logical,  save :: first_call  = .true.
    integer  :: i, Bm, Im, Jm, Km, Fm, Dir, Fi, Fj, Fk
    real(R8) :: r_face(3), r_rel(3), F_wall(3), torque_vec(3)
    real(R8) :: T_axial, p_face, sign_fm
    real(R8) :: N_face(3), A_face

    ! Open the output file on the first call
    if (first_call) then
      open(newunit=unit_torque, file='OUTPUT/torque.dat', &
           status='replace', form='formatted')
      write(unit_torque, '(A)') &
        '#        iter              time   torque_axial (pressure, non-dim)' // &
        '         power (T*omega, non-dim)'
      first_call = .false.
    end if

    T_axial = 0.0_R8

    do i = 1, domain%nbound

      ! Skip non-wall BC types
      if ( .not. Is_Wall(domain%bc(i)%type) ) cycle

      Bm = domain%bc(i)%b
      Fm = domain%bc(i)%f

      ! Skip faces that are stationary in the lab frame (casing / hub shroud)
      if ( RF_Is_Stationary_Face(Bm, Fm) ) cycle

      Im = domain%bc(i)%i
      Jm = domain%bc(i)%j
      Km = domain%bc(i)%k

      ! Direction index: 1 for i-faces (Fm=1,2), 2 for j-faces (3,4), 3 for k-faces (5,6)
      Dir = (Fm + 1) / 2

      ! Index of the boundary face in the precomputed metrics arrays
      select case (Fm)
      case (1,2) ; Fi = Im - mod(Fm,2) ; Fj = Jm             ; Fk = Km
      case (3,4) ; Fi = Im             ; Fj = Jm - mod(Fm,2) ; Fk = Km
      case (5,6) ; Fi = Im             ; Fj = Jm             ; Fk = Km - mod(Fm,2)
      end select

      ! Unit normal and area from precomputed face metrics.
      ! dir(d)%f%N always points in the direction of increasing index d.
      N_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj,Fk)%N
      A_face = domain%blk(Bm)%dir(Dir)%f(Fi,Fj,Fk)%A

      ! Outward-from-fluid sign:  -1 for min faces (Fm odd), +1 for max faces (Fm even)
      sign_fm = real(1 - 2 * mod(Fm, 2), R8)

      ! Pressure at the interior cell adjacent to the wall (zero-order extrapolation)
      p_face = domain%blk(Bm)%P(np, Im, Jm, Km)

      ! Face centre position (4-node average of corner nodes)
      call RF_Face_Center(domain%blk(Bm)%node, Im, Jm, Km, Fm, r_face)

      ! Pressure force on the wall from the fluid:
      !   F_wall = p * (outward-from-fluid normal) * A = p * sign_fm * N * A
      F_wall = p_face * sign_fm * N_face * A_face

      ! Torque contribution about the rotation axis origin:
      !   dT_axial = axis . ( (r_face - origin) x F_wall )
      r_rel        = r_face - obj_rot%origin
      torque_vec(1) = r_rel(2)*F_wall(3) - r_rel(3)*F_wall(2)
      torque_vec(2) = r_rel(3)*F_wall(1) - r_rel(1)*F_wall(3)
      torque_vec(3) = r_rel(1)*F_wall(2) - r_rel(2)*F_wall(1)
      T_axial = T_axial + dot_product(obj_rot%axis, torque_vec)

    end do

    write(unit_torque, '(I9,3E20.10)') iter, time, T_axial, T_axial * obj_rot%omega
    flush(unit_torque)

  end subroutine RF_Torque_Walls


end module MOSE_Lib_RotatingFrame
