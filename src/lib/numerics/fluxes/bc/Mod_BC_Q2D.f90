module MOSE_Mod_BC_Q2D
  use iso_fortran_env, only: R8 => real64
  use MOSE_Types_m
  use MOSE_Global_m

  implicit none

  private
  public :: Setup_Q2D_BC_Data, Q2D_Ghost_Values, q2d_any_active

  !> The faces of type 667 (MOSE_open 410) on the device (2026-10-29): one slot per local record, in the order of the device
  !> boundary arrays (bc_perm). q2d_val(:,k) is the state of the file at the stage time (Q2D_Ghost_Values, every ghost fill),
  !> q2d_has(s,k) = 1 where the file gives the variable s; the ghost kernels (Lib_Extrapolation Ghost_Q2D_Layers/_Edges)
  !> take the boundary cell for the rest, as MOSE_open Ghost_Q2D_Mapped.
  integer, public :: q2d_n = -1                       ! -1: not set up yet
  integer, allocatable, public :: q2d_idx(:), q2d_blk(:), q2d_has(:,:), q2d_dim(:,:)
  integer, allocatable, public :: q2d_all(:)          ! the slots 1..q2d_n (the ghost kernels take a list of slots)
  real(R8), allocatable, public :: q2d_val(:,:)
  integer, allocatable :: q2d_js(:), q2d_grp(:), q2d_ref(:)   ! host: bc(:) index, file group, a record of each group
  integer :: q2d_ng = 0
  integer, allocatable :: q2d_blist(:)   ! host: the distinct blocks of the slots (q2d_any_active)
  ! host: the operands of the time interpolation of every (variable, slot) for the bracket of its group at the last call
  ! (Q2D_Ghost_Values), re-gathered from bc(:)%q2d_map only when a group's bracket or mode changes
  integer,  allocatable :: q2d_gi_c(:), q2d_gm_c(:)
  real(R8), allocatable :: q2d_va(:,:), q2d_vb(:,:), q2d_ta(:,:), q2d_tb(:,:)
  !> The time of the boundary data: the stage time of MOSE_open (Mod_Explicit: t^n + (i_strang - 1 + theta) * dt / N_strang
  !> in a time-accurate step), set by Mod_Explicit around the stages; outside them (setup, steady runs) domain%time.
  real(R8), public :: q2d_time = 0d0
  logical, public  :: q2d_time_set = .false.

contains


  ! ---------------------------------------------------------------------------
  ! Setup Q2D mapped boundary condition data (BC type 667).
  ! Reads the ATLAS Tecplot file(s) specified in bc.txt (one per face),
  ! parses all solution times, maps variables by name to MOSE primitives,
  ! and populates bc(i)%q2d_map for each type 667 boundary cell.
  ! ---------------------------------------------------------------------------
  subroutine Setup_Q2D_BC_Data ( domain )

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, nq2d_count, nfiles, ifile
    character(len=256), allocatable :: file_list(:)
    logical :: found

    ! Count type 667 BCs
    nq2d_count = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type == 667) nq2d_count = nq2d_count + 1
    enddo
    if (nq2d_count == 0) return

    ! Collect unique filenames from bc(:)%q2d_file
    allocate(file_list(nq2d_count))
    nfiles = 0
    do i = 1, domain%nbound
      if (domain%bc(i)%type /= 667) cycle
      if (len_trim(domain%bc(i)%q2d_file) == 0) then
        write(*,*) '[ERROR] BC type 667 entry', i, 'has no Q2D file specified in bc.txt'
        error stop
      endif
      found = .false.
      do ifile = 1, nfiles
        if (trim(domain%bc(i)%q2d_file) == trim(file_list(ifile))) then
          found = .true.; exit
        endif
      enddo
      if (.not. found) then
        nfiles = nfiles + 1
        file_list(nfiles) = domain%bc(i)%q2d_file
      endif
    enddo

    ! Process each unique Q2D file
    do ifile = 1, nfiles
      write(*,'(A,A)') '  Q2D file: ', trim(file_list(ifile))
      call Read_Q2D_Tecplot ( domain, trim(file_list(ifile)) )
    enddo

    deallocate(file_list)

  end subroutine Setup_Q2D_BC_Data


  ! ---------------------------------------------------------------------------
  ! Parse the ATLAS Tecplot ASCII file and distribute data to type 667 BCs.
  !
  ! Two-pass reading:
  !   Pass 1 — scan VARIABLES line + all ZONE headers (no data read)
  !   Pass 2 — rewind, read zone data and distribute to matching BCs
  ! ---------------------------------------------------------------------------
  subroutine Read_Q2D_Tecplot ( domain, filename )
    use strings, only: lowercase

    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    character(len=*), intent(in) :: filename
    ! Local
    integer :: u, ios, iz, nvar_file
    integer :: I_nodes, J_nodes, n1, n2, n_nodes, I_max, J_max
    integer :: b_zone, t_zone, iv, ic, i1, i2, ibc, it
    integer :: ip, jp, ncount
    character(len=1024) :: line
    character(len=32) :: tmp_names(100)
    integer, allocatable :: varmap(:)
    real(R8), allocatable :: node_buf(:), cell_buf(:,:,:)
    integer :: nvar_total
    integer :: nzones_total
    integer, allocatable :: zone_block(:), zone_tidx(:)
    real(R8), allocatable :: zone_time(:), unique_times_tmp(:)
    integer, allocatable :: zone_I(:), zone_J(:)
    integer :: ntimes
    real(R8), allocatable :: unique_times(:)
    integer :: bc_i1, bc_i2
    logical :: found
    integer :: nlist, il
    integer, allocatable :: list(:)

    ! --- Open file ---
    open(newunit=u, file=filename, status='old', iostat=ios, action='read')
    if (ios /= 0) then
      write(*,*) '[ERROR] Cannot open Q2D file: ', trim(filename)
      error stop
    endif

    ! --- Pass 1: VARIABLES line + all ZONE headers ---

    ! Read VARIABLES line — parse quoted names: "x" "y" "z" "var1" ...
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) then
        write(*,*) '[ERROR] Could not find VARIABLES line in Q2D file'
        error stop
      endif
      if (index(lowercase(line), 'variables') > 0) exit
    enddo
    ! Extract quoted variable names (single pass)
    ncount = 0; ip = 1
    do while (ip <= len_trim(line))
      if (line(ip:ip) == '"') then
        jp = index(line(ip+1:), '"')
        if (jp > 0) then
          ncount = ncount + 1
          tmp_names(ncount) = line(ip+1:ip+jp-1)
          ip = ip + jp + 1
        else
          exit
        endif
      else
        ip = ip + 1
      endif
    enddo
    nvar_total = ncount
    nvar_file = nvar_total - 3  ! subtract x, y, z

    if (nvar_file < 1) then
      write(*,*) '[ERROR] Q2D file has no solution variables (only coordinates)'
      error stop
    endif

    ! Build variable name mapping
    allocate(varmap(nvar_file))
    call build_varmap(tmp_names(4:nvar_total), nvar_file, varmap)

    ! Continue scanning from current position (no rewind): count zones
    nzones_total = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(lowercase(line), 'zone') > 0) nzones_total = nzones_total + 1
    enddo

    if (nzones_total == 0) then
      write(*,*) '[ERROR] Q2D file has no zones'
      error stop
    endif

    allocate(zone_block(nzones_total), zone_tidx(nzones_total))
    allocate(zone_time(nzones_total), zone_I(nzones_total), zone_J(nzones_total))
    allocate(unique_times_tmp(nzones_total))

    ! Rewind once: parse zone headers
    rewind(u)
    iz = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(lowercase(line), 'zone') > 0) then
        iz = iz + 1
        call parse_zone_header(line, zone_block(iz), &
                               zone_I(iz), zone_J(iz), zone_time(iz))
      endif
    enddo

    ! Extract unique times and build zone → time index mapping (single pass)
    ntimes = 0
    do iz = 1, nzones_total
      found = .false.
      do it = 1, ntimes
        if (abs(zone_time(iz) - unique_times_tmp(it)) < 1.0d-12) then
          found = .true.; zone_tidx(iz) = it; exit
        endif
      enddo
      if (.not. found) then
        ntimes = ntimes + 1
        unique_times_tmp(ntimes) = zone_time(iz)
        zone_tidx(iz) = ntimes
      endif
    enddo
    allocate(unique_times(ntimes))
    unique_times = unique_times_tmp(1:ntimes)
    deallocate(unique_times_tmp)

    write(*,'(A,I0,A,I0,A)') '  Q2D data: ', ntimes, ' time steps, ', nvar_file, ' variables'

    ! --- Allocate q2d_map for each type 667 BC ---
    do ibc = 1, domain%nbound
      if (domain%bc(ibc)%type == 667 .and. trim(domain%bc(ibc)%q2d_file) == trim(filename)) then
        allocate(domain%bc(ibc)%q2d_map(nprim))
        do iv = 1, nvar_file
          if (varmap(iv) > 0) then
            ic = varmap(iv)
            domain%bc(ibc)%q2d_map(ic)%exists = .true.
            domain%bc(ibc)%q2d_map(ic)%n = ntimes
            allocate(domain%bc(ibc)%q2d_map(ic)%time(ntimes))
            allocate(domain%bc(ibc)%q2d_map(ic)%var(ntimes))
            domain%bc(ibc)%q2d_map(ic)%time = unique_times
            domain%bc(ibc)%q2d_map(ic)%var = 0.0_R8
            domain%bc(ibc)%q2d_map(ic)%periodic = domain%bc(ibc)%q2d_periodic
          endif
        enddo
        ! Extra species not in file: set to zero for all times
        do ic = 1, nsc
          if (.not. domain%bc(ibc)%q2d_map(ic)%exists) then
            domain%bc(ibc)%q2d_map(ic)%exists = .true.
            domain%bc(ibc)%q2d_map(ic)%n = ntimes
            allocate(domain%bc(ibc)%q2d_map(ic)%time(ntimes))
            allocate(domain%bc(ibc)%q2d_map(ic)%var(ntimes))
            domain%bc(ibc)%q2d_map(ic)%time = unique_times
            domain%bc(ibc)%q2d_map(ic)%var = 0.0_R8
            domain%bc(ibc)%q2d_map(ic)%periodic = domain%bc(ibc)%q2d_periodic
          endif
        enddo
      endif
    enddo

    ! --- Pass 2: read zone data and distribute to BCs ---

    ! The records of this file, listed once: every zone is distributed over them only (not over every record of the
    ! domain: 6020 zones x 1.1 M records on the cono wave mesh)
    nlist = 0
    do ibc = 1, domain%nbound
      if (domain%bc(ibc)%type /= 667) cycle
      if (trim(domain%bc(ibc)%q2d_file) /= trim(filename)) cycle
      nlist = nlist + 1
    enddo
    allocate(list(nlist))
    nlist = 0
    do ibc = 1, domain%nbound
      if (domain%bc(ibc)%type /= 667) cycle
      if (trim(domain%bc(ibc)%q2d_file) /= trim(filename)) cycle
      nlist = nlist + 1
      list(nlist) = ibc
    enddo

    ! Allocate buffers once using max dimensions
    I_max = maxval(zone_I(1:nzones_total))
    J_max = maxval(zone_J(1:nzones_total))
    allocate(node_buf(I_max * J_max))
    allocate(cell_buf(nvar_file, I_max - 1, J_max - 1))

    rewind(u)
    ! Skip to first ZONE line
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) then
        write(*,*) '[ERROR] Q2D file: no ZONE found in data pass'
        error stop
      endif
      if (index(lowercase(line), 'zone') > 0) exit
    enddo

    do iz = 1, nzones_total
      ! First zone header already read above; subsequent ones read at end of loop

      I_nodes = zone_I(iz)
      J_nodes = zone_J(iz)
      n1 = I_nodes - 1
      n2 = J_nodes - 1
      n_nodes = I_nodes * J_nodes
      b_zone = zone_block(iz)
      t_zone = zone_tidx(iz)

      ! Read and skip node coordinates (3 * n_nodes values)
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! x
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! y
      read(u, *, iostat=ios) node_buf(1:n_nodes)  ! z

      ! Read cell-centered variables: nvar_file blocks of n_cells values
      do iv = 1, nvar_file
        read(u, *, iostat=ios) ((cell_buf(iv, i1, i2), i1=1,n1), i2=1,n2)
      enddo

      ! Distribute to matching type 667 BCs
      do il = 1, nlist
        ibc = list(il)
        if (domain%bc(ibc)%b /= b_zone) cycle

        ! Face-local indices: (i,j,k,f) -> (i1,i2)
        select case (domain%bc(ibc)%f)
          case (1,2); bc_i1 = domain%bc(ibc)%j; bc_i2 = domain%bc(ibc)%k
          case (3,4); bc_i1 = domain%bc(ibc)%i; bc_i2 = domain%bc(ibc)%k
          case default; bc_i1 = domain%bc(ibc)%i; bc_i2 = domain%bc(ibc)%j
        end select

        ! Check bounds
        if (bc_i1 < 1 .or. bc_i1 > n1 .or. bc_i2 < 1 .or. bc_i2 > n2) then
          write(*,'(A,4I6)') ' [WARNING] Q2D BC index out of range: bc=', ibc, &
            bc_i1, bc_i2, b_zone
          cycle
        endif

        ! Map file variables to MOSE primitives
        do iv = 1, nvar_file
          if (varmap(iv) > 0) then
            domain%bc(ibc)%q2d_map(varmap(iv))%var(t_zone) = cell_buf(iv, bc_i1, bc_i2)
          endif
        enddo
      enddo

      ! Read next zone header (if not the last zone)
      if (iz < nzones_total) read(u, '(A)', iostat=ios) line
    enddo

    close(u)
    deallocate(node_buf, cell_buf, list)
    deallocate(varmap, zone_block, zone_tidx, zone_time, zone_I, zone_J, unique_times)

  end subroutine Read_Q2D_Tecplot


  ! ---------------------------------------------------------------------------
  ! Build mapping from file variable index to MOSE primitive index.
  ! Matching:
  !   (a) "rho(N)" pattern — MOSE canonical (IO_Solution, Mod_Probes)
  !   (b) "u","v","w","p"  — case-insensitive
  ! ---------------------------------------------------------------------------
  subroutine build_varmap(names, nvar, varmap)
    use strings, only: lowercase
    use IR_Precision, only: str

    implicit none
    integer, intent(in) :: nvar
    character(len=*), intent(in) :: names(nvar)
    integer, intent(out) :: varmap(nvar)
    ! Local
    integer :: iv, s
    character(len=64) :: lname

    varmap = 0

    do iv = 1, nvar
      lname = lowercase(trim(names(iv)))

      ! (a) "rho(N)" pattern (MOSE canonical: IO_Solution, Mod_Probes)
      do s = 1, nsc
        if (trim(lname) == 'rho('//trim(str(.true.,s))//')') then
          varmap(iv) = s
          exit
        endif
      enddo
      if (varmap(iv) /= 0) cycle

      ! (b) velocity and pressure (case-insensitive)
      select case (trim(lname))
        case ('u');  varmap(iv) = nu
        case ('v');  varmap(iv) = nv
        case ('w');  varmap(iv) = nw
        case ('p');  varmap(iv) = np
      end select
      if (varmap(iv) /= 0) cycle

      ! (c) the turbulence variables, density-weighted as the solution files hold them, only when the run carries them
      !     (MOSE_open + rans410_main_2026-10-26.patch): "rhoK" "rhoOmega" (or "kappa" "omega") with nrans = 2, "rhoNuTilde"
      !     (or "mi_tilde") with nrans = 1; otherwise the columns stay unmapped and the ghost takes the cell's per unit mass
      select case (trim(lname))
        case ('rhok', 'kappa');           if (nrans == 2) varmap(iv) = nt
        case ('rhoomega', 'omega');       if (nrans == 2) varmap(iv) = nt + 1
        case ('rhonutilde', 'mi_tilde');  if (nrans == 1) varmap(iv) = nt
      end select
    enddo

    ! Log the mapping
    write(*,'(A)') '  Q2D variable mapping:'
    do iv = 1, nvar
      if (varmap(iv) > 0) then
        write(*,'(A,A,A,I3,A)') '    "', trim(names(iv)), '" -> P(', varmap(iv), ')'
      else
        write(*,'(A,A,A)')    '    "', trim(names(iv)), '" -> (unmapped)'
      endif
    enddo

  end subroutine build_varmap


  ! ---------------------------------------------------------------------------
  ! Parse a Tecplot ZONE header line to extract block number,
  ! I (node count), J (node count), and SOLUTIONTIME.
  ! Zone name format: T="B{block}_T{tidx}"
  ! ---------------------------------------------------------------------------
  subroutine parse_zone_header(line, b_num, I_nodes, J_nodes, sol_time)
    use strings, only: parse, lowercase
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: b_num, I_nodes, J_nodes
    real(R8), intent(out) :: sol_time
    ! Local
    integer :: i, pos, pos2, ios
    character(len=100) :: args(20), subargs(4)
    character(len=256) :: zone_title

    b_num = 1; I_nodes = 2; J_nodes = 2; sol_time = 0.0_R8

    ! Split line by commas
    call parse(line, ',', args)

    do i = 1, size(args)
      ! Zone title: T="B{b}_T{t}"
      if (index(args(i), 'T="') > 0 .or. index(args(i), 't="') > 0) then
        pos = index(args(i), '"') + 1
        pos2 = index(args(i)(pos:), '"')
        if (pos2 > 0) then
          zone_title = args(i)(pos:pos+pos2-2)
          ! Parse "B{num}_T{num}"
          pos = index(zone_title, 'B')
          if (pos > 0) then
            pos2 = index(zone_title(pos+1:), '_')
            if (pos2 > 0) read(zone_title(pos+1:pos+pos2-1), *, iostat=ios) b_num
          endif
        endif
      endif

      ! I= (avoid matching STRANDID=)
      if (index(args(i), 'I=') > 0 .and. index(args(i), 'STRANDID') == 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) I_nodes
      endif

      ! J=
      if (index(args(i), 'J=') > 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) J_nodes
      endif

      ! SOLUTIONTIME=
      if (index(lowercase(args(i)), 'solutiontime') > 0) then
        call parse(args(i), '=', subargs)
        read(subargs(2), *, iostat=ios) sol_time
      endif
    enddo

  end subroutine parse_zone_header


  !> A block with faces of type 667 takes part in this ghost fill (always outside a multirate stage; inside it, a block with
  !> cells of the stage's level, Lib_Multirate_Base mr_block_active). True before the slots are built.
  logical function q2d_any_active ()
    use MOSE_Lib_Multirate_Base, only: mr_block_active
    integer :: k, m
    q2d_any_active = .true.
    if ( q2d_n <= 0 ) return
    if ( .not. allocated(q2d_blist) ) then
      allocate ( q2d_blist(0) )
      do k = 1, q2d_n
        if ( .not. any(q2d_blist == q2d_blk(k)) ) q2d_blist = [ q2d_blist, q2d_blk(k) ]
      end do
    end if
    q2d_any_active = .false.
    do m = 1, size(q2d_blist)
      if ( mr_block_active(q2d_blist(m)) ) then
        q2d_any_active = .true.
        return
      end if
    end do
  end function q2d_any_active


  !> The faces of type 667 at the time t, into q2d_val on the device: MOSE_open Ghost_Q2D_Mapped, every variable that the
  !> file gives interpolated as time_bc_type%update (periodic: wrapped into the period; otherwise held at the first and the
  !> last time of the table). Every record of a file holds the times of the file, so the bracket of t is searched once per
  !> file (group) and the interpolation of each record is the line of update, with the same operands. The first call
  !> builds the slots (Q2D_Setup_Slots); with no face of type 667 it returns at once. At the time of the last call q2d_val
  !> (host and device) already holds the state: the call returns (the next RK2 step starts at the time of the last stage, a
  !> multirate macro step has many stages at the same time).
  subroutine Q2D_Ghost_Values ( domain, t )
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    real(R8), intent(in) :: t
    ! Local
    integer :: k, g, s, js, i, n, mode
    real(R8) :: t_, period
    integer, allocatable :: gi(:), gmode(:)
    real(R8), allocatable :: gt(:)
    real(R8), save :: t_last = 0d0
    logical,  save :: have_last = .false.

    if ( q2d_n < 0 ) call Q2D_Setup_Slots ( domain )
    if ( q2d_n == 0 ) return
    if ( have_last .and. t == t_last ) return
    have_last = .true.; t_last = t

    allocate ( gi(q2d_ng), gmode(q2d_ng), gt(q2d_ng) )
    do g = 1, q2d_ng
      associate ( m => domain%bc(q2d_ref(g))%q2d_map(1) )   ! species 1: present in every record (Read_Q2D_Tecplot)
        n = m%n
        t_ = t
        gi(g) = 0
        if ( m%periodic ) then                     ! MOSE_open: update(time), which wraps the time into the period
          period = m%time(n) - m%time(1)
          if ( period > 0d0 ) t_ = m%time(1) + modulo ( t_ - m%time(1), period )
          mode = 1
        else if ( t <= m%time(1) ) then
          mode = 2                                 ! var(1)
        else if ( t >= m%time(n) ) then
          mode = 3                                 ! var(n)
        else
          mode = 1
        end if
        if ( mode == 1 ) then
          if ( t_ == m%time(1) ) then
            mode = 2
          else
            do i = 2, n
              if ( t_ > m%time(i-1) .and. t_ <= m%time(i) ) then
                gi(g) = i
                exit
              end if
            end do
            if ( gi(g) == 0 ) error stop 'Q2D (667): the time falls outside the table of the file'
          end if
        end if
        gmode(g) = mode
        gt(g) = t_
      end associate
    end do

    ! the operands of the slots whose group changed bracket (the nozzle 410 files: samples 123 ns apart, multirate stages
    ! ~0.1 ns: one gather every ~1000 calls), then the interpolation on the dense copies -- the expression of MOSE_open on
    ! the same operands
    if ( .not. allocated(q2d_va) ) then
      allocate ( q2d_va(nprim,q2d_n), q2d_vb(nprim,q2d_n), q2d_ta(nprim,q2d_n), q2d_tb(nprim,q2d_n), &
                 q2d_gi_c(q2d_ng), q2d_gm_c(q2d_ng) )
      q2d_va = 0d0; q2d_vb = 0d0; q2d_ta = 0d0; q2d_tb = 1d0; q2d_gi_c = -1; q2d_gm_c = -1
    end if
    do k = 1, q2d_n
      js = q2d_js(k)
      g  = q2d_grp(k)
      if ( gi(g) == q2d_gi_c(g) .and. gmode(g) == q2d_gm_c(g) ) cycle
      do s = 1, nprim
        if ( q2d_has(s,k) == 0 ) cycle
        associate ( m => domain%bc(js)%q2d_map(s) )
          select case ( gmode(g) )
            case (1)
              i = gi(g)
              q2d_va(s,k) = m%var(i-1); q2d_vb(s,k) = m%var(i); q2d_ta(s,k) = m%time(i-1); q2d_tb(s,k) = m%time(i)
            case (2)
              q2d_va(s,k) = m%var(1)
            case (3)
              q2d_va(s,k) = m%var(m%n)
          end select
        end associate
      end do
    end do
    q2d_gi_c = gi; q2d_gm_c = gmode
    do k = 1, q2d_n
      g  = q2d_grp(k)
      do s = 1, nprim
        if ( q2d_has(s,k) == 0 ) cycle
        if ( gmode(g) == 1 ) then
          q2d_val(s,k) = (q2d_vb(s,k)-q2d_va(s,k))/(q2d_tb(s,k)-q2d_ta(s,k))*(gt(g)-q2d_ta(s,k))+q2d_va(s,k)
        else
          q2d_val(s,k) = q2d_va(s,k)
        end if
      end do
    end do
    deallocate ( gi, gmode, gt )
    !$acc update device(q2d_val)

  end subroutine Q2D_Ghost_Values


  !> One slot per local record of type 667, in the order of the device boundary arrays (bc_perm): its index there, its
  !> block and the dimensions of the block, the variables that its file gives, and its group (file and periodic flag).
  subroutine Q2D_Setup_Slots ( domain )
    use MOSE_Mod_MPI, only: is_local_block
    implicit none
    type(MOSE_domain_type), intent(inout) :: domain
    ! Local
    integer :: i, js, k, g, s, b, nmax
    character(len=256), allocatable :: gfile(:)
    logical, allocatable :: gper(:)

    q2d_n = 0
    do i = 1, domain%nbound
      if ( domain%bc(i)%type == 667 .and. is_local_block ( domain%bc(i)%b ) ) q2d_n = q2d_n + 1
    end do
    if ( q2d_n == 0 ) return
    ! never a run with faces of type 667 and no ghost for them
    if ( .not. allocated ( domain%bc_perm ) ) error stop 'Q2D (667): the device boundary arrays are not built yet'

    nmax = q2d_n
    allocate ( q2d_idx(nmax), q2d_blk(nmax), q2d_js(nmax), q2d_grp(nmax), q2d_has(nprim,nmax), q2d_dim(3,nmax), &
               q2d_val(nprim,nmax), q2d_ref(nmax), gfile(nmax), gper(nmax) )
    q2d_has = 0
    q2d_val = 0d0
    q2d_ng  = 0
    k = 0
    do i = 1, domain%nbound
      js = domain%bc_perm(i)
      if ( domain%bc(js)%type /= 667 ) cycle
      if ( .not. is_local_block ( domain%bc(js)%b ) ) cycle
      k = k + 1
      b = domain%bc(js)%b
      q2d_idx(k)   = i
      q2d_blk(k)   = b
      q2d_js(k)    = js
      q2d_dim(:,k) = domain%blk(b)%dim
      do s = 1, nprim
        if ( domain%bc(js)%q2d_map(s)%exists ) q2d_has(s,k) = 1
      end do
      do g = 1, q2d_ng
        if ( trim(gfile(g)) == trim(domain%bc(js)%q2d_file) .and. ( gper(g) .eqv. domain%bc(js)%q2d_periodic ) ) exit
      end do
      if ( g > q2d_ng ) then
        q2d_ng = q2d_ng + 1
        g = q2d_ng
        gfile(g) = domain%bc(js)%q2d_file
        gper(g)  = domain%bc(js)%q2d_periodic
        q2d_ref(g) = js
      end if
      q2d_grp(k) = g
    end do
    deallocate ( gfile, gper )
    write(*,'(A,I0,A,I0,A)') ' [BC] 667 (MOSE_open 410) on the device: ', q2d_n, ' faces, ', q2d_ng, ' table(s)'
    allocate ( q2d_all(q2d_n) )
    do k = 1, q2d_n
      q2d_all(k) = k
    end do
    !$acc enter data copyin(q2d_idx, q2d_blk, q2d_has, q2d_dim, q2d_val, q2d_all)

  end subroutine Q2D_Setup_Slots


end module MOSE_Mod_BC_Q2D
