program test_eddy_cellmap_standalone
    ! Isolated validation of the proposed droplet-transport scheme.
    !
    ! Deliberately self-contained: the ONLY production code exercised is
    ! triplet_map (plus the kind parameters and the global N that triplet_map's
    ! mod-indexing depends on). The grid, the particle type, the map
    ! composition/inversion and the particle move are all defined locally below,
    ! so this program validates the *scheme* before any source file is changed.
    !
    ! The invariant under test is the one that actually matters physically:
    !
    !     a droplet must end each turbulence step sitting in the fluid it
    !     started the step in.
    !
    ! That is checked directly, by carrying a labelled scalar field through the
    ! same maps and confirming the droplet's final cell holds its own label.

    use globals, only: dp, i4, N, triplet_map
    implicit none

    ! Local grid (independent of the model's z/dz_length)
    real(dp) :: local_H, local_dz
    real(dp), allocatable :: local_z(:)

    ! Local minimal particle
    type :: tracer_particle
        real(dp) :: position
        integer(i4) :: cell
        integer(i4) :: born_in_cell
    end type tracer_particle

    integer :: n_passed, n_failed
    integer(i4) :: rng_state

    n_passed = 0
    n_failed = 0
    rng_state = 20260727_i4

    call test_follows_own_fluid()
    call test_composed_equals_sequential()
    call test_offset_and_domain()
    call test_map_is_bijection()
    call report_gather_vs_forward()

    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    ! -----------------------------------------------------------
    ! Local infrastructure
    ! -----------------------------------------------------------

    subroutine setup_grid(n_cells, domain_height)
        integer(i4), intent(in) :: n_cells
        real(dp), intent(in) :: domain_height
        integer(i4) :: k

        N = n_cells
        local_H = domain_height
        local_dz = local_H / N
        if (allocated(local_z)) deallocate(local_z)
        allocate(local_z(N))
        do k = 1, N
            local_z(k) = local_H * k / N
        end do
    end subroutine setup_grid


    integer(i4) function next_rand(lo, hi) result(r)
        ! Deterministic LCG so failures are reproducible run to run.
        integer(i4), intent(in) :: lo, hi

        rng_state = mod(1103515245_i4 * rng_state + 12345_i4, 2147483647_i4)
        if (rng_state < 0) rng_state = -rng_state
        r = lo + mod(rng_state / 65536_i4, hi - lo + 1_i4)
    end function next_rand


    subroutine random_eddy(eddy_start, eddy_length)
        ! Any eddy length that is a multiple of 3 can occur, from 3 up to the
        ! whole domain; start is unrestricted (wrapping is legal).
        integer(i4), intent(out) :: eddy_start, eddy_length

        eddy_length = 3_i4 * next_rand(1_i4, N / 3_i4)
        eddy_start = next_rand(1_i4, N)
    end subroutine random_eddy


    subroutine build_forward_map(starts, lens, destination_cell)
        ! The scheme under test, written out locally.
        !
        ! Step 1: advect a cell-label tracer through the eddies in order. Because
        ! triplet_map is a gather, the tracer ends up holding origin_cell(j) =
        ! the cell whose fluid now sits at j. Composition order is handled for
        ! free by applying the maps in sequence -- no explicit nesting.
        !
        ! Step 2: invert. destination_cell(c) = where the fluid originally in c
        ! ended up. This is the sender's view a droplet needs.
        integer(i4), intent(in) :: starts(:), lens(:)
        integer(i4), intent(out) :: destination_cell(:)
        integer(i4) :: origin_cell(N)
        integer(i4) :: k, j

        do k = 1, N
            origin_cell(k) = k
        end do
        do k = 1, int(size(starts), i4)
            call triplet_map(lens(k), starts(k), origin_cell)
        end do
        do j = 1, N
            destination_cell(origin_cell(j)) = j
        end do
    end subroutine build_forward_map


    subroutine move_by_cellmap(p, destination_cell)
        ! Displace every particle once by the net cell rearrangement, keeping its
        ! sub-cell offset. The shift is a whole number of cells and
        ! destination_cell is in [1,N], so the result is inside [0,H) with no
        ! modulo needed, wrapping eddies included.
        type(tracer_particle), intent(inout) :: p(:)
        integer(i4), intent(in) :: destination_cell(:)
        integer(i4) :: i, start_cell, end_cell

        do i = 1, int(size(p), i4)
            start_cell = p(i)%cell
            end_cell = destination_cell(start_cell)
            if (end_cell /= start_cell) then
                p(i)%position = local_z(end_cell) + (p(i)%position - local_z(start_cell))
                p(i)%cell = end_cell
            end if
        end do
    end subroutine move_by_cellmap


    integer(i4) function cell_from_position(pos) result(c)
        ! Mirrors the model's update_gridcell: derive the index from position
        ! alone, so the position arithmetic is validated too.
        real(dp), intent(in) :: pos

        c = int(pos / local_dz, i4) + 1_i4
        if (c > N) c = N
        if (c < 1) c = 1
    end function cell_from_position


    subroutine seed_one_per_cell(p)
        type(tracer_particle), intent(out) :: p(:)
        integer(i4) :: c

        do c = 1, N
            p(c)%position = local_z(c) - local_dz / 2   ! cell centre
            p(c)%cell = c
            p(c)%born_in_cell = c
        end do
    end subroutine seed_one_per_cell


    subroutine check_true(name, cond)
        character(*), intent(in) :: name
        logical, intent(in) :: cond

        if (cond) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            n_failed = n_failed + 1
        end if
    end subroutine check_true


    ! -----------------------------------------------------------
    ! Tests
    ! -----------------------------------------------------------

    subroutine test_follows_own_fluid()
        ! THE test. Carry a labelled scalar field (label = origin cell) through
        ! the same eddy sequence the droplets see. Afterwards every droplet must
        ! sit in a cell whose label is its own birth cell -- i.e. it stayed with
        ! its parcel. Repeated over many random sequences of varying length.
        integer(i4), parameter :: n_trials = 300_i4
        integer(i4), parameter :: max_eddies = 6_i4
        type(tracer_particle) :: p(36)
        real(dp) :: label(36)
        integer(i4) :: destination_cell(36)
        integer(i4) :: starts(max_eddies), lens(max_eddies)
        integer(i4) :: trial, k, c, n_eddies
        logical :: ok

        ok = .true.
        do trial = 1, n_trials
            call setup_grid(36_i4, 3.6_dp)
            call seed_one_per_cell(p)
            do c = 1, N
                label(c) = real(c, dp)
            end do

            n_eddies = next_rand(1_i4, max_eddies)
            do k = 1, n_eddies
                call random_eddy(starts(k), lens(k))
                ! scalars get the map...
                call triplet_map(lens(k), starts(k), label)
            end do

            ! ...droplets get the composed forward map, once.
            call build_forward_map(starts(1:n_eddies), lens(1:n_eddies), destination_cell)
            call move_by_cellmap(p, destination_cell)

            do c = 1, N
                ! the cell the droplet now occupies must hold its own fluid
                if (nint(label(p(c)%cell)) /= p(c)%born_in_cell) ok = .false.
                ! and the index must agree with the position it was moved to
                if (cell_from_position(p(c)%position) /= p(c)%cell) ok = .false.
            end do
        end do

        call check_true("droplet stays with its own fluid (300 random sequences)", ok)
    end subroutine test_follows_own_fluid


    subroutine test_composed_equals_sequential()
        ! Composing N maps and moving once must equal applying the eddies one at
        ! a time, each with its own forward map, re-deriving the cell index from
        ! position in between (the ODT-style per-eddy flow).
        integer(i4), parameter :: n_trials = 300_i4
        integer(i4), parameter :: max_eddies = 6_i4
        type(tracer_particle) :: p_seq(36), p_com(36)
        integer(i4) :: destination_cell(36), single_map(36)
        integer(i4) :: starts(max_eddies), lens(max_eddies)
        integer(i4) :: trial, k, c, n_eddies
        logical :: ok

        ok = .true.
        do trial = 1, n_trials
            call setup_grid(36_i4, 3.6_dp)
            call seed_one_per_cell(p_seq)
            call seed_one_per_cell(p_com)

            n_eddies = next_rand(1_i4, max_eddies)
            do k = 1, n_eddies
                call random_eddy(starts(k), lens(k))
            end do

            ! per-eddy: one-element sequence each time, index re-derived between
            do k = 1, n_eddies
                call build_forward_map(starts(k:k), lens(k:k), single_map)
                call move_by_cellmap(p_seq, single_map)
                do c = 1, N
                    p_seq(c)%cell = cell_from_position(p_seq(c)%position)
                end do
            end do

            ! composed: one move for the whole sequence
            call build_forward_map(starts(1:n_eddies), lens(1:n_eddies), destination_cell)
            call move_by_cellmap(p_com, destination_cell)

            do c = 1, N
                if (abs(p_seq(c)%position - p_com(c)%position) > 1.0e-12_dp) ok = .false.
                if (p_seq(c)%cell /= p_com(c)%cell) ok = .false.
            end do
        end do

        call check_true("composed move == per-eddy sequence (300 random sequences)", ok)
    end subroutine test_composed_equals_sequential


    subroutine test_offset_and_domain()
        ! Sub-cell offset is preserved exactly and every particle lands in [0,H)
        ! without a modulo, including full-domain wrapping eddies.
        integer(i4), parameter :: n_trials = 300_i4
        type(tracer_particle) :: p(36)
        integer(i4) :: destination_cell(36)
        integer(i4) :: starts(3), lens(3)
        integer(i4) :: trial, k, c
        real(dp) :: offset_in(36), offset_out
        logical :: ok

        ok = .true.
        do trial = 1, n_trials
            call setup_grid(36_i4, 3.6_dp)
            do c = 1, N
                offset_in(c) = real(next_rand(1_i4, 999_i4), dp) / 1000.0_dp * local_dz
                p(c)%position = local_z(c) - local_dz + offset_in(c)
                p(c)%cell = c
                p(c)%born_in_cell = c
            end do

            ! include a guaranteed full-domain wrapping eddy
            call random_eddy(starts(1), lens(1))
            starts(2) = next_rand(1_i4, N); lens(2) = int(N, i4)
            call random_eddy(starts(3), lens(3))

            call build_forward_map(starts, lens, destination_cell)
            call move_by_cellmap(p, destination_cell)

            do c = 1, N
                if (p(c)%position < 0.0_dp .or. p(c)%position >= local_H) ok = .false.
                offset_out = p(c)%position - (local_z(p(c)%cell) - local_dz)
                if (abs(offset_out - offset_in(c)) > 1.0e-12_dp) ok = .false.
            end do
        end do

        call check_true("offset preserved, position in [0,H), no modulo", ok)
    end subroutine test_offset_and_domain


    subroutine test_map_is_bijection()
        ! destination_cell must always be a permutation of 1..N: droplets are
        ! neither lost nor duplicated, for any sequence.
        integer(i4), parameter :: n_trials = 300_i4
        integer(i4) :: destination_cell(36)
        integer(i4) :: starts(5), lens(5)
        integer(i4) :: trial, k
        logical :: seen(36), ok

        ok = .true.
        do trial = 1, n_trials
            call setup_grid(36_i4, 3.6_dp)
            do k = 1, 5
                call random_eddy(starts(k), lens(k))
            end do
            call build_forward_map(starts, lens, destination_cell)

            seen = .false.
            do k = 1, N
                if (destination_cell(k) < 1 .or. destination_cell(k) > N) then
                    ok = .false.
                else if (seen(destination_cell(k))) then
                    ok = .false.
                else
                    seen(destination_cell(k)) = .true.
                end if
            end do
        end do

        call check_true("forward map is a bijection", ok)
    end subroutine test_map_is_bijection


    subroutine report_gather_vs_forward()
        ! Informational, not an assertion: for each eddy length, report whether
        ! the gather table (what triplet_map leaves in an array) and the forward
        ! table (what a droplet needs) coincide. They coincide exactly when the
        ! permutation is its own inverse.
        integer(i4) :: destination_cell(36), gather(36)
        integer(i4) :: L, k, n_diff

        call setup_grid(36_i4, 3.6_dp)
        write(*,*)
        write(*,'(a)') '  -- gather vs forward table, eddy at start=7 --'
        do L = 3, 36, 3
            do k = 1, N
                gather(k) = k
            end do
            call triplet_map(L, 7_i4, gather)
            call build_forward_map([7_i4], [L], destination_cell)
            n_diff = 0
            do k = 1, N
                if (gather(k) /= destination_cell(k)) n_diff = n_diff + 1
            end do
            write(*,'(a,i3,a,i3,a)') '     L=', L, ':  cells differing = ', n_diff, &
                merge('   (identical)', '   (DIFFERENT)', n_diff == 0)
        end do
        write(*,*)
    end subroutine report_gather_vs_forward

end program test_eddy_cellmap_standalone
