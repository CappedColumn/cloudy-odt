program test_move_particles
    ! Droplet transport through triplet maps, exercised through the unified
    ! eddy-sequence API in globals that both ODT and LEM now use:
    !
    !     begin_eddy_sequence / accumulate_eddy (x n) / finalize_eddy_sequence
    !     -> move_particles_by_cellmap
    !
    ! ODT is the n=1 case, LEM the n>1 case, so these tests cover both.
    use globals, only: dp, i4, N, H, z, dz_length, triplet_map, &
                       begin_eddy_sequence, accumulate_eddy, finalize_eddy_sequence, &
                       destination_cell
    use droplets, only: particle, current_n_particles, move_particles_by_cellmap
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    ! Single-eddy sequences (the ODT case)
    call test_non_wrapping_eddy()
    call test_particle_outside_eddy()
    call test_wrapping_eddy()
    call test_large_wrapping_eddy()

    ! Multi-eddy sequences (the LEM case)
    call test_composed_equals_sequential()
    call test_composed_vs_stale_gridcell()
    call test_composed_offset_preserved()
    call test_composed_in_domain_no_modulo()
    call test_tracer_matches_scalar()

    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    subroutine setup_domain(n_cells, domain_height)
        integer, intent(in) :: n_cells
        real(dp), intent(in) :: domain_height
        integer :: k

        N = n_cells
        H = domain_height
        dz_length = H / N
        if (allocated(z)) deallocate(z)
        allocate(z(N))
        do k = 1, N
            z(k) = H * k / N
        end do
    end subroutine setup_domain


    subroutine build_sequence(starts, lens)
        ! Compose an arbitrary sequence of eddies via the production API.
        integer(i4), intent(in) :: starts(:), lens(:)
        integer :: k

        call begin_eddy_sequence()
        do k = 1, size(starts)
            call accumulate_eddy(lens(k), starts(k))
        end do
        call finalize_eddy_sequence()
    end subroutine build_sequence


    subroutine check_position(name, got, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, expected
        real(dp) :: tol

        tol = 1.0e-12
        if (abs(got - expected) < tol) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,es20.12)') '    expected: ', expected
            write(*,'(a,es20.12)') '    got:      ', got
            n_failed = n_failed + 1
        end if
    end subroutine check_position


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


    logical function is_permutation(arr) result(ok)
        integer(i4), intent(in) :: arr(:)
        logical :: seen(size(arr))
        integer :: k

        seen = .false.
        ok = .true.
        do k = 1, size(arr)
            if (arr(k) < 1 .or. arr(k) > size(arr)) then
                ok = .false.; return
            end if
            if (seen(arr(k))) then
                ok = .false.; return
            end if
            seen(arr(k)) = .true.
        end do
    end function is_permutation


    ! -----------------------------------------------------------
    ! Single-eddy sequences (ODT)
    ! -----------------------------------------------------------

    subroutine test_non_wrapping_eddy()
        ! 12-cell domain, H=1.2. Eddy M=4, L=6 (cells 4-9, no wrapping).
        ! The particle in the eddy follows its fluid to destination_cell(5); the
        ! particle outside is untouched.
        type(particle) :: p(2)
        real(dp) :: offset

        call setup_domain(12, 1.2_dp)
        call build_sequence([4_i4], [6_i4])

        offset = -dz_length / 2

        p(1)%position = z(5) + offset       ! centre of cell 5, inside eddy
        p(1)%gridcell = 5
        p(1)%particle_id = 1

        p(2)%position = z(2) + offset       ! cell 2, outside eddy
        p(2)%gridcell = 2
        p(2)%particle_id = 2

        current_n_particles = 2
        call move_particles_by_cellmap(p, destination_cell)

        call check_position("non-wrapping: particle in eddy follows its fluid", &
            p(1)%position, z(destination_cell(5)) + offset)

        call check_position("non-wrapping: particle outside eddy unchanged", &
            p(2)%position, z(2) + offset)
    end subroutine test_non_wrapping_eddy


    subroutine test_particle_outside_eddy()
        ! A cell completely outside the eddy is a fixed point of the map.
        type(particle) :: p(1)

        call setup_domain(12, 1.2_dp)
        call build_sequence([7_i4], [6_i4])

        p(1)%position = z(1) - dz_length / 2
        p(1)%gridcell = 1
        p(1)%particle_id = 1

        current_n_particles = 1
        call move_particles_by_cellmap(p, destination_cell)

        call check_true("outside eddy: cell is a fixed point", &
            destination_cell(1) == 1)
        call check_position("outside eddy: particle untouched", &
            p(1)%position, z(1) - dz_length / 2)
    end subroutine test_particle_outside_eddy


    subroutine test_wrapping_eddy()
        ! Eddy M=10, L=6 wraps the top boundary: cells 10,11,12,1,2,3.
        type(particle) :: p(1)
        real(dp) :: offset

        call setup_domain(12, 1.2_dp)
        call build_sequence([10_i4], [6_i4])

        offset = -dz_length / 2
        p(1)%position = z(12) + offset
        p(1)%gridcell = 12
        p(1)%particle_id = 1

        current_n_particles = 1
        call move_particles_by_cellmap(p, destination_cell)

        call check_position("wrapping eddy: particle near boundary", &
            p(1)%position, z(destination_cell(12)) + offset)

        call check_true("wrapping eddy: position in [0, H)", &
            p(1)%position >= 0.0_dp .and. p(1)%position < H)
    end subroutine test_wrapping_eddy


    subroutine test_large_wrapping_eddy()
        ! 9-cell domain, eddy M=7, L=9 -- the whole domain, wrapping. L=9 is the
        ! smallest length whose permutation is NOT its own inverse, so this case
        ! distinguishes the forward map from the gather map. Seed every cell so
        ! the check cannot land only on fixed points.
        type(particle) :: p(9)
        real(dp) :: offset
        integer :: c
        logical :: ok

        call setup_domain(9, 0.9_dp)
        call build_sequence([7_i4], [9_i4])

        offset = -dz_length / 2
        do c = 1, N
            p(c)%position = z(c) + offset
            p(c)%gridcell = c
            p(c)%particle_id = c
        end do

        current_n_particles = 9
        call move_particles_by_cellmap(p, destination_cell)

        ok = .true.
        do c = 1, N
            if (abs(p(c)%position - (z(destination_cell(c)) + offset)) > 1.0e-12_dp) ok = .false.
            if (p(c)%position < 0.0_dp .or. p(c)%position >= H) ok = .false.
        end do

        call check_true("full-domain eddy: all particles follow their fluid", ok)
        call check_true("full-domain eddy: map is a bijection", &
            is_permutation(destination_cell))
    end subroutine test_large_wrapping_eddy


    ! -----------------------------------------------------------
    ! Multi-eddy sequences (LEM)
    ! -----------------------------------------------------------

    subroutine test_composed_equals_sequential()
        ! The invariant that justifies composing: applying the eddies one at a
        ! time (each its own single-eddy sequence, gridcell carried forward) must
        ! give exactly the same answer as composing them and moving once. This is
        ! what makes ODT's per-eddy flow and LEM's composed flow the same physics.
        type(particle) :: p_seq(12), p_com(12)
        integer(i4) :: starts(4), lens(4)
        integer :: c, k
        logical :: pos_match, cell_match

        call setup_domain(12, 1.2_dp)

        ! Overlapping, wrapping, and non-self-inverse lengths in one sequence.
        starts = [4_i4, 7_i4, 10_i4, 2_i4]
        lens   = [6_i4, 6_i4, 9_i4, 12_i4]

        do c = 1, N
            p_seq(c)%position = z(c) - dz_length / 2
            p_seq(c)%gridcell = c
            p_seq(c)%particle_id = c
        end do
        p_com = p_seq
        current_n_particles = N

        ! Reference: one eddy at a time.
        do k = 1, size(starts)
            call build_sequence(starts(k:k), lens(k:k))
            call move_particles_by_cellmap(p_seq, destination_cell)
        end do

        ! Under test: compose the same eddies, move once.
        call build_sequence(starts, lens)
        call move_particles_by_cellmap(p_com, destination_cell)

        pos_match = .true.
        cell_match = .true.
        do c = 1, N
            if (abs(p_seq(c)%position - p_com(c)%position) > 1.0e-12_dp) pos_match = .false.
            if (p_seq(c)%gridcell /= p_com(c)%gridcell) cell_match = .false.
        end do

        call check_true("equivalence: composed move matches per-eddy sequence", pos_match)
        call check_true("equivalence: composed gridcell matches per-eddy sequence", cell_match)
    end subroutine test_composed_equals_sequential


    subroutine test_composed_vs_stale_gridcell()
        ! Regression guard for the transport bug this API replaced.
        !
        ! Moving per eddy is NOT wrong in itself -- see the equivalence test
        ! above. It is wrong only if the particle's cell index is not carried
        ! forward between eddies. The old code updated position but never
        ! gridcell, so the second eddy's displacement was looked up against the
        ! cell the particle had already left. That stale-index answer must
        ! differ from the composed one for at least one cell.
        type(particle) :: p(1)
        integer(i4) :: starts(2), lens(2), composed(12)
        real(dp) :: m1z(12), m2z(12), offset, composed_pos, stale_pos
        integer :: c
        logical :: found_diff

        call setup_domain(12, 1.2_dp)
        starts = [4_i4, 7_i4]
        lens = [6_i4, 6_i4]

        call build_sequence(starts, lens)
        composed = destination_cell

        c = 5
        offset = -dz_length / 2
        p(1)%position = z(c) + offset
        p(1)%gridcell = c
        p(1)%particle_id = 1
        current_n_particles = 1
        call move_particles_by_cellmap(p, composed)

        call check_position("composed: particle follows composed map", &
            p(1)%position, z(composed(c)) + offset)

        ! Stale-index behaviour: sum both per-eddy shifts off the frozen start cell.
        m1z = z; call triplet_map(lens(1), starts(1), m1z)
        m2z = z; call triplet_map(lens(2), starts(2), m2z)
        found_diff = .false.
        do c = 1, N
            composed_pos = z(composed(c))
            stale_pos = modulo(z(c) + (m1z(c) - z(c)) + (m2z(c) - z(c)), H)
            if (abs(composed_pos - stale_pos) > 1.0e-9_dp) found_diff = .true.
        end do
        call check_true("composed: differs from stale-gridcell result (bug guard)", &
            found_diff)

        call check_true("composed: destination_cell is a bijection", &
            is_permutation(composed))
    end subroutine test_composed_vs_stale_gridcell


    subroutine test_composed_offset_preserved()
        ! Sub-cell offset is unchanged across a composed sequence of maps.
        type(particle) :: p(1)
        integer(i4) :: starts(3), lens(3)
        real(dp) :: offset_before, offset_after
        integer :: c

        call setup_domain(12, 1.2_dp)
        starts = [2_i4, 5_i4, 9_i4]
        lens = [6_i4, 6_i4, 9_i4]
        call build_sequence(starts, lens)

        c = 6
        offset_before = 0.37_dp * dz_length
        p(1)%position = z(c) - dz_length + offset_before   ! within cell c
        p(1)%gridcell = c
        p(1)%particle_id = 1
        current_n_particles = 1
        call move_particles_by_cellmap(p, destination_cell)

        offset_after = p(1)%position - (z(destination_cell(c)) - dz_length)
        call check_position("composed: sub-cell offset preserved", &
            offset_after, offset_before)
    end subroutine test_composed_offset_preserved


    subroutine test_composed_in_domain_no_modulo()
        ! Every particle lands in [0, H) with no modulo, including a full-domain
        ! wrapping eddy. Seed one particle per cell.
        type(particle) :: p(9)
        integer(i4) :: starts(2), lens(2)
        integer :: c
        logical :: all_in

        call setup_domain(9, 0.9_dp)
        starts = [7_i4, 3_i4]
        lens = [9_i4, 6_i4]
        call build_sequence(starts, lens)

        do c = 1, 9
            p(c)%position = z(c) - dz_length / 2
            p(c)%gridcell = c
            p(c)%particle_id = c
        end do
        current_n_particles = 9
        call move_particles_by_cellmap(p, destination_cell)

        all_in = .true.
        do c = 1, 9
            if (p(c)%position < 0.0_dp .or. p(c)%position >= H) all_in = .false.
            if (p(c)%gridcell /= destination_cell(c)) all_in = .false.
        end do
        call check_true("composed: all positions in [0,H), gridcell consistent", all_in)
    end subroutine test_composed_in_domain_no_modulo


    subroutine test_tracer_matches_scalar()
        ! The integer tracer must undergo exactly the same permutation as a real
        ! scalar field carrying the same cell labels (guards the duplicated
        ! arithmetic in the generic triplet_map interface).
        integer(i4) :: label_int(12), starts(2), lens(2)
        real(dp) :: label_real(12)
        integer :: k
        logical :: match

        call setup_domain(12, 1.2_dp)
        starts = [4_i4, 9_i4]
        lens = [6_i4, 9_i4]
        do k = 1, N
            label_int(k) = k
            label_real(k) = real(k, dp)
        end do
        do k = 1, size(starts)
            call triplet_map(lens(k), starts(k), label_int)
            call triplet_map(lens(k), starts(k), label_real)
        end do

        match = .true.
        do k = 1, N
            if (label_int(k) /= nint(label_real(k))) match = .false.
        end do
        call check_true("tracer: integer and real maps agree", match)
    end subroutine test_tracer_matches_scalar

end program test_move_particles
