program test_move_particles
    use globals, only: dp, i4, N, H, z, dz_length, triplet_map
    use droplets, only: particle, current_n_particles, move_particles_in_eddy
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    call test_non_wrapping_eddy()
    call test_particle_outside_eddy()
    call test_wrapping_eddy()
    call test_large_wrapping_eddy()

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


    subroutine test_non_wrapping_eddy()
        ! 12-cell domain, H=1.2. Eddy M=4, L=6 (cells 4-9, no wrapping).
        ! Place particle at center of cell 5. After triplet map the particle
        ! should move to the rearranged position of cell 5.
        type(particle) :: p(2)
        real(dp) :: z_before(12), z_after(12)
        integer :: k

        call setup_domain(12, 1.2_dp)

        ! Compute expected delta by triplet-mapping the z array
        z_before = z
        z_after = z
        call triplet_map(6, 4, z_after)

        ! Particle in cell 5 (inside eddy)
        p(1)%position = z(5) - dz_length / 2   ! center of cell 5
        p(1)%gridcell = 5
        p(1)%particle_id = 1

        ! Particle in cell 2 (outside eddy)
        p(2)%position = z(2) - dz_length / 2
        p(2)%gridcell = 2
        p(2)%particle_id = 2

        current_n_particles = 2
        call move_particles_in_eddy(p, 4, 6)

        ! Cell 5 delta = z_after(5) - z_before(5)
        call check_position("non-wrapping: particle in eddy moves", &
            p(1)%position, modulo(z(5) - dz_length / 2 + (z_after(5) - z_before(5)), H))

        call check_position("non-wrapping: particle outside eddy unchanged", &
            p(2)%position, z(2) - dz_length / 2)
    end subroutine test_non_wrapping_eddy


    subroutine test_particle_outside_eddy()
        ! Verify that a particle in a cell completely outside the eddy is untouched.
        type(particle) :: p(1)

        call setup_domain(12, 1.2_dp)

        p(1)%position = z(1) - dz_length / 2
        p(1)%gridcell = 1
        p(1)%particle_id = 1

        current_n_particles = 1
        call move_particles_in_eddy(p, 7, 6)

        call check_position("particle outside eddy untouched", &
            p(1)%position, z(1) - dz_length / 2)
    end subroutine test_particle_outside_eddy


    subroutine test_wrapping_eddy()
        ! 12-cell domain, H=1.2. Eddy M=10, L=6 wraps: cells 10,11,12,1,2,3.
        ! Place particle in cell 12 (near top boundary). The triplet map moves
        ! z(12) -> mapped_z(12). The particle keeps its offset within the cell.
        type(particle) :: p(1)
        real(dp) :: mapped_z(12), expected, offset

        call setup_domain(12, 1.2_dp)

        ! Compute mapped z via triplet_map (same as move_particles_in_eddy)
        mapped_z = z
        call triplet_map(6, 10, mapped_z)

        p(1)%position = z(12) - dz_length / 2
        p(1)%gridcell = 12
        p(1)%particle_id = 1
        offset = p(1)%position - z(12)

        current_n_particles = 1
        call move_particles_in_eddy(p, 10, 6)

        expected = modulo(mapped_z(12) + offset, H)
        call check_position("wrapping eddy: particle near boundary", &
            p(1)%position, expected)

        ! Verify position is in valid range [0, H)
        call check_position("wrapping eddy: position in [0, H)", &
            merge(1.0_dp, 0.0_dp, p(1)%position >= 0.0_dp .and. p(1)%position < H), 1.0_dp)
    end subroutine test_wrapping_eddy


    subroutine test_large_wrapping_eddy()
        ! 9-cell domain, H=0.9. Eddy M=7, L=9 (entire domain, wrapping).
        ! Verify the particle position is correct and within [0, H).
        type(particle) :: p(1)
        real(dp) :: mapped_z(9), expected, offset

        call setup_domain(9, 0.9_dp)

        mapped_z = z
        call triplet_map(9, 7, mapped_z)

        p(1)%position = z(7) - dz_length / 2
        p(1)%gridcell = 7
        p(1)%particle_id = 1
        offset = p(1)%position - z(7)

        current_n_particles = 1
        call move_particles_in_eddy(p, 7, 9)

        expected = modulo(mapped_z(7) + offset, H)
        call check_position("large eddy: position correct", &
            p(1)%position, expected)

        ! Verify position is in valid range [0, H)
        call check_position("large eddy: position in [0, H)", &
            merge(1.0_dp, 0.0_dp, p(1)%position >= 0.0_dp .and. p(1)%position < H), 1.0_dp)
    end subroutine test_large_wrapping_eddy

end program test_move_particles
