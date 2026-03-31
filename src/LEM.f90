module LEM
    use globals
    use microphysics, only: virtual_temp, update_supersat
    use droplets, only: particle, particles, current_n_particles
    use writeout, only: write_eddy
    implicit none

    private
    public :: initialize_LEM, lem_diffuse_step, lem_turbulence_step, lem_sync_after_physics
    public :: triplet_map_periodic

    ! Derived LEM parameters (set in initialize_LEM)
    integer(i4) :: maps_per_event        ! Triplet maps per eddy event
    integer(i4) :: steps_between_events  ! Diffusion steps between eddy events
    integer(i4) :: iteration_count       ! Running iteration counter

    real(dp) :: turbulent_diffusivity    ! 0.1 * L_int^(4/3) * epsilon^(1/3) (m^2/s)
    real(dp) :: reynolds_number          ! (L_int / eta)^(4/3)
    real(dp) :: vapor_diffusivity        ! Molecular diffusivity for water vapor (m^2/s)
    real(dp) :: thermal_diffusivity      ! Molecular diffusivity for temperature (m^2/s)

contains

    subroutine initialize_LEM(domain_height)
        real(dp), intent(in) :: domain_height
        real(dp) :: diffusion_timestep, convection_timestep
        real(dp) :: large_eddy_turnover_time, eddy_rate_per_length

        reynolds_number = (integral_length_scale / kolmogorov_length_scale) ** (4./3.)

        turbulent_diffusivity = 0.1 * integral_length_scale**(4./3.) &
                                * dissipation_rate**(1./3.)

        thermal_diffusivity = kT
        vapor_diffusivity   = Dv

        ! Diffusion stability timestep
        diffusion_timestep = 0.2 * dz_length**2 / max(vapor_diffusivity, thermal_diffusivity)

        ! Turbulent convection timestep (Krueger 1993)
        large_eddy_turnover_time = integral_length_scale**2 / turbulent_diffusivity
        eddy_rate_per_length = (54./5. * reynolds_number**1.25) &
                               / (integral_length_scale * large_eddy_turnover_time)
        convection_timestep = 1. / (domain_height * eddy_rate_per_length)

        ! Determine maps_per_event and steps_between_events
        if (diffusion_timestep >= convection_timestep) then
            dt = diffusion_timestep
            maps_per_event = int(diffusion_timestep / convection_timestep)
            steps_between_events = 1
        else
            steps_between_events = int(convection_timestep / diffusion_timestep) + 1
            dt = convection_timestep / real(steps_between_events, dp)
            maps_per_event = 1
        end if

        diffusion_step = dt
        iteration_count = 0

        write(*,*) '--- LEM Configuration ---'
        write(*,*) 'Re:                  ', reynolds_number
        write(*,*) 'turbulent_diffusivity:', turbulent_diffusivity
        write(*,*) 'dt (s):              ', dt
        write(*,*) 'maps_per_event:      ', maps_per_event
        write(*,*) 'steps_between_events:', steps_between_events
        write(*,*) '-------------------------'

    end subroutine initialize_LEM


    ! -----------------------------------------------
    ! Controller subroutines matching abstract interfaces
    ! -----------------------------------------------

    subroutine lem_diffuse_step(ldelta_time, fields_updated)
        real(dp), intent(in) :: ldelta_time
        logical, intent(out) :: fields_updated

        Nd = Nd + 1
        call diffuse_scalar_periodic(T, thermal_diffusivity, ldelta_time)
        call diffuse_scalar_periodic(WV, vapor_diffusivity, ldelta_time)
        call update_virtual_temperature()
        call update_supersat(T, WV, SS, pres)
        fields_updated = .true.

    end subroutine lem_diffuse_step


    subroutine lem_turbulence_step(ldt, ltime, ldelta_time, &
                                   leddy_accepted, eddy_loc, eddy_len)
        real(dp), intent(inout) :: ldt
        real(dp), intent(in) :: ltime, ldelta_time
        logical, intent(out) :: leddy_accepted
        integer(i4), intent(out) :: eddy_loc, eddy_len

        integer(i4) :: j, gridpoints_raw, eddy_gridpoints, eddy_start
        real(dp) :: rand_size, rand_position, sampled_eddy_size

        iteration_count = iteration_count + 1
        leddy_accepted = .false.
        eddy_loc = 0
        eddy_len = 0

        if (mod(iteration_count, steps_between_events) /= 0) return

        do j = 1, maps_per_event
            ! Sample eddy size from -5/3 inertial-subrange spectrum
            call random_number(rand_size)
            sampled_eddy_size = (rand_size &
                * (integral_length_scale**(-5./3.) &
                   - kolmogorov_length_scale**(-5./3.)) &
                + kolmogorov_length_scale**(-5./3.)) ** (-3./5.)

            ! Convert to gridpoints, quantize to nearest multiple of 3
            gridpoints_raw = int(sampled_eddy_size / dz_length)
            eddy_gridpoints = nint(real(gridpoints_raw, dp) / 3.) * 3
            eddy_gridpoints = max(3, min(eddy_gridpoints, N))

            ! Random starting position (periodic wrapping handles boundary)
            call random_number(rand_position)
            eddy_start = int(rand_position * N) + 1

            ! Apply periodic triplet map to dimensional scalars
            call triplet_map_periodic(eddy_gridpoints, eddy_start, T)
            call triplet_map_periodic(eddy_gridpoints, eddy_start, WV)

            ! Move particles
            if (do_microphysics) then
                call move_particles_in_eddy_periodic(particles, eddy_start, eddy_gridpoints)
            end if

            eddy_loc = eddy_start
            eddy_len = eddy_gridpoints
        end do

        call update_virtual_temperature()
        call update_supersat(T, WV, SS, pres)

        if (write_eddies) call write_eddy(eddy_loc, eddy_len, ltime)
        leddy_accepted = .true.

    end subroutine lem_turbulence_step


    subroutine lem_sync_after_physics()
        ! No-op: LEM operates entirely in dimensional space
    end subroutine lem_sync_after_physics


    ! -----------------------------------------------
    ! LEM-specific physics routines
    ! -----------------------------------------------

    subroutine diffuse_scalar_periodic(field, molecular_diffusivity, elapsed_time)
        ! Crank-Nicolson diffusion with periodic BCs via Sherman-Morrison.
        real(dp), intent(inout) :: field(:)
        real(dp), intent(in) :: molecular_diffusivity, elapsed_time
        real(dp) :: De, gamma, correction
        real(dp) :: diag_val, off_diag
        real(dp) :: y_soln(N), q_soln(N)
        real(dp) :: rhs(N)
        real(dp) :: lower(N), diag(N), upper(N)
        real(dp) :: q_rhs(N)
        integer(i4) :: k

        De = (elapsed_time * molecular_diffusivity) / (2.0 * dz_length**2)
        diag_val = 1.0 + 2.0 * De
        off_diag = -De

        ! Build RHS (explicit side) with periodic wrapping
        rhs(1) = (1.0 - 2.0*De)*field(1) + De*(field(2) + field(N))
        do k = 2, N-1
            rhs(k) = (1.0 - 2.0*De)*field(k) + De*(field(k+1) + field(k-1))
        end do
        rhs(N) = (1.0 - 2.0*De)*field(N) + De*(field(1) + field(N-1))

        ! Sherman-Morrison decomposition of cyclic tridiagonal
        gamma = -diag_val

        ! Modified diagonal for standard tridiagonal system
        diag(1) = diag_val - gamma
        do k = 2, N-1
            diag(k) = diag_val
        end do
        diag(N) = diag_val - off_diag * off_diag / gamma

        lower(:) = off_diag
        upper(:) = off_diag

        ! Solve A' * y_soln = rhs
        call tridiagonal_periodic(lower, diag, upper, rhs, y_soln)

        ! Solve A' * q_soln = [gamma, 0, ..., 0, off_diag]
        q_rhs(:) = 0.0
        q_rhs(1) = gamma
        q_rhs(N) = off_diag
        call tridiagonal_periodic(lower, diag, upper, q_rhs, q_soln)

        ! Apply Sherman-Morrison correction
        correction = (y_soln(1) + off_diag * y_soln(N) / gamma) &
                   / (1.0 + q_soln(1) + off_diag * q_soln(N) / gamma)

        do k = 1, N
            field(k) = y_soln(k) - correction * q_soln(k)
        end do

    end subroutine diffuse_scalar_periodic


    subroutine tridiagonal_periodic(l, d, u, rhs, x)
        ! Standard tridiagonal solve (Thomas algorithm) for size N system.
        real(dp), intent(in) :: l(:), d(:), u(:), rhs(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: w(N), b
        integer(i4) :: k

        b = d(1)
        x(1) = rhs(1) / b
        do k = 2, N
            w(k) = u(k-1) / b
            b = d(k) - l(k) * w(k)
            x(k) = (rhs(k) - l(k) * x(k-1)) / b
        end do

        do k = N-1, 1, -1
            x(k) = x(k) - w(k+1) * x(k+1)
        end do

    end subroutine tridiagonal_periodic


    subroutine triplet_map_periodic(eddy_length, eddy_start, field)
        ! Triplet map with periodic wrapping via mod indexing.
        integer(i4), intent(in) :: eddy_length, eddy_start
        real(dp), intent(inout) :: field(:)

        real(dp) :: mapped_values(eddy_length)
        integer(i4) :: j, source_index, dest_index, segment_length

        segment_length = eddy_length / 3

        ! Segment 1: every 3rd element, forward
        do j = 1, segment_length
            source_index = mod(eddy_start + 3*(j-1) - 1, N) + 1
            mapped_values(j) = field(source_index)
        end do

        ! Segment 2: every 3rd element, reversed (block inversion)
        do j = 1, segment_length
            source_index = mod(eddy_start + eddy_length - 3*j, N) + 1
            mapped_values(j + segment_length) = field(source_index)
        end do

        ! Segment 3: every 3rd element, forward offset by 2
        do j = 1, segment_length
            source_index = mod(eddy_start + 3*j - 2, N) + 1
            mapped_values(j + 2*segment_length) = field(source_index)
        end do

        ! Write rearranged values back to the periodic domain
        do j = 1, eddy_length
            dest_index = mod(eddy_start + j - 2, N) + 1
            field(dest_index) = mapped_values(j)
        end do

    end subroutine triplet_map_periodic


    subroutine move_particles_in_eddy_periodic(lparticles, eddy_start, eddy_length)
        ! Move particles within a periodic eddy using position deltas.
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(in) :: eddy_start, eddy_length

        real(dp) :: original_positions(eddy_length), rearranged_positions(eddy_length)
        real(dp) :: position_delta(eddy_length)
        integer(i4) :: eddy_cells(eddy_length)
        integer(i4) :: j, k, particle_cell, segment_length

        ! Build list of gridcells in the eddy (with periodic wrapping)
        do k = 1, eddy_length
            eddy_cells(k) = mod(eddy_start + k - 2, N) + 1
            original_positions(k) = z(eddy_cells(k))
        end do

        ! Apply triplet map to position array (local indices, no wrapping needed)
        segment_length = eddy_length / 3
        rearranged_positions(1:segment_length) = &
            original_positions(1:eddy_length:3)
        rearranged_positions(segment_length+1:2*segment_length) = &
            original_positions(eddy_length:1:-3)
        rearranged_positions(2*segment_length+1:eddy_length) = &
            original_positions(3:eddy_length:3)

        ! Compute position deltas
        position_delta = rearranged_positions - original_positions

        ! Move particles that are in eddy gridcells
        do j = 1, current_n_particles
            particle_cell = lparticles(j)%gridcell
            do k = 1, eddy_length
                if (particle_cell == eddy_cells(k)) then
                    lparticles(j)%position = lparticles(j)%position + position_delta(k)
                    ! Wrap position to [0, H) for periodic domain
                    lparticles(j)%position = modulo(lparticles(j)%position, H)
                    exit
                end if
            end do
        end do

    end subroutine move_particles_in_eddy_periodic


    subroutine update_virtual_temperature()
        integer(i4) :: k
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
    end subroutine update_virtual_temperature

end module LEM
