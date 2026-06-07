! Top-level driver: owns the simulation time loop. All physics lives in the
! modules pulled in below; this module just sequences them each time step. The
! turbulence/diffusion/sync routines are reached through procedure pointers
! (set at init by initialize_simulation) so the same loop serves both chamber
! (ODT) and parcel (LEM) modes.
module CODT
    use write_particle, only: write_trajectory_data
    use globals
    use initialize, only: initialize_simulation, close_simulation
    use writeout, only: write_profiles, write_eddy
    use droplets, only: particles, update_droplets, &
                        total_n_fellout, current_n_particles, n_injected, write_trajectories
    use special_effects, only: run_special_effects
    use parcel, only: do_parcel_ascent, apply_adiabatic_forcing, apply_parcel_entrainment, &
                      pressure_limit_reached
    use radiation, only: compute_radiation
    implicit none

    private
    public :: run_simulation

contains

    subroutine run_simulation()
        ! Initialize, then march from time=0 to tmax. Each step:
        !   1. (parcel) apply adiabatic forcing and optional entrainment
        !   2. write output (before physics, so t=0 state is captured)
        !   3. diffusion backstop: if enough time has elapsed, diffuse then update
        !      droplets -> special effects -> radiation -> sync fields
        !   4. turbulence: if an eddy is accepted, run the same physics chain again
        ! The diffuse -> droplets -> ... -> sync ordering is an invariant: droplets
        ! must see updated fields, and dim/nondim fields must be synced afterward
        ! (see CLAUDE.md). On exit, finalize output and write the _DONE marker.
        real(dp) :: t_start, t_end
        integer :: done_unit
        integer(i4) :: eddy_location, eddy_length
        logical :: eddy_accepted
        character(8) :: date_str
        character(10) :: time_str

        call cpu_time(t_start)

        call initialize_simulation()

        do while (time .le. tmax)

            Nt = Nt + 1
            time = time + dt
            if (do_parcel_ascent) call apply_adiabatic_forcing(dt)
            if (pressure_limit_reached) exit
            if (do_entrainment) call apply_parcel_entrainment()
            delta_time = time - last_time_updated

            ! Output
            call write_profiles(dt)
            if (do_microphysics .and. write_trajectories) call write_trajectory_data(particles, time, dt)

            ! Diffusion backstop: guarantees scalars diffuse at least every
            ! diffusion_step even when few/no eddies are accepted (delta_time is
            ! time accumulated since the last physics update).
            if (delta_time >= diffusion_step) then
                call diffuse_step(delta_time)
                if (do_microphysics) call update_droplets(time, delta_time)
                if (do_special_effects) call run_special_effects(T, WV, delta_time)
                if (do_radiation) call compute_radiation(T, delta_time, time)
                call sync_after_physics()
                last_time_updated = time
            end if

            ! Turbulence
            if (do_turbulence) call turbulence_step(dt, time, delta_time, &
                                                    eddy_accepted, eddy_location, eddy_length)
            if (eddy_accepted) then
                if (write_eddies) call write_eddy(eddy_location, eddy_length, time)
                call diffuse_step(delta_time)
                if (do_microphysics) call update_droplets(time, delta_time)
                if (do_special_effects) call run_special_effects(T, WV, delta_time)
                if (do_radiation) call compute_radiation(T, delta_time, time)
                call sync_after_physics()
                last_time_updated = time
            end if

        end do

        call close_simulation()

        call cpu_time(t_end)

        write(*,*) '--- Run Results ---'
        write(*,*) 'Total Particles: ', current_n_particles
        write(*,*) 'Fallout: ', total_n_fellout
        write(*,*) 'Injected: ', n_injected
        write(*,*) 'Wall-clock time (s): ', t_end - t_start

        call date_and_time(date=date_str, time=time_str)
        open(newunit=done_unit, file=trim(file_prefix)//'_DONE', &
             status='replace', action='write')
        write(done_unit,'(a,a,a,a,a,a,a,a,a)') date_str(1:4), '-', date_str(5:6), '-', date_str(7:8), &
             ' ', time_str(1:2), ':', time_str(3:4)
        close(done_unit)

    end subroutine run_simulation

end module CODT
