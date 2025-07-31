program main
  use write_particle, only: write_trajectory_controller, close_particle_files
  use globals
  use initialize, only: initialize_simulation, close_simulation
  use writeout, only: write_data, add_to_eddy_buffer
  use microphysics
  use ODT
  use droplets!, only: aerosols, particles, move_particles_in_eddy, move_particles_by_gravity, &
               !       total_n_fellout, current_n_particles, droplet_growth_model, update_particle_properties, injection_controller
  use special_effects, only: run_special_effects
  implicit none

  real(dp) :: t_start, t_end

  call cpu_time(t_start)

  ! --- Initialize simulation ---

  call initialize_simulation(output_directory) ! output_dir assigned from namelist
  
  ! -----------------------------

  do while (time_nd .le. tmax_nd)

    ! Update iterators and timing
    Nt = Nt + 1
    time_nd = time_nd + dt_nd
    dt = dt_nd/time_conv_nd
    time = time_nd/time_conv_nd
    write_time_iter = write_time_iter + dt
    delta_time_nd = time_nd - last_time
    delta_time = delta_time_nd/time_conv_nd

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------
    ! Ensure timesteps do not surpass diffusion scheme
    if (delta_time_nd .ge. diffusion_step) then
      Nd = Nd + 1

      call diffusion()

      ! Update particle positions
      if ( do_microphysics ) call update_droplets(time, delta_time)

      ! Implement Special Effects
      if ( do_special_effects ) then
        call run_special_effects(Tdim, WVdim, delta_time)
      end if

      last_time = time_nd ! Remember last "event" call
    end if
    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------
    ! One-Dimensional Turbulence Scheme
    if ( do_turbulence ) then
      ! Set and eddy length, location and acceptance probability
      ! Then either accepts or rejects the tested eddy
      call eddy_acceptance_method(eddy_location, eddy_length, eddy_accepted)
      
      if ( eddy_accepted ) then
        if ( write_eddies ) then
          call add_to_eddy_buffer(eddy_location, eddy_length, time)
        end if

        call diffusion()

        if ( do_microphysics ) then
          call move_particles_in_eddy(particles, eddy_location, eddy_length) ! Eddy movement
          call update_droplets(time, delta_time)
        end if

        ! Implement Special Effects
        if ( do_special_effects ) then
          call run_special_effects(Tdim, WVdim, delta_time)
        end if
        
        last_time = time_nd ! Remember last "event" call
      end if
    end if

    ! Ensure nondim and dim scalars are 
    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    ! Write to netCDF buffer
    if ( write_time_iter >= write_timer ) call write_data()

    if ( write_trajectories ) call write_trajectory_controller(particles, time, dt)


  end do

  write(*,*) 'Total Particles: ', current_n_particles
  write(*,*) 'fallout: ', total_n_fellout
  write(*,*) 'Injected: ', n_injected

  ! Write Out
  call close_simulation()

  call cpu_time(t_end)

  write(*,*) 'Simulation Took: ', t_end-t_start

end program main