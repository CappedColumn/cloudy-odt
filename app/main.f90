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
  integer :: done_unit
  character(8) :: date_str
  character(10) :: time_str

  call cpu_time(t_start)

  ! --- Parse command-line argument ---
  if (command_argument_count() < 1) then
    write(0,*) 'Usage: codt <namelist_path>'
    write(0,*) 'Example: codt /path/to/input/params.nml'
    stop 1
  end if
  call get_command_argument(1, namelist_path)
  if (scan(trim(namelist_path), '/') == 0) then
    write(0,*) 'Error: namelist path must include a directory.'
    write(0,*) 'Use ./params.nml for the current directory.'
    stop 1
  end if
  namelist_dir = parent_directory(namelist_path)

  ! --- Initialize simulation ---

  call initialize_simulation(output_directory) ! output_dir assigned from namelist

  ! Log simulation configuration
  write(*,*) '--- Simulation Configuration ---'
  write(*,*) 'Namelist: ', trim(namelist_path)
  write(*,*) 'N: ', N
  write(*,*) 'tmax (s): ', tmax
  write(*,*) 'Tdiff (K): ', Tdiff
  write(*,*) 'Tref (K): ', Tref
  write(*,*) 'H (m): ', H
  write(*,*) 'volume_scaling: ', volume_scaling
  write(*,*) 'do_turbulence: ', do_turbulence
  write(*,*) 'do_microphysics: ', do_microphysics
  write(*,*) 'do_special_effects: ', do_special_effects
  write(*,*) '--------------------------------'

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

  ! Write Out
  call close_simulation()

  call cpu_time(t_end)

  ! Log run results
  write(*,*) '--- Run Results ---'
  write(*,*) 'Total Particles: ', current_n_particles
  write(*,*) 'Fallout: ', total_n_fellout
  write(*,*) 'Injected: ', n_injected
  write(*,*) 'Wall-clock time (s): ', t_end - t_start

  ! Write DONE marker file to output directory
  call date_and_time(date=date_str, time=time_str)
  open(newunit=done_unit, file=trim(parent_directory(output_directory))//'DONE', &
       status='replace', action='write')
  write(done_unit,'(a,a,a,a,a,a,a,a,a)') date_str(1:4), '-', date_str(5:6), '-', date_str(7:8), &
       ' ', time_str(1:2), ':', time_str(3:4)
  close(done_unit)

end program main