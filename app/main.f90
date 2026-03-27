program main
  use write_particle, only: write_trajectory_data
  use globals
  use initialize, only: initialize_simulation, close_simulation
  use writeout, only: write_profiles, write_eddy
  use microphysics
  use ODT
  use droplets, only: particles, move_particles_in_eddy, update_droplets, &
                      total_n_fellout, current_n_particles, n_injected, write_trajectories
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

  call initialize_simulation()

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

  do while (time .le. tmax)

    ! Update iterators and timing
    Nt = Nt + 1
    time = time + dt
    delta_time = time - last_time_updated

    ! ---------------------------------------------------------
    ! Output (before physics — dt not yet adjusted by eddy method)
    ! ---------------------------------------------------------
    call write_profiles(dt)
    if ( do_microphysics .and. write_trajectories ) call write_trajectory_data(particles, time, dt)

    ! ---------------------------------------------------------
    ! Diffusion event
    ! ---------------------------------------------------------
    if (delta_time .ge. diffusion_step) then
      Nd = Nd + 1

      call diffusion(delta_time)
      call update_dim_scalars(T_nd, WV_nd, Tv_nd, T, WV, Tv)
      call update_supersat(T, WV, SS, pres)

      if ( do_microphysics ) call update_droplets(time, delta_time)
      if ( do_special_effects ) call run_special_effects(T, WV, delta_time)

      ! Sync nondim fields after physics modified dim arrays
      if ( do_microphysics .or. do_special_effects ) then
        call update_nondim_scalars(T, WV, Tv, T_nd, WV_nd, Tv_nd)
      end if

      last_time_updated = time
    end if

    ! ---------------------------------------------------------
    ! ODT eddy event (may adjust dt for next iteration)
    ! ---------------------------------------------------------
    if ( do_turbulence ) then
      call eddy_acceptance_method(dt, eddy_location, eddy_length, eddy_accepted)

      if ( eddy_accepted ) then
        if ( write_eddies ) call write_eddy(eddy_location, eddy_length, time)

        call diffusion(delta_time)
        call update_dim_scalars(T_nd, WV_nd, Tv_nd, T, WV, Tv)
        call update_supersat(T, WV, SS, pres)

        if ( do_microphysics ) then
          call move_particles_in_eddy(particles, eddy_location, eddy_length)
          call update_droplets(time, delta_time)
        end if

        if ( do_special_effects ) call run_special_effects(T, WV, delta_time)

        ! Sync nondim fields after physics modified dim arrays
        if ( do_microphysics .or. do_special_effects ) then
          call update_nondim_scalars(T, WV, Tv, T_nd, WV_nd, Tv_nd)
        end if

        last_time_updated = time
      end if
    end if


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

  ! Write DONE marker file to simulation output directory
  call date_and_time(date=date_str, time=time_str)
  open(newunit=done_unit, file=trim(sim_output_dir)//'DONE', &
       status='replace', action='write')
  write(done_unit,'(a,a,a,a,a,a,a,a,a)') date_str(1:4), '-', date_str(5:6), '-', date_str(7:8), &
       ' ', time_str(1:2), ':', time_str(3:4)
  close(done_unit)

end program main