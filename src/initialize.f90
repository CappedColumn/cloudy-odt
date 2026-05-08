module initialize
    use iso_fortran_env, only: output_unit, error_unit
    use globals
    use microphysics
    use ODT, only: initialize_ODT, close_ODT, odt_init_arrays, &
                   odt_diffuse_step, odt_turbulence_step, &
                   odt_sync_after_physics, Lmin, Lprob, max_accept_prob
    use LEM, only: initialize_LEM, lem_diffuse_step, lem_turbulence_step, lem_sync_after_physics
    use special_effects, only: initialize_special_effects
    use writeout, only: initialize_buffers, deallocate_buffers, create_netcdf, &
                initialize_eddy_file, add_to_profile_buffer, flush_buffer, close_netcdf
    use droplets, only: initialize_microphysics, write_trajectories
    use write_particle, only: initialize_write_particle, close_particle_netcdf
    use collision_coalescence, only: write_collisions, initialize_collision_file, close_collision_file
    use dynamics, only: initialize_dynamics
    implicit none

    private
    public :: initialize_simulation, close_simulation

    integer(i4) :: write_buffer

contains

    subroutine initialize_simulation()

        call read_params()
        call initialize_params()
        call initialize_arrays()
        call initialize_output()

        call create_netcdf(trim(file_prefix)//'.nc', z, ncid, simulation_name, write_buffer)
        if (do_microphysics) then
            call initialize_microphysics()
            if (write_trajectories) call initialize_write_particle(file_prefix)
            if (write_collisions) call initialize_collision_file(file_prefix)
        end if
        call initialize_buffers(write_buffer, N)
        if (do_special_effects) call initialize_special_effects()

        if (write_eddies) call initialize_eddy_file(file_prefix)
        call copy_file(namelist_path, trim(file_prefix)//'.nml')
        if (dynamics_file /= '') then
            call copy_file(resolve_path(namelist_dir, dynamics_file), &
                           trim(sim_output_dir)// &
                           trim(dynamics_file(scan(trim(dynamics_file), '/', back=.true.)+1:)))
        end if
        call add_to_profile_buffer(time, T, WV, Tv, SS)

    end subroutine initialize_simulation


    subroutine close_simulation()

        if (simulation_mode == 'chamber') call close_ODT()
        call update_supersat(T, WV, SS, pres)
        call add_to_profile_buffer(time, T, WV, Tv, SS)
        call flush_buffer()
        call close_netcdf(ncid, Lmin, Lprob, max_accept_prob)
        if (do_microphysics .and. write_trajectories) call close_particle_netcdf()
        if (write_collisions) call close_collision_file()
        call deallocate_buffers()

    end subroutine close_simulation


    subroutine read_params()
        integer     :: ierr, nml_unit
        character(256) :: nml_line, io_emsg

        namelist /PARAMETERS/ N, tmax, Tdiff, Tref, pres, H, volume_scaling, &
        same_random, write_buffer, do_turbulence, do_microphysics, &
        simulation_name, output_directory, write_eddies, do_special_effects, write_timer, &
        overwrite, simulation_mode, dynamics_file

        write(*,*) 'Reading PARAMETERS namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) io_emsg; stop 1
        end if
        read(nml=PARAMETERS, unit=nml_unit, iostat=ierr)
        if (ierr .ne. 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(*,'(a)') 'Invalid Namelist Parameter: '//trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        if (output_directory(1:1) /= '/') then
            write(0,*) 'Error: output_directory must be an absolute path.'
            write(0,*) 'Got: ', trim(output_directory)
            stop 1
        end if

    end subroutine read_params


    subroutine initialize_params()
        integer, allocatable :: rand_seed(:)
        integer :: rand_size

        write(*,*) 'Setting domain variables...'
        Tref = Tref + Tice
        Ttop = Tref - Tdiff
        time = 0.
        last_time_updated = 0.

        dz_length = H/N
        domain_volume = volume_scaling * domain_width**2 * H
        gridcell_volume = domain_volume / N

        WVref = saturation_mixing_ratio(Tref, pres)
        WVtop = saturation_mixing_ratio(Ttop, pres)
        WVdiff = WVref - WVtop
        Tvref = virtual_temp(Tref, WVref)
        Tvtop = virtual_temp(Ttop, WVtop)
        Tvdiff = Tvref - Tvtop

        if (simulation_mode == 'chamber') then
            call initialize_ODT(H)
            diffuse_step       => odt_diffuse_step
            turbulence_step    => odt_turbulence_step
            sync_after_physics => odt_sync_after_physics
        else if (simulation_mode == 'parcel') then
            call initialize_LEM(H)
            diffuse_step       => lem_diffuse_step
            turbulence_step    => lem_turbulence_step
            sync_after_physics => lem_sync_after_physics
            if (dynamics_file /= '') then
                call initialize_dynamics(resolve_path(namelist_dir, dynamics_file))
            end if
        else
            write(0,*) 'Error: unknown simulation_mode: ', trim(simulation_mode)
            stop 1
        end if

        if (same_random) then
            call random_seed(size=rand_size)
            allocate(rand_seed(rand_size))
            rand_seed = 3959
            call random_seed(put=rand_seed)
            deallocate(rand_seed)
        else
            call random_seed()
        end if

    end subroutine initialize_params


    subroutine initialize_arrays()
        integer :: k

        write(*,*) 'Allocating arrays...'
        call allocate_zero_arrays(z, N)
        call allocate_zero_arrays(T, N)
        call allocate_zero_arrays(WV, N)
        call allocate_zero_arrays(Tv, N)
        call allocate_zero_arrays(SS, N)

        do k = 1, N
            z(k) = H*k/N
        end do

        if (simulation_mode == 'chamber') then
            call odt_init_arrays()
        else if (simulation_mode == 'parcel') then
            T(:) = Tref
            WV(:) = WVref
            do k = 1, N
                Tv(k) = virtual_temp(T(k), WV(k))
            end do
        end if
        call update_supersat(T, WV, SS, pres)

    end subroutine initialize_arrays


    subroutine initialize_output()
        integer :: ierr
        logical :: file_exists

        sim_output_dir = trim(output_directory)//'/'//trim(simulation_name)//'/'
        file_prefix = trim(sim_output_dir)//trim(simulation_name)
        call system("mkdir -p "//trim(sim_output_dir))

        inquire(file=trim(file_prefix)//'.nc', exist=file_exists)
        if (file_exists .and. .not. overwrite) then
            write(0,*) 'Error: output file already exists: ', trim(file_prefix)//'.nc'
            write(0,*) 'Set overwrite = .true. in the namelist to allow overwriting.'
            stop 1
        end if

        write(error_unit,*) 'Output: ', trim(sim_output_dir)

        close(output_unit)
        open(output_unit, file=trim(file_prefix)//'.log', &
             status='replace', action='write', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: could not open log file'
            stop 1
        end if


    end subroutine initialize_output


    subroutine allocate_zero_arrays(A, n_array)
        integer(i4), intent(in) ::n_array
        real(dp), intent(inout), allocatable :: A(:)

        allocate(A(n_array))
        A = 0.

    end subroutine allocate_zero_arrays

end module initialize
