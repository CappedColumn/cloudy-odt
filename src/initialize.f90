module initialize
    use iso_fortran_env, only: output_unit, error_unit
    use globals
    use microphysics
    use ODT, only: initialize_ODT, close_ODT, odt_init_arrays, &
                   odt_diffuse_step, odt_turbulence_step, odt_sync_after_physics, &
                   C2, ZC2, Tdiff
    use LEM, only: initialize_LEM, lem_diffuse_step, lem_turbulence_step, lem_sync_after_physics, &
                   integral_length_scale, kolmogorov_length_scale, dissipation_rate, reynolds_number
    use special_effects, only: initialize_special_effects
    use writeout, only: initialize_buffers, deallocate_buffers, create_netcdf, &
                initialize_eddy_file, write_eddy_header_fields, &
                add_to_profile_buffer, flush_buffer, close_netcdf
    use droplets, only: initialize_microphysics, write_trajectories
    use write_particle, only: initialize_write_particle, close_particle_netcdf
    use collision_coalescence, only: write_collisions, initialize_collision_file, close_collision_file
    use parcel, only: initialize_parcel, parcel_file, initial_RH
    use radiation, only: initialize_radiation, finalize_radiation
    implicit none

    private
    public :: initialize_simulation, close_simulation

    integer(i4) :: write_buffer

contains

    subroutine initialize_simulation()

        call read_params()
        call initialize_output()
        call initialize_params()
        call initialize_arrays()
        if (simulation_mode == 'parcel') call initialize_parcel()

        call log_header()
        call create_netcdf(trim(file_prefix)//'.nc', z, ncid, simulation_name, write_buffer)
        if (do_microphysics) then
            call initialize_microphysics()
            if (write_trajectories) call initialize_write_particle(file_prefix)
            if (write_collisions) call initialize_collision_file(file_prefix)
        end if
        call initialize_buffers(write_buffer, N)
        if (do_radiation) call initialize_radiation()
        if (do_special_effects) then
            if (simulation_mode == 'chamber') then
                call initialize_special_effects((g * Tdiff * H**3) / (Tref * nu * kT))
            else
                ! TODO: parcel mode substitutes Re for Ra as a placeholder
                call initialize_special_effects(reynolds_number)
            end if
        end if

        if (write_eddies) then
            call initialize_eddy_file(file_prefix, simulation_mode)
            if (simulation_mode == 'chamber') then
                call write_eddy_header_fields([C2, ZC2, Tdiff, Tref])
            else if (simulation_mode == 'parcel') then
                call write_eddy_header_fields([integral_length_scale, &
                                               kolmogorov_length_scale, dissipation_rate])
            end if
        end if
        call add_to_profile_buffer(time, T, WV, Tv, SS)

    end subroutine initialize_simulation


    subroutine close_simulation()

        if (simulation_mode == 'chamber') call close_ODT()
        call update_supersat(T, WV, SS, pres)
        call add_to_profile_buffer(time, T, WV, Tv, SS)
        call flush_buffer()
        call close_netcdf(ncid)
        if (do_microphysics .and. write_trajectories) call close_particle_netcdf()
        if (write_collisions) call close_collision_file()
        if (do_radiation) call finalize_radiation()
        call deallocate_buffers()

    end subroutine close_simulation


    subroutine read_params()
        integer     :: ierr, nml_unit
        character(256) :: nml_line, io_emsg

        namelist /PARAMETERS/ N, tmax, Tref, pres, H, volume_scaling, &
        same_random, write_buffer, do_turbulence, do_microphysics, &
        simulation_name, output_directory, write_eddies, do_special_effects, write_timer, &
        overwrite, simulation_mode, do_radiation, do_entrainment

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

    end subroutine read_params


    subroutine initialize_params()
        integer, allocatable :: rand_seed(:)
        integer :: rand_size

        write(*,*) 'Setting domain variables...'
        Tref = Tref + Tice
        time = 0.
        last_time_updated = 0.

        dz_length = H/N
        domain_volume = volume_scaling * domain_width**2 * H
        gridcell_volume = domain_volume / N

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
        end if
        call update_supersat(T, WV, SS, pres)

    end subroutine initialize_arrays


    subroutine initialize_output()
        integer :: ierr, dirlen, parent_end
        logical :: file_exists, parent_exists
        character(256) :: cwd

        ! Resolve relative output_directory to absolute
        if (output_directory(1:1) /= '/') then
            call getcwd(cwd)
            output_directory = trim(cwd)//'/'//trim(output_directory)
        end if

        ! Ensure trailing slash
        dirlen = len_trim(output_directory)
        if (output_directory(dirlen:dirlen) /= '/') then
            output_directory = trim(output_directory)//'/'
        end if

        ! Verify parent directory exists
        dirlen = len_trim(output_directory)
        parent_end = scan(output_directory(1:dirlen-1), '/', back=.true.)
        if (parent_end > 0) then
            inquire(file=output_directory(1:parent_end), exist=parent_exists)
            if (.not. parent_exists) then
                write(0,*) 'Error: parent directory does not exist: ', output_directory(1:parent_end)
                write(0,*) 'Full output_directory: ', trim(output_directory)
                stop 1
            end if
        end if

        sim_output_dir = trim(output_directory)
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


    subroutine log_header()
        character(8) :: date_str
        character(10) :: time_str

        call date_and_time(date=date_str, time=time_str)
        write(*,'(a)') ''
        write(*,'(a)') '              .~~.    .~~.'
        write(*,'(a)') '          .~~.    )  (    .~~.'
        write(*,'(a)') '        (    .~~.  ~~  .~~.    )'
        write(*,'(a)') '         )  (                )  ('
        write(*,'(a)') '        (    ####  ###  ####  #####  )'
        write(*,'(a)') '         )  #     #   # #   #   #  ('
        write(*,'(a)') '        (   #     #   # #   #   #    )'
        write(*,'(a)') '         )   ####  ###  ####    #  ('
        write(*,'(a)') '          (                       )'
        write(*,'(a)') '           ~~~~~~~~~~~~~~~~~~~~~~~~'
        write(*,'(a)') '            /  /  /  /  /  /  /  /'
        write(*,'(a)') '              /  /  /  /  /  /  /'
        write(*,'(a)') ''
        write(*,'(a,a,a1,a,a1,a,a,a,a1,a,a1,a)') &
             ' Started: ', date_str(1:4), '-', date_str(5:6), '-', date_str(7:8), &
             ' ', time_str(1:2), ':', time_str(3:4), ':', time_str(5:6)
        write(*,'(a)') ''
        write(*,'(a,a)')    ' Namelist:       ', trim(namelist_path)
        write(*,'(a,a)')    ' Mode:           ', trim(simulation_mode)
        write(*,'(a,i0)')   ' N:              ', N
        write(*,'(a,f0.1)') ' tmax (s):       ', tmax
        write(*,'(a,f0.4)') ' H (m):          ', H
        write(*,'(a,f0.1)') ' volume_scaling: ', volume_scaling

        if (simulation_mode == 'chamber') then
            write(*,'(a,f0.2)') ' Tref (K):       ', Tref
            write(*,'(a,f0.2)') ' Tdiff (K):      ', Tdiff
        else if (simulation_mode == 'parcel') then
            write(*,'(a,f0.2)')  ' Tref (K):       ', Tref
            write(*,'(a,f0.2)')  ' pres (mb):      ', pres / Pa_per_mb
            write(*,'(a,f0.4)')  ' initial_RH:     ', initial_RH
            write(*,'(a,es9.2)') ' L_int (m):      ', integral_length_scale
            write(*,'(a,es9.2)') ' L_kolm (m):     ', kolmogorov_length_scale
            write(*,'(a,es9.2)') ' epsilon (m2/s3):', dissipation_rate
        end if

        write(*,'(a,l1)') ' do_turbulence:       ', do_turbulence
        write(*,'(a,l1)') ' do_microphysics:     ', do_microphysics
        write(*,'(a,l1)') ' do_special_effects:  ', do_special_effects
        write(*,'(a)') ' ============================================'

    end subroutine log_header


    subroutine allocate_zero_arrays(A, n_array)
        integer(i4), intent(in) ::n_array
        real(dp), intent(inout), allocatable :: A(:)

        allocate(A(n_array))
        A = 0.

    end subroutine allocate_zero_arrays

end module initialize
