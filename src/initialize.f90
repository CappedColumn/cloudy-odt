module initialize
    use iso_fortran_env, only: output_unit, error_unit
    use globals
    use microphysics
    use ODT, only: calc_eddy_length_cdf, diffusion
    use special_effects, only: initialize_special_effects
    use writeout, only: initialize_buffers, create_netcdf, initialize_particle_buffers, &
                initialize_eddy_buffer, add_to_profile_buffer, flush_buffer, close_netcdf
    use droplets, only: initialize_microphysics, n_DSD_bins, n_aer_category, size_distribution, &
                write_trajectories
    use write_particle, only: initialize_write_particle, close_particle_netcdf
    implicit none

    integer(i4) :: write_buffer ! Buffer size, n iterations to write to netCDF
    
    private :: allocate_zero_arrays, initialize_velocity_arrays, initialize_params
    public :: initialize_simulation
    
contains

    subroutine initialize_simulation(filename)
        ! I make make main.f90 look approachable by
        ! hiding the monsters under the bed
        character(*), intent(out) :: filename

        ! Life is nothing but a dream, as realistic as it seems
        call initialize_params(filename)

        if ( do_turbulence ) then
            call calc_eddy_length_cdf(prob_eddy_length)
        end if
        
        if ( do_microphysics ) then
            call initialize_microphysics(filename)
            call initialize_particle_buffers(n_aer_category, n_DSD_bins)
            ! Wont work right now if init_drop_each_gridpoint = .false.
            if ( write_trajectories ) call initialize_write_particle(filename)
        end if

        if ( do_special_effects ) then
            call initialize_special_effects(filename)
        end if

        call add_to_profile_buffer(time, Tdim, WVdim, Tvdim, SS, Wdim, size_distribution, statistics)

    end subroutine initialize_simulation

    subroutine close_simulation()

        call diffusion()
        call add_to_profile_buffer(time, Tdim, WVdim, Tvdim, SS, Wdim, size_distribution, statistics)
        call flush_buffer()
        call close_netcdf(ncid)
        if ( do_microphysics .and. write_trajectories ) then
            call close_particle_netcdf()
        end if

    end subroutine close_simulation


    subroutine initialize_params(filename)

        character(256), intent(out) :: filename
        ! I/O Variables
        integer     :: ierr, nml_unit, k
        character(256) :: nml_line, io_emsg
        character(100) :: file_format
        character(100) :: simulation_name
        ! parameters for initializing state of random number generator
        integer, allocatable :: rand_seed(:) ! Some compilers have array of seeds for RNG
        integer :: rand_size


        logical :: file_exists

        namelist /PARAMETERS/ N, Lmin, Lprob, tmax, Tdiff, Tref, pres, H, volume_scaling, &
        max_accept_prob, same_random, write_buffer, do_turbulence, do_microphysics, &
        simulation_name, output_directory, write_eddies, do_special_effects, write_timer, &
        overwrite

        ! Read in namelist
        write(*,*) 'Reading PARAMETERS namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) io_emsg; stop 1
        end if
        read(nml=PARAMETERS, unit=nml_unit, iostat=ierr)
        ! Print value causing namelist read error
        if (ierr .ne. 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(*,'(a)') 'Invalid Namelist Parameter: '//trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        ! Validate output_directory is an absolute path
        if (output_directory(1:1) /= '/') then
            write(0,*) 'Error: output_directory must be an absolute path.'
            write(0,*) 'Got: ', trim(output_directory)
            stop 1
        end if

        ! Calculate additional parameters dependent on namelist variables
        write(*,*) 'Setting domain variables...'
        Tref = Tref + Tice ! Convert to Kelvin
        Ttop = Tref - Tdiff
        time_conv_nd = nu / (H**2)
        tmax_nd = tmax * time_conv_nd
        Lmax = int(N / 3)
        time_nd = 0.
        time = 0.
        last_time = 0.
        dt_nd = 1. / (1. * N * N) ! initial timestep
        diffusion_step = 1. / (1.*N*N) ! Diffusional timestep (fixed)
        
        domain_volume = volume_scaling * domain_width**2 * H
        gridcell_volume = domain_volume / N

        ! Calculate water vapor mixing ratio and virtual temp at boundary conditions
        WVref = saturation_mixing_ratio(Tref, pres) ! kg/kg
        WVtop = saturation_mixing_ratio(Ttop, pres)
        WVdiff = WVref - WVtop
        Tvref = virtual_temp(Tref, WVref)
        Tvtop = virtual_temp(Ttop, WVtop)
        Tvdiff = Tvref - Tvtop

        buoy_nd = (8. * g * alpha * Tvdiff * C2 * H * H * H)/(27. * nu * nu)

        ! ODT parameters
        LpD = 2 * Lprob
        Co = exp(-LpD/(1.*Lmin))
        Cm = exp(-LpD/(1.*Lmax))
        prob_coeff = (exp(-LpD/(1.*Lmax))-exp(-LpD/(1.*Lmin)))*(N/(3.*LpD))
        w_dim_factor = nu / (H * C)

        ! initialize randomness in the model
        if (same_random) then
            call random_seed(size=rand_size)
            allocate(rand_seed(rand_size))
            rand_seed = 3959
            call random_seed(put=rand_seed)
            deallocate(rand_seed)
        else
            call random_seed()
        end if
    
        ! Allocate scalar/vector fields
        write(*,*) 'Allocating arrays...'
        call allocate_zero_arrays(z, N)
        call allocate_nondim_array(W, N)
        call allocate_nondim_array(T, N)
        call allocate_nondim_array(WV, N)
        call allocate_nondim_array(Tv, N)
        call allocate_zero_arrays(Tdim, N)
        call allocate_zero_arrays(WVdim, N)
        call allocate_zero_arrays(Tvdim, N)
        call allocate_zero_arrays(SS, N)
        call allocate_zero_arrays(Wdim, N)
        call allocate_zero_arrays(prob_eddy_length, N)

        ! Initialize scalar/vector fields with some random perturbation
        ! NOTE THIS ONLY HAS POSITIVE PERTURBATIONS - NONDIMENSIONAL???
        !W = W + random_array(W)*2.e-10

        call initialize_linear_array(z)
        call initialize_velocity_arrays(W)
        call update_dim_scalars(W, T, WV, Tv, Wdim, Tdim, WVdim, Tvdim)
        call update_supersat(Tdim, WVdim, SS, pres)

        ! Initialize positional array
        do k = 1, N
            z(k) = H*k/N ! Grid cell position in meters
        end do

        ! Create new subdirectory based on simulation name
        filename = trim(output_directory)//'/'//trim(simulation_name)
        call system("mkdir -p "//filename)

        ! Check for existing output files
        inquire(file=trim(filename)//'/'//trim(simulation_name)//'.nc', exist=file_exists)
        if (file_exists .and. .not. overwrite) then
            write(0,*) 'Error: output file already exists: ', &
                trim(filename)//'/'//trim(simulation_name)//'.nc'
            write(0,*) 'Set overwrite = .true. in the namelist to allow overwriting.'
            stop 1
        end if

        ! Print output location to stderr (visible in terminal/SLURM output)
        write(error_unit,*) 'Output: ', trim(filename)

        ! Redirect stdout to log file in output directory
        close(output_unit)
        open(output_unit, file=trim(filename)//'/'//trim(simulation_name)//'.log', &
             status='replace', action='write', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: could not open log file'
            stop 1
        end if

        ! Create filename based on Tdiff and N
        file_format = '(a, a, a)'
        write(filename, file_format) trim(filename), '/', trim(simulation_name)
        call create_netcdf(trim(filename)//".nc", z, ncid, simulation_name, write_buffer)

        ! Writeout parameters
        write_time_iter = 0.

        ! Initialize buffers for writing to netCDF
        call initialize_buffers(write_buffer, N)
        if ( write_eddies ) call initialize_eddy_buffer(filename)

        ! Copy original namelist file to output directory
        call copy_file(namelist_path, trim(filename)//'.nml')
        
        ! Length of a gridcell
        dz_length = H/N

    end subroutine initialize_params

    subroutine allocate_nondim_array(A, n_array)
        integer(i4), intent(in) ::n_array
        real(dp), intent(inout), allocatable :: A(:)
        integer :: k

        ! Set nondimensional arrays to have a ghost point at top/bottom
        allocate(A(N+1))
        A = 0.

        ! Initialize nondim with lenear profile
        do concurrent (k = 1:N+1)
            A(k) = 1.*k/(N+1)
        end do

    end subroutine allocate_nondim_array

    subroutine allocate_zero_arrays(A, n_array)
        integer(i4), intent(in) ::n_array
        real(dp), intent(inout), allocatable :: A(:)

        ! Dimensional arrays do not have ghost points
        allocate(A(n_array))
        A = 0.

    end subroutine allocate_zero_arrays

    subroutine initialize_linear_array(array)

        real(dp), intent(out) :: array(:)
        integer :: k

        do concurrent (k = 1:N)
            array(k) = 1.*k/N
        end do

    end subroutine initialize_linear_array

    subroutine initialize_velocity_arrays(lw)
        ! Initializes velcoity arrays with some noise. This prevents a numerical issue from
        ! occuring when calculating the kinetic energy (i.e. KE != 0)
        real(dp), intent(inout), allocatable :: lw(:)
        real(dp) :: rand_num
        integer(i4) :: k

        do k = 1, N
            call random_number(rand_num)
            lw(k) = 2.e-10 * (rand_num - 0.5)
        end do
        lw(N+1) = 0.

    end subroutine initialize_velocity_arrays

end module initialize