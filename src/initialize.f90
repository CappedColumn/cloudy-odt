module initialize
    use globals
    use microphysics
    use ODT, only: calc_eddy_length_cdf, diffusion
    use special_effects, only: initialize_special_effects
    use writeout, only: initialize_buffers, create_netcdf, initialize_particle_buffers, &
    initialize_eddy_buffer, add_to_profile_buffer, flush_buffer, close_netcdf
    use droplets, only: initialize_microphysics, n_DSD_bins, size_distribution, write_trajectories, close_droplets
    use write_particle, only: initialize_write_particle, close_particle_files
    implicit none

    integer(i4) :: write_buffer ! Buffer size, n iterations to write to netCDF
    
    private :: allocate_zero_arrays, initialize_velocity_arrays, initialize_scalar_arrays, initialize_params
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
            call initialize_particle_buffers(n_DSD_bins)
            ! Wont work right now if init_drop_each_gridpoint = .false.
            if ( write_trajectories ) call initialize_write_particle(output_directory)
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
        if ( write_trajectories ) then
            call close_droplets()
            call close_particle_files()
        end if

    end subroutine close_simulation


    subroutine initialize_params(filename)

        character(100), intent(out) :: filename
        ! I/O Variables
        integer     :: ierr, nml_unit, k
        character(100) :: nml_line, io_emsg
        character(100) :: file_format
        character(100) :: simulation_name
        ! parameters for initializing state of random number generator
        integer, allocatable :: rand_seed(:) ! Some compilers have array of seeds for RNG
        integer :: rand_size


        namelist /PARAMETERS/ N, Lmin, Lprob, tmax, Tdiff, Tref, pres, H, volume_scaling, &
        max_accept_prob, same_random, write_buffer, do_turbulence, do_microphysics, &
        simulation_name, output_directory, write_eddies, do_special_effects, write_timer

        ! Read in namelist
        write(*,*) 'Reading PARAMETERS namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) io_emsg; stop
        end if
        read(nml=PARAMETERS, unit=nml_unit, iostat=ierr)
        ! Print value causing namelist read error
        if (ierr .ne. 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(*,'(a)') 'Invalid Namelist Parameter: '//trim(nml_line)
            stop
        end if
        close(nml_unit)

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
        call allocate_zero_arrays(W, N)
        call allocate_zero_arrays(T, N)
        call allocate_zero_arrays(WV, N)
        call allocate_zero_arrays(Tv, N)
        call allocate_zero_arrays(Tdim, N)
        call allocate_zero_arrays(WVdim, N)
        call allocate_zero_arrays(Tvdim, N)
        call allocate_zero_arrays(SS, N)
        call allocate_zero_arrays(Wdim, N)
        call allocate_zero_arrays(prob_eddy_length, N)

        ! Initialize scalar/vector fields with some random perturbation
        ! NOTE THIS ONLY HAS POSITIVE PERTURBATIONS - NONDIMENSIONAL???
        !W = W + random_array(W)*2.e-10

        call initialize_scalar_arrays(T, WV)
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

        ! Create filename based on Tdiff and N
        file_format = '(a, a, a)'
        write(filename, file_format) trim(filename), '/', trim(simulation_name)
        call create_netcdf(trim(filename)//".nc", z, ncid)

        ! Writeout parameters
        write_time_iter = 0.

        ! Initialize buffers for writing to netCDF
        call initialize_buffers(write_buffer, N)
        if ( write_eddies ) call initialize_eddy_buffer(filename)

        ! Write out namelist file
        open(newunit=nml_unit, file=trim(filename)//'_nml.txt', action='write')
        write(nml_unit, nml=PARAMETERS)
        close(nml_unit)
        
        ! Length of a gridcell
        dz_length = H/N

    end subroutine initialize_params



    subroutine allocate_zero_arrays(A, n_array)
        integer(i4), intent(in) ::n_array
        real(dp), intent(inout), allocatable :: A(:)

        allocate(A(n_array))
        A = 0.

    end subroutine allocate_zero_arrays

    subroutine initialize_velocity_arrays(lw)
        ! Initializes velcoity arrays with some noise. This prevents a numerical issue from
        ! occuring when calculating the kinetic energy (i.e. KE != 0)
        real(dp), intent(inout), allocatable :: lw(:)
        real(dp) :: rand_num
        integer(i4) :: k

        do k = 1, N-1
            call random_number(rand_num)
            lw(k) = 2.e-10 * (rand_num - 0.5)
        end do

    end subroutine initialize_velocity_arrays


    subroutine initialize_scalar_arrays(lT, lWV)
        real(dp), intent(out) :: lT(:), lWV(:)
        integer(i4) :: k
        real(dp) :: lin_prof

        ! Non-Dimensional Values
        do concurrent (k = 1:N)
            lin_prof = 1.*k/N
            lT(k) = lin_prof
            lWV(k) = lin_prof
        end do

    end subroutine initialize_scalar_arrays

end module initialize