module writeout
    use globals
    use netcdf
    use droplets, only: particles, calculate_droplet_statistics, bin_droplet_radii, particle_bin_edges, &
                        size_distribution, n_aer_category, init_drop_each_gridpoint, &
                        expected_Ndrops_per_gridpoint, write_trajectories, trajectory_start, &
                        trajectory_end, trajectory_timer, initial_wet_radius, &
                        dsd_varid, aerDSD_varids
    use special_effects, only: do_sidewalls, do_random_fallout, area_sw, area_bot, C_sw, T_sw, &
                               RH_sw, P_sw, sw_nudging_time, random_fallout_rate
    implicit none

    ! Buffer variables and arrays for writing to netCDF
    integer(i4) :: buffer_size, buffer_count, nc_write_iter
    real(dp), allocatable :: buffer_T(:,:), buffer_WV(:,:), buffer_Tv(:,:) ! dims (buffer_size, N_grid)
    real(dp), allocatable :: buffer_SS(:,:), buffer_time(:), buffer_stats(:,:)
    real(dp), allocatable :: buffer_budgets(:,:) ! (n_budgets, buffer_size)
    integer(i4), allocatable :: buffer_DSD(:,:,:) ! n_DSDs, rbins, buffer_size

    ! Eddy output file unit (unformatted stream binary)
    integer(i4) :: eddy_unit


    ! Namelist metadata for global attributes (set by create_netcdf)
    character(100) :: nc_simulation_name
    integer(i4) :: nc_write_buffer

    ! Cached NetCDF variable IDs (populated by create_netcdf)
    integer :: varid_time, varid_T, varid_QV, varid_Tv, varid_S
    integer :: varid_stats(7)  ! Np, Nact, Nun, Ravg, LWC, N_collisions, N_coalescences
    integer :: varid_budgets(n_budgets)

    public  :: create_netcdf, initialize_buffers, initialize_particle_buffers, &
               add_to_profile_buffer, flush_buffer, close_netcdf, &
               write_profiles, write_eddy, initialize_eddy_file
    private :: write_netcdf_profiles
    
contains

    subroutine write_profiles(ldt)
        ! Accumulates the write timer and writes profile data when the interval is reached.
        real(dp), intent(in) :: ldt
        real(dp), save :: write_time_iter = 0.

        write_time_iter = write_time_iter + ldt
        if (write_time_iter >= write_timer) then
            write(*,*) 'Writing time: ', time
            if (do_microphysics) then
                call calculate_droplet_statistics(particles, statistics)
                call bin_droplet_radii(particles, particle_bin_edges, size_distribution)
            end if
            call add_to_profile_buffer(time, T, WV, Tv, SS, size_distribution, statistics)
            call reset_budgets()
            write_time_iter = mod(write_time_iter, write_timer)
        end if

    end subroutine write_profiles

    subroutine initialize_buffers(buff_len, N_grid)
        ! Initializes buffer size and arrays for writing to netCDF
        integer(i4), intent(in) :: buff_len, N_grid

        ! Set variables in module scope
        buffer_size = buff_len
        buffer_count = 0
        
        ! Allocate and initialize buffer arrays
        allocate(buffer_T(N_grid, buff_len))
        allocate(buffer_WV(N_grid, buff_len))
        allocate(buffer_Tv(N_grid, buff_len))
        allocate(buffer_SS(N_grid, buff_len))
        allocate(buffer_time(buff_len))
        allocate(buffer_stats(7, buff_len))
        allocate(buffer_budgets(n_budgets, buff_len))

        buffer_T = 0.
        buffer_WV = 0.
        buffer_Tv = 0.
        buffer_SS = 0.
        buffer_time = 0.
        buffer_stats = 0.
        buffer_budgets = 0.

    end subroutine initialize_buffers
    
    subroutine initialize_particle_buffers(n_cat, n_bins)

        integer(i4), intent(in) :: n_cat, n_bins
        integer(i4) :: n_DSDs

        ! Calculate number of DSDs (n_cat + 1 for multiple categories, 1 otherwise)
        
        if ( n_cat > 1 ) then
            n_DSDs = n_cat + 1
        else
            n_DSDs = 1
        end if
        allocate(buffer_DSD(n_DSDs, n_bins, buffer_size))
        buffer_DSD = 0

    end subroutine initialize_particle_buffers

    recursive subroutine add_to_profile_buffer(ltime, lT, lWV, lTv, lSS, lDSD, lstats)
        ! Writes the profile arrays into the profile buffers, will flush
        ! the buffer to write netCDF if buffer is full
        real(dp), intent(in) :: ltime, lT(:), lWV(:), lTv(:), lSS(:), lstats(:)
        integer(i4), intent(in) :: lDSD(:,:)

        ! Write profile data into buffers
        if (buffer_count < buffer_size) then
            buffer_count = buffer_count + 1 ! recall buffer_count resets to 0
            buffer_time(buffer_count) = ltime
            buffer_T(:, buffer_count) = lT - Tice ! (C)
            buffer_WV(:, buffer_count) = lWV
            buffer_Tv(:, buffer_count) = lTv - Tice ! (C)
            buffer_SS(:, buffer_count) = lSS
            buffer_stats(:, buffer_count) = lstats
            buffer_DSD(:, :, buffer_count) = lDSD
            buffer_budgets(1,  buffer_count) = budget_inject_solute_mass
            buffer_budgets(2,  buffer_count) = budget_inject_liquid_mass
            buffer_budgets(3,  buffer_count) = budget_fallout_liquid_mass
            buffer_budgets(4,  buffer_count) = budget_fallout_solute_mass
            buffer_budgets(5,  buffer_count) = budget_condensation
            buffer_budgets(6,  buffer_count) = budget_dgm_delta_T
            buffer_budgets(7,  buffer_count) = budget_diffusion_delta_T
            buffer_budgets(8,  buffer_count) = budget_diffusion_delta_WV
            buffer_budgets(9,  buffer_count) = budget_sidewall_delta_T
            buffer_budgets(10, buffer_count) = budget_sidewall_delta_WV
            buffer_budgets(11, buffer_count) = real(budget_n_injected, dp)
            buffer_budgets(12, buffer_count) = real(budget_n_fellout, dp)
        else
            ! Flush buffer and start new buffer
            call flush_buffer()
            call add_to_profile_buffer(ltime, lT, lWV, lTv, lSS, lDSD, lstats)
        end if

    end subroutine add_to_profile_buffer

    subroutine initialize_eddy_file(filename)
        character(*), intent(in) :: filename
        integer(i4) :: ierr

        open(newunit=eddy_unit, file=trim(filename)//'_eddies.bin', &
             form='unformatted', access='stream', status='replace', iostat=ierr)
        if (ierr /= 0) then
            print *, "Error opening eddy data file."
            stop 1
        end if

        ! Header: grid size, domain height, ODT constants, temperature BCs
        write(eddy_unit) N, H, C2, ZC2, Tdiff, Tref

    end subroutine initialize_eddy_file

    subroutine write_eddy(loc, len, ltime)
        integer(i4), intent(in) :: loc, len
        real(dp), intent(in) :: ltime

        write(eddy_unit) loc, len, ltime

    end subroutine write_eddy

    subroutine define_budget_var(lncid, t_dimid, idx, varname, long_name, units)
        integer, intent(in) :: lncid, t_dimid, idx
        character(*), intent(in) :: varname, long_name, units

        call nc_verify( nf90_def_var(lncid, varname, NF90_DOUBLE, t_dimid, varid_budgets(idx), &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: "//varname )
        call nc_verify( nf90_put_att(lncid, varid_budgets(idx), "long_name", long_name), "nf90_put_att: "//varname//", name" )
        call nc_verify( nf90_put_att(lncid, varid_budgets(idx), "units", units), "nf90_put_att: "//varname//", units" )

    end subroutine define_budget_var

    subroutine flush_buffer()
        ! writes the buffer to the netCDF file and resets buffer position
        call write_netcdf_profiles(ncid, buffer_time(1:buffer_count), &
                                        buffer_T(:, 1:buffer_count), &
                                        buffer_WV(:, 1:buffer_count), &
                                        buffer_Tv(:, 1:buffer_count), &
                                        buffer_SS(:, 1:buffer_count), &
                                        buffer_DSD(:, :, 1:buffer_count), &
                                        buffer_stats(:, 1:buffer_count), &
                                        buffer_budgets(:, 1:buffer_count))
        call nc_verify( nf90_sync(ncid) )
        buffer_count = 0

    end subroutine flush_buffer


    subroutine create_netcdf(file_name, z_m, lncid, sim_name, lwrite_buffer)
        character(*), intent(in):: file_name, sim_name
        real(dp), intent(in) :: z_m(:) !, scalar_vars(:) ! Establish z-dimension
        integer, intent(out) :: lncid
        integer(i4), intent(in) :: lwrite_buffer

        ! time/height dimensions
        logical :: file_exists
        integer :: old_nc, stat
        integer :: t_dimid, z_dimid, dimids(2)
        integer :: z_varid
        integer :: j, k, nz, nbins

        nz = size(z_m)
        nc_write_iter = 1

        ! Store namelist metadata for writing as global attributes at close time
        nc_simulation_name = sim_name
        nc_write_buffer = lwrite_buffer
        
        ! This subroutine creates/overwrites a netcdf file, defines all
        ! dimensions and variables to be written out from the simulation
        ! and exits define mode. Note, the netCDF is still open for writing.
    
        ! Test for existance, overwrite if present
        inquire(file=trim(file_name), exist=file_exists)
        if (file_exists) then
            write(*,'(a)') "NetCDF file already detected..."
            open(newunit=old_nc, file=file_name, status='old', iostat=stat)
            if (stat == 0) close(old_nc, status='delete')
            write(*,'(a)') "Deleted old netCDF file."
        end if
    
        ! Create initial netCDF file
        write(*,*) file_name
        call nc_verify( nf90_create(file_name, NF90_NETCDF4, lncid), "nf90_create" )
    
        ! Establish dimensions
        call nc_verify( nf90_def_dim(lncid, "time", NF90_UNLIMITED, t_dimid), "nf90_def_dim: time"  )
        call nc_verify( nf90_def_dim(lncid, "z", nz, z_dimid), "nf90_def_dim: z" )
    


        ! Establish variables and attributes for coordinate(dimension) variables
        call nc_verify( nf90_def_var(lncid, "z", NF90_FLOAT, z_dimid, z_varid, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: z" )
        call nc_verify( nf90_put_att(lncid, z_varid, "long_name", "Height"), "nf90_put_att: z, name" )
        call nc_verify( nf90_put_att(lncid, z_varid, "units", "meters"), "nf90_put_att: z, units" )
        call nc_verify( nf90_def_var(lncid, "time", NF90_FLOAT, t_dimid, varid_time, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: time" )
        call nc_verify( nf90_put_att(lncid, varid_time, "long_name", "Time"), "nf90_put_att: time, name" )
        call nc_verify( nf90_put_att(lncid, varid_time, "units", "seconds"), "nf90_put_att: time, units" )


        ! Note, netCDF will write out variables in (time, height) order
        dimids = (/ z_dimid, t_dimid /) ! Standard time/height dimensions
        call nc_verify( nf90_def_var(lncid, "T", NF90_FLOAT, dimids, varid_T, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: T" )
        call nc_verify( nf90_put_att(lncid, varid_T, "long_name", "Temperature"), "nf90_put_att: T, name" )
        call nc_verify( nf90_put_att(lncid, varid_T, "units", "celsius"), "nf90_put_att: T, units" )

        call nc_verify( nf90_def_var(lncid, "QV", NF90_FLOAT, dimids, varid_QV, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: QV" )
        call nc_verify( nf90_put_att(lncid, varid_QV, "long_name", "Water Vapor Mixing Ratio"), "nf90_put_att: QV, name" )
        call nc_verify( nf90_put_att(lncid, varid_QV, "units", "g/kg"), "nf90_put_att: QV, units")

        call nc_verify( nf90_def_var(lncid, "Tv", NF90_FLOAT, dimids, varid_Tv, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: Tv" )
        call nc_verify( nf90_put_att(lncid, varid_Tv, "long_name", "Virtual Temperature"), "nf90_put_att: Tv, name" )
        call nc_verify( nf90_put_att(lncid, varid_Tv, "units", "celsius"), "nf90_put_att: Tv, units")

        call nc_verify( nf90_def_var(lncid, "S", NF90_FLOAT, dimids, varid_S, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: S" )
        call nc_verify( nf90_put_att(lncid, varid_S, "long_name", "Supersaturation"), "nf90_put_att: S, name" )
        call nc_verify( nf90_put_att(lncid, varid_S, "units", "%"), "nf90_put_att: S, units")



        ! Particle statistic variables (only when microphysics is enabled)
        if ( do_microphysics ) then
            call nc_verify( nf90_def_var(lncid, "Np", NF90_INT, t_dimid, varid_stats(1), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Np" )
            call nc_verify( nf90_put_att(lncid, varid_stats(1), "long_name", "Number of Particles"), "nf90_put_att: Np, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(1), "units", "#"), "nf90_put_att: Np units")

            call nc_verify( nf90_def_var(lncid, "Nact", NF90_INT, t_dimid, varid_stats(2), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Nact" )
            call nc_verify( nf90_put_att(lncid, varid_stats(2), "long_name", "Number of Activated Particles"), "nf90_put_att: Nact, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(2), "units", "#"), "nf90_put_att: Nact units")

            call nc_verify( nf90_def_var(lncid, "Nun", NF90_INT, t_dimid, varid_stats(3), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Nun" )
            call nc_verify( nf90_put_att(lncid, varid_stats(3), "long_name", "Number of Unactivated Particles"), "nf90_put_att: Nun, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(3), "units", "#"), "nf90_put_att: Nun units")

            call nc_verify( nf90_def_var(lncid, "Ravg", NF90_FLOAT, t_dimid, varid_stats(4), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Ravg" )
            call nc_verify( nf90_put_att(lncid, varid_stats(4), "long_name", "Average Particle Radius (wet)"), "nf90_put_att: Ravg, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(4), "units", "um"), "nf90_put_att: Ravg units")

            call nc_verify( nf90_def_var(lncid, "LWC", NF90_FLOAT, t_dimid, varid_stats(5), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: LWC" )
            call nc_verify( nf90_put_att(lncid, varid_stats(5), "long_name", "Liquid Water Content"), "nf90_put_att: LWC, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(5), "units", "g/m3"), "nf90_put_att: LWC units")

            call nc_verify( nf90_def_var(lncid, "N_collisions", NF90_INT, t_dimid, varid_stats(6), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: N_collisions" )
            call nc_verify( nf90_put_att(lncid, varid_stats(6), "long_name", "Number of Collisions Since Last Write"), &
                            "nf90_put_att: N_collisions, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(6), "units", "#"), "nf90_put_att: N_collisions units")

            call nc_verify( nf90_def_var(lncid, "N_coalescences", NF90_INT, t_dimid, varid_stats(7), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: N_coalescences" )
            call nc_verify( nf90_put_att(lncid, varid_stats(7), "long_name", "Number of Coalescences Since Last Write"), &
                            "nf90_put_att: N_coalescences, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(7), "units", "#"), "nf90_put_att: N_coalescences units")
        end if

        ! Budget variables (time dimension only, defined unconditionally)
        call define_budget_var(lncid, t_dimid, 1,  "budget_inject_solute_mass",  "Accumulated injected solute mass",          "kg")
        call define_budget_var(lncid, t_dimid, 2,  "budget_inject_liquid_mass",  "Accumulated injected liquid water mass",    "kg")
        call define_budget_var(lncid, t_dimid, 3,  "budget_fallout_liquid_mass", "Accumulated liquid water removed by fallout","kg")
        call define_budget_var(lncid, t_dimid, 4,  "budget_fallout_solute_mass", "Accumulated solute mass removed by fallout","kg")
        call define_budget_var(lncid, t_dimid, 5,  "budget_condensation",        "Net liquid water change from cond/evap",    "kg")
        call define_budget_var(lncid, t_dimid, 6,  "budget_dgm_delta_T",         "Sum of per-droplet temperature changes from DGM", "K")
        call define_budget_var(lncid, t_dimid, 7,  "budget_diffusion_delta_T",   "Domain-sum T change from diffusion",        "K")
        call define_budget_var(lncid, t_dimid, 8,  "budget_diffusion_delta_WV",  "Domain-sum WV change from diffusion",       "kg/kg")
        call define_budget_var(lncid, t_dimid, 9,  "budget_sidewall_delta_T",    "Domain-sum T change from sidewall nudging", "K")
        call define_budget_var(lncid, t_dimid, 10, "budget_sidewall_delta_WV",   "Domain-sum WV change from sidewall nudging","kg/kg")
        call define_budget_var(lncid, t_dimid, 11, "budget_n_injected",          "Number of particles injected",              "#")
        call define_budget_var(lncid, t_dimid, 12, "budget_n_fellout",           "Number of particles fallen out",            "#")

        ! Exit define mode, however netCDF is still open
        call nc_verify( nf90_enddef(lncid), "nf90_enddef" )
    
        ! Fill variables for known dimensions
        call nc_verify( nf90_put_var(lncid, z_varid, z_m), "nf90_put_var: z" )

        write(*,'(a)') "NetCDF file fully created."
    
    end subroutine create_netcdf


    subroutine write_netcdf_profiles(lncid, ltime, lT, lWV, lTv, lSS, lDSD, lstats, lbudgets)
        integer(i4), intent(in) :: lncid
        real(dp), intent(in) :: ltime(:), lT(:,:), lWV(:,:), lTv(:,:), lSS(:,:), lstats(:,:)
        real(dp), intent(in) :: lbudgets(:,:)
        integer(i4), intent(in) :: lDSD(:,:,:)

        integer :: i, time_len, z_len, bin_len
        integer :: count_dim(2), start_dim(2)

        ! netCDF profile variables are (time, height) dimensions
        time_len = size(ltime,1)
        z_len = size(lT,1)
        bin_len = size(lDSD,2)
        if (size(lT,2) /= time_len) then
            write(*,'(a)') "Error: Time and profile arrays are not the same length."
            stop 1
        end if
        ! set the starting/ending positions for netCDF file
        count_dim = (/z_len, time_len/)
        start_dim = (/1, nc_write_iter/)

        ! Write profile variables using cached varids
        call nc_verify( nf90_put_var(lncid, varid_time, ltime, start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(lncid, varid_T,  lT,  start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_QV, lWV, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_Tv, lTv, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_S,  lSS, start=start_dim, count=count_dim) )

        if ( do_microphysics ) then
            ! Particle statistics
            do i = 1, 7
                call nc_verify( nf90_put_var(lncid, varid_stats(i), lstats(i,:), start=(/nc_write_iter/)) )
            end do

            ! Write total DSD
            count_dim = (/ bin_len, time_len /)
            call nc_verify( nf90_put_var(lncid, dsd_varid, lDSD(1,:,:), start=start_dim, count=count_dim) )

            ! Write the subcategory DSDs
            if ( n_aer_category > 1 ) then
                do i = 1, n_aer_category
                    call nc_verify( nf90_put_var(lncid, aerDSD_varids(i), lDSD(i+1,:,:), &
                                    start=start_dim, count=count_dim) )
                end do
            end if
        end if

        ! Write budget time series using cached varids
        do i = 1, n_budgets
            call nc_verify( nf90_put_var(lncid, varid_budgets(i), lbudgets(i,:), &
                            start=(/nc_write_iter/)) )
        end do

        ! Move 'start' time location to end of buffer for next write
        nc_write_iter = nc_write_iter + buffer_count

    end subroutine write_netcdf_profiles


    subroutine close_netcdf(lncid)
        integer, intent(in) :: lncid

        call write_namelist_attributes(lncid, nc_simulation_name, nc_write_buffer)
        call nc_verify( nf90_close(lncid), 'nf90_close')
        if ( write_eddies ) close(eddy_unit)

    end subroutine




    

    subroutine write_namelist_attributes(lncid, sim_name, lwrite_buffer)
        ! Writes all namelist parameters as global attributes to the netCDF file.
        ! Re-enters define mode, writes attributes, then exits define mode.
        integer, intent(in) :: lncid
        character(*), intent(in) :: sim_name
        integer(i4), intent(in) :: lwrite_buffer

        call nc_verify( nf90_redef(lncid), "nf90_redef: namelist attributes" )

        ! PARAMETERS namelist (19 attributes)
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.simulation_name", trim(sim_name)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.N", N) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Lmin", Lmin) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Lprob", Lprob) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.tmax", tmax) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Tdiff", Tdiff) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Tref", Tref) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.pres", pres) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.H", H) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.volume_scaling", volume_scaling) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.max_accept_prob", max_accept_prob) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_timer", write_timer) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_buffer", lwrite_buffer) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.same_random", merge(1, 0, same_random)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_turbulence", merge(1, 0, do_turbulence)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_microphysics", merge(1, 0, do_microphysics)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_eddies", merge(1, 0, write_eddies)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_special_effects", merge(1, 0, do_special_effects)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.overwrite", merge(1, 0, overwrite)) )

        ! MICROPHYSICS namelist (7 attributes, only when microphysics is enabled)
        if ( do_microphysics ) then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.init_drop_each_gridpoint", &
                            merge(1, 0, init_drop_each_gridpoint)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.expected_Ndrops_per_gridpoint", &
                            expected_Ndrops_per_gridpoint) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.write_trajectories", &
                            merge(1, 0, write_trajectories)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_start", trajectory_start) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_end", trajectory_end) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_timer", trajectory_timer) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.initial_wet_radius", initial_wet_radius) )
        end if

        ! SPECIALEFFECTS namelist (10 attributes, only when special effects is enabled)
        if ( do_special_effects ) then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.do_sidewalls", &
                            merge(1, 0, do_sidewalls)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.do_random_fallout", &
                            merge(1, 0, do_random_fallout)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.area_sw", area_sw) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.area_bot", area_bot) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.C_sw", C_sw) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.T_sw", T_sw) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.RH_sw", RH_sw) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.P_sw", P_sw) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.sw_nudging_time", sw_nudging_time) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.random_fallout_rate", random_fallout_rate) )
        end if

        call nc_verify( nf90_enddef(lncid), "nf90_enddef: namelist attributes" )

    end subroutine write_namelist_attributes

end module writeout