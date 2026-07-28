! Buffered NetCDF output of the main results file ({name}.nc): vertical profiles
! and time series (and parcel diagnostics in parcel mode). Accumulates writes in
! a buffer and flushes periodically to limit I/O overhead. Embeds code_version
! and git_commit as global attributes; the conventions string is the format
! contract with the CODT_tools readers.
module writeout
    use globals
    use version, only: code_version, git_commit
    use netcdf
    use droplets, only: particles, calculate_droplet_statistics, bin_droplet_radii, particle_bin_edges, &
                        size_distribution, n_aer_category, n_DSD_bins, init_drop_each_gridpoint, &
                        expected_Ndrops_per_gridpoint, write_trajectories, trajectory_start, &
                        trajectory_end, trajectory_timer, initial_wet_radius, &
                        dsd_varid, aerDSD_varids, aerosol_file, aerosol_concentration, &
                        do_collisions, do_coalescence, wmax_collision, write_collisions, &
                        coalescence_kernel, do_seeding, seed_hydration, seed_growth_time
    use special_effects, only: do_sidewalls, do_random_fallout, area_sw, area_bot, C_sw, T_sw, &
                               RH_sw, P_sw, sw_nudging_time, random_fallout_rate
    use parcel, only: do_parcel_ascent, parcel_height, parcel_velocity, &
                      parcel_file, initial_RH, pressure_limit, &
                      pressure_mode, parcel_height_env, write_height_env, &
                      initial_height, vertical_axis
    use entrainment, only: ent_rate, n_blob, psigma, random_entrainment, &
                           time_varying_entrainment
    use radiation, only: rad_F_net, rad_heating_rate, radiation_method, mie_data_file, &
                         eps_top, eps_bot, sky_temp, sky_cooling_flag, rad_call_interval, &
                         nPhotons, nBins, Lx_rad, Ly_rad, T_side, max_droplets_per_cell
    use ODT, only: Tdiff, Lmin, Lprob, max_accept_prob, C2, ZC2
    use LEM, only: integral_length_scale, dissipation_rate, &
                   actual_kolmogorov_scale, grid_eddy_scale, diffusivity_length_scale, &
                   smallest_eddy_scale, smallest_eddy_gridpoints, diffusivity_enhancement
    implicit none

    private

    ! Buffer variables and arrays for writing to netCDF
    integer(i4) :: buffer_size, buffer_count, nc_write_iter
    real(dp), allocatable :: buffer_T(:,:), buffer_WV(:,:), buffer_Tv(:,:) ! dims (buffer_size, N_grid)
    real(dp), allocatable :: buffer_SS(:,:), buffer_time(:), buffer_stats(:,:)
    real(dp), allocatable :: buffer_field_budgets(:,:)  ! (n_field_budgets, buffer_size)
    real(dp), allocatable :: buffer_micro_budgets(:,:)  ! (n_micro_budgets, buffer_size)
    integer(i4), allocatable :: buffer_DSD(:,:,:) ! n_DSDs, rbins, buffer_size

    ! Eddy output file unit (unformatted stream binary)
    integer(i4) :: eddy_unit



    ! Namelist metadata for global attributes (set by create_netcdf)
    character(100) :: nc_simulation_name
    integer(i4) :: nc_write_buffer

    ! Cached NetCDF variable IDs (populated by create_netcdf)
    integer :: varid_time, varid_T, varid_QV, varid_Tv, varid_S
    integer :: varid_stats(7)  ! Np, Nact, Nun, Ravg, LWC, N_collisions, N_coalescences
    integer :: varid_field_budgets(n_field_budgets)
    integer :: varid_micro_budgets(n_micro_budgets)
    integer :: varid_entrain_budgets(n_entrain_budgets)
    integer :: varid_parcel_height, varid_parcel_pressure, varid_parcel_velocity
    integer :: varid_parcel_height_env
    real(dp), allocatable :: buffer_parcel_height_env(:)

    ! Parcel ascent buffers
    real(dp), allocatable :: buffer_parcel_height(:), buffer_parcel_pressure(:), buffer_parcel_velocity(:)

    ! Entrainment budget varids and buffers (only when do_entrainment = .true.)
    real(dp), allocatable :: buffer_entrain_budgets(:,:)  ! (n_entrain_budgets, buffer_size)

    ! Time-varying entrainment parameter series (only when the schedule is active)
    integer :: varid_ent_rate, varid_n_blob, varid_psigma
    real(dp), allocatable :: buffer_ent_rate(:), buffer_psigma(:)
    integer(i4), allocatable :: buffer_n_blob(:)

    ! Radiation varids and buffers (only when do_radiation = .true.)
    integer :: varid_rad_F_net, varid_rad_heating_rate, varid_rad_budget
    real(dp), allocatable :: buffer_rad_F_net(:,:), buffer_rad_heating_rate(:,:)
    real(dp), allocatable :: buffer_rad_budget(:)

    public  :: create_netcdf, initialize_buffers, deallocate_buffers, &
               add_to_profile_buffer, flush_buffer, close_netcdf, &
               write_profiles, write_eddy, initialize_eddy_file, write_eddy_header_fields
    private :: write_netcdf_profiles
    
contains

    subroutine write_profiles(ldt)
        ! Accumulates the write timer and writes profile data when the interval is reached.
        real(dp), intent(in) :: ldt
        real(dp), save :: write_time_iter = 0.

        write_time_iter = write_time_iter + ldt
        if (write_time_iter >= write_timer) then
            if (do_parcel_ascent) then
                write(*,'(a,f10.2,a,f8.1,a)') ' Writing time: ', time, '  P: ', pres / 100.0, ' mb'
            else
                write(*,*) 'Writing time: ', time
            end if
            if (do_microphysics) then
                call calculate_droplet_statistics(particles, statistics)
                call bin_droplet_radii(particles, particle_bin_edges, size_distribution)
            end if
            call add_to_profile_buffer(time, T, WV, Tv, SS)
            call reset_budgets()
            write_time_iter = mod(write_time_iter, write_timer)
        end if

    end subroutine write_profiles

    subroutine initialize_buffers(buff_len, N_grid)
        integer(i4), intent(in) :: buff_len, N_grid
        integer(i4) :: n_DSDs

        buffer_size = buff_len
        buffer_count = 0

        allocate(buffer_T(N_grid, buff_len))
        allocate(buffer_WV(N_grid, buff_len))
        allocate(buffer_Tv(N_grid, buff_len))
        allocate(buffer_SS(N_grid, buff_len))
        allocate(buffer_time(buff_len))
        allocate(buffer_field_budgets(n_field_budgets, buff_len))

        buffer_T = 0.
        buffer_WV = 0.
        buffer_Tv = 0.
        buffer_SS = 0.
        buffer_time = 0.
        buffer_field_budgets = 0.

        if (do_microphysics) then
            allocate(buffer_stats(7, buff_len))
            buffer_stats = 0.
            allocate(buffer_micro_budgets(n_micro_budgets, buff_len))
            buffer_micro_budgets = 0.
        end if

        if (do_microphysics) then
            if (n_aer_category > 1) then
                n_DSDs = n_aer_category + 1
            else
                n_DSDs = 1
            end if
            allocate(buffer_DSD(n_DSDs, n_DSD_bins, buff_len))
            buffer_DSD = 0
        end if

        if (do_parcel_ascent) then
            allocate(buffer_parcel_height(buff_len))
            allocate(buffer_parcel_pressure(buff_len))
            allocate(buffer_parcel_velocity(buff_len))
            buffer_parcel_height = 0.
            buffer_parcel_pressure = 0.
            buffer_parcel_velocity = 0.
            if (write_height_env) then
                allocate(buffer_parcel_height_env(buff_len))
                buffer_parcel_height_env = 0.
            end if
        end if

        if (do_entrainment) then
            allocate(buffer_entrain_budgets(n_entrain_budgets, buff_len))
            buffer_entrain_budgets = 0.
        end if

        if (do_entrainment .and. time_varying_entrainment) then
            allocate(buffer_ent_rate(buff_len))
            allocate(buffer_n_blob(buff_len))
            allocate(buffer_psigma(buff_len))
            buffer_ent_rate = 0.
            buffer_n_blob = 0
            buffer_psigma = 0.
        end if

        if (do_radiation) then
            allocate(buffer_rad_F_net(N_grid, buff_len))
            allocate(buffer_rad_heating_rate(N_grid, buff_len))
            allocate(buffer_rad_budget(buff_len))
            buffer_rad_F_net = 0.
            buffer_rad_heating_rate = 0.
            buffer_rad_budget = 0.
        end if

    end subroutine initialize_buffers


    subroutine deallocate_buffers()
        deallocate(buffer_T, buffer_WV, buffer_Tv, buffer_SS, buffer_time, buffer_field_budgets)
        if (allocated(buffer_stats)) deallocate(buffer_stats)
        if (allocated(buffer_micro_budgets)) deallocate(buffer_micro_budgets)
        if (allocated(buffer_DSD)) deallocate(buffer_DSD)
        if (allocated(buffer_parcel_height)) deallocate(buffer_parcel_height)
        if (allocated(buffer_parcel_pressure)) deallocate(buffer_parcel_pressure)
        if (allocated(buffer_parcel_velocity)) deallocate(buffer_parcel_velocity)
        if (allocated(buffer_rad_F_net)) deallocate(buffer_rad_F_net)
        if (allocated(buffer_rad_heating_rate)) deallocate(buffer_rad_heating_rate)
        if (allocated(buffer_rad_budget)) deallocate(buffer_rad_budget)
        if (allocated(buffer_entrain_budgets)) deallocate(buffer_entrain_budgets)
        if (allocated(buffer_parcel_height_env)) deallocate(buffer_parcel_height_env)
        if (allocated(buffer_ent_rate)) deallocate(buffer_ent_rate)
        if (allocated(buffer_n_blob)) deallocate(buffer_n_blob)
        if (allocated(buffer_psigma)) deallocate(buffer_psigma)
    end subroutine deallocate_buffers


    recursive subroutine add_to_profile_buffer(ltime, lT, lWV, lTv, lSS)
        real(dp), intent(in) :: ltime, lT(:), lWV(:), lTv(:), lSS(:)

        if (buffer_count < buffer_size) then
            buffer_count = buffer_count + 1
            buffer_time(buffer_count) = ltime
            buffer_T(:, buffer_count) = lT - Tice ! (C)
            buffer_WV(:, buffer_count) = lWV
            buffer_Tv(:, buffer_count) = lTv - Tice ! (C)
            buffer_SS(:, buffer_count) = lSS
            buffer_field_budgets(1, buffer_count) = budget_diffusion_delta_T
            buffer_field_budgets(2, buffer_count) = budget_diffusion_delta_WV
            buffer_field_budgets(3, buffer_count) = budget_sidewall_delta_T
            buffer_field_budgets(4, buffer_count) = budget_sidewall_delta_WV
            if (do_microphysics) then
                buffer_stats(:, buffer_count) = statistics
                buffer_DSD(:, :, buffer_count) = size_distribution
                buffer_micro_budgets(1, buffer_count) = budget_inject_solute_mass
                buffer_micro_budgets(2, buffer_count) = budget_inject_liquid_mass
                buffer_micro_budgets(3, buffer_count) = budget_fallout_liquid_mass
                buffer_micro_budgets(4, buffer_count) = budget_fallout_solute_mass
                buffer_micro_budgets(5, buffer_count) = budget_condensation
                buffer_micro_budgets(6, buffer_count) = budget_dgm_delta_T
                buffer_micro_budgets(7, buffer_count) = real(budget_n_injected, dp)
                buffer_micro_budgets(8, buffer_count) = real(budget_n_fellout, dp)
                buffer_micro_budgets(9, buffer_count) = real(budget_n_coalesced, dp)
            end if
            if (do_parcel_ascent) then
                buffer_parcel_height(buffer_count) = parcel_height
                buffer_parcel_pressure(buffer_count) = pres / Pa_per_mb
                buffer_parcel_velocity(buffer_count) = parcel_velocity
                if (write_height_env) then
                    buffer_parcel_height_env(buffer_count) = parcel_height_env
                end if
            end if
            if (do_entrainment) then
                buffer_entrain_budgets(1, buffer_count) = budget_detrain_liquid_mass
                buffer_entrain_budgets(2, buffer_count) = budget_detrain_solute_mass
                buffer_entrain_budgets(3, buffer_count) = budget_entrain_liquid_mass
                buffer_entrain_budgets(4, buffer_count) = budget_entrain_solute_mass
                buffer_entrain_budgets(5, buffer_count) = real(budget_n_detrained, dp)
                buffer_entrain_budgets(6, buffer_count) = real(budget_n_entrained, dp)
            end if
            if (do_entrainment .and. time_varying_entrainment) then
                buffer_ent_rate(buffer_count) = ent_rate * m_per_km   ! 1/m -> 1/km
                buffer_n_blob(buffer_count) = n_blob
                buffer_psigma(buffer_count) = psigma
            end if
            if (do_radiation) then
                buffer_rad_F_net(:, buffer_count) = rad_F_net
                buffer_rad_heating_rate(:, buffer_count) = rad_heating_rate
                buffer_rad_budget(buffer_count) = budget_radiation_delta_T
            end if
        else
            ! Flush buffer and start new buffer
            call flush_buffer()
            call add_to_profile_buffer(ltime, lT, lWV, lTv, lSS)
        end if

    end subroutine add_to_profile_buffer

    subroutine initialize_eddy_file(filename, mode)
        character(*), intent(in) :: filename, mode
        integer(i4) :: ierr
        integer(i1) :: mode_flag

        open(newunit=eddy_unit, file=trim(filename)//'_eddies.bin', &
             form='unformatted', access='stream', status='replace', iostat=ierr)
        if (ierr /= 0) then
            write(error_unit,*) 'Error opening eddy data file.'
            call exit(1)
        end if

        if (trim(mode) == 'chamber') then
            mode_flag = 0_i1
        else
            mode_flag = 1_i1
        end if
        write(eddy_unit) mode_flag
        write(eddy_unit) N, H

    end subroutine initialize_eddy_file


    subroutine write_eddy_header_fields(fields)
        real(dp), intent(in) :: fields(:)

        write(eddy_unit) fields

    end subroutine write_eddy_header_fields

    subroutine write_eddy(loc, len, ltime)
        integer(i4), intent(in) :: loc, len
        real(dp), intent(in) :: ltime

        write(eddy_unit) loc, len, ltime

    end subroutine write_eddy

    subroutine define_budget_var(lncid, t_dimid, varids, idx, varname, long_name, units)
        integer, intent(in) :: lncid, t_dimid, idx
        integer, intent(inout) :: varids(:)
        character(*), intent(in) :: varname, long_name, units

        call nc_verify( nf90_def_var(lncid, varname, NF90_DOUBLE, t_dimid, varids(idx), &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: "//varname )
        call nc_verify( nf90_put_att(lncid, varids(idx), "long_name", long_name), "nf90_put_att: "//varname//", name" )
        call nc_verify( nf90_put_att(lncid, varids(idx), "units", units), "nf90_put_att: "//varname//", units" )

    end subroutine define_budget_var

    subroutine flush_buffer()
        call write_netcdf_profiles(ncid, buffer_time(1:buffer_count), &
                                        buffer_T(:, 1:buffer_count), &
                                        buffer_WV(:, 1:buffer_count), &
                                        buffer_Tv(:, 1:buffer_count), &
                                        buffer_SS(:, 1:buffer_count))
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
        integer :: j, k, nz

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
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "conventions", "CODT_output_v1"), &
                        "nf90_put_att: conventions" )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "code_version", code_version), &
                        "nf90_put_att: code_version" )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "git_commit", git_commit), &
                        "nf90_put_att: git_commit" )

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
        call nc_verify( nf90_put_att(lncid, varid_QV, "units", "kg/kg"), "nf90_put_att: QV, units")

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
            call nc_verify( nf90_put_att(lncid, varid_stats(2), "long_name", &
                            "Number of Activated Particles"), "nf90_put_att: Nact, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(2), "units", "#"), "nf90_put_att: Nact units")

            call nc_verify( nf90_def_var(lncid, "Nun", NF90_INT, t_dimid, varid_stats(3), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Nun" )
            call nc_verify( nf90_put_att(lncid, varid_stats(3), "long_name", &
                            "Number of Unactivated Particles"), "nf90_put_att: Nun, name" )
            call nc_verify( nf90_put_att(lncid, varid_stats(3), "units", "#"), "nf90_put_att: Nun units")

            call nc_verify( nf90_def_var(lncid, "Ravg", NF90_FLOAT, t_dimid, varid_stats(4), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: Ravg" )
            call nc_verify( nf90_put_att(lncid, varid_stats(4), "long_name", &
                            "Average Particle Radius (wet)"), "nf90_put_att: Ravg, name" )
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

        ! Field budget variables (always active)
        call define_budget_var(lncid, t_dimid, varid_field_budgets, 1, "budget_diffusion_delta_T", &
                              "Domain-sum T change from diffusion", "K")
        call define_budget_var(lncid, t_dimid, varid_field_budgets, 2, "budget_diffusion_delta_WV", &
                              "Domain-sum WV change from diffusion", "kg/kg")
        call define_budget_var(lncid, t_dimid, varid_field_budgets, 3, "budget_sidewall_delta_T", &
                              "Domain-sum T change from sidewall nudging", "K")
        call define_budget_var(lncid, t_dimid, varid_field_budgets, 4, "budget_sidewall_delta_WV", &
                              "Domain-sum WV change from sidewall nudging", "kg/kg")

        ! Microphysics budget variables
        if (do_microphysics) then
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 1, &
                                  "budget_inject_solute_mass", "Accumulated injected solute mass", "kg")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 2, &
                                  "budget_inject_liquid_mass", "Accumulated injected liquid water mass", "kg")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 3, &
                                  "budget_fallout_liquid_mass", "Accumulated liquid water removed by fallout", "kg")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 4, &
                                  "budget_fallout_solute_mass", "Accumulated solute mass removed by fallout", "kg")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 5, &
                                  "budget_condensation", "Net liquid water change from cond/evap", "kg")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 6, &
                                  "budget_dgm_delta_T", "Sum of per-droplet temperature changes from DGM", "K")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 7, &
                                  "budget_n_injected", "Number of particles injected", "#")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 8, &
                                  "budget_n_fellout", "Number of particles fallen out", "#")
            call define_budget_var(lncid, t_dimid, varid_micro_budgets, 9, &
                                  "budget_n_coalesced", "Number of particles removed by coalescence", "#")
        end if

        ! Entrainment budget variables
        if (do_entrainment) then
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 1, &
                                  "budget_detrain_liquid_mass", "Liquid water removed by detrainment", "kg")
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 2, &
                                  "budget_detrain_solute_mass", "Solute mass removed by detrainment", "kg")
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 3, &
                                  "budget_entrain_liquid_mass", "Liquid water added by entrainment", "kg")
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 4, &
                                  "budget_entrain_solute_mass", "Solute mass added by entrainment", "kg")
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 5, &
                                  "budget_n_detrained", "Number of particles removed by detrainment", "#")
            call define_budget_var(lncid, t_dimid, varid_entrain_budgets, 6, &
                                  "budget_n_entrained", "Number of particles added by entrainment", "#")
        end if

        ! Parcel ascent variables (parcel mode with dynamics only)
        if (do_parcel_ascent) then
            call nc_verify( nf90_def_var(lncid, "parcel_height", NF90_FLOAT, t_dimid, &
                            varid_parcel_height, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: parcel_height" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_height, "long_name", "Parcel Height"), &
                            "nf90_put_att: parcel_height, name" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_height, "units", "m"), &
                            "nf90_put_att: parcel_height, units" )

            call nc_verify( nf90_def_var(lncid, "parcel_pressure", NF90_FLOAT, t_dimid, &
                            varid_parcel_pressure, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: parcel_pressure" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_pressure, "long_name", "Parcel Pressure"), &
                            "nf90_put_att: parcel_pressure, name" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_pressure, "units", "mb"), &
                            "nf90_put_att: parcel_pressure, units" )

            call nc_verify( nf90_def_var(lncid, "parcel_velocity", NF90_FLOAT, t_dimid, &
                            varid_parcel_velocity, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: parcel_velocity" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_velocity, "long_name", "Parcel Vertical Velocity"), &
                            "nf90_put_att: parcel_velocity, name" )
            call nc_verify( nf90_put_att(lncid, varid_parcel_velocity, "units", "m/s"), &
                            "nf90_put_att: parcel_velocity, units" )

            ! Diagnostic (v3 + hydrostatic pressure_mode): environment height at
            ! the parcel's pressure; drift from parcel_height measures how far the
            ! self-integrated pressure has left the sounding's p(z) curve.
            if (write_height_env) then
                call nc_verify( nf90_def_var(lncid, "parcel_height_env", NF90_FLOAT, t_dimid, &
                                varid_parcel_height_env, deflate_level=1, shuffle=.true.), &
                                "nf90_def_var: parcel_height_env" )
                call nc_verify( nf90_put_att(lncid, varid_parcel_height_env, "long_name", &
                                "Environment Height at Parcel Pressure"), &
                                "nf90_put_att: parcel_height_env, name" )
                call nc_verify( nf90_put_att(lncid, varid_parcel_height_env, "units", "m"), &
                                "nf90_put_att: parcel_height_env, units" )
            end if
        end if

        ! Time-varying entrainment parameter series (active value each output step)
        if (do_entrainment .and. time_varying_entrainment) then
            call nc_verify( nf90_def_var(lncid, "ent_rate", NF90_FLOAT, t_dimid, &
                            varid_ent_rate, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: ent_rate" )
            call nc_verify( nf90_put_att(lncid, varid_ent_rate, "long_name", &
                            "Entrainment Rate"), "nf90_put_att: ent_rate, name" )
            call nc_verify( nf90_put_att(lncid, varid_ent_rate, "units", "1/km"), &
                            "nf90_put_att: ent_rate, units" )

            call nc_verify( nf90_def_var(lncid, "n_blob", NF90_INT, t_dimid, &
                            varid_n_blob, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: n_blob" )
            call nc_verify( nf90_put_att(lncid, varid_n_blob, "long_name", &
                            "Number of Blobs per Entrainment Event"), "nf90_put_att: n_blob, name" )
            call nc_verify( nf90_put_att(lncid, varid_n_blob, "units", "1"), &
                            "nf90_put_att: n_blob, units" )

            call nc_verify( nf90_def_var(lncid, "psigma", NF90_FLOAT, t_dimid, &
                            varid_psigma, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: psigma" )
            call nc_verify( nf90_put_att(lncid, varid_psigma, "long_name", &
                            "Blob Fraction of Domain"), "nf90_put_att: psigma, name" )
            call nc_verify( nf90_put_att(lncid, varid_psigma, "units", "1"), &
                            "nf90_put_att: psigma, units" )
        end if

        ! Radiation variables (chamber mode with do_radiation only)
        if (do_radiation) then
            call nc_verify( nf90_def_var(lncid, "rad_F_net", NF90_FLOAT, dimids, varid_rad_F_net, &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: rad_F_net" )
            call nc_verify( nf90_put_att(lncid, varid_rad_F_net, "long_name", "Net Radiative Flux"), &
                            "nf90_put_att: rad_F_net, name" )
            call nc_verify( nf90_put_att(lncid, varid_rad_F_net, "units", "W/m2"), &
                            "nf90_put_att: rad_F_net, units" )

            call nc_verify( nf90_def_var(lncid, "rad_heating_rate", NF90_FLOAT, dimids, &
                            varid_rad_heating_rate, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: rad_heating_rate" )
            call nc_verify( nf90_put_att(lncid, varid_rad_heating_rate, "long_name", &
                            "Radiative Heating Rate"), "nf90_put_att: rad_heating_rate, name" )
            call nc_verify( nf90_put_att(lncid, varid_rad_heating_rate, "units", "K/s"), &
                            "nf90_put_att: rad_heating_rate, units" )

            call nc_verify( nf90_def_var(lncid, "budget_radiation_delta_T", NF90_DOUBLE, &
                            t_dimid, varid_rad_budget, deflate_level=1, shuffle=.true.), &
                            "nf90_def_var: budget_radiation_delta_T" )
            call nc_verify( nf90_put_att(lncid, varid_rad_budget, "long_name", &
                            "Domain-sum T change from radiation"), &
                            "nf90_put_att: budget_radiation_delta_T, name" )
            call nc_verify( nf90_put_att(lncid, varid_rad_budget, "units", "K"), &
                            "nf90_put_att: budget_radiation_delta_T, units" )
        end if

        ! Exit define mode, however netCDF is still open
        call nc_verify( nf90_enddef(lncid), "nf90_enddef" )
    
        ! Fill variables for known dimensions
        call nc_verify( nf90_put_var(lncid, z_varid, z_m), "nf90_put_var: z" )

        write(*,'(a)') "NetCDF file fully created."
    
    end subroutine create_netcdf


    subroutine write_netcdf_profiles(lncid, ltime, lT, lWV, lTv, lSS)
        integer(i4), intent(in) :: lncid
        real(dp), intent(in) :: ltime(:), lT(:,:), lWV(:,:), lTv(:,:), lSS(:,:)

        integer :: i, time_len, z_len, bin_len
        integer :: count_dim(2), start_dim(2)

        time_len = size(ltime,1)
        z_len = size(lT,1)

        count_dim = (/z_len, time_len/)
        start_dim = (/1, nc_write_iter/)

        call nc_verify( nf90_put_var(lncid, varid_time, ltime, start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(lncid, varid_T,  lT,  start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_QV, lWV, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_Tv, lTv, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(lncid, varid_S,  lSS, start=start_dim, count=count_dim) )

        if ( do_microphysics ) then
            do i = 1, 7
                call nc_verify( nf90_put_var(lncid, varid_stats(i), buffer_stats(i, 1:buffer_count), &
                                start=(/nc_write_iter/)) )
            end do

            bin_len = size(buffer_DSD, 2)
            count_dim = (/ bin_len, time_len /)
            call nc_verify( nf90_put_var(lncid, dsd_varid, buffer_DSD(1, :, 1:buffer_count), &
                            start=start_dim, count=count_dim) )

            if ( n_aer_category > 1 ) then
                do i = 1, n_aer_category
                    call nc_verify( nf90_put_var(lncid, aerDSD_varids(i), &
                                    buffer_DSD(i+1, :, 1:buffer_count), &
                                    start=start_dim, count=count_dim) )
                end do
            end if
        end if

        do i = 1, n_field_budgets
            call nc_verify( nf90_put_var(lncid, varid_field_budgets(i), &
                            buffer_field_budgets(i, 1:buffer_count), &
                            start=(/nc_write_iter/)) )
        end do

        if (do_microphysics) then
            do i = 1, n_micro_budgets
                call nc_verify( nf90_put_var(lncid, varid_micro_budgets(i), &
                                buffer_micro_budgets(i, 1:buffer_count), &
                                start=(/nc_write_iter/)) )
            end do
        end if

        if (do_entrainment) then
            do i = 1, n_entrain_budgets
                call nc_verify( nf90_put_var(lncid, varid_entrain_budgets(i), &
                                buffer_entrain_budgets(i, 1:buffer_count), &
                                start=(/nc_write_iter/)) )
            end do
        end if

        if (do_parcel_ascent) then
            call nc_verify( nf90_put_var(lncid, varid_parcel_height, &
                            buffer_parcel_height(1:buffer_count), start=(/nc_write_iter/)) )
            call nc_verify( nf90_put_var(lncid, varid_parcel_pressure, &
                            buffer_parcel_pressure(1:buffer_count), start=(/nc_write_iter/)) )
            call nc_verify( nf90_put_var(lncid, varid_parcel_velocity, &
                            buffer_parcel_velocity(1:buffer_count), start=(/nc_write_iter/)) )
            if (write_height_env) then
                call nc_verify( nf90_put_var(lncid, varid_parcel_height_env, &
                                buffer_parcel_height_env(1:buffer_count), start=(/nc_write_iter/)) )
            end if
        end if

        if (do_entrainment .and. time_varying_entrainment) then
            call nc_verify( nf90_put_var(lncid, varid_ent_rate, &
                            buffer_ent_rate(1:buffer_count), start=(/nc_write_iter/)) )
            call nc_verify( nf90_put_var(lncid, varid_n_blob, &
                            buffer_n_blob(1:buffer_count), start=(/nc_write_iter/)) )
            call nc_verify( nf90_put_var(lncid, varid_psigma, &
                            buffer_psigma(1:buffer_count), start=(/nc_write_iter/)) )
        end if

        if (do_radiation) then
            call nc_verify( nf90_put_var(lncid, varid_rad_F_net, &
                            buffer_rad_F_net(:, 1:buffer_count), &
                            start=start_dim, count=count_dim) )
            call nc_verify( nf90_put_var(lncid, varid_rad_heating_rate, &
                            buffer_rad_heating_rate(:, 1:buffer_count), &
                            start=start_dim, count=count_dim) )
            call nc_verify( nf90_put_var(lncid, varid_rad_budget, &
                            buffer_rad_budget(1:buffer_count), start=(/nc_write_iter/)) )
        end if

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
        integer, intent(in) :: lncid
        character(*), intent(in) :: sim_name
        integer(i4), intent(in) :: lwrite_buffer

        call nc_verify( nf90_redef(lncid), "nf90_redef: namelist attributes" )

        ! PARAMETERS — shared
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.simulation_name", trim(sim_name)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.simulation_mode", trim(simulation_mode)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.N", N) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.tmax", tmax) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Tref", Tref) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.pres", pres) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.H", H) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.volume_scaling", volume_scaling) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_timer", write_timer) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_buffer", lwrite_buffer) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.same_random", merge(1, 0, same_random)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_turbulence", merge(1, 0, do_turbulence)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_microphysics", merge(1, 0, do_microphysics)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.write_eddies", merge(1, 0, write_eddies)) )
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.overwrite", merge(1, 0, overwrite)) )

        ! PARAMETERS — chamber only
        if (simulation_mode == 'chamber') then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.Tdiff", Tdiff) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_special_effects", &
                            merge(1, 0, do_special_effects)) )
        end if

        ! PARCEL namelist — parcel only
        if (simulation_mode == 'parcel') then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.parcel_file", trim(parcel_file)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.initial_RH", initial_RH) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.pressure_mode", trim(pressure_mode)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.vertical_axis", trim(vertical_axis)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.initial_height", initial_height) )
            if (pressure_limit > 0.0) &
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.pressure_limit", pressure_limit) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_entrainment", &
                            merge(1, 0, do_entrainment)) )
            if (do_entrainment) then
                ! ent_rate attribute in 1/km, matching the namelist/file interface
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.ent_rate", ent_rate * m_per_km) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.n_blob", n_blob) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.psigma", psigma) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARCEL.random_entrainment", &
                                merge(1, 0, random_entrainment)) )
            end if
        end if

        ! TURBULENCE_ODT — chamber only
        if (simulation_mode == 'chamber') then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_ODT.Lmin", Lmin) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_ODT.Lprob", Lprob) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_ODT.max_accept_prob", max_accept_prob) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_ODT.C2", C2) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_ODT.ZC2", ZC2) )
        end if

        ! TURBULENCE_LEM — parcel only
        if (simulation_mode == 'parcel') then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_LEM.integral_length_scale", &
                            integral_length_scale) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "TURBULENCE_LEM.dissipation_rate", &
                            dissipation_rate) )
            ! Derived turbulence scales (not namelist inputs). smallest_eddy_scale
            ! governs; actual_kolmogorov_scale is reported for reference only.
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.actual_kolmogorov_scale", &
                            actual_kolmogorov_scale) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.grid_eddy_scale", &
                            grid_eddy_scale) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.diffusivity_length_scale", &
                            diffusivity_length_scale) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.smallest_eddy_scale", &
                            smallest_eddy_scale) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.smallest_eddy_gridpoints", &
                            smallest_eddy_gridpoints) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "LEM.diffusivity_enhancement", &
                            diffusivity_enhancement) )
        end if

        ! MICROPHYSICS — shared
        if (do_microphysics) then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.aerosol_file", trim(aerosol_file)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.initial_wet_radius", initial_wet_radius) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.write_trajectories", &
                            merge(1, 0, write_trajectories)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_start", trajectory_start) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_end", trajectory_end) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.trajectory_timer", trajectory_timer) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.do_collisions", &
                            merge(1, 0, do_collisions)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.do_coalescence", &
                            merge(1, 0, do_coalescence)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.wmax_collision", wmax_collision) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.write_collisions", &
                            merge(1, 0, write_collisions)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.coalescence_kernel", &
                            trim(coalescence_kernel)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.do_seeding", &
                            merge(1, 0, do_seeding)) )
            if (do_seeding) then
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.seed_hydration", &
                                trim(seed_hydration)) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.seed_growth_time", &
                                seed_growth_time) )
            end if

            ! MICROPHYSICS — chamber only
            if (simulation_mode == 'chamber') then
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.init_drop_each_gridpoint", &
                                merge(1, 0, init_drop_each_gridpoint)) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.expected_Ndrops_per_gridpoint", &
                                expected_Ndrops_per_gridpoint) )
            end if

            ! MICROPHYSICS — parcel only
            if (simulation_mode == 'parcel') then
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "MICROPHYSICS.aerosol_concentration", &
                                aerosol_concentration) )
            end if
        end if

        ! RADIATION
        if (do_radiation) then
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_radiation", 1) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.radiation_method", &
                            trim(radiation_method)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.mie_data_file", &
                            trim(mie_data_file)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.eps_top", eps_top) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.eps_bot", eps_bot) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.sky_temp", sky_temp) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.sky_cooling_flag", &
                            merge(1, 0, sky_cooling_flag)) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.rad_call_interval", &
                            rad_call_interval) )
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.max_droplets_per_cell", &
                            max_droplets_per_cell) )
            if (radiation_method == '3d') then
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.nPhotons", nPhotons) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.nBins", nBins) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.Lx_rad", Lx_rad) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.Ly_rad", Ly_rad) )
                call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "RADIATION.T_side", T_side) )
            end if
        else
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "PARAMETERS.do_radiation", 0) )
        end if

        ! SPECIALEFFECTS — chamber only
        if (simulation_mode == 'chamber' .and. do_special_effects) then
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
            call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "SPECIALEFFECTS.random_fallout_rate", &
                            random_fallout_rate) )
        end if

        call nc_verify( nf90_enddef(lncid), "nf90_enddef: namelist attributes" )

    end subroutine write_namelist_attributes

end module writeout