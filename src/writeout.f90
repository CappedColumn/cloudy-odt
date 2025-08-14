module writeout
    use globals
    use netcdf
    use droplets, only: particles, calculate_droplet_statistics, bin_droplet_radii, particle_bin_edges, &
                        size_distribution, n_aer_category
    implicit none

    ! Buffer variables and arrays for writing to netCDF
    integer(i4) :: buffer_size, buffer_count, nc_write_iter
    real(dp), allocatable :: buffer_T(:,:), buffer_WV(:,:), buffer_Tv(:,:) ! dims (buffer_size, N_grid)
    real(dp), allocatable :: buffer_SS(:,:), buffer_W(:,:), buffer_time(:), buffer_stats(:,:)
    integer(i4), allocatable :: buffer_DSD(:,:,:) ! n_DSDs, rbins, buffer_size

    ! Buffer arrays for writting eddy data
    integer(i4) :: buffer_eddy_size = 50000 ! Will write to disk every 50000 eddies
    integer(i4) :: eddy_counter = 0
    integer(i4) :: eddy_unit
    real(dp), allocatable :: buffer_eddy(:,:)

    !private :: nc_verify
    public  :: create_netcdf, initialize_buffers, add_to_profile_buffer, write_netcdf_profiles, flush_buffer, close_netcdf
    
contains

    subroutine write_data()

        write(*,*) 'Writing time: ', time
        call calculate_droplet_statistics(particles, statistics)
        call bin_droplet_radii(particles, particle_bin_edges, size_distribution)
        call add_to_profile_buffer(time, Tdim, WVdim, Tvdim, SS, Wdim, size_distribution, statistics)
        ! Reset write timer, note "extra" time from eddy method is
        ! considered to prevent time compounding
        write_time_iter = mod(write_time_iter, write_timer)

    end subroutine write_data

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
        allocate(buffer_W(N_grid, buff_len))
        allocate(buffer_time(buff_len))
        allocate(buffer_stats(5, buff_len))

        buffer_T = 0.
        buffer_WV = 0.
        buffer_Tv = 0.
        buffer_SS = 0.
        buffer_W = 0.
        buffer_time = 0.
        buffer_stats = 0.

    end subroutine initialize_buffers
    
    subroutine initialize_particle_buffers(n_cat, n_bins)

        integer(i4), intent(in) :: n_cat, n_bins

        if ( n_cat > 1 ) then
            ! Account for subcategorical DSDs
            allocate(buffer_DSD(n_cat+1, n_bins, buffer_size))
        else
            allocate(buffer_DSD(1, n_bins, buffer_size))
        ! Calculate number of DSDs (n_cat + 1 for multiple categories, 1 otherwise)
        integer(i4) :: n_DSDs
        if ( n_cat > 1 ) then
            n_DSDs = n_cat + 1
        else
            n_DSDs = 1
        end if
        allocate(buffer_DSD(n_DSDs, n_bins, buffer_size))
        buffer_DSD = 0

    end subroutine initialize_particle_buffers

    recursive subroutine add_to_profile_buffer(ltime, lT, lWV, lTv, lSS, lW, lDSD, lstats)
        ! Writes the profile arrays into the profile buffers, will flush
        ! the buffer to write netCDF if buffer is full
        real(dp), intent(in) :: ltime, lT(:), lWV(:), lTv(:), lSS(:), lW(:), lstats(:)
        integer(i4), intent(in) :: lDSD(:,:)

        ! Write profile data into buffers
        if (buffer_count < buffer_size) then
            buffer_count = buffer_count + 1 ! recall buffer_count resets to 0
            buffer_time(buffer_count) = ltime
            buffer_T(:, buffer_count) = lT - Tice ! (C)
            buffer_WV(:, buffer_count) = lWV
            buffer_Tv(:, buffer_count) = lTv - Tice ! (C)
            buffer_SS(:, buffer_count) = lSS
            buffer_W(:, buffer_count) = lW
            buffer_stats(:, buffer_count) = lstats
            buffer_DSD(:, :, buffer_count) = lDSD
        else
            ! Flush buffer and start new buffer
            call flush_buffer()
            call add_to_profile_buffer(ltime, lT, lWV, lTv, lSS, lW, lDSD, lstats)
        end if

    end subroutine add_to_profile_buffer

    subroutine initialize_eddy_buffer(filename)

        character(*), intent(in) :: filename
        integer(i4) :: ierr

        allocate(buffer_eddy(buffer_eddy_size,3))

        open(newunit=eddy_unit, file=trim(filename)//'_eddies.txt', form='formatted', status='replace', iostat=ierr)
        if (ierr /= 0) then
            print *, "Error opening eddy data file. "
            stop
        end if

    end subroutine initialize_eddy_buffer

    recursive subroutine add_to_eddy_buffer(loc, len, ltime)

        integer(i4), intent(in) :: loc, len
        real(dp), intent(in) :: ltime

        if ( eddy_counter < buffer_eddy_size ) then
            eddy_counter = eddy_counter + 1
            buffer_eddy(eddy_counter, 1) = (1.*loc+(len/2.))/N ! Midpoint (z)
            buffer_eddy(eddy_counter, 2) = 1.*len/(2.*N)       ! 1/2 length (real)
            buffer_eddy(eddy_counter, 3) = ltime               ! time (s)
        else
            call write_eddy_buffer() ! Will reset eddy_counter to 0 for recursive call
            call add_to_eddy_buffer(loc, len, ltime)
        end if

    end subroutine

    subroutine flush_buffer()
        ! writes the buffer to the netCDF file and resets buffer position
        call write_netcdf_profiles(ncid, buffer_time(1:buffer_count), &
                                        buffer_T(:, 1:buffer_count), &
                                        buffer_WV(:, 1:buffer_count), &
                                        buffer_Tv(:, 1:buffer_count), &
                                        buffer_SS(:, 1:buffer_count), &
                                        buffer_W(:, 1:buffer_count), &
                                        buffer_DSD(:, :, 1:buffer_count), &
                                        buffer_stats(:, 1:buffer_count))
        buffer_count = 0

    end subroutine flush_buffer

    subroutine write_eddy_buffer()

        integer(i4) :: i

        do i = 1, eddy_counter
            write(eddy_unit, '(f10.3, 1X, f10.3, 1X, f10.6)') buffer_eddy(i, 1), buffer_eddy(i, 2), buffer_eddy(i, 3)
        end do
        eddy_counter = 0

    end subroutine write_eddy_buffer

    subroutine create_netcdf(file_name, z_m, lncid)
        character(*), intent(in):: file_name
        real(dp), intent(in) :: z_m(:) !, scalar_vars(:) ! Establish z-dimension
        integer, intent(out) :: lncid

        ! time/height dimensions
        logical :: file_exists
        integer :: old_nc, stat
        integer :: t_dimid, z_dimid, dimids(2)
        ! Variable initialization
        integer :: t_varid, z_varid
        integer :: tc_varid, qv_varid, s_varid, tv_varid, w_varid
        integer :: statids(5)
        integer :: r_prtcl_varid
        integer :: j, k, nz, nbins

        nz = size(z_m)
        nc_write_iter = 1
        
        ! Static Variables
        ! integer :: i, nvars
        ! nvars = size(scalar_vars)
        ! character(*) :: name(nvars), long_name(nvars), units(nvars)
        ! integer :: varids(nvars)

        
        ! allocate(name(nvars))
        ! allocate(long_name(nvars))
        ! allocate(units(nvars))
        ! allocate(varids(nvars))

        ! name = ['Pr', 'Sc', 'Nu', 'kT', 'Dv', 'g', 'C2', 'ZC2']
        ! long_name = ['Prandtl Number', 'Schmidt Number', 'Kinematic Viscosity', &
        !             'Thermal Diffusivity', 'Water Vapor Diffusivity', 'Gravity', &
        !             'Turbulent Strength', 'Viscous Cut-Off Parameter']
        ! units = ['non-dim', 'non-dim', 'kg/m/s', 'm2/s', 'm2/s', 'm/s2', 'non-dim', &
        !         'non-dim']
    
    
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
    
        ! Simulation global attributes NF90_GLOBAL
        call nc_verify( nf90_put_att(lncid, NF90_GLOBAL, "Attribute", "Global"), "nf90_put_att: global")

        ! Create variables for scalar properties of simulation
        ! do i = 1,nvars
        !     call nc_verify( nf90_def_var(lncid, name(i), NF90_FLOAT, varids(i)) )
        ! end do
    
        ! Establish variables and attributes for coordinate(dimension) variables
        call nc_verify( nf90_def_var(lncid, "z", NF90_FLOAT, z_dimid, z_varid), "nf90_def_var: z" )
        call nc_verify( nf90_put_att(lncid, z_varid, "units", "meters"), "nf90_put_att: z, units" )
        call nc_verify( nf90_def_var(lncid, "time", NF90_FLOAT, t_dimid, t_varid), "nf90_def_var: time" )
        call nc_verify( nf90_put_att(lncid, t_varid, "units", "seconds"), "nf90_put_att: time, units" )

    
        ! Note, netCDF will write out variables in (time, height) order
        dimids = (/ z_dimid, t_dimid /) ! Standard time/height dimensions
        call nc_verify( nf90_def_var(lncid, "T", NF90_FLOAT, dimids, tc_varid), "nf90_def_var: T" )
        call nc_verify( nf90_put_att(lncid, tc_varid, "long name", "Temperature"), "nf90_put_att: T, name" )
        call nc_verify( nf90_put_att(lncid, tc_varid, "units", "celsius"), "nf90_put_att: T, units" )
        
        call nc_verify( nf90_def_var(lncid, "QV", NF90_FLOAT, dimids, qv_varid), "nf90_def_var: QV" )
        call nc_verify( nf90_put_att(lncid, qv_varid, "long name", "Water Vapor Mixing Ratio"), "nf90_put_att: QV, name" )
        call nc_verify( nf90_put_att(lncid, qv_varid, "units", "g/kg"), "nf90_put_att: QV, units")

        call nc_verify( nf90_def_var(lncid, "Tv", NF90_FLOAT, dimids, tv_varid), "nf90_def_var: Tv" )
        call nc_verify( nf90_put_att(lncid, tv_varid, "long name", "Virtual Temperature"), "nf90_put_att: Tv, name" )
        call nc_verify( nf90_put_att(lncid, tv_varid, "units", "celcius"), "nf90_put_att: Tv, units")
        
        call nc_verify( nf90_def_var(lncid, "S", NF90_FLOAT, dimids, s_varid), "nf90_def_var: S" )
        call nc_verify( nf90_put_att(lncid, s_varid, "long name", "Supersaturation"), "nf90_put_att: S, name" )
        call nc_verify( nf90_put_att(lncid, s_varid, "units", "%"), "nf90_put_att: S, units")

        ! Vertical Velocity
        call nc_verify( nf90_def_var(lncid, "W", NF90_FLOAT, dimids, w_varid), "nf90_def_var: W" )
        call nc_verify( nf90_put_att(lncid, w_varid, "long name", "W-Velocity"), "nf90_put_att: W, name" )
        call nc_verify( nf90_put_att(lncid, w_varid, "units", "m/s"), "nf90_put_att: W, units")



        ! Statistic Variable creation
        call nc_verify( nf90_def_var(lncid, "Np", NF90_INT, t_dimid, statids(1)), "nf90_def_var: Np" )
        call nc_verify( nf90_put_att(lncid, statids(1), "long name", "Number of Particles"), "nf90_put_att: Np, name" )
        call nc_verify( nf90_put_att(lncid, statids(1), "units", "#"), "nf90_put_att: Np units")

        call nc_verify( nf90_def_var(lncid, "Nact", NF90_INT, t_dimid, statids(2)), "nf90_def_var: Nact" )
        call nc_verify( nf90_put_att(lncid, statids(2), "long name", "Number of Activated Particles"), "nf90_put_att: Nact, name" )
        call nc_verify( nf90_put_att(lncid, statids(2), "units", "#"), "nf90_put_att: Nact units")

        call nc_verify( nf90_def_var(lncid, "Nun", NF90_INT, t_dimid, statids(3)), "nf90_def_var: Nun" )
        call nc_verify( nf90_put_att(lncid, statids(3), "long name", "Number of Unactivated Particles"), "nf90_put_att: Nun, name" )
        call nc_verify( nf90_put_att(lncid, statids(3), "units", "#"), "nf90_put_att: Nun units")

        call nc_verify( nf90_def_var(lncid, "Ravg", NF90_FLOAT, t_dimid, statids(4)), "nf90_def_var: Ravg" )
        call nc_verify( nf90_put_att(lncid, statids(4), "long name", "Average Particle Radius (wet)"), "nf90_put_att: Ravg, name" )
        call nc_verify( nf90_put_att(lncid, statids(4), "units", "um"), "nf90_put_att: Ravg units")

        call nc_verify( nf90_def_var(lncid, "LWC", NF90_FLOAT, t_dimid, statids(5)), "nf90_def_var: LWC" )
        call nc_verify( nf90_put_att(lncid, statids(5), "long name", "Liquid Water Content"), "nf90_put_att: LWC, name" )
        call nc_verify( nf90_put_att(lncid, statids(5), "units", "g/m3"), "nf90_put_att: LWC units")

    
        ! Variable length particle variables
        ! call nc_verify( nf90_def_vlen(lncid, "r", NF90_FLOAT, r_prtcl_varid) )
        ! call nc_verify( nf90_put_att(lncid, r_prtcl_varid, "long name", "Particle Radius"), "nf90_put_att: r, name" )
        ! call nc_verify( nf90_put_att(lncid, r_prtcl_varid, "units", "microns"), "nf90_put_att: r, units")
    
        ! Exit define mode, however netCDF is still open
        call nc_verify( nf90_enddef(lncid), "nf90_enddef" )
    
        ! Fill variables for known dimensions
        call nc_verify( nf90_put_var(lncid, z_varid, z_m), "nf90_put_var: z" )

        ! Fill variables for already determined scalar variables
        ! do i = 1,nvars
        !     call nc_verify( nf90_put_var(lncid, varids(i), scalar_vars(i)), "nf90_put_var: scalar" )
        !     call nc_verify( nf90_put_att(lncid, varids(i), "long name", long_name(i)), "nf90_put_att: scalar" )
        !     call nc_verify( nf90_put_att(lncid, varids(i), "units", units(i)), "nf90_put_att: scalar" )
        ! end do
    
        write(*,'(a)') "NetCDF file fully created."
    
    end subroutine create_netcdf


    subroutine write_netcdf_profiles(lncid, ltime, lT, lWV, lTv, lSS, lw, lDSD, lstats)
        integer(i4), intent(in) :: lncid
        real(dp), intent(in) :: ltime(:), lT(:,:), lWV(:,:), lTv(:,:), lSS(:,:), lw(:,:), lstats(:,:)
        integer(i4), intent(in) :: lDSD(:,:,:)

        integer :: varids(6), DSDid, i ! time, T, WV, Tv, SS, w, DSD are coordinates
        integer :: statids(5), time_len, z_len, bin_len
        integer :: count_dim(2), start_dim(2)
        character(100) :: varname, strint
    
        ! netCDF profile variables are (time, height) dimensions
        time_len = size(ltime,1)
        z_len = size(lT,1)
        bin_len = size(lDSD,2)
        if (size(lT,2) /= time_len) then
            write(*,'(a)') "Error: Time and profile arrays are not the same length."
            stop
        end if
        ! set the starting/ending positions for netCDF file
        count_dim = (/z_len, time_len/)
        start_dim = (/1, nc_write_iter/)
        
        ! Get variable IDs
        call nc_verify( nf90_inq_varid(lncid, "time", varids(1)) )   
        call nc_verify( nf90_inq_varid(lncid, "T", varids(2)) )
        call nc_verify( nf90_inq_varid(lncid, "QV", varids(3)) )
        call nc_verify( nf90_inq_varid(lncid, "Tv", varids(4)) )
        call nc_verify( nf90_inq_varid(lncid, "S", varids(5)) )
        call nc_verify( nf90_inq_varid(lncid, "W", varids(6)) )
        
        call nc_verify( nf90_inq_varid(ncid, "Np", statids(1)) )
        call nc_verify( nf90_inq_varid(ncid, "Nact", statids(2)) )
        call nc_verify( nf90_inq_varid(ncid, "Nun", statids(3)) )
        call nc_verify( nf90_inq_varid(ncid, "Ravg", statids(4)) )
        call nc_verify( nf90_inq_varid(ncid, "LWC", statids(5)) )

        ! Write to variable slots
        call nc_verify( nf90_put_var(ncid, varids(1), ltime, start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(ncid, varids(2), lT, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(ncid, varids(3), lWV, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(ncid, varids(4), lTv, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(ncid, varids(5), lSS, start=start_dim, count=count_dim) )
        call nc_verify( nf90_put_var(ncid, varids(6), lw, start=start_dim, count=count_dim) )
        
        if ( do_microphysics ) then

            ! Write total DSD
            call nc_verify( nf90_inq_varid(lncid, "DSD", DSDid ))
            count_dim = (/ bin_len, time_len /)
            call nc_verify( nf90_put_var(ncid, DSDid, lDSD(1,:,:), start=start_dim, count=count_dim) )

            ! Write the subcategory DSDs
            if ( n_aer_category > 1 ) then
                do i = 1, n_aer_category
                    write(strint,*) i
                    varname = "DSD_" // adjustl(strint)
                    call nc_verify( nf90_inq_varid(lncid, trim(varname), DSDid ))
                    call nc_verify( nf90_put_var(ncid, DSDid, lDSD(i+1,:,:), start=start_dim, count=count_dim) )
                end do
            end if

        end if

        call nc_verify( nf90_put_var(ncid, statids(1), lstats(1,:), start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(ncid, statids(2), lstats(2,:), start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(ncid, statids(3), lstats(3,:), start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(ncid, statids(4), lstats(4,:), start=(/nc_write_iter/)) )
        call nc_verify( nf90_put_var(ncid, statids(5), lstats(5,:), start=(/nc_write_iter/)) )

        ! Move 'start' time location to end of buffer for next write
        nc_write_iter = nc_write_iter + buffer_count

    end subroutine write_netcdf_profiles

    subroutine close_netcdf(lncid)
        integer, intent(in) :: lncid

        call nc_verify( nf90_close(lncid), 'nf90_close')
        if ( write_eddies ) then
            if ( eddy_counter > 0 ) call write_eddy_buffer()
            close(eddy_unit)
        end if

    end subroutine




    

end module writeout