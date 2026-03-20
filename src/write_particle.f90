module write_particle
    use netcdf
    use droplets, only: particle, particles, current_n_particles, write_trajectories, trajectory_start, trajectory_end, trajectory_timer
    use globals
    use writeout, only: nc_verify
    implicit none

    ! tracking state management
    logical :: tracking_activated = .false.
    real(dp) :: trajectory_time_iter = 0.

    ! NetCDF particle output state
    integer :: pnc_id
    integer :: record_count = 0
    integer :: time_step_count = 0

    ! Cached NetCDF variable IDs (set once in create_particle_netcdf)
    integer :: vid_particle_id, vid_aerosol_id, vid_gridcell, vid_position
    integer :: vid_temperature, vid_water_vapor, vid_supersaturation, vid_radius
    integer :: vid_solute_radius, vid_activated, vid_aerosol_category
    integer :: vid_time, vid_row_sizes

contains

    subroutine write_trajectory_controller(lparticles, ltime, ldt)

        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(in) :: ltime, ldt

        if ( .not. tracking_activated ) then
            if (ltime >= trajectory_start .and. ltime < trajectory_end) then
                tracking_activated = .true.
            end if
        end if

        trajectory_time_iter = trajectory_time_iter + ldt
        if ( trajectory_start <= ltime .and. trajectory_end > ltime ) then
            if ( trajectory_time_iter >= trajectory_timer ) then
                call write_particle_data_nc(lparticles, ltime)
                trajectory_time_iter = mod(trajectory_time_iter, trajectory_timer)
            end if
        end if

    end subroutine write_trajectory_controller

    subroutine initialize_write_particle(filename)

        character(*), intent(in) :: filename

        call create_particle_netcdf(filename)

    end subroutine initialize_write_particle

    subroutine create_particle_netcdf(filename)
        ! Creates the _particles.nc file using CF contiguous ragged array convention.
        ! Defines dimensions (record, time_step) and all per-record / per-time-step
        ! variables with compression and CF attributes.
        character(*), intent(in) :: filename

        integer :: rec_dimid, ts_dimid

        ! Reset counters
        record_count = 0
        time_step_count = 0

        ! Create file
        call nc_verify( nf90_create(trim(filename)//'_particles.nc', NF90_NETCDF4, pnc_id), &
                         "create_particle_netcdf: nf90_create" )

        ! Dimensions (both unlimited)
        call nc_verify( nf90_def_dim(pnc_id, "record", NF90_UNLIMITED, rec_dimid) )
        call nc_verify( nf90_def_dim(pnc_id, "time_step", NF90_UNLIMITED, ts_dimid) )

        ! --- Per-record variables (on record dimension) ---

        call nc_verify( nf90_def_var(pnc_id, "particle_id", NF90_INT, rec_dimid, vid_particle_id, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_particle_id, "long_name", "Particle ID") )

        call nc_verify( nf90_def_var(pnc_id, "aerosol_id", NF90_INT, rec_dimid, vid_aerosol_id, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_aerosol_id, "long_name", "Aerosol Type ID") )

        call nc_verify( nf90_def_var(pnc_id, "gridcell", NF90_INT, rec_dimid, vid_gridcell, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_gridcell, "long_name", "Grid Cell Index") )

        call nc_verify( nf90_def_var(pnc_id, "position", NF90_FLOAT, rec_dimid, vid_position, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_position, "long_name", "Vertical Position") )
        call nc_verify( nf90_put_att(pnc_id, vid_position, "units", "m") )

        call nc_verify( nf90_def_var(pnc_id, "temperature", NF90_FLOAT, rec_dimid, vid_temperature, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_temperature, "long_name", "Temperature") )
        call nc_verify( nf90_put_att(pnc_id, vid_temperature, "units", "celsius") )

        call nc_verify( nf90_def_var(pnc_id, "water_vapor", NF90_FLOAT, rec_dimid, vid_water_vapor, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_water_vapor, "long_name", "Water Vapor Mixing Ratio") )
        call nc_verify( nf90_put_att(pnc_id, vid_water_vapor, "units", "g/kg") )

        call nc_verify( nf90_def_var(pnc_id, "supersaturation", NF90_FLOAT, rec_dimid, vid_supersaturation, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_supersaturation, "long_name", "Supersaturation") )
        call nc_verify( nf90_put_att(pnc_id, vid_supersaturation, "units", "%") )

        call nc_verify( nf90_def_var(pnc_id, "radius", NF90_FLOAT, rec_dimid, vid_radius, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_radius, "long_name", "Droplet Radius") )
        call nc_verify( nf90_put_att(pnc_id, vid_radius, "units", "um") )

        call nc_verify( nf90_def_var(pnc_id, "solute_radius", NF90_FLOAT, rec_dimid, vid_solute_radius, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_solute_radius, "long_name", "Dry Solute Radius") )
        call nc_verify( nf90_put_att(pnc_id, vid_solute_radius, "units", "m") )

        call nc_verify( nf90_def_var(pnc_id, "activated", NF90_INT, rec_dimid, vid_activated, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_activated, "long_name", "Activation Flag") )
        call nc_verify( nf90_put_att(pnc_id, vid_activated, "flag_values", (/0, 1/)) )
        call nc_verify( nf90_put_att(pnc_id, vid_activated, "flag_meanings", "unactivated activated") )

        call nc_verify( nf90_def_var(pnc_id, "aerosol_category", NF90_INT, rec_dimid, vid_aerosol_category, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_aerosol_category, "long_name", "Aerosol Category") )

        ! --- Per-time-step variables (on time_step dimension) ---

        call nc_verify( nf90_def_var(pnc_id, "time", NF90_DOUBLE, ts_dimid, vid_time, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_time, "long_name", "Time") )
        call nc_verify( nf90_put_att(pnc_id, vid_time, "units", "seconds") )

        call nc_verify( nf90_def_var(pnc_id, "row_sizes", NF90_INT, ts_dimid, vid_row_sizes, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, vid_row_sizes, "long_name", "Number of particles per time step") )
        call nc_verify( nf90_put_att(pnc_id, vid_row_sizes, "cf_role", "ragged_row_sizes") )
        call nc_verify( nf90_put_att(pnc_id, vid_row_sizes, "sample_dimension", "record") )

        ! Exit define mode
        call nc_verify( nf90_enddef(pnc_id), "create_particle_netcdf: nf90_enddef" )

    end subroutine create_particle_netcdf

    subroutine write_particle_data_nc(lparticles, ltime)
        ! Writes one time step of particle data to the _particles.nc file.
        ! Extracts fields into temporary 1D arrays and appends to the
        ! record and time_step dimensions.
        type(particle), intent(in) :: lparticles(:)
        real(dp), intent(in) :: ltime

        integer :: np, i
        integer(i4), allocatable :: int_buf(:)
        real, allocatable :: float_buf(:)

        np = current_n_particles
        if (np == 0) return

        ! --- Per-time-step variables ---
        time_step_count = time_step_count + 1

        call nc_verify( nf90_put_var(pnc_id, vid_time, ltime, start=(/time_step_count/)) )
        call nc_verify( nf90_put_var(pnc_id, vid_row_sizes, np, start=(/time_step_count/)) )

        ! --- Per-record variables ---
        allocate(int_buf(np), float_buf(np))

        ! particle_id
        do i = 1, np; int_buf(i) = lparticles(i)%particle_id; end do
        call nc_verify( nf90_put_var(pnc_id, vid_particle_id, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! aerosol_id
        do i = 1, np; int_buf(i) = lparticles(i)%solute_type%aerosol_id; end do
        call nc_verify( nf90_put_var(pnc_id, vid_aerosol_id, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! gridcell
        do i = 1, np; int_buf(i) = lparticles(i)%gridcell; end do
        call nc_verify( nf90_put_var(pnc_id, vid_gridcell, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! position
        do i = 1, np; float_buf(i) = real(lparticles(i)%position); end do
        call nc_verify( nf90_put_var(pnc_id, vid_position, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! temperature (K -> C)
        do i = 1, np; float_buf(i) = real(lparticles(i)%temperature - Tice); end do
        call nc_verify( nf90_put_var(pnc_id, vid_temperature, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! water_vapor (kg/kg -> g/kg)
        do i = 1, np; float_buf(i) = real(lparticles(i)%water_vapor * g_per_kg); end do
        call nc_verify( nf90_put_var(pnc_id, vid_water_vapor, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! supersaturation
        do i = 1, np; float_buf(i) = real(lparticles(i)%supersaturation); end do
        call nc_verify( nf90_put_var(pnc_id, vid_supersaturation, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! radius (m -> um)
        do i = 1, np; float_buf(i) = real(lparticles(i)%radius * um_per_m); end do
        call nc_verify( nf90_put_var(pnc_id, vid_radius, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! solute_radius
        do i = 1, np; float_buf(i) = real(lparticles(i)%solute_radius); end do
        call nc_verify( nf90_put_var(pnc_id, vid_solute_radius, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! activated (logical -> 0/1)
        do i = 1, np
            if (lparticles(i)%activated) then
                int_buf(i) = 1
            else
                int_buf(i) = 0
            end if
        end do
        call nc_verify( nf90_put_var(pnc_id, vid_activated, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! aerosol_category
        do i = 1, np; int_buf(i) = lparticles(i)%aerosol_category; end do
        call nc_verify( nf90_put_var(pnc_id, vid_aerosol_category, int_buf, start=(/record_count+1/), count=(/np/)) )

        deallocate(int_buf, float_buf)

        ! Update record counter
        record_count = record_count + np

    end subroutine write_particle_data_nc

    subroutine close_particle_netcdf()

        call nc_verify( nf90_close(pnc_id), "close_particle_netcdf: nf90_close" )

    end subroutine close_particle_netcdf


end module write_particle
