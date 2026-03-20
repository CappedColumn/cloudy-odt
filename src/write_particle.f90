module write_particle
    use, intrinsic :: iso_c_binding
    use netcdf
    use droplets, only: particle, particles, current_n_particles, write_trajectories, trajectory_start, trajectory_end, write_particle_id, trajectory_timer
    use globals
    use writeout, only: nc_verify
    implicit none

    type :: w_particle
        sequence ! Parsing scripts expect values written to be in this order
        ! "(W)rite-particle" type using only fields which will be written to netCDF
        ! Combines particle and aerosol types for simplicity
        integer(i4) :: id
        integer(i4) :: gridcell
        real(dp) :: position
        real(dp) :: temperature
        real(dp) :: water_vapor
        real(dp) :: supersaturation
        real(dp) :: radius
        logical(1) :: activated ! Aligning with 1-bit value from parser

    end type

    !type(w_particle), allocatable :: wparticle_array(:)
    type(w_particle) :: mydrop

    ! file unit  handling
    integer :: data_unit, meta_unit

    ! tracking state management
    logical :: tracking_activated = .false.
    real(dp) :: trajectory_time_iter = 0.

    ! NetCDF particle output state
    integer :: pnc_id
    integer :: record_count = 0
    integer :: time_step_count = 0


    ! Main - write_trajectory_controller

    contains

    subroutine write_trajectory_controller(lparticles, ltime, ldt)

        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(in) :: ltime, ldt
        integer :: i

        if ( .not. tracking_activated ) then
            if (ltime >= trajectory_start .and. ltime < trajectory_end) then
            ! We should be tracking - turn on if not already active
                tracking_activated = .true.
                do i = 1, current_n_particles
                    call write_particle_id(lparticles(i))
                end do
            end if
        end if

        trajectory_time_iter = trajectory_time_iter + ldt
        if ( trajectory_start <= time .and. trajectory_end > ltime ) then
            if ( trajectory_time_iter >= trajectory_timer ) then
                call write_particle_data(lparticles, ltime)
                trajectory_time_iter = mod(trajectory_time_iter, trajectory_timer)
            end if
        end if

    end subroutine write_trajectory_controller

    subroutine initialize_write_particle(filename)

        character(*), intent(in) :: filename

        call wparticle_initialize(mydrop, particles(1))
        call open_particle_files(filename)
        call write_particle_meta_data()

    end subroutine initialize_write_particle

    subroutine open_particle_files(filename)

        character(*), intent(in) :: filename
        integer :: ios

        ! Open the file for binary writing
        open(newunit=data_unit, file=trim(filename)//'_particle.bin', form='unformatted', &
        access='stream', status='replace', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: "
            stop 1
        end if

        open(newunit=meta_unit, file=trim(filename)//'_particle_meta.txt', status='replace', &
        iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: "
            stop 1
        end if

    end subroutine open_particle_files

    subroutine close_particle_files()

        close(data_unit)
        close(meta_unit)

    end subroutine

    subroutine write_particle_data(lparticles, ltime)

        type(particle), intent(in) :: lparticles(:)
        real(dp), intent(in) :: ltime
        type(w_particle), allocatable :: w_particles(:)
        integer :: i

        allocate(w_particles(current_n_particles))

        do i = 1, current_n_particles
            call wparticle_initialize(w_particles(i), lparticles(i))
        end do

        ! Write the time and number of particles
        write(meta_unit, *) time, current_n_particles

        ! Write each particle to the file
        do i = 1, current_n_particles
            write(data_unit) w_particles(i), time
        end do

        deallocate(w_particles)
        
    end subroutine write_particle_data

    subroutine write_particle_meta_data()

        character(11) :: format = '(a4,a8,i10)'

        ! Write properties of aerosol
        write(meta_unit, *)

        ! Write format strucutre of particle data
        write(meta_unit, *) 'Data format for Particle Stream File (sequential order)'
        write(meta_unit, '(a4,a8,a14)') 'Name', 'Type', 'No. of Bits'
        write(meta_unit, format) 'PID', 'int', sizeof(mydrop%id)*8
        write(meta_unit, format) 'GRD', 'int', sizeof(mydrop%gridcell)*8
        write(meta_unit, format) 'POS', 'float', sizeof(mydrop%position)*8
        write(meta_unit, format) 'TPC', 'float', sizeof(mydrop%temperature)*8
        write(meta_unit, format) 'WVM', 'float', sizeof(mydrop%water_vapor)*8
        write(meta_unit, format) 'SSN', 'float', sizeof(mydrop%supersaturation)*8
        write(meta_unit, format) 'RAD', 'float', sizeof(mydrop%radius)*8
        write(meta_unit, format) 'ACT', 'bool', sizeof(mydrop%activated)
        write(meta_unit, *)

        write(meta_unit, *) 'Time - N-Particles for each particle write command'

    end subroutine write_particle_meta_data



    subroutine wparticle_initialize(this, droplet)
        ! Copies over properties of particle type into the
        ! write-out version
        type(w_particle), intent(out) :: this
        type(particle), intent(in) :: droplet

        this%id = droplet%particle_id
        this%gridcell = droplet%gridcell
        this%position = droplet%position
        this%temperature = droplet%temperature
        this%water_vapor = droplet%water_vapor
        this%supersaturation = droplet%supersaturation
        this%radius = droplet%radius * um_per_m
        this%activated = droplet%activated

    end subroutine wparticle_initialize

    subroutine create_particle_netcdf(filename)
        ! Creates the _particles.nc file using CF contiguous ragged array convention.
        ! Defines dimensions (record, time_step) and all per-record / per-time-step
        ! variables with compression and CF attributes.
        character(*), intent(in) :: filename

        integer :: rec_dimid, ts_dimid
        integer :: varid

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

        call nc_verify( nf90_def_var(pnc_id, "particle_id", NF90_INT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Particle ID") )

        call nc_verify( nf90_def_var(pnc_id, "aerosol_id", NF90_INT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Aerosol Type ID") )

        call nc_verify( nf90_def_var(pnc_id, "gridcell", NF90_INT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Grid Cell Index") )

        call nc_verify( nf90_def_var(pnc_id, "position", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Vertical Position") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "m") )

        call nc_verify( nf90_def_var(pnc_id, "temperature", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Temperature") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "K") )

        call nc_verify( nf90_def_var(pnc_id, "water_vapor", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Water Vapor Mixing Ratio") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "kg/kg") )

        call nc_verify( nf90_def_var(pnc_id, "supersaturation", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Supersaturation") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "%") )

        call nc_verify( nf90_def_var(pnc_id, "radius", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Droplet Radius") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "um") )

        call nc_verify( nf90_def_var(pnc_id, "solute_radius", NF90_FLOAT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Dry Solute Radius") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "m") )

        call nc_verify( nf90_def_var(pnc_id, "activated", NF90_INT, rec_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Activation Flag") )
        call nc_verify( nf90_put_att(pnc_id, varid, "flag_values", (/0, 1/)) )
        call nc_verify( nf90_put_att(pnc_id, varid, "flag_meanings", "unactivated activated") )

        ! --- Per-time-step variables (on time_step dimension) ---

        call nc_verify( nf90_def_var(pnc_id, "time", NF90_DOUBLE, ts_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Time") )
        call nc_verify( nf90_put_att(pnc_id, varid, "units", "seconds") )

        call nc_verify( nf90_def_var(pnc_id, "row_sizes", NF90_INT, ts_dimid, varid, &
                         deflate_level=1, shuffle=.true.) )
        call nc_verify( nf90_put_att(pnc_id, varid, "long_name", "Number of particles per time step") )
        call nc_verify( nf90_put_att(pnc_id, varid, "cf_role", "ragged_row_sizes") )
        call nc_verify( nf90_put_att(pnc_id, varid, "sample_dimension", "record") )

        ! Exit define mode
        call nc_verify( nf90_enddef(pnc_id), "create_particle_netcdf: nf90_enddef" )

    end subroutine create_particle_netcdf

    subroutine write_particle_data_nc(lparticles, ltime)
        ! Writes one time step of particle data to the _particles.nc file.
        ! Extracts fields into temporary 1D arrays and appends to the
        ! record and time_step dimensions.
        type(particle), intent(in) :: lparticles(:)
        real(dp), intent(in) :: ltime

        integer :: np, i, varid
        integer(i4), allocatable :: int_buf(:)
        real, allocatable :: float_buf(:)

        np = current_n_particles
        if (np == 0) return

        ! --- Per-time-step variables ---
        time_step_count = time_step_count + 1

        call nc_verify( nf90_inq_varid(pnc_id, "time", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, ltime, start=(/time_step_count/)) )

        call nc_verify( nf90_inq_varid(pnc_id, "row_sizes", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, np, start=(/time_step_count/)) )

        ! --- Per-record variables ---
        allocate(int_buf(np), float_buf(np))

        ! particle_id
        do i = 1, np; int_buf(i) = lparticles(i)%particle_id; end do
        call nc_verify( nf90_inq_varid(pnc_id, "particle_id", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! aerosol_id
        do i = 1, np; int_buf(i) = lparticles(i)%solute_type%aerosol_id; end do
        call nc_verify( nf90_inq_varid(pnc_id, "aerosol_id", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! gridcell
        do i = 1, np; int_buf(i) = lparticles(i)%gridcell; end do
        call nc_verify( nf90_inq_varid(pnc_id, "gridcell", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, int_buf, start=(/record_count+1/), count=(/np/)) )

        ! position
        do i = 1, np; float_buf(i) = real(lparticles(i)%position); end do
        call nc_verify( nf90_inq_varid(pnc_id, "position", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! temperature
        do i = 1, np; float_buf(i) = real(lparticles(i)%temperature); end do
        call nc_verify( nf90_inq_varid(pnc_id, "temperature", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! water_vapor
        do i = 1, np; float_buf(i) = real(lparticles(i)%water_vapor); end do
        call nc_verify( nf90_inq_varid(pnc_id, "water_vapor", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! supersaturation
        do i = 1, np; float_buf(i) = real(lparticles(i)%supersaturation); end do
        call nc_verify( nf90_inq_varid(pnc_id, "supersaturation", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! radius (m -> um)
        do i = 1, np; float_buf(i) = real(lparticles(i)%radius * um_per_m); end do
        call nc_verify( nf90_inq_varid(pnc_id, "radius", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! solute_radius
        do i = 1, np; float_buf(i) = real(lparticles(i)%solute_radius); end do
        call nc_verify( nf90_inq_varid(pnc_id, "solute_radius", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, float_buf, start=(/record_count+1/), count=(/np/)) )

        ! activated (logical -> 0/1)
        do i = 1, np
            if (lparticles(i)%activated) then
                int_buf(i) = 1
            else
                int_buf(i) = 0
            end if
        end do
        call nc_verify( nf90_inq_varid(pnc_id, "activated", varid) )
        call nc_verify( nf90_put_var(pnc_id, varid, int_buf, start=(/record_count+1/), count=(/np/)) )

        deallocate(int_buf, float_buf)

        ! Update record counter
        record_count = record_count + np

    end subroutine write_particle_data_nc

    subroutine close_particle_netcdf()

        call nc_verify( nf90_close(pnc_id), "close_particle_netcdf: nf90_close" )

    end subroutine close_particle_netcdf


end module write_particle