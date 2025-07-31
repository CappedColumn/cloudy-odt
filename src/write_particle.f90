module write_particle
    use, intrinsic :: iso_c_binding
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
            stop
        end if

        open(newunit=meta_unit, file=trim(filename)//'_particle_meta.txt', status='replace', &
        iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: "
            stop
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


end module write_particle