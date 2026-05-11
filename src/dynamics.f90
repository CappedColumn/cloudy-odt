module dynamics
    use globals
    use netcdf
    use microphysics, only: virtual_temp, update_supersat
    implicit none

    private
    public :: initialize_dynamics, apply_adiabatic_forcing, &
              do_parcel_ascent, parcel_height, parcel_velocity

    integer(i4) :: n_segments
    real(dp), allocatable :: segment_times(:)
    real(dp), allocatable :: segment_velocity(:)

    real(dp) :: parcel_height   = 0.0
    real(dp) :: parcel_velocity = 0.0
    logical  :: do_parcel_ascent = .false.

contains

    subroutine initialize_dynamics(filepath)
        character(*), intent(in) :: filepath
        integer :: dyn_ncid, varid, dimid, i
        character(64) :: conventions

        write(*,*) 'Reading dynamics data from: ', trim(filepath)

        call nc_verify(nf90_open(trim(filepath), NF90_NOWRITE, dyn_ncid), &
                       'opening dynamics file')
        call nc_verify(nf90_get_att(dyn_ncid, NF90_GLOBAL, 'conventions', conventions), &
                       'reading conventions attribute')
        if (trim(conventions) /= 'CODT_dynamics_input_v1') then
            write(0,*) 'Error: expected CODT_dynamics_input_v1, got: ', trim(conventions)
            stop 1
        end if

        call nc_verify(nf90_inq_dimid(dyn_ncid, 'segment', dimid), 'finding segment dim')
        call nc_verify(nf90_inquire_dimension(dyn_ncid, dimid, len=n_segments), 'reading segment dim')

        allocate(segment_times(n_segments))
        call nc_verify(nf90_inq_varid(dyn_ncid, 'time', varid), 'finding time')
        call nc_verify(nf90_get_var(dyn_ncid, varid, segment_times), 'reading time')

        allocate(segment_velocity(n_segments))
        call nc_verify(nf90_inq_varid(dyn_ncid, 'velocity', varid), 'finding velocity')
        call nc_verify(nf90_get_var(dyn_ncid, varid, segment_velocity), 'reading velocity')

        call nc_verify(nf90_close(dyn_ncid), 'closing dynamics file')

        if (abs(segment_times(1)) > 1.0e-10) then
            write(0,*) 'Error: first segment time must be 0, got: ', segment_times(1)
            stop 1
        end if
        do i = 2, n_segments
            if (segment_times(i) <= segment_times(i-1)) then
                write(0,*) 'Error: segment times must be monotonically increasing'
                stop 1
            end if
        end do

        do_parcel_ascent = .true.
        parcel_velocity = segment_velocity(1)

        write(*,*) '--- Dynamics Configuration ---'
        write(*,*) 'n_segments:        ', n_segments
        write(*,*) 'initial velocity:  ', segment_velocity(1), ' m/s'
        write(*,*) 'initial pressure:  ', pres / Pa_per_mb, ' mb'
        write(*,*) '------------------------------'

    end subroutine initialize_dynamics


    pure function get_velocity(query_time) result(vel)
        real(dp), intent(in) :: query_time
        real(dp) :: vel
        integer :: i

        vel = segment_velocity(1)
        do i = 2, n_segments
            if (query_time < segment_times(i)) exit
            vel = segment_velocity(i)
        end do
    end function get_velocity


    subroutine apply_adiabatic_forcing(ldt)
        real(dp), intent(in) :: ldt
        real(dp) :: rho_air, qv_mean, cp_m, dT_adi
        integer :: k

        parcel_velocity = get_velocity(time)
        if (abs(parcel_velocity) < 1.0e-30) return

        parcel_height = parcel_height + parcel_velocity * ldt

        rho_air = pres / (Rd * sum(Tv) / N)
        pres = pres - rho_air * g * parcel_velocity * ldt

        qv_mean = sum(WV) / N
        cp_m = cp * (1.0 + cp_wv / cp * qv_mean) / (1.0 + qv_mean)
        dT_adi = -(g / cp_m) * parcel_velocity * ldt

        T(:) = T(:) + dT_adi
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

    end subroutine apply_adiabatic_forcing

end module dynamics
