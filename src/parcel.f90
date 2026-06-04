module parcel
    use globals
    use netcdf
    use microphysics, only: virtual_temp, update_supersat, saturation_mixing_ratio
    use entrainment, only: initialize_entrainment, apply_entrainment
    implicit none

    private
    public :: initialize_parcel, apply_adiabatic_forcing, apply_parcel_entrainment, &
              do_parcel_ascent, parcel_height, parcel_velocity, &
              parcel_file, initial_RH, pressure_limit, pressure_limit_reached

    ! --- PARCEL namelist variables ---
    character(512) :: parcel_file = ''
    real(dp) :: initial_RH = 1.0
    real(dp) :: pressure_limit = 0.0

    ! --- Velocity segments ---
    integer(i4) :: n_segments
    real(dp), allocatable :: segment_times(:)
    real(dp), allocatable :: segment_velocity(:)

    ! --- Parcel state ---
    real(dp) :: parcel_height   = 0.0
    real(dp) :: parcel_velocity = 0.0
    logical  :: do_parcel_ascent = .false.
    logical  :: pressure_limit_reached = .false.

    ! --- Environmental profile for entrainment ---
    integer(i4) :: n_env_levels
    real(dp), allocatable :: env_pressure(:)
    real(dp), allocatable :: env_temperature(:)
    real(dp), allocatable :: env_RH(:)

contains

    subroutine initialize_parcel()
        integer :: nml_unit, ierr, i
        character(256) :: nml_line, io_emsg

        namelist /PARCEL/ parcel_file, initial_RH, pressure_limit

        ! --- Read PARCEL namelist ---
        write(*,*) 'Reading PARCEL namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, &
             action='read', status='old')
        if (ierr /= 0) then
            write(error_unit,*) io_emsg; call exit(1)
        end if
        read(nml=PARCEL, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) call namelist_read_error(nml_unit, 'PARCEL')
        close(nml_unit)

        ! --- Validate ---
        if (initial_RH < 0.0 .or. initial_RH > 1.0) then
            write(error_unit,*) 'Error: initial_RH must be between 0 and 1, got: ', initial_RH
            call exit(1)
        end if

        ! --- Initialize parcel arrays ---
        T(:) = Tref
        WV(:) = initial_RH * saturation_mixing_ratio(Tref, pres)
        do i = 1, N
            Tv(i) = virtual_temp(T(i), WV(i))
        end do
        call update_supersat(T, WV, SS, pres)

        ! --- Read parcel input file ---
        if (parcel_file /= '') then
            call read_parcel_file(resolve_path(namelist_dir, parcel_file))
        end if

    end subroutine initialize_parcel


    subroutine read_parcel_file(filepath)
        character(*), intent(in) :: filepath
        integer :: dyn_ncid, varid, dimid, i
        character(64) :: conventions

        write(*,*) 'Reading parcel data from: ', trim(filepath)

        call nc_verify(nf90_open(trim(filepath), NF90_NOWRITE, dyn_ncid), &
                       'opening parcel file')
        call nc_verify(nf90_get_att(dyn_ncid, NF90_GLOBAL, 'conventions', conventions), &
                       'reading conventions attribute')

        if (trim(conventions) /= 'CODT_parcel_input_v1' .and. &
            trim(conventions) /= 'CODT_parcel_input_v2') then
            write(error_unit,*) 'Error: expected CODT_parcel_input_v1 or v2, got: ', trim(conventions)
            call exit(1)
        end if

        ! --- Velocity segments ---
        call nc_verify(nf90_inq_dimid(dyn_ncid, 'segment', dimid), 'finding segment dim')
        call nc_verify(nf90_inquire_dimension(dyn_ncid, dimid, len=n_segments), 'reading segment dim')

        allocate(segment_times(n_segments))
        call nc_verify(nf90_inq_varid(dyn_ncid, 'time', varid), 'finding time')
        call nc_verify(nf90_get_var(dyn_ncid, varid, segment_times), 'reading time')

        allocate(segment_velocity(n_segments))
        call nc_verify(nf90_inq_varid(dyn_ncid, 'velocity', varid), 'finding velocity')
        call nc_verify(nf90_get_var(dyn_ncid, varid, segment_velocity), 'reading velocity')

        if (abs(segment_times(1)) > 1.0e-10) then
            write(error_unit,*) 'Error: first segment time must be 0, got: ', segment_times(1)
            call exit(1)
        end if
        do i = 2, n_segments
            if (segment_times(i) <= segment_times(i-1)) then
                write(error_unit,*) 'Error: segment times must be monotonically increasing'
                call exit(1)
            end if
        end do

        do_parcel_ascent = .true.
        parcel_velocity = segment_velocity(1)

        ! --- Environmental profile (entrainment) ---
        if (do_entrainment) then
            if (trim(conventions) /= 'CODT_parcel_input_v2') then
                write(error_unit,*) 'Error: do_entrainment requires CODT_parcel_input_v2'
                call exit(1)
            end if
            call load_env_profile(dyn_ncid)
            call initialize_entrainment(parcel_velocity)
        end if

        call nc_verify(nf90_close(dyn_ncid), 'closing parcel file')

        ! --- Log ---
        write(*,*) '--- Parcel Configuration ---'
        write(*,*) 'n_segments:        ', n_segments
        write(*,'(a,f8.2,a)')  '  initial velocity:  ', segment_velocity(1), ' m/s'
        write(*,'(a,f8.1,a)')  '  initial pressure:  ', pres / Pa_per_mb, ' mb'
        write(*,'(a,f8.3)')    '  initial RH:        ', initial_RH
        if (do_entrainment) then
            write(*,*) '  env levels:      ', n_env_levels
        end if
        write(*,*) '----------------------------'

    end subroutine read_parcel_file


    subroutine load_env_profile(ncid)
        integer, intent(in) :: ncid
        integer :: varid, dimid, i

        call nc_verify(nf90_inq_dimid(ncid, 'level', dimid), 'finding level dim')
        call nc_verify(nf90_inquire_dimension(ncid, dimid, len=n_env_levels), 'reading level dim')

        allocate(env_pressure(n_env_levels))
        allocate(env_temperature(n_env_levels))
        allocate(env_RH(n_env_levels))

        call nc_verify(nf90_inq_varid(ncid, 'env_pressure', varid), 'finding env_pressure')
        call nc_verify(nf90_get_var(ncid, varid, env_pressure), 'reading env_pressure')

        call nc_verify(nf90_inq_varid(ncid, 'env_temperature', varid), 'finding env_temperature')
        call nc_verify(nf90_get_var(ncid, varid, env_temperature), 'reading env_temperature')

        call nc_verify(nf90_inq_varid(ncid, 'env_RH', varid), 'finding env_RH')
        call nc_verify(nf90_get_var(ncid, varid, env_RH), 'reading env_RH')

        do i = 2, n_env_levels
            if (env_pressure(i) >= env_pressure(i-1)) then
                write(error_unit,*) 'Error: env_pressure must be monotonically decreasing'
                call exit(1)
            end if
        end do

    end subroutine load_env_profile


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

        if (pressure_limit > 0.0 .and. pres <= pressure_limit) then
            pres = pressure_limit
            write(*,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            write(error_unit,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            pressure_limit_reached = .true.
            return
        end if

        qv_mean = sum(WV) / N
        cp_m = cp * (1.0 + cp_wv / cp * qv_mean) / (1.0 + qv_mean)
        dT_adi = -(g / cp_m) * parcel_velocity * ldt

        T(:) = T(:) + dT_adi
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

    end subroutine apply_adiabatic_forcing


    subroutine apply_parcel_entrainment()
        real(dp) :: T_env, qv_env

        call interp_env(pres, T_env, qv_env)
        call apply_entrainment(T_env, qv_env, parcel_velocity)

    end subroutine apply_parcel_entrainment


    subroutine interp_env(p_current, T_env, qv_env)
        real(dp), intent(in) :: p_current
        real(dp), intent(out) :: T_env, qv_env
        real(dp) :: RH_env, frac
        integer :: k

        if (p_current >= env_pressure(1)) then
            T_env = env_temperature(1)
            RH_env = env_RH(1)
        else if (p_current <= env_pressure(n_env_levels)) then
            T_env = env_temperature(n_env_levels)
            RH_env = env_RH(n_env_levels)
        else
            do k = 2, n_env_levels
                if (p_current > env_pressure(k)) exit
            end do
            frac = (p_current - env_pressure(k)) / (env_pressure(k-1) - env_pressure(k))
            T_env = env_temperature(k) + frac * (env_temperature(k-1) - env_temperature(k))
            RH_env = env_RH(k) + frac * (env_RH(k-1) - env_RH(k))
        end if

        qv_env = RH_env * saturation_mixing_ratio(T_env, p_current)
    end subroutine interp_env

end module parcel
