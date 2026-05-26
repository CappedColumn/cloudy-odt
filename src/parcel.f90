module parcel
    use globals
    use netcdf
    use microphysics, only: virtual_temp, update_supersat, saturation_mixing_ratio
    implicit none

    private
    public :: initialize_parcel, apply_adiabatic_forcing, apply_entrainment, &
              do_parcel_ascent, parcel_height, parcel_velocity, &
              parcel_file, initial_RH, pressure_limit, pressure_limit_reached, &
              do_entrainment, ent_rate, n_blob, psigma, random_entrainment

    ! --- PARCEL namelist variables ---
    character(512) :: parcel_file = ''
    real(dp) :: initial_RH = 1.0
    logical  :: do_entrainment = .false.
    real(dp) :: ent_rate = 2.0
    integer(i4) :: n_blob = 1
    real(dp) :: psigma = 0.1
    logical  :: random_entrainment = .true.
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
    real(dp) :: t_next_entrain = 0.0

contains

    subroutine initialize_parcel()
        integer :: nml_unit, ierr, i
        character(256) :: nml_line, io_emsg

        namelist /PARCEL/ parcel_file, initial_RH, pressure_limit, &
            do_entrainment, ent_rate, n_blob, psigma, random_entrainment

        ! --- Read PARCEL namelist ---
        write(*,*) 'Reading PARCEL namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, &
             action='read', status='old')
        if (ierr /= 0) then
            write(0,*) io_emsg; stop 1
        end if
        read(nml=PARCEL, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(0,'(a)') 'Invalid PARCEL namelist parameter: '//trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        ! --- Validate ---
        if (initial_RH < 0.0 .or. initial_RH > 1.0) then
            write(0,*) 'Error: initial_RH must be between 0 and 1, got: ', initial_RH
            stop 1
        end if
        if (do_entrainment) then
            if (psigma <= 0.0 .or. psigma >= 1.0) then
                write(0,*) 'Error: psigma must be in (0, 1), got: ', psigma
                stop 1
            end if
            if (psigma * n_blob >= 1.0) then
                write(0,*) 'Error: psigma * n_blob must be < 1'
                stop 1
            end if
            if (ent_rate <= 0.0) then
                write(0,*) 'Error: ent_rate must be > 0, got: ', ent_rate
                stop 1
            end if
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
            write(0,*) 'Error: expected CODT_parcel_input_v1 or v2, got: ', trim(conventions)
            stop 1
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

        ! --- Environmental profile (entrainment) ---
        if (do_entrainment) then
            if (trim(conventions) /= 'CODT_parcel_input_v2') then
                write(0,*) 'Error: do_entrainment requires CODT_parcel_input_v2'
                stop 1
            end if
            call nc_verify(nf90_inq_dimid(dyn_ncid, 'level', dimid), 'finding level dim')
            call nc_verify(nf90_inquire_dimension(dyn_ncid, dimid, len=n_env_levels), 'reading level dim')

            allocate(env_pressure(n_env_levels))
            allocate(env_temperature(n_env_levels))
            allocate(env_RH(n_env_levels))

            call nc_verify(nf90_inq_varid(dyn_ncid, 'env_pressure', varid), 'finding env_pressure')
            call nc_verify(nf90_get_var(dyn_ncid, varid, env_pressure), 'reading env_pressure')

            call nc_verify(nf90_inq_varid(dyn_ncid, 'env_temperature', varid), 'finding env_temperature')
            call nc_verify(nf90_get_var(dyn_ncid, varid, env_temperature), 'reading env_temperature')

            call nc_verify(nf90_inq_varid(dyn_ncid, 'env_RH', varid), 'finding env_RH')
            call nc_verify(nf90_get_var(dyn_ncid, varid, env_RH), 'reading env_RH')

            do i = 2, n_env_levels
                if (env_pressure(i) >= env_pressure(i-1)) then
                    write(0,*) 'Error: env_pressure must be monotonically decreasing'
                    stop 1
                end if
            end do

            t_next_entrain = compute_dt_entm()
        end if

        call nc_verify(nf90_close(dyn_ncid), 'closing parcel file')

        ! --- Log ---
        write(*,*) '--- Parcel Configuration ---'
        write(*,*) 'n_segments:        ', n_segments
        write(*,'(a,f8.2,a)')  '  initial velocity:  ', segment_velocity(1), ' m/s'
        write(*,'(a,f8.1,a)')  '  initial pressure:  ', pres / Pa_per_mb, ' mb'
        write(*,'(a,f8.3)')    '  initial RH:        ', initial_RH
        if (do_entrainment) then
            write(*,*) 'entrainment:       ON'
            write(*,*) '  ent_rate:        ', ent_rate
            write(*,*) '  n_blob:          ', n_blob
            write(*,*) '  psigma:          ', psigma
            write(*,*) '  random_entrain:  ', random_entrainment
            write(*,*) '  env levels:      ', n_env_levels
            write(*,*) '  first dt_entm:   ', t_next_entrain, ' s'
        end if
        write(*,*) '----------------------------'

    end subroutine read_parcel_file


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
            write(0,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
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


    subroutine apply_entrainment()
        integer, parameter :: max_blobs = 10
        integer :: blob_start(max_blobs), blob_end(max_blobs), n_final
        real(dp) :: T_env, qv_env
        integer :: i, k

        if (abs(parcel_velocity) < 1.0e-30) return
        if (time < t_next_entrain) return

        call interp_env(pres, T_env, qv_env)
        call place_blobs(blob_start, blob_end, n_final)

        do i = 1, n_final
            do k = blob_start(i), blob_end(i)
                T(k) = T_env
                WV(k) = qv_env
            end do
        end do

        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

        t_next_entrain = time + compute_dt_entm()

    end subroutine apply_entrainment


    function compute_dt_entm() result(dt_entm)
        real(dp) :: dt_entm, u

        dt_entm = (real(n_blob, dp) / ent_rate) * (psigma / (1.0 - psigma)) &
                  / abs(parcel_velocity)

        if (random_entrainment) then
            call random_number(u)
            dt_entm = dt_entm * (-log(1.0 - u))
        end if
    end function compute_dt_entm


    subroutine place_blobs(starts, ends, n_final)
        integer, intent(out) :: starts(:), ends(:), n_final
        integer :: xn, blob_size, s, i
        real(dp) :: u

        blob_size = int(psigma * N)
        xn = N - blob_size * n_blob

        call random_number(u)
        s = int(u * xn) + 1

        n_final = 0
        do i = 1, n_blob
            if (i > 1) s = ends(n_final) + int(real(xn, dp) / n_blob) + 1

            if (s <= N .and. s + blob_size - 1 > N) then
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = N
                n_final = n_final + 1
                starts(n_final) = 1
                ends(n_final) = blob_size - (N - s + 1)
            else if (s > N) then
                n_final = n_final + 1
                starts(n_final) = mod(s - 1, N) + 1
                ends(n_final) = mod(s - 1 + blob_size - 1, N) + 1
                if (ends(n_final) < starts(n_final)) then
                    ends(n_final) = N
                    n_final = n_final + 1
                    starts(n_final) = 1
                    ends(n_final) = blob_size - (N - starts(n_final - 1) + 1)
                end if
            else
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = s + blob_size - 1
            end if
        end do
    end subroutine place_blobs


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
