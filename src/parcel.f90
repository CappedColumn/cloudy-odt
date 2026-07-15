! Parcel-mode forcing: drives an adiabatically ascending air parcel (LEM mode).
! Reads a piecewise-constant vertical-velocity profile from a NetCDF parcel file,
! and each step lifts the parcel, drops its pressure hydrostatically, and cools it
! at the (moist-weighted) dry adiabatic rate. Optionally stops at pressure_limit.
! When do_entrainment is set, also owns the environmental sounding (env_*) and
! feeds interpolated environmental air to the entrainment module.
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
        ! Reads &PARCEL, sets the uniform initial parcel state (T = Tref, WV from
        ! initial_RH at saturation), then loads the velocity profile (and, if
        ! entraining, the environmental sounding) from the parcel NetCDF file.
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
            trim(conventions) /= 'CODT_parcel_input_v2' .and. &
            trim(conventions) /= 'CODT_parcel_input_v3') then
            write(error_unit,*) 'Error: expected CODT_parcel_input_v1, v2 or v3, got: ', trim(conventions)
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
            if (trim(conventions) /= 'CODT_parcel_input_v2' .and. &
                trim(conventions) /= 'CODT_parcel_input_v3') then
                write(error_unit,*) 'Error: do_entrainment requires CODT_parcel_input_v2 or v3'
                call exit(1)
            end if
            call load_env_profile(dyn_ncid)
            ! v3 adds per-segment entrainment parameters (time-varying schedule);
            ! v2 uses the constant &ENTRAINMENT namelist values.
            if (trim(conventions) == 'CODT_parcel_input_v3') then
                call load_entrainment_schedule(dyn_ncid)
            else
                call initialize_entrainment(parcel_velocity)
            end if
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


    subroutine load_entrainment_schedule(ncid)
        ! Reads the per-segment entrainment parameters (parcel input v3) on the
        ! shared velocity-segment axis and hands them to the entrainment module,
        ! which validates them and drives a time-varying schedule.
        integer, intent(in) :: ncid
        integer :: varid
        real(dp), allocatable :: sched_ent_rate(:), sched_psigma(:)
        integer(i4), allocatable :: sched_n_blob(:)

        allocate(sched_ent_rate(n_segments), sched_n_blob(n_segments), &
                 sched_psigma(n_segments))

        call nc_verify(nf90_inq_varid(ncid, 'ent_rate', varid), 'finding ent_rate')
        call nc_verify(nf90_get_var(ncid, varid, sched_ent_rate), 'reading ent_rate')

        call nc_verify(nf90_inq_varid(ncid, 'n_blob', varid), 'finding n_blob')
        call nc_verify(nf90_get_var(ncid, varid, sched_n_blob), 'reading n_blob')

        call nc_verify(nf90_inq_varid(ncid, 'psigma', varid), 'finding psigma')
        call nc_verify(nf90_get_var(ncid, varid, sched_psigma), 'reading psigma')

        call initialize_entrainment(parcel_velocity, segment_times, sched_ent_rate, &
                                    sched_n_blob, sched_psigma)

    end subroutine load_entrainment_schedule


    pure function get_velocity(query_time) result(vel)
        ! Piecewise-constant ascent velocity: returns the velocity of the segment
        ! whose start time most recently preceded query_time (segments are sorted).
        real(dp), intent(in) :: query_time   ! simulation time (s)
        real(dp) :: vel                      ! ascent velocity (m/s)
        integer :: i

        vel = segment_velocity(1)
        do i = 2, n_segments
            if (query_time < segment_times(i)) exit
            vel = segment_velocity(i)
        end do
    end function get_velocity


    subroutine apply_adiabatic_forcing(ldt)
        ! Advances the parcel one step: raise height, drop pressure hydrostatically
        ! (dp = -rho*g*w*dt), and cool dry-adiabatically (dT = -(g/cp_moist)*w*dt).
        ! WV is unchanged here (condensation is handled by droplet growth); only T,
        ! Tv, SS and pres update. Sets pressure_limit_reached and returns early if
        ! the parcel reaches the target pressure.
        real(dp), intent(in) :: ldt   ! time step (s)
        real(dp) :: rho_air, qv_mean, cp_m, dT_adi   ! air density; mean qv; moist cp; adiabatic dT
        integer :: k

        parcel_velocity = get_velocity(time)
        if (abs(parcel_velocity) < 1.0e-30) return

        parcel_height = parcel_height + parcel_velocity * ldt

        ! Hydrostatic pressure change over the height ascended this step
        rho_air = pres / (Rd * sum(Tv) / N)
        pres = pres - rho_air * g * parcel_velocity * ldt

        if (pressure_limit > 0.0 .and. pres <= pressure_limit) then
            pres = pressure_limit
            write(*,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            write(error_unit,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            pressure_limit_reached = .true.
            return
        end if

        ! Dry adiabatic cooling at the moist-weighted heat capacity cp_m
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
        ! Linearly interpolates the environmental sounding (T, RH vs pressure) to
        ! the parcel's current pressure, clamping at the profile ends, then converts
        ! RH to a vapor mixing ratio. Supplies the "environment" for entrainment.
        real(dp), intent(in) :: p_current    ! current parcel pressure (Pa)
        real(dp), intent(out) :: T_env       ! interpolated environmental temperature (K)
        real(dp), intent(out) :: qv_env      ! environmental vapor mixing ratio (kg/kg)
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
