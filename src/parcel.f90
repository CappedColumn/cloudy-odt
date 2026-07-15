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
              parcel_file, initial_RH, pressure_limit, pressure_limit_reached, &
              pressure_mode, parcel_height_env, write_height_env

    ! --- PARCEL namelist variables ---
    character(512) :: parcel_file = ''
    real(dp) :: initial_RH = 1.0
    real(dp) :: pressure_limit = 0.0
    ! Pressure evolution: 'hydrostatic' integrates dp = -rho_parcel*g*w*dt (v1/v2
    ! behavior); 'environment' (v3 only) sets pres = p_env(parcel_height) from the
    ! sounding, so the parcel stays on the environment's pressure-height curve.
    character(16) :: pressure_mode = 'hydrostatic'
    ! v3 only: which vertical coordinate the file's segment_coord values are on
    ! ('height' [m] or 'pressure' [Pa]); governs velocity/entrainment segment
    ! changes. Ignored for v1/v2 (their segments are in time).
    character(16) :: vertical_axis = 'height'

    ! --- Velocity segments ---
    ! Piecewise-constant velocity on segment_coords, whose meaning is set by
    ! segment_axis: time [s] for input v1/v2, height [m] or pressure [Pa] for v3.
    integer(i4) :: n_segments
    integer(i4) :: segment_axis = AXIS_TIME
    real(dp), allocatable :: segment_coords(:)
    real(dp), allocatable :: segment_velocity(:)

    ! --- Parcel state ---
    real(dp) :: parcel_height   = 0.0
    real(dp) :: parcel_velocity = 0.0
    logical  :: do_parcel_ascent = .false.
    logical  :: pressure_limit_reached = .false.
    ! Diagnostic: environment height at the parcel's current pressure (v3 +
    ! hydrostatic mode). Its drift from parcel_height quantifies how far the
    ! self-integrated pressure has left the sounding's p(z) curve.
    real(dp) :: parcel_height_env = 0.0
    logical  :: write_height_env = .false.

    ! --- Environmental profile (v3: always required; v2: entrainment only) ---
    integer(i4) :: n_env_levels
    real(dp), allocatable :: env_pressure(:)
    real(dp), allocatable :: env_temperature(:)
    real(dp), allocatable :: env_RH(:)
    real(dp), allocatable :: env_height(:)   ! v3 only

contains

    subroutine initialize_parcel()
        ! Reads &PARCEL, sets the uniform initial parcel state (T = Tref, WV from
        ! initial_RH at saturation), then loads the velocity profile (and, if
        ! entraining, the environmental sounding) from the parcel NetCDF file.
        integer :: nml_unit, ierr, i
        character(256) :: nml_line, io_emsg

        namelist /PARCEL/ parcel_file, initial_RH, pressure_limit, pressure_mode, &
                          vertical_axis

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
        if (trim(pressure_mode) /= 'hydrostatic' .and. trim(pressure_mode) /= 'environment') then
            write(error_unit,*) "Error: pressure_mode must be 'hydrostatic' or 'environment', got: ", &
                trim(pressure_mode)
            call exit(1)
        end if
        if (trim(vertical_axis) /= 'height' .and. trim(vertical_axis) /= 'pressure') then
            write(error_unit,*) "Error: vertical_axis must be 'height' or 'pressure', got: ", &
                trim(vertical_axis)
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

        ! --- Segment coordinate + velocity ---
        ! v1/v2: piecewise-constant velocity in time. v3: velocity (and the
        ! optional entrainment schedule) on a vertical coordinate — exactly one
        ! of 'height' [m, ascending] or 'pressure' [Pa, descending].
        call nc_verify(nf90_inq_dimid(dyn_ncid, 'segment', dimid), 'finding segment dim')
        call nc_verify(nf90_inquire_dimension(dyn_ncid, dimid, len=n_segments), 'reading segment dim')

        ! Allow re-initialization (unit tests exercise multiple configurations)
        if (allocated(segment_coords)) deallocate(segment_coords, segment_velocity)

        allocate(segment_coords(n_segments))
        if (trim(conventions) == 'CODT_parcel_input_v3') then
            ! The &PARCEL vertical_axis says whether segment_coord holds heights
            ! [m] or pressures [Pa]; the sounding always carries both env_height
            ! and env_pressure.
            if (trim(vertical_axis) == 'pressure') then
                segment_axis = AXIS_PRESSURE
            else
                segment_axis = AXIS_HEIGHT
            end if
            call nc_verify(nf90_inq_varid(dyn_ncid, 'segment_coord', varid), &
                           'finding segment_coord')
            call nc_verify(nf90_get_var(dyn_ncid, varid, segment_coords), &
                           'reading segment_coord')
        else
            segment_axis = AXIS_TIME
            call nc_verify(nf90_inq_varid(dyn_ncid, 'time', varid), 'finding time')
            call nc_verify(nf90_get_var(dyn_ncid, varid, segment_coords), 'reading time')
        end if

        allocate(segment_velocity(n_segments))
        call nc_verify(nf90_inq_varid(dyn_ncid, 'velocity', varid), 'finding velocity')
        call nc_verify(nf90_get_var(dyn_ncid, varid, segment_velocity), 'reading velocity')

        call validate_segment_coords()

        if (trim(pressure_mode) == 'environment' .and. &
            trim(conventions) /= 'CODT_parcel_input_v3') then
            write(error_unit,*) "Error: pressure_mode = 'environment' requires a v3 parcel file"
            call exit(1)
        end if

        do_parcel_ascent = .true.
        parcel_velocity = segment_velocity(1)

        ! --- Environmental profile ---
        ! v3 always carries the sounding (env_height/pressure/temperature/RH):
        ! pressure_mode = 'environment' needs p(z) even without entrainment, and
        ! the parcel_height_env diagnostic needs it in hydrostatic mode.
        ! v2 carries it only for entrainment.
        if (trim(conventions) == 'CODT_parcel_input_v3') then
            call load_env_profile(dyn_ncid, require_height=.true.)
            write_height_env = (trim(pressure_mode) == 'hydrostatic')
            if (do_entrainment) call load_entrainment_schedule(dyn_ncid)
        else if (do_entrainment) then
            if (trim(conventions) /= 'CODT_parcel_input_v2') then
                write(error_unit,*) 'Error: do_entrainment requires CODT_parcel_input_v2 or v3'
                call exit(1)
            end if
            call load_env_profile(dyn_ncid, require_height=.false.)
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


    subroutine validate_segment_coords()
        ! Time starts at 0; all axes must be strictly monotonic (time/height
        ! increasing, pressure decreasing). Height and pressure need not start
        ! at the parcel's launch value — lookups clamp, so the first segment
        ! covers everything on its side of the profile.
        integer :: i

        select case (segment_axis)
        case (AXIS_PRESSURE)
            do i = 2, n_segments
                if (segment_coords(i) >= segment_coords(i-1)) then
                    write(error_unit,*) 'Error: segment pressures must be strictly decreasing'
                    call exit(1)
                end if
            end do
        case default   ! AXIS_TIME, AXIS_HEIGHT
            if (segment_axis == AXIS_TIME .and. abs(segment_coords(1)) > 1.0e-10) then
                write(error_unit,*) 'Error: first segment time must be 0, got: ', &
                    segment_coords(1)
                call exit(1)
            end if
            do i = 2, n_segments
                if (segment_coords(i) <= segment_coords(i-1)) then
                    write(error_unit,*) 'Error: segment coordinates must be strictly increasing'
                    call exit(1)
                end if
            end do
        end select

    end subroutine validate_segment_coords


    subroutine load_env_profile(lncid, require_height)
        integer, intent(in) :: lncid
        logical, intent(in) :: require_height   ! v3: sounding must carry env_height
        integer :: varid, dimid, i

        call nc_verify(nf90_inq_dimid(lncid, 'level', dimid), 'finding level dim')
        call nc_verify(nf90_inquire_dimension(lncid, dimid, len=n_env_levels), 'reading level dim')

        ! Allow re-initialization (unit tests exercise multiple configurations)
        if (allocated(env_pressure)) deallocate(env_pressure, env_temperature, env_RH)
        if (allocated(env_height)) deallocate(env_height)

        allocate(env_pressure(n_env_levels))
        allocate(env_temperature(n_env_levels))
        allocate(env_RH(n_env_levels))

        call nc_verify(nf90_inq_varid(lncid, 'env_pressure', varid), 'finding env_pressure')
        call nc_verify(nf90_get_var(lncid, varid, env_pressure), 'reading env_pressure')

        call nc_verify(nf90_inq_varid(lncid, 'env_temperature', varid), 'finding env_temperature')
        call nc_verify(nf90_get_var(lncid, varid, env_temperature), 'reading env_temperature')

        call nc_verify(nf90_inq_varid(lncid, 'env_RH', varid), 'finding env_RH')
        call nc_verify(nf90_get_var(lncid, varid, env_RH), 'reading env_RH')

        do i = 2, n_env_levels
            if (env_pressure(i) >= env_pressure(i-1)) then
                write(error_unit,*) 'Error: env_pressure must be monotonically decreasing'
                call exit(1)
            end if
        end do

        if (require_height) then
            allocate(env_height(n_env_levels))
            call nc_verify(nf90_inq_varid(lncid, 'env_height', varid), 'finding env_height')
            call nc_verify(nf90_get_var(lncid, varid, env_height), 'reading env_height')
            do i = 2, n_env_levels
                if (env_height(i) <= env_height(i-1)) then
                    write(error_unit,*) 'Error: env_height must be monotonically increasing'
                    call exit(1)
                end if
            end do
        end if

    end subroutine load_env_profile


    subroutine load_entrainment_schedule(lncid)
        ! Reads the per-segment entrainment parameters (parcel input v3) on the
        ! shared vertical segment coordinate and hands them to the entrainment
        ! module. The schedule variables are optional: when absent, the constant
        ! &ENTRAINMENT namelist values apply for the whole run.
        ! File ent_rate is in 1/km; internal physics uses 1/m.
        integer, intent(in) :: lncid
        integer :: varid, ierr
        real(dp), allocatable :: sched_ent_rate(:), sched_psigma(:)
        integer(i4), allocatable :: sched_n_blob(:)

        ierr = nf90_inq_varid(lncid, 'ent_rate', varid)
        if (ierr /= NF90_NOERR) then
            call initialize_entrainment(parcel_velocity)
            return
        end if

        allocate(sched_ent_rate(n_segments), sched_n_blob(n_segments), &
                 sched_psigma(n_segments))

        call nc_verify(nf90_get_var(lncid, varid, sched_ent_rate), 'reading ent_rate')
        sched_ent_rate = sched_ent_rate / m_per_km   ! 1/km -> 1/m

        call nc_verify(nf90_inq_varid(lncid, 'n_blob', varid), 'finding n_blob')
        call nc_verify(nf90_get_var(lncid, varid, sched_n_blob), 'reading n_blob')

        call nc_verify(nf90_inq_varid(lncid, 'psigma', varid), 'finding psigma')
        call nc_verify(nf90_get_var(lncid, varid, sched_psigma), 'reading psigma')

        call initialize_entrainment(parcel_velocity, segment_axis, segment_coords, &
                                    sched_ent_rate, sched_n_blob, sched_psigma)

    end subroutine load_entrainment_schedule


    ! Current value of the segment coordinate for schedule lookups: simulation
    ! time (v1/v2) or the parcel's height/pressure (v3).
    pure function parcel_axis_query() result(query)
        real(dp) :: query

        select case (segment_axis)
        case (AXIS_HEIGHT)
            query = parcel_height
        case (AXIS_PRESSURE)
            query = pres
        case default
            query = time
        end select
    end function parcel_axis_query


    pure function get_velocity(query) result(vel)
        ! Piecewise-constant ascent velocity: returns the velocity of the segment
        ! containing query on the segment axis. Time/height segments are sorted
        ! ascending, pressure descending; either way segment i spans from
        ! segment_coords(i) toward segment_coords(i+1).
        real(dp), intent(in) :: query   ! time [s], height [m], or pressure [Pa]
        real(dp) :: vel                 ! ascent velocity (m/s)
        integer :: i

        vel = segment_velocity(1)
        do i = 2, n_segments
            if (segment_axis == AXIS_PRESSURE) then
                if (query > segment_coords(i)) exit
            else
                if (query < segment_coords(i)) exit
            end if
            vel = segment_velocity(i)
        end do
    end function get_velocity


    pure function interp_profile(coords, values, query) result(v)
        ! Clamped piecewise-linear lookup on a strictly monotonic coordinate
        ! array (ascending, e.g. env_height, or descending, e.g. env_pressure).
        ! Used both ways on the v3 sounding: p_env(z) = interp_profile(env_height,
        ! env_pressure, z) and z_env(p) = interp_profile(env_pressure, env_height, p).
        real(dp), intent(in) :: coords(:)   ! monotonic coordinate array
        real(dp), intent(in) :: values(:)   ! values on the same levels
        real(dp), intent(in) :: query       ! coordinate to interpolate at
        real(dp) :: v
        real(dp) :: frac, direction
        integer :: k, n_levels

        n_levels = size(coords)
        direction = merge(1.0, -1.0, coords(n_levels) > coords(1))   ! +1 ascending, -1 descending

        if (direction * (query - coords(1)) <= 0.0) then
            v = values(1)
        else if (direction * (query - coords(n_levels)) >= 0.0) then
            v = values(n_levels)
        else
            do k = 2, n_levels
                if (direction * (query - coords(k)) < 0.0) exit
            end do
            frac = (query - coords(k-1)) / (coords(k) - coords(k-1))
            v = values(k-1) + frac * (values(k) - values(k-1))
        end if
    end function interp_profile


    subroutine apply_adiabatic_forcing(ldt)
        ! Advances the parcel one step: raise height, update pressure, and change
        ! temperature adiabatically. Pressure evolves per pressure_mode:
        !   'hydrostatic'  — self-integration dp = -rho_parcel*g*w*dt with cooling
        !                    dT = -(g/cp_moist)*w*dt (v1/v2 behavior, default);
        !   'environment'  — pres = p_env(parcel_height) from the sounding (v3),
        !                    with dT = (Rd*Tv/(cp_m*p))*dp from the actual dp
        !                    (identical to the hydrostatic form when the sounding
        !                    is hydrostatic in the parcel's Tv, correct otherwise;
        !                    reversible under descent since p is a function of z).
        ! WV is unchanged here (condensation is handled by droplet growth); only T,
        ! Tv, SS and pres update. Sets pressure_limit_reached and returns early if
        ! the parcel reaches the target pressure.
        real(dp), intent(in) :: ldt   ! time step (s)
        real(dp) :: rho_air, qv_mean, cp_m, dT_adi   ! air density; mean qv; moist cp; adiabatic dT
        real(dp) :: pres_old
        integer :: k

        parcel_velocity = get_velocity(parcel_axis_query())
        if (abs(parcel_velocity) < 1.0e-30) return

        parcel_height = parcel_height + parcel_velocity * ldt

        qv_mean = sum(WV) / N
        cp_m = cp * (1.0 + cp_wv / cp * qv_mean) / (1.0 + qv_mean)

        if (trim(pressure_mode) == 'environment') then
            ! Follow the environment's pressure-height curve
            pres_old = pres
            pres = interp_profile(env_height, env_pressure, parcel_height)
            dT_adi = (Rd * (sum(Tv) / N) / (cp_m * pres_old)) * (pres - pres_old)
        else
            ! Hydrostatic pressure change over the height ascended this step
            rho_air = pres / (Rd * sum(Tv) / N)
            pres = pres - rho_air * g * parcel_velocity * ldt
            dT_adi = -(g / cp_m) * parcel_velocity * ldt
        end if

        if (pressure_limit > 0.0 .and. pres <= pressure_limit) then
            pres = pressure_limit
            write(*,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            write(error_unit,'(a,f8.1,a)') ' Pressure limit reached: ', pres / Pa_per_mb, ' mb. Stopping.'
            pressure_limit_reached = .true.
            return
        end if

        T(:) = T(:) + dT_adi
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

        ! Diagnostic: where this pressure sits in the environment (v3 hydrostatic
        ! mode). Drift from parcel_height measures the departure of the
        ! self-integrated pressure from the sounding's p(z).
        if (write_height_env) then
            parcel_height_env = interp_profile(env_pressure, env_height, pres)
        end if

    end subroutine apply_adiabatic_forcing


    subroutine apply_parcel_entrainment()
        real(dp) :: T_env, qv_env

        call interp_env(pres, T_env, qv_env)
        call apply_entrainment(T_env, qv_env, parcel_velocity, parcel_axis_query())

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
