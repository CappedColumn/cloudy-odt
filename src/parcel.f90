! Parcel-mode forcing: drives an adiabatic air parcel (LEM mode) along a
! waypoint trajectory read from a v3 parcel NetCDF file. The trajectory is an
! ordered sequence of legs — "proceed to this level at this signed velocity" —
! so ascent histories that revisit levels (up, down, up again) are expressible.
! Completing the final leg ends the simulation (trajectory_complete). Pressure
! evolves per pressure_mode (hydrostatic self-integration, or following the
! sounding's p(z)). When present, the module also owns the environmental
! sounding (env_*) and feeds interpolated environmental air to the entrainment
! module. Support for the time-based v1/v2 parcel inputs was removed.
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
              pressure_mode, parcel_height_env, write_height_env, &
              initial_height, trajectory_complete, vertical_axis

    ! --- PARCEL namelist variables ---
    character(512) :: parcel_file = ''
    real(dp) :: initial_RH = 1.0
    real(dp) :: pressure_limit = 0.0
    ! Pressure evolution: 'hydrostatic' integrates dp = -rho_parcel*g*w*dt;
    ! 'environment' sets pres = p_env(parcel_height) from the sounding, so the
    ! parcel stays on the environment's pressure-height curve (requires the
    ! sounding; initial pres is then taken from p_env(initial_height)).
    character(16) :: pressure_mode = 'hydrostatic'
    ! Which vertical coordinate the file's segment_coord (leg target) values are
    ! on: 'height' [m] or 'pressure' [Pa].
    character(16) :: vertical_axis = 'height'
    ! Launch height of the parcel [m] (leg-1 direction is validated against it;
    ! on the pressure axis the launch level is the initial pres instead).
    real(dp) :: initial_height = 0.0

    ! --- Trajectory legs ---
    ! Leg i: move toward leg_target(i) at signed leg_velocity(i). The active leg
    ! advances when its target is reached; completing the last leg sets
    ! trajectory_complete, which ends the run (main-loop exit, like
    ! pressure_limit_reached).
    integer(i4) :: n_legs
    integer(i4) :: leg_axis = AXIS_HEIGHT
    integer(i4) :: current_leg = 1
    real(dp), allocatable :: leg_target(:)
    real(dp), allocatable :: leg_velocity(:)
    logical :: trajectory_complete = .false.

    ! --- Parcel state ---
    real(dp) :: parcel_height   = 0.0
    real(dp) :: parcel_velocity = 0.0
    logical  :: do_parcel_ascent = .false.
    logical  :: pressure_limit_reached = .false.
    ! Diagnostic: environment height at the parcel's current pressure
    ! (hydrostatic mode with a sounding). Its drift from parcel_height
    ! quantifies how far the self-integrated pressure has left the sounding's
    ! p(z) curve.
    real(dp) :: parcel_height_env = 0.0
    logical  :: write_height_env = .false.

    ! --- Environmental sounding (optional) ---
    ! Required when do_entrainment or pressure_mode = 'environment'; a bare
    ! adiabatic hydrostatic parcel runs without one.
    logical :: have_sounding = .false.
    integer(i4) :: n_env_levels
    real(dp), allocatable :: env_pressure(:)
    real(dp), allocatable :: env_temperature(:)
    real(dp), allocatable :: env_RH(:)
    real(dp), allocatable :: env_height(:)

contains

    subroutine initialize_parcel()
        ! Reads &PARCEL, loads the trajectory (and sounding, if present) from the
        ! parcel NetCDF file, then sets the uniform initial parcel state
        ! (T = Tref, WV from initial_RH at saturation). The file is read before
        ! the state init so that pressure_mode = 'environment' can place the
        ! initial pressure on the sounding at initial_height.
        integer :: nml_unit, ierr, i
        character(256) :: nml_line, io_emsg

        namelist /PARCEL/ parcel_file, initial_RH, pressure_limit, pressure_mode, &
                          vertical_axis, initial_height

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

        parcel_height = initial_height

        ! --- Read parcel input file (legs, sounding, entrainment schedule) ---
        if (parcel_file /= '') then
            call read_parcel_file(resolve_path(namelist_dir, parcel_file))
        end if

        ! --- Initialize parcel state (pres is final here in both modes) ---
        T(:) = Tref
        WV(:) = initial_RH * saturation_mixing_ratio(Tref, pres)
        do i = 1, N
            Tv(i) = virtual_temp(T(i), WV(i))
        end do
        call update_supersat(T, WV, SS, pres)

    end subroutine initialize_parcel


    subroutine read_parcel_file(filepath)
        character(*), intent(in) :: filepath
        integer :: dyn_ncid, varid, ierr, dimid
        character(64) :: conventions
        logical :: need_sounding

        write(*,*) 'Reading parcel data from: ', trim(filepath)

        call nc_verify(nf90_open(trim(filepath), NF90_NOWRITE, dyn_ncid), &
                       'opening parcel file')
        call nc_verify(nf90_get_att(dyn_ncid, NF90_GLOBAL, 'conventions', conventions), &
                       'reading conventions attribute')

        if (trim(conventions) /= 'CODT_parcel_input_v3') then
            write(error_unit,*) 'Error: expected CODT_parcel_input_v3, got: ', &
                trim(conventions)
            write(error_unit,*) '  (v1/v2 parcel inputs are no longer supported; ' // &
                'regenerate the file in the v3 waypoint format)'
            call exit(1)
        end if

        ! --- Trajectory legs ---
        call nc_verify(nf90_inq_dimid(dyn_ncid, 'segment', dimid), 'finding segment dim')
        call nc_verify(nf90_inquire_dimension(dyn_ncid, dimid, len=n_legs), 'reading segment dim')

        ! Allow re-initialization (unit tests exercise multiple configurations)
        if (allocated(leg_target)) deallocate(leg_target, leg_velocity)

        allocate(leg_target(n_legs), leg_velocity(n_legs))

        if (trim(vertical_axis) == 'pressure') then
            leg_axis = AXIS_PRESSURE
        else
            leg_axis = AXIS_HEIGHT
        end if
        call nc_verify(nf90_inq_varid(dyn_ncid, 'segment_coord', varid), &
                       'finding segment_coord')
        call nc_verify(nf90_get_var(dyn_ncid, varid, leg_target), &
                       'reading segment_coord')
        call nc_verify(nf90_inq_varid(dyn_ncid, 'velocity', varid), 'finding velocity')
        call nc_verify(nf90_get_var(dyn_ncid, varid, leg_velocity), 'reading velocity')

        ! --- Optional environmental sounding ---
        ierr = nf90_inq_varid(dyn_ncid, 'env_pressure', varid)
        have_sounding = (ierr == NF90_NOERR)
        need_sounding = do_entrainment .or. trim(pressure_mode) == 'environment'
        if (need_sounding .and. .not. have_sounding) then
            write(error_unit,*) 'Error: the parcel file has no environmental sounding ' // &
                '(env_height/env_pressure/env_temperature/env_RH), which is required ' // &
                "when do_entrainment = .true. or pressure_mode = 'environment'"
            call exit(1)
        end if
        if (have_sounding) call load_env_profile(dyn_ncid)

        ! In environment mode the parcel starts on the sounding's p(z) curve.
        if (trim(pressure_mode) == 'environment') then
            pres = interp_profile(env_height, env_pressure, initial_height)
        end if
        write_height_env = have_sounding .and. trim(pressure_mode) == 'hydrostatic'

        call validate_legs()

        current_leg = 1
        trajectory_complete = .false.
        do_parcel_ascent = .true.
        parcel_velocity = leg_velocity(1)

        if (do_entrainment) call load_entrainment_schedule(dyn_ncid)

        call nc_verify(nf90_close(dyn_ncid), 'closing parcel file')

        ! --- Log ---
        write(*,*) '--- Parcel Configuration ---'
        write(*,*) 'trajectory legs:   ', n_legs
        write(*,'(a,f8.2,a)')  '  initial velocity:  ', leg_velocity(1), ' m/s'
        write(*,'(a,f8.1,a)')  '  initial height:    ', parcel_height, ' m'
        write(*,'(a,f8.1,a)')  '  initial pressure:  ', pres / Pa_per_mb, ' mb'
        write(*,'(a,f8.3)')    '  initial RH:        ', initial_RH
        if (have_sounding) then
            write(*,*) '  env levels:      ', n_env_levels
        end if
        write(*,*) '----------------------------'

    end subroutine read_parcel_file


    subroutine validate_legs()
        ! Each leg's signed velocity must point from the previous level (the
        ! launch level for leg 1) toward its target: on the height axis "up"
        ! means target > previous and requires velocity > 0; on the pressure
        ! axis "up" means target < previous (pressure falls with ascent) and
        ! also requires velocity > 0. Zero velocity or a target equal to the
        ! previous level can never complete and is rejected.
        real(dp) :: prev, toward
        integer :: i

        if (leg_axis == AXIS_PRESSURE) then
            prev = pres
        else
            prev = initial_height
        end if

        do i = 1, n_legs
            if (leg_velocity(i) == 0.0) then
                write(error_unit,*) 'Error: leg ', i, ' has zero velocity ' // &
                    '(the leg could never complete)'
                call exit(1)
            end if
            if (leg_target(i) == prev) then
                write(error_unit,*) 'Error: leg ', i, ' target equals the previous level: ', prev
                call exit(1)
            end if
            ! Upward displacement is positive on the height axis, negative on
            ! the pressure axis.
            toward = leg_target(i) - prev
            if (leg_axis == AXIS_PRESSURE) toward = -toward
            if (toward * leg_velocity(i) < 0.0) then
                write(error_unit,*) 'Error: leg ', i, ' velocity ', leg_velocity(i), &
                    ' points away from its target ', leg_target(i), ' (previous level ', prev, ')'
                call exit(1)
            end if
            prev = leg_target(i)
        end do

    end subroutine validate_legs


    subroutine load_env_profile(lncid)
        ! Loads the full sounding: env_height, env_pressure, env_temperature,
        ! env_RH (all required together).
        integer, intent(in) :: lncid
        integer :: varid, dimid, i

        call nc_verify(nf90_inq_dimid(lncid, 'level', dimid), 'finding level dim')
        call nc_verify(nf90_inquire_dimension(lncid, dimid, len=n_env_levels), 'reading level dim')

        ! Allow re-initialization (unit tests exercise multiple configurations)
        if (allocated(env_pressure)) deallocate(env_pressure, env_temperature, env_RH)
        if (allocated(env_height)) deallocate(env_height)

        allocate(env_pressure(n_env_levels))
        allocate(env_temperature(n_env_levels))
        allocate(env_RH(n_env_levels))
        allocate(env_height(n_env_levels))

        call nc_verify(nf90_inq_varid(lncid, 'env_pressure', varid), 'finding env_pressure')
        call nc_verify(nf90_get_var(lncid, varid, env_pressure), 'reading env_pressure')

        call nc_verify(nf90_inq_varid(lncid, 'env_temperature', varid), 'finding env_temperature')
        call nc_verify(nf90_get_var(lncid, varid, env_temperature), 'reading env_temperature')

        call nc_verify(nf90_inq_varid(lncid, 'env_RH', varid), 'finding env_RH')
        call nc_verify(nf90_get_var(lncid, varid, env_RH), 'reading env_RH')

        call nc_verify(nf90_inq_varid(lncid, 'env_height', varid), 'finding env_height')
        call nc_verify(nf90_get_var(lncid, varid, env_height), 'reading env_height')

        do i = 2, n_env_levels
            if (env_pressure(i) >= env_pressure(i-1)) then
                write(error_unit,*) 'Error: env_pressure must be monotonically decreasing'
                call exit(1)
            end if
            if (env_height(i) <= env_height(i-1)) then
                write(error_unit,*) 'Error: env_height must be monotonically increasing'
                call exit(1)
            end if
        end do

    end subroutine load_env_profile


    subroutine load_entrainment_schedule(lncid)
        ! Reads the optional per-leg entrainment parameters and hands them to
        ! the entrainment module. When absent, the constant &ENTRAINMENT
        ! namelist values apply for the whole run.
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

        allocate(sched_ent_rate(n_legs), sched_n_blob(n_legs), &
                 sched_psigma(n_legs))

        call nc_verify(nf90_get_var(lncid, varid, sched_ent_rate), 'reading ent_rate')
        sched_ent_rate = sched_ent_rate / m_per_km   ! 1/km -> 1/m

        call nc_verify(nf90_inq_varid(lncid, 'n_blob', varid), 'finding n_blob')
        call nc_verify(nf90_get_var(lncid, varid, sched_n_blob), 'reading n_blob')

        call nc_verify(nf90_inq_varid(lncid, 'psigma', varid), 'finding psigma')
        call nc_verify(nf90_get_var(lncid, varid, sched_psigma), 'reading psigma')

        call initialize_entrainment(parcel_velocity, sched_ent_rate, sched_n_blob, &
                                    sched_psigma)

    end subroutine load_entrainment_schedule


    pure function interp_profile(coords, values, query) result(v)
        ! Clamped piecewise-linear lookup on a strictly monotonic coordinate
        ! array (ascending, e.g. env_height, or descending, e.g. env_pressure).
        ! Used both ways on the sounding: p_env(z) = interp_profile(env_height,
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
        ! Advances the parcel one step along the active trajectory leg: move at
        ! the leg's signed velocity, update pressure, change temperature
        ! adiabatically, then advance the leg if its target was reached
        ! (completing the last leg sets trajectory_complete, which ends the
        ! run). Pressure evolves per pressure_mode:
        !   'hydrostatic'  — self-integration dp = -rho_parcel*g*w*dt with
        !                    cooling dT = -(g/cp_moist)*w*dt;
        !   'environment'  — pres = p_env(parcel_height) from the sounding, with
        !                    dT = (Rd*Tv/(cp_m*p))*dp from the actual dp
        !                    (identical to the hydrostatic form when the
        !                    sounding is hydrostatic in the parcel's Tv, correct
        !                    otherwise; reversible under descent since p is a
        !                    function of z).
        ! WV is unchanged here (condensation is handled by droplet growth); only
        ! T, Tv, SS and pres update. Sets pressure_limit_reached and returns
        ! early if the parcel reaches the target pressure.
        real(dp), intent(in) :: ldt   ! time step (s)
        real(dp) :: rho_air, qv_mean, cp_m, dT_adi   ! air density; mean qv; moist cp; adiabatic dT
        real(dp) :: pres_old
        integer :: k

        if (trajectory_complete) return

        parcel_velocity = leg_velocity(current_leg)

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

        ! Diagnostic: where this pressure sits in the environment (hydrostatic
        ! mode with a sounding). Drift from parcel_height measures the departure
        ! of the self-integrated pressure from the sounding's p(z).
        if (write_height_env) then
            parcel_height_env = interp_profile(env_pressure, env_height, pres)
        end if

        call advance_leg()

    end subroutine apply_adiabatic_forcing


    subroutine advance_leg()
        ! Advance past every leg whose target the parcel has reached (a single
        ! step can overshoot more than one target). A reversing leg is never
        ! "already reached", so the loop terminates. Completing the last leg
        ! sets trajectory_complete, which exits the main loop.
        logical :: reached

        do
            if (leg_axis == AXIS_PRESSURE) then
                reached = (leg_velocity(current_leg) > 0.0 .and. pres <= leg_target(current_leg)) &
                     .or. (leg_velocity(current_leg) < 0.0 .and. pres >= leg_target(current_leg))
            else
                reached = (leg_velocity(current_leg) > 0.0 .and. parcel_height >= leg_target(current_leg)) &
                     .or. (leg_velocity(current_leg) < 0.0 .and. parcel_height <= leg_target(current_leg))
            end if
            if (.not. reached) return

            if (current_leg == n_legs) then
                trajectory_complete = .true.
                write(*,'(a,f8.1,a)') ' Trajectory complete at height ', parcel_height, ' m. Stopping.'
                write(error_unit,'(a,f8.1,a)') ' Trajectory complete at height ', parcel_height, ' m. Stopping.'
                return
            end if
            current_leg = current_leg + 1
        end do

    end subroutine advance_leg


    subroutine apply_parcel_entrainment()
        real(dp) :: T_env, qv_env

        call interp_env(pres, T_env, qv_env)
        call apply_entrainment(T_env, qv_env, parcel_velocity, current_leg)

    end subroutine apply_parcel_entrainment


    subroutine interp_env(p_current, T_env, qv_env)
        ! Linearly interpolates the environmental sounding (T, RH vs pressure) to
        ! the parcel's current pressure, clamping at the profile ends, then converts
        ! RH to a vapor mixing ratio. Supplies the "environment" for entrainment.
        real(dp), intent(in) :: p_current    ! current parcel pressure (Pa)
        real(dp), intent(out) :: T_env       ! interpolated environmental temperature (K)
        real(dp), intent(out) :: qv_env      ! environmental vapor mixing ratio (kg/kg)
        real(dp) :: RH_env

        T_env = interp_profile(env_pressure, env_temperature, p_current)
        RH_env = interp_profile(env_pressure, env_RH, p_current)

        qv_env = RH_env * saturation_mixing_ratio(T_env, p_current)
    end subroutine interp_env

end module parcel
