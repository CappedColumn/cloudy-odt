! Unit test for v3 waypoint-leg parcel trajectories.
!
! Builds v3 parcel input NetCDFs (with and without a sounding) and matching
! namelists, then drives initialize_parcel + apply_adiabatic_forcing:
!
!   1. environment mode, up-down-up trajectory (0 -> 400 m at +2, -> 100 m at
!      -1, -> 300 m at +1): heights are revisited with different velocities,
!      legs advance at their targets, pres lies exactly on the sounding's p(z)
!      every step, and trajectory_complete is set when the last target is
!      reached (and the parcel stops moving afterwards).
!   2. hydrostatic mode with a sounding: the parcel_height_env diagnostic
!      drifts away from the integrated parcel_height (deliberately
!      non-matching sounding), and the trajectory is not yet complete.
!   3. hydrostatic mode without a sounding: a bare adiabatic parcel runs, with
!      the diagnostic disabled.
!
! Leg-direction validation errors call exit(1), so negative cases are not
! exercised here (they are covered by the reader's validation messages).
program test_parcel_v3_ascent
    use globals, only: dp, i4, N, T, WV, Tv, SS, pres, time, Tref, &
                       namelist_path, namelist_dir
    use netcdf
    use parcel, only: initialize_parcel, apply_adiabatic_forcing, &
                      parcel_height, parcel_velocity, parcel_height_env, &
                      write_height_env, trajectory_complete
    implicit none

    integer, parameter :: n_leg = 3, n_lev = 4
    real(dp), parameter :: leg_target(n_leg) = [400.0, 100.0, 300.0]
    real(dp), parameter :: leg_w(n_leg) = [2.0, -1.0, 1.0]
    real(dp), parameter :: env_z(n_lev) = [0.0, 1000.0, 2000.0, 3000.0]
    ! Deliberately steep (non-hydrostatic for the parcel's Tv) pressure profile,
    ! so the hydrostatic-mode diagnostic drift is unmistakable.
    real(dp), parameter :: env_p(n_lev) = [95000.0, 80000.0, 67000.0, 56000.0]
    real(dp), parameter :: env_T(n_lev) = [285.0, 278.0, 271.0, 264.0]
    real(dp), parameter :: env_rh(n_lev) = [0.5, 0.5, 0.5, 0.5]
    real(dp), parameter :: dt = 1.0

    integer :: n_passed, n_failed
    real(dp) :: z_expect, p_expect, w200_leg1, w200_leg2, w200_leg3
    integer :: step, leg_expect
    logical :: ok

    n_passed = 0
    n_failed = 0

    ! --- Minimal global state (normally done by base initialization) ---
    N = 8
    allocate(T(N), WV(N), Tv(N), SS(N))
    Tref = 290.0        ! K (test sets the post-conversion value directly)

    call write_parcel_file('test_parcel_v3.nc', with_sounding=.true.)

    ! =======================================================================
    ! Subtest 1: environment mode, up-down-up waypoint trajectory
    ! =======================================================================
    call write_namelist('test_parcel_v3.nml', 'test_parcel_v3.nc', 'environment')
    time = 0.0
    call initialize_parcel()

    call check_true("env: initial pressure on sounding at launch", &
                    abs(pres - env_p(1)) < 1.0e-9)
    call check_true("env: no height_env diagnostic", .not. write_height_env)

    ! Mirror the leg state machine step by step. Timeline (dt = 1 s):
    !   leg 1 (+2): reaches 400 m at t = 200
    !   leg 2 (-1): reaches 100 m at t = 500
    !   leg 3 (+1): reaches 300 m at t = 700 -> trajectory_complete
    z_expect = 0.0
    leg_expect = 1
    w200_leg1 = 0.0; w200_leg2 = 0.0; w200_leg3 = 0.0
    ok = .true.
    do step = 1, 700
        z_expect = z_expect + leg_w(leg_expect) * dt
        p_expect = interp_linear(env_z, env_p, z_expect)

        call apply_adiabatic_forcing(dt)
        time = time + dt

        if (abs(parcel_height - z_expect) > 1.0e-9 .or. &
            abs(pres - p_expect) > 1.0e-6) then
            ok = .false.
            write(*,'(a,i0,4f14.5)') 'mismatch at step ', step, parcel_height, &
                z_expect, pres, p_expect
            exit
        end if

        ! Record the velocity used when passing (about) 200 m on each leg.
        if (abs(z_expect - 200.0) < 0.5) then
            select case (leg_expect)
            case (1); w200_leg1 = parcel_velocity
            case (2); w200_leg2 = parcel_velocity
            case (3); w200_leg3 = parcel_velocity
            end select
        end if

        ! Advance the mirror when the leg target is reached.
        if ((leg_w(leg_expect) > 0.0 .and. z_expect >= leg_target(leg_expect)) .or. &
            (leg_w(leg_expect) < 0.0 .and. z_expect <= leg_target(leg_expect))) then
            leg_expect = leg_expect + 1
            if (leg_expect > n_leg) exit
        end if
    end do
    call check_true("env: height and pressure track the waypoint path", ok)
    call check_true("env: trajectory completes at the last target", &
                    trajectory_complete .and. abs(parcel_height - 300.0) < 1.0e-9)
    call check_true("env: 200 m visited at +2, -1, +1 on the three legs", &
                    w200_leg1 == 2.0 .and. w200_leg2 == -1.0 .and. w200_leg3 == 1.0)

    ! After completion the parcel must not move further.
    z_expect = parcel_height
    call apply_adiabatic_forcing(dt)
    call check_true("env: parcel frozen after completion", parcel_height == z_expect)

    ! =======================================================================
    ! Subtest 2: hydrostatic mode with sounding — diagnostic drift
    ! =======================================================================
    call write_namelist('test_parcel_v3.nml', 'test_parcel_v3.nc', 'hydrostatic')
    pres = env_p(1)
    time = 0.0
    T(:) = Tref
    call initialize_parcel()

    call check_true("hydro: height_env diagnostic enabled", write_height_env)

    do step = 1, 150   ! 150 s at +2 -> 300 m, still inside leg 1
        call apply_adiabatic_forcing(dt)
        time = time + dt
    end do
    call check_true("hydro: parcel ascended on leg 1", parcel_height == 300.0)
    call check_true("hydro: height_env drifts from integrated height", &
                    abs(parcel_height_env - parcel_height) > 1.0)
    call check_true("hydro: trajectory not yet complete", .not. trajectory_complete)

    ! =======================================================================
    ! Subtest 3: hydrostatic mode without a sounding
    ! =======================================================================
    call write_parcel_file('test_parcel_v3_bare.nc', with_sounding=.false.)
    call write_namelist('test_parcel_v3.nml', 'test_parcel_v3_bare.nc', 'hydrostatic')
    pres = 95000.0
    time = 0.0
    T(:) = Tref
    call initialize_parcel()

    call check_true("bare: no height_env diagnostic without a sounding", &
                    .not. write_height_env)
    do step = 1, 50
        call apply_adiabatic_forcing(dt)
        time = time + dt
    end do
    call check_true("bare: parcel ascends hydrostatically", &
                    parcel_height == 100.0 .and. pres < 95000.0)

    ! --- Cleanup ---
    call delete_file('test_parcel_v3.nml')
    call delete_file('test_parcel_v3.nc')
    call delete_file('test_parcel_v3_bare.nc')

    write(*,'(a,i0,a,i0,a)') 'test_parcel_v3_ascent: ', n_passed, ' passed, ', &
                             n_failed, ' failed'
    if (n_failed > 0) call exit(1)

contains

    ! Clamped linear interpolation on an ascending coordinate (test-side mirror
    ! of the model's lookup).
    pure function interp_linear(xs, ys, x) result(y)
        real(dp), intent(in) :: xs(:), ys(:), x
        real(dp) :: y, frac
        integer :: k, nl

        nl = size(xs)
        if (x <= xs(1)) then
            y = ys(1)
        else if (x >= xs(nl)) then
            y = ys(nl)
        else
            do k = 2, nl
                if (x < xs(k)) exit
            end do
            frac = (x - xs(k-1)) / (xs(k) - xs(k-1))
            y = ys(k-1) + frac * (ys(k) - ys(k-1))
        end if
    end function interp_linear


    subroutine write_parcel_file(fname, with_sounding)
        character(*), intent(in) :: fname
        logical, intent(in) :: with_sounding
        integer :: ncid, dim_seg, dim_lev, vid

        call nc_check(nf90_create(fname, NF90_CLOBBER, ncid))
        call nc_check(nf90_def_dim(ncid, 'segment', n_leg, dim_seg))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'conventions', &
                                   'CODT_parcel_input_v3'))
        call nc_check(nf90_def_var(ncid, 'segment_coord', NF90_DOUBLE, dim_seg, vid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, vid, leg_target))
        call def_put(ncid, 'velocity', dim_seg, leg_w)

        if (with_sounding) then
            call nc_check(nf90_redef(ncid))
            call nc_check(nf90_def_dim(ncid, 'level', n_lev, dim_lev))
            call nc_check(nf90_enddef(ncid))
            call def_put(ncid, 'env_height', dim_lev, env_z)
            call def_put(ncid, 'env_pressure', dim_lev, env_p)
            call def_put(ncid, 'env_temperature', dim_lev, env_T)
            call def_put(ncid, 'env_RH', dim_lev, env_rh)
        end if

        call nc_check(nf90_close(ncid))
    end subroutine write_parcel_file


    subroutine def_put(ncid, vname, dimid, vals)
        integer, intent(in) :: ncid, dimid
        character(*), intent(in) :: vname
        real(dp), intent(in) :: vals(:)
        integer :: vid

        call nc_check(nf90_redef(ncid))
        call nc_check(nf90_def_var(ncid, vname, NF90_DOUBLE, dimid, vid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, vid, vals))
    end subroutine def_put


    subroutine nc_check(status)
        integer, intent(in) :: status

        if (status /= NF90_NOERR) then
            write(*,'(a)') 'NetCDF error: '//trim(nf90_strerror(status))
            call exit(1)
        end if
    end subroutine nc_check


    subroutine write_namelist(fname, parcel_nc, mode)
        character(*), intent(in) :: fname, parcel_nc, mode
        integer :: u

        namelist_path = fname
        namelist_dir = './'
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(a)') '&PARCEL'
        write(u, '(a)') "  parcel_file = '"//parcel_nc//"'"
        write(u, '(a)') "  vertical_axis = 'height'"
        write(u, '(a)') "  pressure_mode = '"//mode//"'"
        write(u, '(a)') '  initial_height = 0.0'
        write(u, '(a)') '/'
        close(u)
    end subroutine write_namelist


    subroutine delete_file(fname)
        character(*), intent(in) :: fname
        integer :: u

        open(newunit=u, file=fname, status='old')
        close(u, status='delete')
    end subroutine delete_file


    subroutine check_true(label, cond)
        character(*), intent(in) :: label
        logical, intent(in) :: cond

        if (cond) then
            n_passed = n_passed + 1
            write(*,'(a)') '  PASS  '//label
        else
            n_failed = n_failed + 1
            write(*,'(a)') '  FAIL  '//label
        end if
    end subroutine check_true

end program test_parcel_v3_ascent
