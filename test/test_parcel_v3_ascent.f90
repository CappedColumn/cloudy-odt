! Unit test for v3 parcel ascent on a height segment coordinate.
!
! Builds a small v3 parcel input NetCDF (height-axis segments including a
! negative-velocity segment, plus a full sounding with env_height) and a
! matching namelist, then drives initialize_parcel + apply_adiabatic_forcing:
!
!   1. pressure_mode = 'environment': the parcel's pressure must lie exactly on
!      the sounding's p(z) at the (kinematic) parcel height every step, and the
!      velocity must switch segments as the parcel crosses the segment levels.
!      The w>0 / w<0 boundary at the top segment acts as a stagnation level:
!      the parcel must stay bounded near it.
!   2. pressure_mode = 'hydrostatic': the self-integrated pressure leaves the
!      sounding's (deliberately non-matching) p(z) curve, so the
!      parcel_height_env diagnostic must drift away from parcel_height.
program test_parcel_v3_ascent
    use globals, only: dp, i4, N, T, WV, Tv, SS, pres, time, Tref, &
                       namelist_path, namelist_dir
    use netcdf
    use parcel, only: initialize_parcel, apply_adiabatic_forcing, &
                      parcel_height, parcel_velocity, parcel_height_env, &
                      write_height_env
    implicit none

    integer, parameter :: n_seg = 3, n_lev = 4
    real(dp), parameter :: seg_z(n_seg) = [0.0, 500.0, 1000.0]
    real(dp), parameter :: seg_w(n_seg) = [1.0, 2.0, -1.0]
    real(dp), parameter :: env_z(n_lev) = [0.0, 1000.0, 2000.0, 3000.0]
    ! Deliberately steep (non-hydrostatic for the parcel's Tv) pressure profile,
    ! so the hydrostatic-mode diagnostic drift is unmistakable.
    real(dp), parameter :: env_p(n_lev) = [95000.0, 80000.0, 67000.0, 56000.0]
    real(dp), parameter :: env_T(n_lev) = [285.0, 278.0, 271.0, 264.0]
    real(dp), parameter :: env_rh(n_lev) = [0.5, 0.5, 0.5, 0.5]
    real(dp), parameter :: dt = 1.0

    integer :: n_passed, n_failed
    real(dp) :: z_expect, w_expect, p_expect
    integer :: step
    logical :: ok

    n_passed = 0
    n_failed = 0

    ! --- Minimal global state (normally done by base initialization) ---
    N = 8
    allocate(T(N), WV(N), Tv(N), SS(N))
    Tref = 290.0        ! K (test sets the post-conversion value directly)

    call write_parcel_file('test_parcel_v3.nc')

    ! =======================================================================
    ! Subtest 1: environment mode — pressure rides the sounding's p(z)
    ! =======================================================================
    call write_namelist('test_parcel_v3.nml', 'environment')
    pres = env_p(1)
    time = 0.0
    parcel_height = 0.0
    call initialize_parcel()

    call check_true("env mode: no height_env diagnostic", .not. write_height_env)

    z_expect = 0.0
    ok = .true.
    do step = 1, 700
        ! Expected kinematics (mirror of the model's forward-Euler stepping)
        w_expect = expected_velocity(z_expect)
        z_expect = z_expect + w_expect * dt
        p_expect = interp_linear(env_z, env_p, z_expect)

        call apply_adiabatic_forcing(dt)
        time = time + dt

        if (abs(parcel_height - z_expect) > 1.0e-9) then
            ok = .false.
            write(*,'(a,i0,2f12.4)') 'height mismatch at step ', step, &
                parcel_height, z_expect
            exit
        end if
        if (abs(pres - p_expect) > 1.0e-6) then
            ok = .false.
            write(*,'(a,i0,2f14.4)') 'pressure off p_env(z) at step ', step, &
                pres, p_expect
            exit
        end if
    end do
    call check_true("env mode: pres == p_env(parcel_height) every step", ok)

    ! After 700 s: ascent reaches z=500 at t=500 (w=1), then w=2 up to z=1000
    ! at t=750... but the w>0/w<0 boundary at 1000 m is reached at
    ! t = 500 + 250 = 750 > 700, so at t=700: z = 500 + 200*2 = 900, segment 2.
    call check_true("env mode: velocity in segment 2 after crossing z=500", &
                    parcel_velocity == seg_w(2))

    ! Keep stepping past the stagnation boundary at z=1000 (w=2 below, w=-1
    ! above): the parcel must stay bounded near it, not run away.
    do step = 1, 300
        call apply_adiabatic_forcing(dt)
        time = time + dt
    end do
    call check_true("env mode: parcel settles at the stagnation level", &
                    parcel_height >= 1000.0 - 2.0*dt .and. &
                    parcel_height <= 1000.0 + 2.0*dt)

    ! =======================================================================
    ! Subtest 2: hydrostatic mode — parcel_height_env drifts from parcel_height
    ! =======================================================================
    call write_namelist('test_parcel_v3.nml', 'hydrostatic')
    pres = env_p(1)
    time = 0.0
    parcel_height = 0.0
    T(:) = Tref
    call initialize_parcel()

    call check_true("hydro mode: height_env diagnostic enabled", write_height_env)

    do step = 1, 400   ! 400 s at w=1 -> 400 m, stays in segment 1
        call apply_adiabatic_forcing(dt)
        time = time + dt
    end do

    call check_true("hydro mode: parcel ascended", parcel_height == 400.0)
    call check_true("hydro mode: height_env drifts from integrated height", &
                    abs(parcel_height_env - parcel_height) > 1.0)

    ! --- Cleanup ---
    call delete_file('test_parcel_v3.nml')
    call delete_file('test_parcel_v3.nc')

    write(*,'(a,i0,a,i0,a)') 'test_parcel_v3_ascent: ', n_passed, ' passed, ', &
                             n_failed, ' failed'
    if (n_failed > 0) call exit(1)

contains

    ! Piecewise-constant expected velocity on the height segments.
    pure function expected_velocity(z) result(w)
        real(dp), intent(in) :: z
        real(dp) :: w
        integer :: i

        w = seg_w(1)
        do i = 2, n_seg
            if (z < seg_z(i)) exit
            w = seg_w(i)
        end do
    end function expected_velocity


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


    subroutine write_parcel_file(fname)
        character(*), intent(in) :: fname
        integer :: ncid, dim_seg, dim_lev, vid

        call nc_check(nf90_create(fname, NF90_CLOBBER, ncid))
        call nc_check(nf90_def_dim(ncid, 'segment', n_seg, dim_seg))
        call nc_check(nf90_def_dim(ncid, 'level', n_lev, dim_lev))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'conventions', &
                                   'CODT_parcel_input_v3'))

        call nc_check(nf90_def_var(ncid, 'segment_coord', NF90_DOUBLE, dim_seg, vid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, vid, seg_z))

        call def_put(ncid, 'velocity', dim_seg, seg_w)
        call def_put(ncid, 'env_height', dim_lev, env_z)
        call def_put(ncid, 'env_pressure', dim_lev, env_p)
        call def_put(ncid, 'env_temperature', dim_lev, env_T)
        call def_put(ncid, 'env_RH', dim_lev, env_rh)

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


    subroutine write_namelist(fname, mode)
        character(*), intent(in) :: fname, mode
        integer :: u

        namelist_path = fname
        namelist_dir = './'
        open(newunit=u, file=fname, status='replace', action='write')
        write(u, '(a)') '&PARCEL'
        write(u, '(a)') "  parcel_file = 'test_parcel_v3.nc'"
        write(u, '(a)') "  vertical_axis = 'height'"
        write(u, '(a)') "  pressure_mode = '"//mode//"'"
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
