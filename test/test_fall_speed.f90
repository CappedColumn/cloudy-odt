program test_fall_speed
    ! Unit test for the droplet fall-speed formulas in particle_types:
    ! beard_fall_speed (Beard 1976) and the fall_speed_scheme switch behind
    ! calculate_terminal_velocity.
    !
    ! Reference speeds marked "transcription" come from an independent Python
    ! copy of the same published formula (same coefficients, same viscosity
    ! fit), so they catch coding errors, not errors in the formula. The
    ! Gunn & Kinzer (1949) comparison is the check against measurements.
    use globals, only: dp, g, rho_l, nu, Rd, pres
    use particle_types, only: particle, calculate_terminal_velocity, beard_fall_speed, &
                              fall_speed_scheme, set_fall_speed_selector
    implicit none

    real(dp), parameter :: um = 1.0e-6_dp
    real(dp), parameter :: mm = 1.0e-3_dp

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    call test_beard_transcription()
    call test_beard_against_measurements()
    call test_beard_shape()
    call test_scheme_switch()

    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' checks failed'
        stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' checks passed'

contains

    ! -----------------------------------------------------------------------
    ! Subtest 1: the Fortran reproduces the reference implementation
    ! -----------------------------------------------------------------------
    subroutine test_beard_transcription()
        ! Fall speeds at 283 K and 1000 hPa, dry air (density = p / (Rd T)).
        ! The source constants are single-precision literals, so agreement is
        ! expected to about 1e-7, not to round-off.
        real(dp), parameter :: ltemp = 283.0_dp, lpres = 1.0e5_dp
        real(dp), parameter :: rtol = 1.0e-5_dp
        real(dp) :: lrho_air

        lrho_air = lpres / (Rd * ltemp)

        ! Slip-corrected Stokes range (diameter <= 19 um)
        call check_close("transcription: r =   5 um", &
            beard_fall_speed(5.0_dp * um, ltemp, lpres, lrho_air), 3.13133770e-3_dp, rtol)
        ! Reynolds-number fit (19 um < diameter <= 1.07 mm)
        call check_close("transcription: r =  10 um", &
            beard_fall_speed(10.0_dp * um, ltemp, lpres, lrho_air), 1.23894675e-2_dp, rtol)
        call check_close("transcription: r =  20 um", &
            beard_fall_speed(20.0_dp * um, ltemp, lpres, lrho_air), 4.84260754e-2_dp, rtol)
        call check_close("transcription: r =  30 um", &
            beard_fall_speed(30.0_dp * um, ltemp, lpres, lrho_air), 1.04548232e-1_dp, rtol)
        call check_close("transcription: r =  40 um", &
            beard_fall_speed(40.0_dp * um, ltemp, lpres, lrho_air), 1.74967636e-1_dp, rtol)
        call check_close("transcription: r = 100 um", &
            beard_fall_speed(100.0_dp * um, ltemp, lpres, lrho_air), 7.01886963e-1_dp, rtol)
        call check_close("transcription: r = 400 um", &
            beard_fall_speed(400.0_dp * um, ltemp, lpres, lrho_air), 3.26371230_dp, rtol)

        ! Below 0 C the viscosity fit gains a quadratic term.
        lrho_air = 7.0e4_dp / (Rd * 263.15_dp)
        call check_close("transcription: r =  20 um at 263.15 K, 700 hPa", &
            beard_fall_speed(20.0_dp * um, 263.15_dp, 7.0e4_dp, lrho_air), 5.15094176e-2_dp, rtol)
    end subroutine test_beard_transcription


    ! -----------------------------------------------------------------------
    ! Subtest 2: agreement with measured fall speeds
    ! -----------------------------------------------------------------------
    subroutine test_beard_against_measurements()
        ! Gunn & Kinzer (1949), J. Meteor., 6, 243-248: water drops in stagnant
        ! air at 20 C and 1013 hPa. Beard's fit should lie within a few percent.
        real(dp), parameter :: ltemp = 293.15_dp, lpres = 101325.0_dp
        real(dp), parameter :: rtol = 0.04_dp
        real(dp) :: lrho_air

        lrho_air = lpres / (Rd * ltemp)

        call check_close("Gunn-Kinzer: d = 0.2 mm (0.72 m/s)", &
            beard_fall_speed(0.1_dp * mm, ltemp, lpres, lrho_air), 0.72_dp, rtol)
        call check_close("Gunn-Kinzer: d = 0.5 mm (2.06 m/s)", &
            beard_fall_speed(0.25_dp * mm, ltemp, lpres, lrho_air), 2.06_dp, rtol)
        call check_close("Gunn-Kinzer: d = 1.0 mm (4.03 m/s)", &
            beard_fall_speed(0.5_dp * mm, ltemp, lpres, lrho_air), 4.03_dp, rtol)
    end subroutine test_beard_against_measurements


    ! -----------------------------------------------------------------------
    ! Subtest 3: shape of the curve
    ! -----------------------------------------------------------------------
    subroutine test_beard_shape()
        ! The step where the two size ranges meet, growth with size away from
        ! it, the 1.07 mm cap, and non-positive radii.
        !
        ! Beard's two fits do not join exactly: at d = 19 um the speed steps
        ! DOWN by about 0.2%, so in a band about 0.1% wide in radius a smaller
        ! drop falls faster than a larger one. That is a property of the
        ! formula (docs/known_issues.md), and the test pins it down so that a
        ! change in it is noticed.
        real(dp), parameter :: ltemp = 283.0_dp, lpres = 1.0e5_dp
        real(dp), parameter :: r_join = 9.5_dp * um    ! diameter 19 um
        real(dp), parameter :: r_cap = 0.535_dp * mm   ! diameter 1.07 mm
        real(dp) :: lrho_air, v_below, v_above, v_prev, v_now, lradius
        logical :: increasing
        integer :: k, n_down_steps

        lrho_air = lpres / (Rd * ltemp)

        v_below = beard_fall_speed(r_join * (1.0_dp - 1.0e-6_dp), ltemp, lpres, lrho_air)
        v_above = beard_fall_speed(r_join * (1.0_dp + 1.0e-6_dp), ltemp, lpres, lrho_air)
        call check_close("shape: ranges agree to 0.5% at d = 19 um", v_above, v_below, 5.0e-3_dp)
        call check_true("shape: the step at d = 19 um is downward", v_above < v_below)

        ! Fine scan across the join, 9.4 to 9.6 um in steps of 1e-5 um: the
        ! speed must fall at exactly one step (the join) and rise at all others.
        n_down_steps = 0
        v_prev = beard_fall_speed(9.4_dp * um, ltemp, lpres, lrho_air)
        do k = 1, 20000
            lradius = (9.4_dp + 1.0e-5_dp * real(k, dp)) * um
            v_now = beard_fall_speed(lradius, ltemp, lpres, lrho_air)
            if (v_now <= v_prev) n_down_steps = n_down_steps + 1
            v_prev = v_now
        end do
        call check_true("shape: one downward step between 9.4 and 9.6 um", n_down_steps == 1)

        ! Coarse scan, 1 um to 500 um radius in 2% steps (too coarse to see the join)
        increasing = .true.
        lradius = 1.0_dp * um
        v_prev = beard_fall_speed(lradius, ltemp, lpres, lrho_air)
        do k = 1, 314
            lradius = lradius * 1.02_dp
            v_now = beard_fall_speed(lradius, ltemp, lpres, lrho_air)
            if (v_now <= v_prev) increasing = .false.
            v_prev = v_now
        end do
        call check_true("shape: speed increases with radius in 2% steps, 1 to 500 um", increasing)

        call check_close("shape: d = 3 mm gets the d = 1.07 mm speed", &
            beard_fall_speed(1.5_dp * mm, ltemp, lpres, lrho_air), &
            beard_fall_speed(r_cap, ltemp, lpres, lrho_air), 1.0e-14_dp)

        call check_true("shape: zero radius gives zero speed", &
            beard_fall_speed(0.0_dp, ltemp, lpres, lrho_air) == 0.0_dp)
        call check_true("shape: negative radius gives zero speed", &
            beard_fall_speed(-1.0_dp * um, ltemp, lpres, lrho_air) == 0.0_dp)
    end subroutine test_beard_shape


    ! -----------------------------------------------------------------------
    ! Subtest 4: the fall_speed_scheme switch
    ! -----------------------------------------------------------------------
    subroutine test_scheme_switch()
        ! calculate_terminal_velocity must give Stokes by default, Beard when
        ! selected, a negative (downward) velocity in both, and must reject an
        ! unknown name without changing the active formula.
        type(particle) :: p
        real(dp) :: lrho_air, v_stokes, v_beard
        logical :: lvalid

        pres = 1.0e5_dp
        p%radius = 30.0_dp * um
        p%temperature = 283.0_dp
        p%virt_temp = 284.0_dp
        lrho_air = pres / (p%virt_temp * Rd)

        ! Default, before any selector call
        v_stokes = -(2.0_dp * g * rho_l / (9.0_dp * nu * lrho_air)) * p%radius**2
        call check_close("switch: default is Stokes", &
            calculate_terminal_velocity(p), v_stokes, 1.0e-14_dp)

        fall_speed_scheme = 'beard'
        call set_fall_speed_selector(lvalid)
        v_beard = -beard_fall_speed(p%radius, p%temperature, pres, lrho_air)
        call check_true("switch: 'beard' accepted", lvalid)
        call check_close("switch: 'beard' selects beard_fall_speed", &
            calculate_terminal_velocity(p), v_beard, 1.0e-14_dp)
        call check_true("switch: Beard velocity is downward", calculate_terminal_velocity(p) < 0.0_dp)
        call check_true("switch: Beard is slower than Stokes at r = 30 um", &
            abs(v_beard) < abs(v_stokes))

        fall_speed_scheme = 'bread'
        call set_fall_speed_selector(lvalid)
        call check_true("switch: unknown name rejected", .not. lvalid)
        call check_close("switch: unknown name leaves Beard active", &
            calculate_terminal_velocity(p), v_beard, 1.0e-14_dp)

        fall_speed_scheme = 'stokes'
        call set_fall_speed_selector(lvalid)
        call check_true("switch: 'stokes' accepted", lvalid)
        call check_close("switch: 'stokes' selects Stokes again", &
            calculate_terminal_velocity(p), v_stokes, 1.0e-14_dp)
    end subroutine test_scheme_switch


    ! -----------------------------------------------------------------------
    ! Helpers
    ! -----------------------------------------------------------------------
    subroutine check_true(label, cond)
        character(*), intent(in) :: label
        logical, intent(in) :: cond

        if (cond) then
            n_passed = n_passed + 1
            print '(a,a)', '  PASS  ', label
        else
            n_failed = n_failed + 1
            print '(a,a)', '  FAIL  ', label
        end if
    end subroutine check_true


    subroutine check_close(label, got, expected, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, rtol
        real(dp) :: relerr

        relerr = abs(got - expected) / max(abs(expected), 1.0e-30_dp)
        if (relerr <= rtol) then
            n_passed = n_passed + 1
            print '(a,a,a,es12.5,a,es10.3,a)', '  PASS  ', label, '  (got ', got, ', rel err ', relerr, ')'
        else
            n_failed = n_failed + 1
            print '(a,a,a,es12.5,a,es12.5,a,es10.3,a)', '  FAIL  ', label, &
                '  (got ', got, ', expected ', expected, ', rel err ', relerr, ')'
        end if
    end subroutine check_close

end program test_fall_speed
