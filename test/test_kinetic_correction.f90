program test_kinetic_correction
    ! Unit test for the do_kinetic_correction switch in DGM.
    !
    ! Checks that:
    !   - the switch defaults to on
    !   - with the correction off, |dr/dt| is larger (the jump lengths only
    !     ever reduce K and Dv), for growth and for evaporation
    !   - the on/off difference shrinks as the radius grows past the jump
    !     lengths, and from 0.06 um (just above the 50.5 nm radius floor of a
    !     50 nm dry particle) to 200 um the off/on ratio equals the expected
    !     D_on / D_off. That check uses the same property formulas as DGM, so
    !     it catches coding errors in the switch, not errors in the formulas.
    !   - in both modes the analytic Jacobians (growth_jacobian and
    !     growth_rhs_jac) match central finite differences of growth_rhs in
    !     the radius and temperature columns
    use globals, only: dp, i4, pi, pi_43, pres, domain_width, H, N, Rd, Rv, cv, Tice
    use microphysics, only: saturation_vapor_pressure
    use DGM, only: set_aerosol_properties, growth_rhs, growth_jacobian, growth_rhs_jac, &
                   do_kinetic_correction
    implicit none

    real(dp), parameter :: um = 1.0e-6
    real(dp), parameter :: t0 = 0.0               ! time argument of the DGM routines (unused)
    real(dp), parameter :: temp_K = 280.0         ! air temperature [K]
    real(dp), parameter :: rho_solute = 1769.0    ! ammonium sulfate [kg/m3]
    real(dp), parameter :: r_solute = 0.05 * um   ! dry radius [m]
    real(dp), parameter :: jac_rtol = 1.0e-4      ! Jacobian vs finite differences
    ! off/on ratio vs expected. Not 1e-15: es here comes from
    ! microphysics::saturation_vapor_pressure, whose Flatau coefficients are
    ! single-precision literals, while DGM's private esat uses d0 literals; the
    ! two differ by ~1e-8.
    real(dp), parameter :: ratio_rtol = 1.0e-7

    integer :: n_passed, n_failed
    real(dp) :: inv_grid_mass

    n_passed = 0
    n_failed = 0

    inv_grid_mass = 1.0 / (50.0 * domain_width**2 * (H / N) * pres / (Rd * temp_K))

    call check(do_kinetic_correction, 'default: do_kinetic_correction is .true.')

    call test_growth_ratio()
    call test_jacobians(.true.)
    call test_jacobians(.false.)
    do_kinetic_correction = .true.

    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' checks failed'
        error stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' checks passed'

contains

    ! dr/dt from growth_rhs with the switch set to lkinetic.
    function drdt_at(lradius, lss, lkinetic) result(ldrdt)
        real(dp), intent(in) :: lradius   ! droplet radius [m]
        real(dp), intent(in) :: lss       ! supersaturation [fraction]
        logical,  intent(in) :: lkinetic
        real(dp) :: ldrdt

        real(dp) :: ly(3), ldydt(3)
        integer(i4) :: lierr

        do_kinetic_correction = lkinetic
        call set_state(lradius, lss, ly)
        call growth_rhs(t0, ly, ldydt, lierr)
        ldrdt = ldydt(1)
    end function drdt_at


    ! Sets the DGM aerosol state and fills the state vector y.
    subroutine set_state(lradius, lss, ly)
        real(dp), intent(in)  :: lradius, lss
        real(dp), intent(out) :: ly(3)

        real(dp) :: lmass, les, lqv_sat

        lmass = pi_43 * r_solute**3 * rho_solute
        call set_aerosol_properties(2_i4, lmass, r_solute, inv_grid_mass, lss)
        les = saturation_vapor_pressure(temp_K)
        lqv_sat = 0.622 * les / (pres - les)
        ly = [lradius, lqv_sat * (1.0 + lss), temp_K]
    end subroutine set_state


    ! Expected off/on dr/dt ratio, D_on / D_off. The Kohler terms and the
    ! solution density cancel, so it depends only on r, T and p:
    !   D_on / D_off = [A (r + l_v)/r + B (r + l_T)/r] / (A + B)
    ! with A = Rv T / (Dv es) (vapor) and B = Lv^2 / (K Rv T^2) (heat).
    ! Same property formulas and coefficients (alpha = 0.04, alpha_T = 1) as DGM.
    function expected_ratio(lradius, ljump_v, ljump_T) result(lratio)
        real(dp), intent(in)  :: lradius            ! droplet radius [m]
        real(dp), intent(out) :: ljump_v, ljump_T   ! vapor, thermal jump lengths [m]
        real(dp) :: lratio

        real(dp) :: lK, lDv, lLv, les, lA, lB

        lK  = 7.7e-5 * (temp_K - Tice) + 0.02399
        lDv = (1.57e-7 * (temp_K - Tice) + 2.211e-5) * 1.0e5 / pres
        lLv = (2.501 - 0.00237 * (temp_K - Tice)) * 1.0e6
        les = saturation_vapor_pressure(temp_K)

        ljump_T = lK * sqrt(2.0 * pi * Rd * temp_K) / (1.0 * pres * (cv + Rd / 2.0))
        ljump_v = sqrt(2.0 * pi / (Rv * temp_K)) * lDv / 0.04

        lA = Rv * temp_K / (lDv * les)
        lB = lLv**2 / (lK * Rv * temp_K**2)
        lratio = (lA * (lradius + ljump_v) / lradius + lB * (lradius + ljump_T) / lradius) &
               / (lA + lB)
    end function expected_ratio


    subroutine test_growth_ratio()
        ! The first four radii are below 1 um
        real(dp), parameter :: radii(8) = [0.06, 0.1, 0.2, 0.5, 1.0, 5.0, 20.0, 200.0] * um
        real(dp), parameter :: ss_values(2) = [0.005, -0.05]   ! growth, evaporation
        real(dp) :: ratio(size(radii)), lexpect, ljump_v, ljump_T, lerr_max
        real(dp) :: lon, loff
        integer :: i, j
        character(80) :: label

        lexpect = expected_ratio(1.0 * um, ljump_v, ljump_T)

        print *
        print '(a)', ' dr/dt with the correction off / on (ammonium sulfate, r_dry = 50 nm, 280 K)'
        print '(a,f6.3,a,f6.3,a)', '   jump lengths: vapor ', ljump_v / um, ' um, thermal ', &
            ljump_T / um, ' um'
        do j = 1, size(ss_values)
            lerr_max = 0.0
            do i = 1, size(radii)
                lon  = drdt_at(radii(i), ss_values(j), .true.)
                loff = drdt_at(radii(i), ss_values(j), .false.)
                ratio(i) = loff / lon
                lexpect = expected_ratio(radii(i), ljump_v, ljump_T)
                lerr_max = max(lerr_max, abs(ratio(i) / lexpect - 1.0))
                print '(a,f6.3,a,f7.2,a,es11.3,a,es11.3,a,f8.4,a,f8.4)', '   S = ', &
                    100.0 * ss_values(j), ' %  r = ', radii(i) / um, ' um  on ', lon, &
                    '  off ', loff, '  ratio ', ratio(i), '  expected ', lexpect
            end do

            write(label, '(a,f6.3,a)') 'S = ', 100.0 * ss_values(j), &
                ' %: off/on equals D_on/D_off, 0.06-200 um'
            call check(lerr_max < ratio_rtol, trim(label), lerr_max)
            write(label, '(a,f6.3,a)') 'S = ', 100.0 * ss_values(j), ' %: off/on > 3 below 1 um'
            call check(all(ratio(:4) > 3.0), trim(label))
            write(label, '(a,f6.3,a)') 'S = ', 100.0 * ss_values(j), ' %: off/on > 1 at every radius'
            call check(all(ratio > 1.0), trim(label))
            write(label, '(a,f6.3,a)') 'S = ', 100.0 * ss_values(j), ' %: off/on falls with radius'
            call check(all(ratio(2:) < ratio(:size(radii) - 1)), trim(label))
            write(label, '(a,f6.3,a)') 'S = ', 100.0 * ss_values(j), ' %: off/on < 1.05 at 200 um'
            call check(ratio(size(radii)) < 1.05, trim(label))
        end do
    end subroutine test_growth_ratio


    ! Compare both analytic Jacobians with central differences of growth_rhs.
    subroutine test_jacobians(lkinetic)
        logical, intent(in) :: lkinetic

        real(dp), parameter :: radii(5) = [0.1, 0.3, 1.0, 5.0, 30.0] * um
        real(dp), parameter :: ss_values(2) = [0.005, -0.05]
        real(dp) :: ly(3), lyp(3), lym(3), lfp(3), lfm(3), ldydt(3)
        real(dp) :: ljac_a(3, 3), ljac_f(3, 3), lfd(3, 3), lstep, lerr_a, lerr_f
        integer(i4) :: lierr
        integer :: i, j, col
        character(100) :: label

        do j = 1, size(ss_values)
            do i = 1, size(radii)
                do_kinetic_correction = lkinetic
                call set_state(radii(i), ss_values(j), ly)

                lfd = 0.0
                do col = 1, 3, 2   ! radius and temperature columns
                    if (col == 1) then
                        lstep = 1.0e-6 * ly(1)
                    else
                        lstep = 1.0e-4
                    end if
                    lyp = ly;  lyp(col) = ly(col) + lstep
                    lym = ly;  lym(col) = ly(col) - lstep
                    call growth_rhs(t0, lyp, lfp, lierr)
                    call growth_rhs(t0, lym, lfm, lierr)
                    lfd(:, col) = (lfp - lfm) / (2.0 * lstep)
                end do

                call growth_jacobian(t0, ly, ljac_a)
                call growth_rhs_jac(t0, ly, ldydt, ljac_f, lierr)

                lerr_a = max_rel_err(ljac_a, lfd)
                lerr_f = max_rel_err(ljac_f, lfd)

                write(label, '(a,l1,a,f5.1,a,f6.3,a)') 'kinetic=', lkinetic, ' r=', radii(i) / um, &
                    ' um S=', 100.0 * ss_values(j), ' %: growth_jacobian matches FD'
                call check(lerr_a < jac_rtol, trim(label), lerr_a)
                write(label, '(a,l1,a,f5.1,a,f6.3,a)') 'kinetic=', lkinetic, ' r=', radii(i) / um, &
                    ' um S=', 100.0 * ss_values(j), ' %: growth_rhs_jac matches FD'
                call check(lerr_f < jac_rtol, trim(label), lerr_f)
            end do
        end do
    end subroutine test_jacobians


    ! Largest relative error over the r and T columns, each entry scaled by
    ! the largest finite-difference entry in its row.
    function max_rel_err(ljac, lfd) result(lerr)
        real(dp), intent(in) :: ljac(3, 3), lfd(3, 3)
        real(dp) :: lerr
        integer :: row, col
        real(dp) :: lscale

        lerr = 0.0
        do row = 1, 3
            lscale = max(abs(lfd(row, 1)), abs(lfd(row, 3)))
            if (lscale == 0.0) cycle
            do col = 1, 3, 2
                lerr = max(lerr, abs(ljac(row, col) - lfd(row, col)) / lscale)
            end do
        end do
    end function max_rel_err


    subroutine check(lcondition, label, lvalue)
        logical, intent(in) :: lcondition
        character(*), intent(in) :: label
        real(dp), intent(in), optional :: lvalue

        if (lcondition) then
            n_passed = n_passed + 1
            if (present(lvalue)) then
                print '(a,a,a,es10.3,a)', '  PASS  ', label, '  (rel err ', lvalue, ')'
            else
                print '(a,a)', '  PASS  ', label
            end if
        else
            n_failed = n_failed + 1
            if (present(lvalue)) then
                print '(a,a,a,es10.3,a)', '  FAIL  ', label, '  (rel err ', lvalue, ')'
            else
                print '(a,a)', '  FAIL  ', label
            end if
        end if
    end subroutine check

end program test_kinetic_correction
