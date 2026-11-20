program test_dgm_composition
    ! Unit test: droplet growth uses the droplet's own solute composition.
    !
    ! Before the fix, growth used NaCl constants for every droplet whatever its
    ! aerosol type. Now set_aerosol_properties takes n_ions, molar mass and
    ! solute density (from the aerosol file), and derives the solution-density
    ! factor c7 = 1 - rho_w/rho_s (solute and water volumes additive).
    !
    ! Checks:
    !   - equilibrium: for each composition and radius, the Kohler equilibrium
    !     S_eq = Kelvin - Raoult computed here from the physical picture (water
    !     mass rho_w*(V - m_s/rho_s), rho_sol = (m_w + m_s)/V; not DGM's c7
    !     algebra) gives dr/dt ~ 0 in growth_rhs; S_eq + 0.1% does not
    !   - distinct materials: (NH4)2SO4 and NaCl have different S_eq, i.e. the
    !     composition actually reaches the growth equation
    !   - the old table is reproduced: with densities 2165 (NaCl) and 1774
    !     ((NH4)2SO4), S_eq built from the former hard-coded c7 values is still
    !     an equilibrium of growth_rhs to within their single-precision rounding
    use globals, only: dp, i4, pi_43, pres, Rv, Mw, rho_l
    use DGM, only: set_aerosol_properties, growth_rhs
    implicit none

    real(dp), parameter :: um = 1.0e-6
    real(dp), parameter :: t0 = 0.0
    real(dp), parameter :: temp_K = 280.0
    real(dp), parameter :: sigma_w = 7.392730e-2     ! as in DGM [N/m]
    real(dp), parameter :: r_dry = 0.05 * um          ! dry radius [m]
    real(dp), parameter :: radii(4) = [0.08, 0.2, 1.0, 5.0] * um
    real(dp), parameter :: offset = 1.0e-3            ! S offset for the reference rate
    real(dp), parameter :: eq_tol = 1.0e-8            ! |dr/dt(S_eq)| / |dr/dt(S_eq + offset)|
    ! Former hard-coded c7 (single-precision literals) and its tolerance. The
    ! constants were given to 7 digits (0.5381062 vs 1 - 1000/2165 =
    ! 0.53810624, ~6e-8 relative), which at 0.08 um moves S_eq by ~1e-8, i.e.
    ! ~1e-5 of offset; a mismatch in the formula itself would be O(1).
    real(dp), parameter :: c7_nacl_old = 0.5381062, c7_as_old = 0.4363021
    real(dp), parameter :: old_tol = 1.0e-4

    ! Compositions: n_ions, molar mass [kg/mol], density [kg/m3]
    integer(i4), parameter :: n_ions(3) = [2, 3, 4]
    real(dp), parameter :: molar_mass(3) = [58.4428e-3, 132.1395e-3, 0.2]
    real(dp), parameter :: density(3) = [2163.0, 1769.0, 1300.0]
    character(9), parameter :: names(3) = ['NaCl     ', '(NH4)2SO4', 'made-up  ']

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    call test_equilibrium()
    call test_distinct_materials()
    call test_old_table()

    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' checks failed'
        error stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' checks passed'

contains

    ! Kohler equilibrium S (fraction) at radius lr, written from the physical
    ! picture rather than DGM's algebra: solute and water volumes add, so the
    ! water occupies V - m_s/rho_s, its mass is rho_w*(V - m_s/rho_s), and the
    ! solution density is (m_w + m_s)/V. Raoult = n_ions*(Mw/M_s)*m_s/m_w.
    ! A wrong c7 = 1 - rho_w/rho_s derivation in DGM would not match this.
    pure function s_equilibrium_additive(lr, ln_ions, lmolar_mass, lrho_s) result(ls_eq)
        real(dp), intent(in) :: lr, lmolar_mass, lrho_s
        integer(i4), intent(in) :: ln_ions
        real(dp) :: ls_eq
        real(dp) :: lm_s, lvol, lm_w, lrho_sol

        lm_s = pi_43 * r_dry**3 * lrho_s
        lvol = pi_43 * lr**3
        lm_w = rho_l * (lvol - lm_s / lrho_s)
        lrho_sol = (lm_w + lm_s) / lvol
        ls_eq = 2.0 * sigma_w / (Rv * temp_K * lrho_sol * lr) &
              - real(ln_ions, dp) * (Mw / lmolar_mass) * lm_s / lm_w
    end function s_equilibrium_additive


    ! The same equilibrium written with a given solution-density factor lc7
    ! (rho_sol = rho_w + lc7*m_s/V), used only to reproduce the former
    ! hard-coded c7 constants.
    pure function s_equilibrium(lr, ln_ions, lmolar_mass, lrho_s, lc7) result(ls_eq)
        real(dp), intent(in) :: lr, lmolar_mass, lrho_s, lc7
        integer(i4), intent(in) :: ln_ions
        real(dp) :: ls_eq
        real(dp) :: lm_s, lvol, lrho_sol

        lm_s = pi_43 * r_dry**3 * lrho_s
        lvol = pi_43 * lr**3
        lrho_sol = rho_l + lc7 * lm_s / lvol
        ls_eq = 2.0 * sigma_w / (Rv * temp_K * lrho_sol * lr) &
              - real(ln_ions, dp) * (Mw / lmolar_mass) * lm_s / (lvol * lrho_sol - lm_s)
    end function s_equilibrium


    ! dr/dt from growth_rhs at radius lr and supersaturation lss.
    function drdt(lr, ln_ions, lmolar_mass, lrho_s, lss) result(lrate)
        real(dp), intent(in) :: lr, lmolar_mass, lrho_s, lss
        integer(i4), intent(in) :: ln_ions
        real(dp) :: lrate
        real(dp) :: ly(3), ldydt(3)
        integer(i4) :: lierr

        call set_aerosol_properties(ln_ions, lmolar_mass, lrho_s, pi_43 * r_dry**3 * lrho_s, &
                                    r_dry, 1.0d-6, lss)
        ly = [lr, 6.0d-3, temp_K]
        call growth_rhs(t0, ly, ldydt, lierr)
        lrate = ldydt(1)
    end function drdt


    subroutine test_equilibrium()
        real(dp) :: ls_eq, lratio, lworst
        integer :: k, i
        character(80) :: label

        do k = 1, size(n_ions)
            lworst = 0.0
            do i = 1, size(radii)
                ls_eq = s_equilibrium_additive(radii(i), n_ions(k), molar_mass(k), density(k))
                lratio = abs(drdt(radii(i), n_ions(k), molar_mass(k), density(k), ls_eq)) &
                       / abs(drdt(radii(i), n_ions(k), molar_mass(k), density(k), ls_eq + offset))
                lworst = max(lworst, lratio)
            end do
            write(label, '(a,a,a)') 'equilibrium: dr/dt = 0 at the Kohler S_eq, ', trim(names(k)), &
                ', 4 radii'
            call check(lworst < eq_tol, trim(label), lworst)
        end do
    end subroutine test_equilibrium


    subroutine test_distinct_materials()
        real(dp) :: ls_nacl, ls_as, lrate_as_at_nacl, lrate_ref

        ls_nacl = s_equilibrium_additive(0.2 * um, n_ions(1), molar_mass(1), density(1))
        ls_as   = s_equilibrium_additive(0.2 * um, n_ions(2), molar_mass(2), density(2))
        print '(a,es12.4,a,es12.4)', '   S_eq at 0.2 um: NaCl ', ls_nacl, '   (NH4)2SO4 ', ls_as
        call check(abs(ls_as - ls_nacl) > 0.01 * abs(ls_nacl), &
                   'distinct: (NH4)2SO4 and NaCl equilibria differ by > 1%')

        ! An (NH4)2SO4 droplet held at NaCl's equilibrium S is not in equilibrium
        lrate_as_at_nacl = drdt(0.2 * um, n_ions(2), molar_mass(2), density(2), ls_nacl)
        lrate_ref = drdt(0.2 * um, n_ions(2), molar_mass(2), density(2), ls_as + offset)
        call check(abs(lrate_as_at_nacl) > 1.0e-3 * abs(lrate_ref), &
                   'distinct: (NH4)2SO4 droplet at the NaCl S_eq still grows or evaporates')
    end subroutine test_distinct_materials


    subroutine test_old_table()
        real(dp) :: ls_eq, lratio, lworst
        integer :: i

        ! NaCl, density 2165 (what the former c7 implied)
        lworst = 0.0
        do i = 1, size(radii)
            ls_eq = s_equilibrium(radii(i), 2_i4, 58.4428d-3, 2165.0d0, c7_nacl_old)
            lratio = abs(drdt(radii(i), 2_i4, 58.4428d-3, 2165.0d0, ls_eq)) &
                   / abs(drdt(radii(i), 2_i4, 58.4428d-3, 2165.0d0, ls_eq + offset))
            lworst = max(lworst, lratio)
        end do
        call check(lworst < old_tol, 'old table: density 2165 reproduces the former NaCl c7', lworst)

        ! (NH4)2SO4, density 1774
        lworst = 0.0
        do i = 1, size(radii)
            ls_eq = s_equilibrium(radii(i), 3_i4, 132.1395d-3, 1774.0d0, c7_as_old)
            lratio = abs(drdt(radii(i), 3_i4, 132.1395d-3, 1774.0d0, ls_eq)) &
                   / abs(drdt(radii(i), 3_i4, 132.1395d-3, 1774.0d0, ls_eq + offset))
            lworst = max(lworst, lratio)
        end do
        call check(lworst < old_tol, 'old table: density 1774 reproduces the former (NH4)2SO4 c7', lworst)
    end subroutine test_old_table


    subroutine check(lcondition, label, lvalue)
        logical, intent(in) :: lcondition
        character(*), intent(in) :: label
        real(dp), intent(in), optional :: lvalue

        if (lcondition) then
            n_passed = n_passed + 1
            if (present(lvalue)) then
                print '(a,a,a,es10.3,a)', '  PASS  ', label, '  (', lvalue, ')'
            else
                print '(a,a)', '  PASS  ', label
            end if
        else
            n_failed = n_failed + 1
            if (present(lvalue)) then
                print '(a,a,a,es10.3,a)', '  FAIL  ', label, '  (', lvalue, ')'
            else
                print '(a,a)', '  FAIL  ', label
            end if
        end if
    end subroutine check

end program test_dgm_composition
