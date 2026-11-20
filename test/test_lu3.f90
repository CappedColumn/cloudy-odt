program test_lu3
    ! Unit test for lu_factor3 / lu_solve3 in ode_integrators.
    !
    ! ros3_integrate used to call solve3 three times per step with the same
    ! matrix. The factorization is now done once (lu_factor3) and reused for
    ! the three stage solves (lu_solve3). This test checks that factor + solve
    ! gives bit-for-bit the result of the former solve3, kept below as a
    ! verbatim copy (solve3_reference), for:
    !   - random 3x3 matrices with entries over 24 orders of magnitude and both
    !     signs, counting how often each of the six pivot paths is taken
    !   - the real ROS3 iteration matrices W = I/(gamma h) - J, with J from
    !     growth_rhs_jac, over droplet states and step sizes
    ! and that reusing one factorization for several right-hand sides gives
    ! the same bits as factoring each time.
    use iso_fortran_env, only: int64
    use globals, only: dp, i4, pi_43, pres, domain_width, H, N, Rd
    use ode_integrators, only: lu_factor3, lu_solve3
    use microphysics, only: saturation_vapor_pressure
    use DGM, only: set_aerosol_properties, growth_rhs_jac
    implicit none

    integer, parameter :: n_random = 200000
    real(dp), parameter :: gamma_ros3 = 0.43586652150845899941601945119356  ! ROS3 gamma

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    call test_random_matrices()
    call test_ros3_matrices()

    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' checks failed'
        error stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' checks passed'

contains

    ! Former ode_integrators::solve3, copied verbatim (only renamed).
    subroutine solve3_reference(A_in, b)
      real(dp), intent(in)    :: A_in(3,3)
      real(dp), intent(inout) :: b(3)

      real(dp) :: A(3,3), tmp
      integer  :: piv, i, j

      A = A_in

      ! Forward elimination — column 1
      piv = 1
      if (abs(A(2,1)) > abs(A(piv,1))) piv = 2
      if (abs(A(3,1)) > abs(A(piv,1))) piv = 3
      if (piv /= 1) then
        do j = 1, 3; tmp = A(1,j); A(1,j) = A(piv,j); A(piv,j) = tmp; end do
        tmp = b(1); b(1) = b(piv); b(piv) = tmp
      end if
      do i = 2, 3
        A(i,1) = A(i,1) / A(1,1)
        do j = 2, 3
          A(i,j) = A(i,j) - A(i,1) * A(1,j)
        end do
        b(i) = b(i) - A(i,1) * b(1)
      end do

      ! Forward elimination — column 2
      if (abs(A(3,2)) > abs(A(2,2))) then
        do j = 2, 3; tmp = A(2,j); A(2,j) = A(3,j); A(3,j) = tmp; end do
        tmp = b(2); b(2) = b(3); b(3) = tmp
      end if
      A(3,2) = A(3,2) / A(2,2)
      A(3,3) = A(3,3) - A(3,2) * A(2,3)
      b(3)   = b(3) - A(3,2) * b(2)

      ! Back-substitution
      b(3) = b(3) / A(3,3)
      b(2) = (b(2) - A(2,3) * b(3)) / A(2,2)
      b(1) = (b(1) - A(1,2) * b(2) - A(1,3) * b(3)) / A(1,1)

    end subroutine solve3_reference


    ! True when every element has the same bit pattern.
    pure logical function same_bits(la, lb)
        real(dp), intent(in) :: la(:), lb(:)
        same_bits = all(transfer(la, [0_int64]) == transfer(lb, [0_int64]))
    end function same_bits


    ! Solve A x = b three ways for three right-hand sides: reference solve3
    ! each time, factor once + three solves, and factor each time. Returns
    ! .true. when all three agree bit for bit. Also returns the pivot path.
    logical function agrees(lA, lrhs, lpiv1, lswap2)
        real(dp), intent(in) :: lA(3,3), lrhs(3,3)   ! columns are right-hand sides
        integer, intent(out) :: lpiv1
        logical, intent(out) :: lswap2

        real(dp) :: lLU(3,3), lLU2(3,3), lx_ref(3), lx_once(3), lx_each(3)
        integer :: k, lp2
        logical :: ls2

        call lu_factor3(lA, lLU, lpiv1, lswap2)
        agrees = .true.
        do k = 1, 3
            lx_ref = lrhs(:, k)
            call solve3_reference(lA, lx_ref)

            lx_once = lrhs(:, k)
            call lu_solve3(lLU, lpiv1, lswap2, lx_once)

            call lu_factor3(lA, lLU2, lp2, ls2)
            lx_each = lrhs(:, k)
            call lu_solve3(lLU2, lp2, ls2, lx_each)

            agrees = agrees .and. same_bits(lx_ref, lx_once) .and. same_bits(lx_ref, lx_each)
        end do
    end function agrees


    subroutine test_random_matrices()
        real(dp) :: la(3,3), lrhs(3,3), lu_exp(3,3), lu_sign(3,3), lu_rhs(3,3)
        integer :: k, lpiv1, n_bad, path_count(3, 0:1)
        logical :: lswap2
        character(80) :: label

        n_bad = 0
        path_count = 0
        do k = 1, n_random
            call random_number(la)
            call random_number(lu_exp)
            call random_number(lu_sign)
            ! entries 1e-12 .. 1e12, random sign
            la = (0.5 + la) * 10.0**(24.0 * lu_exp - 12.0) * merge(1.0, -1.0, lu_sign > 0.5)
            call random_number(lu_rhs)
            call random_number(lu_exp)
            lrhs = (lu_rhs - 0.5) * 10.0**(12.0 * lu_exp - 6.0)

            if (.not. agrees(la, lrhs, lpiv1, lswap2)) n_bad = n_bad + 1
            path_count(lpiv1, merge(1, 0, lswap2)) = path_count(lpiv1, merge(1, 0, lswap2)) + 1
        end do

        print '(a)', ' random matrices, pivot path counts (piv1 x swap2):'
        do k = 1, 3
            print '(a,i1,a,i8,a,i8)', '   piv1 = ', k, ':  no swap2 ', path_count(k, 0), &
                '   swap2 ', path_count(k, 1)
        end do
        write(label, '(a,i0,a)') 'random: factor+solve == former solve3, ', n_random, ' matrices x 3 rhs'
        call check(n_bad == 0, trim(label))
        call check(all(path_count > 0), 'random: all six pivot paths exercised')
    end subroutine test_random_matrices


    subroutine test_ros3_matrices()
        real(dp), parameter :: temp_K = 280.0
        real(dp), parameter :: radii(6) = [0.06, 0.3, 1.0, 5.0, 30.0, 200.0] * 1.0e-6
        real(dp), parameter :: ss_values(3) = [0.005, -0.05, -0.5]
        real(dp), parameter :: steps(5) = [1.0e-6, 1.0e-4, 1.0e-2, 0.1, 1.0]
        real(dp), parameter :: r_solute = 0.05e-6, rho_solute = 1769.0
        real(dp), parameter :: molar_mass_as = 132.1395e-3   ! ammonium sulfate [kg/mol]
        real(dp) :: ly(3), ldydt(3), ljac(3,3), lW(3,3), lrhs(3,3), lu_rhs(3,3)
        real(dp) :: lmass, les, lqv_sat, inv_grid_mass
        integer(i4) :: lierr
        integer :: i, j, k, m, lpiv1, n_bad, n_cases
        logical :: lswap2

        inv_grid_mass = 1.0 / (50.0 * domain_width**2 * (H / N) * pres / (Rd * temp_K))
        lmass = pi_43 * r_solute**3 * rho_solute
        les = saturation_vapor_pressure(temp_K)
        lqv_sat = 0.622 * les / (pres - les)

        n_bad = 0
        n_cases = 0
        do i = 1, size(radii)
            do j = 1, size(ss_values)
                call set_aerosol_properties(3_i4, molar_mass_as, rho_solute, lmass, r_solute, &
                                            inv_grid_mass, ss_values(j))
                ly = [radii(i), lqv_sat * (1.0 + ss_values(j)), temp_K]
                call growth_rhs_jac(0.0d0, ly, ldydt, ljac, lierr)
                do k = 1, size(steps)
                    lW = -ljac
                    do m = 1, 3
                        lW(m, m) = lW(m, m) + 1.0 / (gamma_ros3 * steps(k))
                    end do
                    call random_number(lu_rhs)
                    lrhs = spread(ldydt, 2, 3) * (0.5 + lu_rhs)
                    if (.not. agrees(lW, lrhs, lpiv1, lswap2)) n_bad = n_bad + 1
                    n_cases = n_cases + 1
                end do
            end do
        end do
        call check(n_bad == 0, 'ROS3 W = I/(gamma h) - J from growth_rhs_jac: factor+solve == former solve3')
        print '(a,i0,a)', '   (', n_cases, ' matrices: 6 radii x 3 S x 5 step sizes, 3 rhs each)'
    end subroutine test_ros3_matrices


    subroutine check(lcondition, label)
        logical, intent(in) :: lcondition
        character(*), intent(in) :: label

        if (lcondition) then
            n_passed = n_passed + 1
            print '(a,a)', '  PASS  ', label
        else
            n_failed = n_failed + 1
            print '(a,a)', '  FAIL  ', label
        end if
    end subroutine check

end program test_lu3
