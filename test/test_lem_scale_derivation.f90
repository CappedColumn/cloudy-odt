program test_lem_scale_derivation
    ! Exercises LEM::derive_turbulence_scales over a sweep of (N, H, eps) that
    ! covers both regimes and the transition between them.
    !
    ! The invariant under test is the one the scheme exists to guarantee: the
    ! LEM diffusivities are NEVER below molecular. That held for neither the
    ! pre-v3.0.0 code (which sampled eddies below the grid) nor the first slice-2
    ! attempt (which replaced kT/Dv with 0.1*eps^(1/3)*(6dz)^(4/3), a value below
    ! molecular for typical grids).
    use globals, only: dp, i4, N, H, dz_length, kT, Dv, nu, min_eddy_gridpoints
    use LEM, only: derive_turbulence_scales, integral_length_scale, dissipation_rate, &
                   actual_kolmogorov_scale, grid_eddy_scale, diffusivity_length_scale, &
                   smallest_eddy_scale, smallest_eddy_gridpoints, diffusivity_enhancement, &
                   thermal_diffusivity, vapor_diffusivity, reynolds_number
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    call test_regime_sweep()
    call test_quantization_bounds()
    call test_reference_values()

    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    subroutine check_true(name, cond)
        character(*), intent(in) :: name
        logical, intent(in) :: cond
        if (cond) then
            write(*,'(a,a)') '   PASS  ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '   FAIL  ', name
            n_failed = n_failed + 1
        end if
    end subroutine check_true

    subroutine check_close(name, got, expected, rtol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, expected, rtol
        logical :: ok
        ok = abs(got - expected) <= rtol * max(abs(expected), tiny(1.0_dp))
        if (ok) then
            write(*,'(a,a)') '   PASS  ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '   FAIL  ', name
            write(*,'(a,es14.6,a,es14.6)') '         got ', got, '  expected ', expected
            n_failed = n_failed + 1
        end if
    end subroutine check_close

    subroutine configure(grid_cells, domain_height, eps, l_int)
        ! Sets the globals derive_turbulence_scales reads. dz_length is normally
        ! set by initialize.f90:187 as H/N; replicated here.
        integer(i4), intent(in) :: grid_cells
        real(dp), intent(in) :: domain_height, eps, l_int
        N = grid_cells
        H = domain_height
        dz_length = H / real(N, dp)
        dissipation_rate = eps
        integral_length_scale = l_int
        call derive_turbulence_scales()
    end subroutine configure

    subroutine test_regime_sweep()
        ! Every invariant, at every grid resolution and dissipation rate that
        ! keeps the config valid (validate_turbulence_scales exits on the rest).
        integer(i4), parameter :: grids(5) = [700, 1000, 2000, 4000, 8000]
        real(dp), parameter :: epsilons(3) = [0.001, 0.01, 0.1]
        integer :: ig, ie
        logical :: ok_floor, ok_mult3, ok_ge_ld, ok_f, ok_heat, ok_vapor, ok_ratio, ok_scale
        real(dp) :: molecular_ratio

        ok_floor = .true.; ok_mult3 = .true.; ok_ge_ld = .true.; ok_f = .true.
        ok_heat = .true.;  ok_vapor = .true.; ok_ratio = .true.; ok_scale = .true.
        molecular_ratio = kT / Dv

        do ig = 1, size(grids)
            do ie = 1, size(epsilons)
                call configure(grids(ig), 1.0_dp, epsilons(ie), 0.1_dp)

                ok_floor = ok_floor .and. (smallest_eddy_gridpoints >= min_eddy_gridpoints)
                ok_mult3 = ok_mult3 .and. (mod(smallest_eddy_gridpoints, 3) == 0)
                ! The `ceiling` guarantee.
                ok_ge_ld = ok_ge_ld .and. (smallest_eddy_scale >= diffusivity_length_scale)
                ok_f     = ok_f     .and. (diffusivity_enhancement >= 1.0_dp)
                ! The point of the whole scheme.
                ok_heat  = ok_heat  .and. (thermal_diffusivity >= kT)
                ok_vapor = ok_vapor .and. (vapor_diffusivity   >= Dv)
                ! One common factor => Pr and Sc preserved.
                ok_ratio = ok_ratio .and. &
                    (abs(thermal_diffusivity/vapor_diffusivity - molecular_ratio) <= 1.0e-12_dp)
                ok_scale = ok_scale .and. &
                    (abs(smallest_eddy_scale - smallest_eddy_gridpoints*dz_length) <= 1.0e-15_dp)
            end do
        end do

        write(*,'(a)') ' Regime sweep (5 grids x 3 dissipation rates):'
        call check_true("smallest_eddy_gridpoints >= min_eddy_gridpoints", ok_floor)
        call check_true("smallest_eddy_gridpoints is a multiple of 3", ok_mult3)
        call check_true("smallest_eddy_scale >= diffusivity_length_scale (ceiling)", ok_ge_ld)
        call check_true("diffusivity_enhancement >= 1", ok_f)
        call check_true("thermal_diffusivity >= kT (never sub-molecular)", ok_heat)
        call check_true("vapor_diffusivity >= Dv (never sub-molecular)", ok_vapor)
        call check_true("kT:Dv ratio preserved (Pr, Sc unchanged)", ok_ratio)
        call check_true("smallest_eddy_scale == gridpoints * dz_length", ok_scale)
    end subroutine test_regime_sweep

    subroutine test_quantization_bounds()
        ! The enhancement is a STEP function of dz, because smallest_eddy_scale
        ! is quantized to multiples of 3 cells. The two regimes behave
        ! differently and the difference is worth pinning:
        !
        !   grid-limited (6dz >= l_D)  6dz is exactly 2 quanta, so `ceiling`
        !                              never overshoots: n == 6 and
        !                              f == (6dz/l_D)^(4/3) exactly. Continuous
        !                              in dz, and unbounded above as the grid
        !                              coarsens (which is the intent).
        !   diffusivity-limited        the quantum 3dz can be a large fraction
        !   (6dz < l_D)                of l_D -- up to l_D/2 right at the
        !                              crossover -- so `ceiling` can overshoot
        !                              by up to 50%, giving f up to
        !                              (3/2)^(4/3) = 1.7171. f -> 1 only as
        !                              dz -> 0.
        !
        ! In both cases f >= 1, which is the invariant that matters. What is
        ! ruled out is the first slice-2 attempt's ~15x jump at the crossover.
        real(dp), parameter :: eps = 0.01
        real(dp), parameter :: overshoot_bound = 1.5_dp ** (4./3.)
        integer(i4), parameter :: grids(9) = [1005, 1025, 1045, 1200, 1500, 2000, 3000, 4000, 8000]
        integer :: ig
        logical :: ok_bound, ok_grid_exact, ok_diff_bounded
        real(dp) :: f_expected, f_coarsest, f_finest

        ok_bound = .true.; ok_grid_exact = .true.; ok_diff_bounded = .true.

        do ig = 1, size(grids)
            call configure(grids(ig), 1.0_dp, eps, 0.1_dp)

            ! Universal: overshoot is at most one quantum.
            ok_bound = ok_bound .and. &
                (smallest_eddy_scale <= max(grid_eddy_scale, diffusivity_length_scale) &
                                        + 3.0_dp * dz_length + 1.0e-15_dp)

            if (grid_eddy_scale >= diffusivity_length_scale) then
                ! Grid-limited: no overshoot at all, f is exact.
                f_expected = (grid_eddy_scale / diffusivity_length_scale) ** (4./3.)
                ok_grid_exact = ok_grid_exact .and. &
                    (smallest_eddy_gridpoints == min_eddy_gridpoints) .and. &
                    (abs(diffusivity_enhancement - f_expected) <= 1.0e-12_dp)
            else
                ! Diffusivity-limited: bounded by the 50% worst case.
                ok_diff_bounded = ok_diff_bounded .and. &
                    (diffusivity_enhancement < overshoot_bound + 1.0e-12_dp)
            end if
        end do

        write(*,'(a)') ' Quantization behaviour:'
        call check_true("smallest_eddy_scale overshoots by at most one 3-cell quantum", ok_bound)
        call check_true("grid-limited: n == 6 and f == (6dz/l_D)^(4/3) exactly", ok_grid_exact)
        call check_true("diffusivity-limited: f < (3/2)^(4/3) = 1.717", ok_diff_bounded)

        ! Convergence: refining the grid drives the overshoot, and so f, to 1.
        call configure(1045_i4, 1.0_dp, eps, 0.1_dp)
        f_coarsest = diffusivity_enhancement
        call configure(200000_i4, 1.0_dp, eps, 0.1_dp)
        f_finest = diffusivity_enhancement
        call check_true("worst-case overshoot near the crossover is real (f > 1.6)", &
                        f_coarsest > 1.6_dp)
        call check_close("f -> 1 as the grid refines", f_finest, 1.0_dp, 1.0e-3_dp)
    end subroutine test_quantization_bounds

    subroutine test_reference_values()
        ! Pins the bundled input/params.nml configuration and the analytic
        ! identity 0.1*eps^(1/3)*eta^(4/3) == 0.1*nu, which is why
        ! actual_kolmogorov_scale can never govern.
        real(dp) :: turbulent_d_at_eta

        write(*,'(a)') ' Reference configuration (N=2000, H=1, eps=0.01, L=0.01):'
        call configure(2000_i4, 1.0_dp, 0.01_dp, 0.01_dp)

        call check_close("actual_kolmogorov_scale", actual_kolmogorov_scale, 7.5762e-4_dp, 1.0e-4_dp)
        call check_close("grid_eddy_scale (6*dz)",  grid_eddy_scale,         3.0e-3_dp,    1.0e-12_dp)
        call check_close("diffusivity_length_scale", diffusivity_length_scale, 5.8491e-3_dp, 1.0e-4_dp)
        call check_close("smallest_eddy_scale",     smallest_eddy_scale,     6.0e-3_dp,    1.0e-12_dp)
        call check_true ("smallest_eddy_gridpoints == 12", smallest_eddy_gridpoints == 12)
        call check_close("diffusivity_enhancement", diffusivity_enhancement, 1.03454_dp,   1.0e-4_dp)
        call check_close("reynolds_number",         reynolds_number,         1.97605_dp,   1.0e-4_dp)
        call check_true ("this configuration is diffusivity-limited", &
                         grid_eddy_scale < diffusivity_length_scale)

        ! eps cancels: 0.1*eps^(1/3) * ((nu^3/eps)^(1/4))^(4/3) == 0.1*nu.
        turbulent_d_at_eta = 0.1 * dissipation_rate**(1./3.) * actual_kolmogorov_scale**(4./3.)
        call check_close("turbulent D at eta == 0.1*nu (why eta never governs)", &
                         turbulent_d_at_eta, 0.1_dp * nu, 1.0e-6_dp)
        call check_true ("actual_kolmogorov_scale is below diffusivity_length_scale", &
                         actual_kolmogorov_scale < diffusivity_length_scale)
    end subroutine test_reference_values

end program test_lem_scale_derivation
