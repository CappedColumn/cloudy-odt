program test_collection_efficiency
    use globals, only: dp
    use collection_efficiency, only: coalescence_kernel, set_kernel_selector, &
                                     collection_efficiency_E
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    ! ---- Unity kernel ----
    coalescence_kernel = 'unity'
    call set_kernel_selector()

    call check_range("unity: small equal droplets", &
        collection_efficiency_E(10e-6_dp, 10e-6_dp), 1.0_dp, 1.0_dp)
    call check_range("unity: large collector, small collectee", &
        collection_efficiency_E(200e-6_dp, 5e-6_dp), 1.0_dp, 1.0_dp)

    ! ---- Long kernel (stand-in formula, tests disabled) ----
    ! coalescence_kernel = 'long'
    ! call set_kernel_selector()
    ! TODO: enable once Long (1974) coefficients are verified

    ! ---- Hall kernel ----
    coalescence_kernel = 'hall'
    call set_kernel_selector()

    ! Table corner: R=300um, p=0 -> E=0.97
    call check_range("hall: R=300um, p~0 -> E~0.97", &
        collection_efficiency_E(300e-6_dp, 1e-9_dp), 0.90_dp, 1.0_dp)

    ! Small collector: R=10um, p=0.5 -> E=0.035
    call check_range("hall: R=10um, r=5um -> E~0.035", &
        collection_efficiency_E(10e-6_dp, 5e-6_dp), 0.01_dp, 0.08_dp)

    ! Mid-range: R=60um, p=0.5 -> E~0.91
    call check_range("hall: R=60um, r=30um -> E~0.91", &
        collection_efficiency_E(60e-6_dp, 30e-6_dp), 0.80_dp, 1.0_dp)

    ! Interpolation: R=75um (between 70 and 100), p=0.5
    call check_range("hall: R=75um, r=37.5um -> interpolated", &
        collection_efficiency_E(75e-6_dp, 37.5e-6_dp), 0.85_dp, 1.0_dp)

    ! Below table minimum: R=3um should clamp to R=5um row
    call check_range("hall: R=3um clamps to table min", &
        collection_efficiency_E(3e-6_dp, 1.5e-6_dp), 0.0_dp, 0.1_dp)

    ! E must be non-negative
    call check_range("hall: E >= 0 everywhere", &
        collection_efficiency_E(20e-6_dp, 10e-6_dp), 0.0_dp, 1.0_dp)

    ! ---- Monotonicity: larger collector -> higher E (Hall) ----
    coalescence_kernel = 'hall'
    call set_kernel_selector()
    ! Varying collector size (fixed collectee)
    call check_monotone("hall: R=20->100um, r=10um fixed", &
        collection_efficiency_E(20e-6_dp, 10e-6_dp), &
        collection_efficiency_E(100e-6_dp, 10e-6_dp))
    call check_monotone("hall: R=50->200um, r=20um fixed", &
        collection_efficiency_E(50e-6_dp, 20e-6_dp), &
        collection_efficiency_E(200e-6_dp, 20e-6_dp))

    ! Varying collectee size (fixed collector)
    call check_monotone("hall: R=100um, r=5->30um", &
        collection_efficiency_E(100e-6_dp, 5e-6_dp), &
        collection_efficiency_E(100e-6_dp, 30e-6_dp))
    call check_monotone("hall: R=200um, r=10->50um", &
        collection_efficiency_E(200e-6_dp, 10e-6_dp), &
        collection_efficiency_E(200e-6_dp, 50e-6_dp))
    call check_monotone("hall: R=60um, r=5->20um", &
        collection_efficiency_E(60e-6_dp, 5e-6_dp), &
        collection_efficiency_E(60e-6_dp, 20e-6_dp))

    ! ---- Summary ----
    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' tests failed'
        stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' tests passed'

contains

    subroutine check_range(label, val, lo, hi)
        character(*), intent(in) :: label
        real(dp), intent(in) :: val, lo, hi

        if (val >= lo .and. val <= hi) then
            n_passed = n_passed + 1
            print '(a,a,a,es10.3,a,es10.3,a,es10.3,a)', &
                '  PASS  ', label, '  (E=', val, ', range [', lo, ', ', hi, '])'
        else
            n_failed = n_failed + 1
            print '(a,a,a,es10.3,a,es10.3,a,es10.3,a)', &
                '  FAIL  ', label, '  (E=', val, ', expected [', lo, ', ', hi, '])'
        end if
    end subroutine check_range

    subroutine check_monotone(label, e_small, e_large)
        character(*), intent(in) :: label
        real(dp), intent(in) :: e_small, e_large

        if (e_large > e_small) then
            n_passed = n_passed + 1
            print '(a,a,a,es10.3,a,es10.3,a)', &
                '  PASS  ', label, '  (', e_small, ' < ', e_large, ')'
        else
            n_failed = n_failed + 1
            print '(a,a,a,es10.3,a,es10.3,a)', &
                '  FAIL  ', label, '  (', e_small, ' >= ', e_large, ')'
        end if
    end subroutine check_monotone

end program test_collection_efficiency
