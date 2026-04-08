program test_triplet_map
    use globals, only: dp, i4, N, triplet_map
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    ! --- Equivalence tests: periodic == non-periodic for non-wrapping eddies ---

    call test_equivalence_start()
    call test_equivalence_middle()
    call test_equivalence_end()

    ! --- Periodic wrapping edge-case tests ---

    call test_wrap_around_boundary()
    call test_eddy_starts_at_last_cell()
    call test_full_domain_eddy()

    ! --- Summary ---

    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    subroutine check_array(name, got, expected, n_elem)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got(:), expected(:)
        integer, intent(in) :: n_elem
        integer :: k
        logical :: match

        match = .true.
        do k = 1, n_elem
            if (abs(got(k) - expected(k)) > 1.0e-12) then
                match = .false.
                exit
            end if
        end do

        if (match) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,i0,a)') '    first mismatch at k=', k, ':'
            write(*,'(a,es15.7)') '      expected: ', expected(k)
            write(*,'(a,es15.7)') '      got:      ', got(k)
            n_failed = n_failed + 1
        end if
    end subroutine check_array


    subroutine test_equivalence_start()
        ! Eddy at start of domain: L=9, M=1, N=18
        real(dp) :: field_periodic(18), field_standard(18)
        integer(i4) :: k

        N = 18
        do k = 1, N
            field_periodic(k) = real(k, dp) * 10.0
            field_standard(k) = real(k, dp) * 10.0
        end do

        call triplet_map(9, 1, field_standard)
        call triplet_map(9, 1, field_periodic)

        call check_array("equivalence: eddy at start (L=9, M=1)", &
                          field_periodic, field_standard, N)
    end subroutine test_equivalence_start


    subroutine test_equivalence_middle()
        ! Eddy in the middle: L=6, M=4, N=18
        real(dp) :: field_periodic(18), field_standard(18)
        integer(i4) :: k

        N = 18
        do k = 1, N
            field_periodic(k) = real(k, dp) * 10.0
            field_standard(k) = real(k, dp) * 10.0
        end do

        call triplet_map(6, 4, field_standard)
        call triplet_map(6, 4, field_periodic)

        call check_array("equivalence: eddy in middle (L=6, M=4)", &
                          field_periodic, field_standard, N)
    end subroutine test_equivalence_middle


    subroutine test_equivalence_end()
        ! Eddy at end of domain (no wrapping): L=6, M=13, N=18
        real(dp) :: field_periodic(18), field_standard(18)
        integer(i4) :: k

        N = 18
        do k = 1, N
            field_periodic(k) = real(k, dp) * 10.0
            field_standard(k) = real(k, dp) * 10.0
        end do

        call triplet_map(6, 13, field_standard)
        call triplet_map(6, 13, field_periodic)

        call check_array("equivalence: eddy at end (L=6, M=13)", &
                          field_periodic, field_standard, N)
    end subroutine test_equivalence_end


    subroutine test_wrap_around_boundary()
        ! Eddy wraps: N=12, eddy_start=11, eddy_length=6
        ! Eddy covers cells: 11, 12, 1, 2, 3, 4
        ! Permutation [0,3,4,1,2,5] on values [110,120,10,20,30,40]
        ! Result: cells 11=110, 12=20, 1=30, 2=120, 3=10, 4=40
        real(dp) :: field(12), expected(12)
        integer(i4) :: k

        N = 12
        do k = 1, N
            field(k) = real(k, dp) * 10.0
        end do

        expected = [30.0_dp, 120.0_dp, 10.0_dp, 40.0_dp, 50.0_dp, 60.0_dp, &
                    70.0_dp, 80.0_dp, 90.0_dp, 100.0_dp, 110.0_dp, 20.0_dp]

        call triplet_map(6, 11, field)

        call check_array("periodic wrap: eddy across boundary (start=11, L=6, N=12)", &
                          field, expected, N)
    end subroutine test_wrap_around_boundary


    subroutine test_eddy_starts_at_last_cell()
        ! Eddy starts at last cell: N=12, eddy_start=12, eddy_length=6
        ! Eddy covers cells: 12, 1, 2, 3, 4, 5
        ! Values at those cells: [120, 10, 20, 30, 40, 50]
        ! Permutation [0,3,4,1,2,5]: [120, 30, 40, 10, 20, 50]
        ! Written to cells 12,1,2,3,4,5
        real(dp) :: field(12), expected(12)
        integer(i4) :: k

        N = 12
        do k = 1, N
            field(k) = real(k, dp) * 10.0
        end do

        expected = [30.0_dp, 40.0_dp, 10.0_dp, 20.0_dp, 50.0_dp, 60.0_dp, &
                    70.0_dp, 80.0_dp, 90.0_dp, 100.0_dp, 110.0_dp, 120.0_dp]

        call triplet_map(6, 12, field)

        call check_array("periodic wrap: eddy starts at last cell (start=12, L=6, N=12)", &
                          field, expected, N)
    end subroutine test_eddy_starts_at_last_cell


    subroutine test_full_domain_eddy()
        ! Full-domain eddy: N=9, eddy_start=1, eddy_length=9
        ! Should match non-periodic triplet_map exactly
        real(dp) :: field_periodic(9), field_standard(9)
        integer(i4) :: k

        N = 9
        do k = 1, N
            field_periodic(k) = real(k, dp) * 10.0
            field_standard(k) = real(k, dp) * 10.0
        end do

        call triplet_map(9, 1, field_standard)
        call triplet_map(9, 1, field_periodic)

        call check_array("periodic wrap: full domain eddy (L=N=9)", &
                          field_periodic, field_standard, N)
    end subroutine test_full_domain_eddy

end program test_triplet_map
