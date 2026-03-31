program demo_triplet_map
    use globals, only: dp, i4, N
    use LEM, only: triplet_map_periodic
    implicit none

    integer(i4) :: n_domain, eddy_len, eddy_start
    integer :: u, k

    open(newunit=u, file='triplet_map_demo.dat', status='replace', action='write')

    ! --- Case 1: non-wrapping eddy ---
    n_domain = 12
    eddy_len = 6
    eddy_start = 4
    call write_case(u, n_domain, eddy_len, eddy_start)

    ! --- Case 2: wrapping eddy ---
    eddy_start = 11
    call write_case(u, n_domain, eddy_len, eddy_start)

    ! --- Case 3: larger domain, wrapping ---
    n_domain = 21
    eddy_len = 9
    eddy_start = 18
    call write_case(u, n_domain, eddy_len, eddy_start)

    close(u)
    write(*,*) 'Wrote triplet_map_demo.dat'

contains

    subroutine write_case(unit_num, nd, el, es)
        integer, intent(in) :: unit_num
        integer(i4), intent(in) :: nd, el, es
        real(dp) :: before(nd), after(nd)
        integer(i4) :: k

        N = nd
        do k = 1, nd
            before(k) = real(k, dp)
        end do
        after = before
        call triplet_map_periodic(el, es, after)

        ! Header: n_domain eddy_length eddy_start
        write(unit_num, '(3i6)') nd, el, es
        ! Data: cell_index before after
        do k = 1, nd
            write(unit_num, '(i6, 2f12.4)') k, before(k), after(k)
        end do

    end subroutine write_case

end program demo_triplet_map
