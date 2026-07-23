program test_remap_compare
    ! Side-by-side dump of the two eddy-remapping methods.
    !
    ! OLD: materialize the whole mapped grid with triplet_map, then read
    !      mapped_z(gc) for each cell.
    ! NEW: query one cell at a time with triplet_map_cell, i.e. z(src).
    !
    ! For each cell we print both mapped values and their difference. Any
    ! nonzero difference is a real divergence between the methods. Two eddies
    ! are shown: one entirely inside the domain, one that wraps the periodic
    ! boundary.
    use globals, only: dp, i4, N, H, z, dz_length, triplet_map, triplet_map_cell
    implicit none

    integer :: n_mismatch

    n_mismatch = 0

    call setup_domain(24, 2.4_dp)

    ! Mid-domain eddy: M=10, L=6 -> cells 10..15, no wrapping.
    call dump_eddy("MID-DOMAIN eddy  (M=10, L=6, cells 10-15)", 6, 10)

    ! Boundary-wrapping eddy: M=22, L=6 -> cells 22,23,24,1,2,3.
    call dump_eddy("WRAPPING eddy    (M=22, L=6, cells 22,23,24,1,2,3)", 6, 22)

    write(*,*)
    if (n_mismatch == 0) then
        write(*,'(a)') ' Results: methods agree on every cell (0 mismatches)'
    else
        write(*,'(a,i0,a)') ' Results: ', n_mismatch, ' mismatches'
        stop 1
    end if

contains

    subroutine setup_domain(n_cells, domain_height)
        integer, intent(in) :: n_cells
        real(dp), intent(in) :: domain_height
        integer :: k

        N = n_cells
        H = domain_height
        dz_length = H / N
        if (allocated(z)) deallocate(z)
        allocate(z(N))
        do k = 1, N
            z(k) = H * k / N
        end do
    end subroutine setup_domain


    subroutine dump_eddy(label, L, M)
        character(*), intent(in) :: label
        integer(i4), intent(in) :: L, M
        real(dp) :: mapped_z(N)
        real(dp) :: old_val, new_val, diff
        integer(i4) :: gc, src

        ! OLD method: build the full mapped grid.
        mapped_z = z
        call triplet_map(L, M, mapped_z)

        write(*,*)
        write(*,'(a)') ' '//label
        write(*,'(a)') '   cell   src   old mapped_z(gc)   new z(src)         diff       moved'
        write(*,'(a)') '   ----   ---   ----------------   ----------------   --------   -----'

        do gc = 1, N
            ! NEW method: one-cell lookup.
            src = triplet_map_cell(L, M, gc)
            old_val = mapped_z(gc)
            new_val = z(src)
            diff = new_val - old_val

            write(*,'(i6,i6,es20.10,es19.10,es13.4,a6)') &
                gc, src, old_val, new_val, diff, &
                trim(merge('  yes', '   no', src /= gc))

            if (diff /= 0.0_dp) n_mismatch = n_mismatch + 1
        end do
    end subroutine dump_eddy

end program test_remap_compare
