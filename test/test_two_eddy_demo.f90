program test_two_eddy_demo
    ! A deliberately simple, printed walkthrough of two triplet maps on a
    ! 20-cell domain, and of the maps the droplet transport is built from.
    !
    !   eddy 1: location  4, length 12
    !   eddy 2: location  2, length 15   (applied to the result of eddy 1)
    !
    ! The scalar array starts as 1..20, so each entry is the label of the parcel
    ! that started in that cell. Reading the array after the maps therefore says
    ! "cell k currently holds parcel <value>".
    use globals, only: dp, i4, N, triplet_map, &
                       begin_eddy_sequence, accumulate_eddy, finalize_eddy_sequence, &
                       origin_cell, destination_cell
    implicit none

    integer(i4), parameter :: n_cells = 20_i4
    integer(i4), parameter :: loc1 = 4_i4,  len1 = 12_i4
    integer(i4), parameter :: loc2 = 2_i4,  len2 = 15_i4

    integer(i4) :: field(n_cells)
    integer(i4) :: after_map1(n_cells), after_map2(n_cells)
    integer(i4) :: cells(n_cells)
    integer(i4) :: k
    logical :: tracer_matches, is_inverse

    N = n_cells

    do k = 1, n_cells
        cells(k) = k
        field(k) = k
    end do

    write(*,'(a)') ' Two triplet maps on a 20-cell domain'
    write(*,'(a)') ' ------------------------------------'
    write(*,'(a,i0,a,i0)') '   eddy 1:  location ', loc1, ',  length ', len1
    write(*,'(a,i0,a,i0)') '   eddy 2:  location ', loc2, ',  length ', len2
    write(*,*)

    call show('cell index      ', cells)
    call show('initial field   ', field)
    write(*,*)

    ! --- eddy 1 -------------------------------------------------------------
    call triplet_map(len1, loc1, field)
    after_map1 = field
    call show('after eddy 1    ', after_map1)

    ! --- eddy 2, applied to the already-remapped array -----------------------
    call triplet_map(len2, loc2, field)
    after_map2 = field
    call show('after eddy 2    ', after_map2)
    write(*,*)

    ! --- the same two eddies through the production API ----------------------
    call begin_eddy_sequence()
    call accumulate_eddy(len1, loc1)
    call accumulate_eddy(len2, loc2)

    call show('origin_cell     ', origin_cell)

    call finalize_eddy_sequence()
    call show('destination_cell', destination_cell)
    write(*,*)

    ! --- what the two rows mean ---------------------------------------------
    write(*,'(a)') ' Reading them:'
    write(*,'(a)') '   origin_cell(k)      = parcel now sitting in cell k   ("what is here?")'
    write(*,'(a)') '   destination_cell(c) = cell the parcel from c ended in ("where did I go?")'
    write(*,*)
    write(*,'(a,i0,a,i0,a)') '   e.g. cell 1 now holds parcel ', origin_cell(1), &
        ', while parcel 1 ended up in cell ', destination_cell(1), '.'
    write(*,*)

    ! --- checks --------------------------------------------------------------
    tracer_matches = all(origin_cell == after_map2)
    is_inverse = .true.
    do k = 1, n_cells
        if (destination_cell(origin_cell(k)) /= k) is_inverse = .false.
    end do

    write(*,'(a,l1)') ' origin_cell equals the twice-mapped scalar field : ', tracer_matches
    write(*,'(a,l1)') ' destination_cell(origin_cell(j)) == j for all j   : ', is_inverse

    if (.not. tracer_matches .or. .not. is_inverse) stop 1

contains

    subroutine show(label, arr)
        character(*), intent(in) :: label
        integer(i4), intent(in) :: arr(:)

        write(*,'(a,a,20i4)') '  ', label, arr
    end subroutine show

end program test_two_eddy_demo
