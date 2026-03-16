program test_paths
    use globals, only: parent_directory, resolve_path, copy_file
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    ! --- parent_directory tests ---

    call check("parent_directory: absolute path", &
        trim(parent_directory("/home/user/input/params.nml")), &
        "/home/user/input/")

    call check("parent_directory: nested path", &
        trim(parent_directory("/a/b/c/file.txt")), &
        "/a/b/c/")

    call check("parent_directory: root file", &
        trim(parent_directory("/params.nml")), &
        "/")

    call check("parent_directory: no directory (bare filename)", &
        trim(parent_directory("params.nml")), &
        "./")

    ! --- resolve_path tests ---

    call check("resolve_path: relative file", &
        trim(resolve_path("/home/user/input/", "injection_data.txt")), &
        "/home/user/input/injection_data.txt")

    call check("resolve_path: absolute file unchanged", &
        trim(resolve_path("/home/user/input/", "/shared/bin_data.txt")), &
        "/shared/bin_data.txt")

    call check("resolve_path: relative with subdirectory", &
        trim(resolve_path("/home/user/input/", "data/file.txt")), &
        "/home/user/input/data/file.txt")

    call check("resolve_path: relative with parent traversal", &
        trim(resolve_path("/home/user/input/", "../shared/bin_data.txt")), &
        "/home/user/input/../shared/bin_data.txt")

    ! --- copy_file tests ---

    call test_copy_file()

    ! --- Summary ---

    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) then
        stop 1
    end if

contains

    subroutine check(name, got, expected)
        character(*), intent(in) :: name, got, expected

        if (got == expected) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,a)') '    expected: ', expected
            write(*,'(a,a)') '    got:      ', got
            n_failed = n_failed + 1
        end if
    end subroutine check

    subroutine test_copy_file()
        integer :: u, ierr
        character(256) :: line1, line2
        character(*), parameter :: src = '/tmp/codt_test_src.txt'
        character(*), parameter :: dst = '/tmp/codt_test_dst.txt'

        ! Write a test source file
        open(newunit=u, file=src, status='replace', action='write')
        write(u, '(a)') 'Line one of test file'
        write(u, '(a)') 'Line two with numbers 12345'
        close(u)

        ! Copy it
        call copy_file(src, dst)

        ! Read back and verify
        open(newunit=u, file=dst, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            call check("copy_file: destination exists", "missing", "exists")
            return
        end if

        read(u, '(a)') line1
        read(u, '(a)') line2
        close(u)

        call check("copy_file: line 1 matches", trim(line1), "Line one of test file")
        call check("copy_file: line 2 matches", trim(line2), "Line two with numbers 12345")

        ! Clean up
        open(newunit=u, file=src, status='old'); close(u, status='delete')
        open(newunit=u, file=dst, status='old'); close(u, status='delete')

    end subroutine test_copy_file

end program test_paths
