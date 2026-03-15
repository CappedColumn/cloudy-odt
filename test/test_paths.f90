program test_paths
    use globals, only: parent_directory, resolve_path
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

end program test_paths
