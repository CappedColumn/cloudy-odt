program main
    use version, only: code_version, git_commit, print_usage
    implicit none

    character(256) :: arg

    if (command_argument_count() < 1) then
        call print_usage()
        call exit(1)
    end if

    call get_command_argument(1, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
        call print_usage()
        call exit(0)
    else if (trim(arg) == '--version' .or. trim(arg) == '-v') then
        write(*,'(a,a,a,a,a)') 'CODT ', trim(code_version), ' (', trim(git_commit), ')'
        call exit(0)
    else if (arg(1:2) == '--') then
        write(0,'(a,a)') 'Error: unknown option: ', trim(arg)
        write(0,'(a)') 'Run "codt --help" for usage information.'
        call exit(1)
    end if

    call run_codt(trim(arg))

contains

    subroutine run_codt(namelist_path_arg)
        use globals, only: namelist_path, namelist_dir, parent_directory
        use CODT, only: run_simulation
        character(*), intent(in) :: namelist_path_arg

        if (scan(namelist_path_arg, '/') == 0) then
            write(0,*) 'Error: namelist path must include a directory.'
            write(0,*) 'Use ./params.nml for the current directory.'
            call exit(1)
        end if

        namelist_path = namelist_path_arg
        namelist_dir = parent_directory(namelist_path)
        call run_simulation()

    end subroutine run_codt

end program main
