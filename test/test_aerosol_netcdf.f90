program test_aerosol_netcdf
    use globals, only: dp, i4
    use droplets, only: read_aerosol_netcdf, aerosol_size_edges, injection_times, &
                        injection_rates, aerosol_partition, aerosol_bin_freq, &
                        aerosol_radii, aerosols, particle_bin_edges, n_DSD_bins
    implicit none

    integer :: n_passed, n_failed
    character(*), parameter :: test_file = 'input/aerosol_input.nc'

    n_passed = 0
    n_failed = 0

    call read_aerosol_netcdf(test_file)

    ! --- Dimension checks ---
    call check_int("n edges", size(aerosol_size_edges), 3)
    call check_int("n bins", size(aerosol_radii), 2)
    call check_int("n times (injection_times)", size(injection_times), 1)
    call check_int("n times (injection_rates)", size(injection_rates), 1)
    call check_int("n partition", size(aerosol_partition), 2)
    call check_int("bin_freq dim 1 (bins)", size(aerosol_bin_freq, 1), 2)
    call check_int("bin_freq dim 2 (times)", size(aerosol_bin_freq, 2), 1)

    ! --- Value checks ---
    call check_real("edge_radii(1)", aerosol_size_edges(1), 60.0_dp)
    call check_real("edge_radii(2)", aerosol_size_edges(2), 70.0_dp)
    call check_real("edge_radii(3)", aerosol_size_edges(3), 4930.0_dp)

    call check_real("midpoint radii(1)", aerosol_radii(1), 65.0_dp)
    call check_real("midpoint radii(2)", aerosol_radii(2), 2500.0_dp)

    call check_real("injection_time(1)", injection_times(1), 0.0_dp)
    call check_real("injection_rate(1)", injection_rates(1), 5.5e5_dp)

    call check_int("category(1)", aerosol_partition(1), 1)
    call check_int("category(2)", aerosol_partition(2), 2)

    ! CDF values: bin_freq(bin, time)
    call check_real("cdf(1,1)", aerosol_bin_freq(1,1), 1.0_dp)
    call check_real("cdf(2,1)", aerosol_bin_freq(2,1), 1.0_dp)

    ! DSD bin edges
    call check_int("n_DSD_bins", n_DSD_bins, 200)
    call check_int("n dsd_bin_edges", size(particle_bin_edges), 201)
    call check_real("dsd_bin_edges(1)", particle_bin_edges(1), 0.049_dp)
    call check_real("dsd_bin_edges(201)", particle_bin_edges(201), 65.262_dp)

    ! Aerosol properties
    call check_int("n_ions", aerosols(1)%n_ions, 2)
    call check_real("molar_mass", aerosols(1)%solute_molar_mass, 58.4428e-3_dp)
    call check_real("solute_density", aerosols(1)%solute_density, 2163.0e0_dp)

    ! --- Summary ---
    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    subroutine check_int(name, got, expected)
        character(*), intent(in) :: name
        integer, intent(in) :: got, expected

        if (got == expected) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,i0)') '    expected: ', expected
            write(*,'(a,i0)') '    got:      ', got
            n_failed = n_failed + 1
        end if
    end subroutine check_int

    subroutine check_real(name, got, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, expected

        if (abs(got - expected) < 1.0e-10_dp) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,es20.12)') '    expected: ', expected
            write(*,'(a,es20.12)') '    got:      ', got
            n_failed = n_failed + 1
        end if
    end subroutine check_real

end program test_aerosol_netcdf
