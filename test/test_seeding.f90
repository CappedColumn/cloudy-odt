! Unit test for aerosol seeding: the seed group reader, the bin -> type mapping,
! the two-population sampler, and the Kohler equilibrium radius solve.
!
! Writes its own aerosol file with a seed group so the checks are against known
! values rather than whatever the bundled input happens to contain. The file has
! two background bins (categories 1-2, type 1) and two seed bins (category 3,
! type 2), which is the multi-type case the bundled input cannot exercise.
program test_seeding
    use netcdf
    use globals, only: dp, i4, m_per_nm, pi_43
    use particle_types, only: particle, aerosol
    use droplets, only: read_aerosol_netcdf, sample_radius, aerosols, &
                        aerosol_radii, aerosol_partition, aerosol_bin_type, &
                        aerosol_bin_freq, n_seed_bins, seed_radii, seed_partition, &
                        seed_bin_type, seed_bin_freq, seed_event_coord, seed_event_conc
    implicit none

    integer :: n_passed, n_failed
    character(*), parameter :: test_file = 'test_seed_aerosol.nc'

    n_passed = 0
    n_failed = 0

    call write_seeded_aerosol_file(test_file)
    call read_aerosol_netcdf(test_file)

    call check_reader()
    call check_type_mapping()
    call check_sampler()
    call check_kohler_equilibrium()

    ! --- Summary ---
    write(*,*)
    write(*,'(a,i0,a,i0,a)') ' test_seeding: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) stop 1

contains

    subroutine check_reader()
        ! The seed group is read as its own population, separate from background.

        call check_int("n seed bins", n_seed_bins, 2)
        call check_int("n background bins", size(aerosol_radii), 2)
        call check_int("n aerosol types", size(aerosols), 2)

        call check_real("seed midpoint radius(1)", seed_radii(1), 750.0_dp)
        call check_real("seed midpoint radius(2)", seed_radii(2), 1500.0_dp)

        call check_int("n seed events", size(seed_event_coord), 2)
        call check_real("seed event coord(1)", seed_event_coord(1), 200.0_dp)
        call check_real("seed event coord(2)", seed_event_coord(2), 400.0_dp)
        call check_real("seed concentration(1)", seed_event_conc(1), 5.0_dp)

        ! Per-event CDF rows: the two events release different size mixes
        call check_real("seed cdf(1, event 1)", seed_bin_freq(1,1), 0.6_dp)
        call check_real("seed cdf(2, event 1)", seed_bin_freq(2,1), 1.0_dp)
        call check_real("seed cdf(1, event 2)", seed_bin_freq(1,2), 0.25_dp)
        call check_real("seed cdf(2, event 2)", seed_bin_freq(2,2), 1.0_dp)

    end subroutine check_reader

    subroutine check_type_mapping()
        ! A type is seed material iff a seed bin maps to it. Nothing in the file
        ! states this directly, so it is worth pinning down.

        call check_int("background bin_type(1)", aerosol_bin_type(1), 1)
        call check_int("background bin_type(2)", aerosol_bin_type(2), 1)
        call check_int("seed_bin_type(1)", seed_bin_type(1), 2)
        call check_int("seed_bin_type(2)", seed_bin_type(2), 2)

        call check_logical("type 1 is not seed", aerosols(1)%is_seed, .false.)
        call check_logical("type 2 is seed", aerosols(2)%is_seed, .true.)

        ! The two types carry genuinely different composition, and the seed row
        ! is not just a copy of row 1.
        call check_int("type 1 n_ions", aerosols(1)%n_ions, 2)
        call check_int("type 2 n_ions", aerosols(2)%n_ions, 3)
        call check_real("type 2 solute density", aerosols(2)%solute_density, 1725.0_dp)

    end subroutine check_type_mapping

    subroutine check_sampler()
        ! Each population is drawn from its own arrays. The CDFs here are
        ! degenerate on purpose (all mass in one bin) so the draw is deterministic
        ! and the test does not depend on the RNG.
        real(dp) :: radius
        integer(i4) :: category, type_idx
        real(dp) :: forced_cdf(2)

        ! All mass in seed bin 1 -> must sample bin 1: category 3, type 2
        forced_cdf = [1.0_dp, 1.0_dp]
        call sample_radius(forced_cdf, seed_radii, seed_partition, seed_bin_type, &
                           radius, category, type_idx)
        call check_real("seed sample takes bin 1 radius", radius, 750.0_dp * m_per_nm)
        call check_int("seed sample carries its category", category, 3)
        call check_int("seed sample carries its type", type_idx, 2)

        ! Same for the background: category 1, type 1
        call sample_radius(forced_cdf, aerosol_radii, aerosol_partition, aerosol_bin_type, &
                           radius, category, type_idx)
        call check_int("background sample carries its category", category, 1)
        call check_int("background sample carries its type", type_idx, 1)

        ! Zero mass in bin 1 -> must fall through to bin 2
        forced_cdf = [0.0_dp, 1.0_dp]
        call sample_radius(forced_cdf, seed_radii, seed_partition, seed_bin_type, &
                           radius, category, type_idx)
        call check_real("seed sample falls through to bin 2", radius, 1500.0_dp * m_per_nm)

    end subroutine check_sampler

    subroutine check_kohler_equilibrium()
        ! The equilibrium radius must sit on the stable branch and actually satisfy
        ! the Kohler equation at the requested humidity.
        type(particle) :: p
        real(dp) :: r_eq, RH, s_eq

        ! A 100 nm NaCl haze particle at 293 K
        p%temperature = 293.0
        p%solute_type = aerosols(1)
        p%solute_radius = 100.0 * m_per_nm
        p%solute_gross_mass = pi_43 * p%solute_type%solute_density * p%solute_radius**3
        call p%critical_kohler()

        RH = 0.95
        r_eq = p%kohler_equilibrium_radius(RH)

        call check_logical("equilibrium radius exceeds the dry radius", &
                           r_eq > p%solute_radius, .true.)
        call check_logical("equilibrium radius stays below the critical radius", &
                           r_eq < p%critical_radius, .true.)

        ! Residual of the Kohler equation at the solved radius: S(r) - 1 should
        ! equal RH - 1. Tolerance is loose because a and b are approximations, but
        ! a wrong branch or a bad bracket would miss by orders of magnitude.
        s_eq = kohler_supersaturation(p, r_eq)
        call check_logical("solved radius satisfies the Kohler equation at RH=0.95", &
                           abs(s_eq - (RH - 1.0)) < 1.0e-6, .true.)

        ! Drier air must hold a smaller haze droplet
        call check_logical("lower RH gives a smaller equilibrium radius", &
                           p%kohler_equilibrium_radius(0.80_dp) < r_eq, .true.)

    end subroutine check_kohler_equilibrium

    real(dp) function kohler_supersaturation(p, radius) result(s)
        ! S(r) - 1 = a/r - b/r^3, independently of the solver, to check its answer.
        ! Mirrors the coefficients in particle_types (cm units).
        use globals, only: a_RY, m_per_cm
        type(particle), intent(in) :: p
        real(dp), intent(in) :: radius
        real(dp) :: a, b, r_cm

        a = a_RY / p%temperature
        b = 4.3 * p%solute_gross_mass * p%solute_type%n_ions / p%solute_type%solute_molar_mass
        r_cm = radius / m_per_cm
        s = a / r_cm - b / r_cm**3

    end function kohler_supersaturation

    subroutine write_seeded_aerosol_file(path)
        ! Minimal CODT_aerosol_input_v1 file carrying a seed group. Doubles as a
        ! worked example of the schema.
        character(*), intent(in) :: path
        integer :: ncid, varid
        integer :: d_type, d_edge, d_bin, d_time, d_dsd
        integer :: d_sbin, d_sedge, d_sevent
        integer :: i
        real(dp) :: dsd_edges(3)

        call nc_check(nf90_create(path, NF90_CLOBBER, ncid))

        call nc_check(nf90_def_dim(ncid, 'aerosol_type', 2, d_type))
        call nc_check(nf90_def_dim(ncid, 'edge', 3, d_edge))
        call nc_check(nf90_def_dim(ncid, 'bin', 2, d_bin))
        call nc_check(nf90_def_dim(ncid, 'time', 1, d_time))
        call nc_check(nf90_def_dim(ncid, 'dsd_edge', 3, d_dsd))
        call nc_check(nf90_def_dim(ncid, 'seed_bin', 2, d_sbin))
        call nc_check(nf90_def_dim(ncid, 'seed_edge', 3, d_sedge))
        call nc_check(nf90_def_dim(ncid, 'seed_event', 2, d_sevent))

        call def_and_close(ncid, 'n_ions', NF90_INT, [d_type])
        call def_and_close(ncid, 'molar_mass', NF90_DOUBLE, [d_type])
        call def_and_close(ncid, 'solute_density', NF90_DOUBLE, [d_type])
        call def_and_close(ncid, 'edge_radii', NF90_DOUBLE, [d_edge])
        call def_and_close(ncid, 'category', NF90_INT, [d_bin])
        call def_and_close(ncid, 'bin_type', NF90_INT, [d_bin])
        call def_and_close(ncid, 'cumulative_frequency', NF90_DOUBLE, [d_bin, d_time])
        call def_and_close(ncid, 'injection_time', NF90_DOUBLE, [d_time])
        call def_and_close(ncid, 'injection_rate', NF90_DOUBLE, [d_time])
        call def_and_close(ncid, 'dsd_bin_edges', NF90_DOUBLE, [d_dsd])
        call def_and_close(ncid, 'seed_edge_radii', NF90_DOUBLE, [d_sedge])
        call def_and_close(ncid, 'seed_category', NF90_INT, [d_sbin])
        call def_and_close(ncid, 'seed_bin_type', NF90_INT, [d_sbin])
        call def_and_close(ncid, 'seed_frequency', NF90_DOUBLE, [d_sbin, d_sevent])
        call def_and_close(ncid, 'seed_coord', NF90_DOUBLE, [d_sevent])
        call def_and_close(ncid, 'seed_concentration', NF90_DOUBLE, [d_sevent])

        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'conventions', 'CODT_aerosol_input_v1'))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'aerosol_name', 'NaCl'))
        call nc_check(nf90_enddef(ncid))

        ! Type 1 = NaCl background, type 2 = a distinct seed material
        call put_int(ncid, 'n_ions', [2, 3])
        call put_real(ncid, 'molar_mass', [58.4428e-3_dp, 132.14e-3_dp])
        call put_real(ncid, 'solute_density', [2163.0_dp, 1725.0_dp])

        call put_real(ncid, 'edge_radii', [60.0_dp, 70.0_dp, 4930.0_dp])
        call put_int(ncid, 'category', [1, 2])
        call put_int(ncid, 'bin_type', [1, 1])
        call put_real_2d(ncid, 'cumulative_frequency', reshape([1.0_dp, 1.0_dp], [2, 1]))
        call put_real(ncid, 'injection_time', [0.0_dp])
        call put_real(ncid, 'injection_rate', [5.5e5_dp])
        dsd_edges = [(0.1_dp * i, i = 1, 3)]
        call put_real(ncid, 'dsd_bin_edges', dsd_edges)

        ! Seed: coarse bins, one category, two events with different size mixes
        call put_real(ncid, 'seed_edge_radii', [500.0_dp, 1000.0_dp, 2000.0_dp])
        call put_int(ncid, 'seed_category', [3, 3])
        call put_int(ncid, 'seed_bin_type', [2, 2])
        ! (seed_bin, seed_event) in Fortran order: event 1 is 60/40, event 2 is 25/75
        call put_real_2d(ncid, 'seed_frequency', &
                         reshape([0.6_dp, 1.0_dp, 0.25_dp, 1.0_dp], [2, 2]))
        call put_real(ncid, 'seed_coord', [200.0_dp, 400.0_dp])
        call put_real(ncid, 'seed_concentration', [5.0_dp, 10.0_dp])

        call nc_check(nf90_close(ncid))

    end subroutine write_seeded_aerosol_file

    subroutine def_and_close(ncid, name, xtype, dims)
        integer, intent(in) :: ncid, xtype, dims(:)
        character(*), intent(in) :: name
        integer :: varid

        call nc_check(nf90_def_var(ncid, name, xtype, dims, varid))

    end subroutine def_and_close

    subroutine put_int(ncid, name, values)
        integer, intent(in) :: ncid, values(:)
        character(*), intent(in) :: name
        integer :: varid

        call nc_check(nf90_inq_varid(ncid, name, varid))
        call nc_check(nf90_put_var(ncid, varid, values))

    end subroutine put_int

    subroutine put_real(ncid, name, values)
        integer, intent(in) :: ncid
        real(dp), intent(in) :: values(:)
        character(*), intent(in) :: name
        integer :: varid

        call nc_check(nf90_inq_varid(ncid, name, varid))
        call nc_check(nf90_put_var(ncid, varid, values))

    end subroutine put_real

    subroutine put_real_2d(ncid, name, values)
        integer, intent(in) :: ncid
        real(dp), intent(in) :: values(:,:)
        character(*), intent(in) :: name
        integer :: varid

        call nc_check(nf90_inq_varid(ncid, name, varid))
        call nc_check(nf90_put_var(ncid, varid, values))

    end subroutine put_real_2d

    subroutine nc_check(status)
        integer, intent(in) :: status

        if (status /= NF90_NOERR) then
            write(*,*) 'NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if

    end subroutine nc_check

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

        if (abs(got - expected) < 1.0e-10_dp * max(1.0_dp, abs(expected))) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,es20.12)') '    expected: ', expected
            write(*,'(a,es20.12)') '    got:      ', got
            n_failed = n_failed + 1
        end if
    end subroutine check_real

    subroutine check_logical(name, got, expected)
        character(*), intent(in) :: name
        logical, intent(in) :: got, expected

        if (got .eqv. expected) then
            write(*,'(a,a)') '  PASS: ', name
            n_passed = n_passed + 1
        else
            write(*,'(a,a)') '  FAIL: ', name
            write(*,'(a,l1,a,l1)') '    expected: ', expected, '  got: ', got
            n_failed = n_failed + 1
        end if
    end subroutine check_logical

end program test_seeding
