! Longwave radiative heating for chamber-mode CODT.
!
! Two solvers available via radiation_method namelist parameter:
!   '1d' — 1D two-stream (diffusivity-factor approximation)
!   '3d' — 3D Monte Carlo photon path-length method
!
! Ported from Suryadev Singh's radiation module (Oct 2025).
module radiation
    use globals, only: dp, i4, N, H, z, T, gridcell_volume, dz_length, &
                       do_radiation, budget_radiation_delta_T, &
                       nc_verify, resolve_path, namelist_path, namelist_dir, &
                       pi, pi_43, rho_l, cp, c_l
    use droplets, only: particles, current_n_particles
    implicit none

    private
    public :: initialize_radiation, compute_radiation, finalize_radiation
    public :: rad_F_net, rad_heating_rate
    public :: radiation_method, mie_data_file, eps_top, eps_bot, sky_temp, &
              sky_cooling_flag, max_droplets_per_cell, rad_call_interval, &
              nPhotons, nBins, Lx_rad, Ly_rad, T_side

    ! --- RADIATION namelist parameters ---
    character(4)  :: radiation_method = '1d'
    character(256) :: mie_data_file = ''
    real(dp) :: eps_top = 1.0
    real(dp) :: eps_bot = 1.0
    real(dp) :: sky_temp = 263.15
    logical :: sky_cooling_flag = .false.
    integer(i4) :: max_droplets_per_cell = 20
    real(dp) :: rad_call_interval = 0.0
    integer(i4) :: nPhotons = 700000
    integer(i4) :: nBins = 30
    real(dp) :: Lx_rad = 2.0
    real(dp) :: Ly_rad = 2.0
    real(dp) :: T_side = 293.15

    ! --- Public output arrays (written to NetCDF by writeout) ---
    real(dp), allocatable :: rad_F_net(:)
    real(dp), allocatable :: rad_heating_rate(:)

    ! --- Shared physical constants ---
    real(dp), parameter :: SIGMA_SB = 5.670374e-8   ! Stefan-Boltzmann [W m^-2 K^-4]
    ! Fixed air density — valid for Pi-Chamber near STP.
    ! Must be computed from equation of state if radiation is
    ! extended to parcel mode or chamber at lower pressures.
    real(dp), parameter :: RHO_A = 1.2              ! air density [kg m^-3]

    ! --- Mie table (read once at init) ---
    integer(i4), parameter :: NROWS_MIE = 501, NCOLS_MIE = 112
    real(dp), allocatable :: mie_table(:,:)
    real(dp), allocatable :: mie_wavelength(:), mie_delta_lambda(:), mie_radius(:)

    ! --- 3D MC geometry (computed once at init when method='3d') ---
    real(dp), allocatable :: mc_z_edges(:), mc_z_centers(:)
    real(dp) :: mc_A_bottom, mc_A_top, mc_A_sides

    ! --- Interval-based calling accumulators ---
    real(dp), allocatable :: kappa_sum(:), T_sum(:)
    real(dp) :: dt_accumulated = 0.0
    real(dp) :: next_rad_time = 0.0

contains

! =========================================================================
! Shared physics subroutines
! =========================================================================

    ! Planck spectral radiance B(lambda, T) [W m^-3 sr^-1] for an array of wavelengths.
    subroutine planck_lambda(n_wave, wavelength_m, temp, B_lambda)
        integer(i4), intent(in) :: n_wave
        real(dp), intent(in) :: wavelength_m(n_wave), temp
        real(dp), intent(out) :: B_lambda(n_wave)

        real(dp), parameter :: H_PLANCK = 6.62607015e-34
        real(dp), parameter :: C_LIGHT = 2.99792458e8
        real(dp), parameter :: K_BOLTZ = 1.380649e-23

        B_lambda = (2.0 * H_PLANCK * C_LIGHT**2 / wavelength_m**5) / &
                   (exp(H_PLANCK * C_LIGHT / (wavelength_m * K_BOLTZ * temp)) - 1.0)

    end subroutine planck_lambda


    ! Planck-mean absorption coefficient profile from droplet radii.
    ! Uses module-level Mie table arrays (loaded once at init).
    subroutine compute_kappa_prof(nrows, ncols_max, rad_box, nums, dv, T_profile, kappa_prof)
        integer(i4), intent(in) :: nrows, ncols_max
        real(dp), intent(in) :: rad_box(nrows, ncols_max)
        integer(i4), intent(in) :: nums(nrows)
        real(dp), intent(in) :: dv(nrows), T_profile(nrows)
        real(dp), intent(out) :: kappa_prof(nrows)

        integer(i4), parameter :: N_WAVE = NCOLS_MIE - 1
        integer(i4), parameter :: N_RADIUS = NROWS_MIE - 1
        real(dp) :: C_abs(N_WAVE), kappa(N_WAVE), kappa_mean
        real(dp) :: B_planck(N_WAVE), qabs(N_WAVE)
        integer :: i, j, idx

        kappa_prof = 0.0

        do i = 1, nrows
            if (nums(i) == 0) cycle

            C_abs = 0.0
            do j = 1, nums(i)
                idx = 0
                do while (idx < N_RADIUS .and. mie_radius(idx+1) <= rad_box(i,j))
                    idx = idx + 1
                end do
                qabs = mie_table(idx+1, 2:N_WAVE+1)
                C_abs = C_abs + pi * rad_box(i,j)**2 * qabs
            end do

            kappa = C_abs / dv(i)
            call planck_lambda(N_WAVE, mie_wavelength, T_profile(i), B_planck)

            kappa_mean = sum(kappa * B_planck * mie_delta_lambda) / &
                         sum(B_planck * mie_delta_lambda)
            kappa_prof(i) = kappa_mean
        end do

    end subroutine compute_kappa_prof


! =========================================================================
! 1D Two-Stream Solver
! =========================================================================

    ! Optical depth profile from absorption coefficient.
    ! Applies diffusivity factor (1.66) to account for cosine-weighted
    ! angular integration in the two-stream approximation.
    subroutine compute_tau(nrows, kappa_prof, dz_arr, tau, kappa_d)
        integer(i4), intent(in) :: nrows
        real(dp), intent(in) :: kappa_prof(nrows), dz_arr(nrows)
        real(dp), intent(out) :: tau(nrows), kappa_d(nrows)
        integer :: i

        real(dp), parameter :: DIFF_FAC = 1.66

        kappa_d = kappa_prof * DIFF_FAC
        tau(1) = 0.0
        do i = 2, nrows
            tau(i) = tau(i-1) + kappa_d(i-1) * dz_arr(i-1)
        end do

    end subroutine compute_tau


    ! Upward media-emitted radiation at z-level z_ind.
    subroutine upward_media_emitted(z_ind, B_emm, kappa_d, tau, dz_arr, sum_val)
        integer(i4), intent(in) :: z_ind
        real(dp), intent(in) :: B_emm(:), kappa_d(:), tau(:), dz_arr(:)
        real(dp), intent(out) :: sum_val
        integer :: i

        sum_val = 0.0
        if (z_ind > 1) then
            do i = 1, z_ind - 1
                sum_val = sum_val + B_emm(i) * kappa_d(i) * dz_arr(i) * &
                          exp(-(tau(z_ind) - tau(i)))
            end do
        end if

    end subroutine upward_media_emitted


    ! Downward media-emitted radiation at z-level z_ind.
    subroutine downward_media_emitted(z_ind, B_emm, kappa_d, tau, dz_arr, sum_val)
        integer(i4), intent(in) :: z_ind
        real(dp), intent(in) :: B_emm(:), kappa_d(:), tau(:), dz_arr(:)
        real(dp), intent(out) :: sum_val
        integer :: i, n_lev

        n_lev = size(B_emm)
        sum_val = 0.0
        if (z_ind <= n_lev) then
            do i = z_ind, n_lev
                sum_val = sum_val + B_emm(i) * kappa_d(i) * dz_arr(i) * &
                          exp(-(tau(i) - tau(z_ind)))
            end do
        end if

    end subroutine downward_media_emitted


    ! Compute upward/downward intensities, net flux, and flux divergence.
    subroutine compute_fluxes(nrows, dz_arr, z_arr, tau, kappa_d, T_profile, &
                              l_eps_top, l_eps_bot, T_bot, T_top, &
                              I_plus, I_minus, F_net, dF_dz)
        integer(i4), intent(in) :: nrows
        real(dp), intent(in) :: dz_arr(nrows), z_arr(nrows), tau(nrows)
        real(dp), intent(in) :: kappa_d(nrows), T_profile(nrows)
        real(dp), intent(in) :: l_eps_bot, l_eps_top, T_bot, T_top
        real(dp), intent(out) :: I_plus(nrows), I_minus(nrows)
        real(dp), intent(out) :: F_net(nrows), dF_dz(nrows)

        real(dp) :: B_emm(nrows)
        real(dp) :: J_t, J_b
        real(dp) :: t1, t2_1, t2_2, t2_3, t_3
        real(dp) :: tmp_up, tmp_down
        integer :: i

        B_emm = SIGMA_SB * T_profile**4

        t1 = l_eps_top * SIGMA_SB * T_top**4
        t2_1 = l_eps_bot * SIGMA_SB * T_bot**4 * exp(-tau(nrows))

        call downward_media_emitted(1, B_emm, kappa_d, tau, dz_arr, tmp_down)
        t2_2 = (1.0 - l_eps_bot) * exp(-tau(nrows)) * tmp_down

        call upward_media_emitted(nrows, B_emm, kappa_d, tau, dz_arr, tmp_up)
        t2_3 = tmp_up

        t_3 = 1.0 - ((1.0 - l_eps_bot) * (1.0 - l_eps_top) * exp(-tau(nrows))**2)
        J_t = (t1 + (1.0 - l_eps_top) * (t2_1 + t2_2 + t2_3)) / t_3

        J_b = l_eps_bot * SIGMA_SB * T_bot**4 + &
              (1.0 - l_eps_bot) * (J_t * exp(-tau(nrows)) + tmp_down)

        I_plus = 0.0
        I_minus = 0.0

        do i = 1, nrows
            call upward_media_emitted(i, B_emm, kappa_d, tau, dz_arr, tmp_up)
            call downward_media_emitted(i, B_emm, kappa_d, tau, dz_arr, tmp_down)
            I_plus(i) = J_b * exp(-tau(i)) + tmp_up
            I_minus(i) = J_t * exp(-(tau(nrows) - tau(i))) + tmp_down
        end do

        F_net = I_plus - I_minus

        dF_dz(1) = (F_net(2) - F_net(1)) / (z_arr(2) - z_arr(1))
        do i = 2, nrows - 1
            dF_dz(i) = (F_net(i+1) - F_net(i-1)) / (z_arr(i+1) - z_arr(i-1))
        end do
        dF_dz(nrows) = (F_net(nrows) - F_net(nrows-1)) / &
                        (z_arr(nrows) - z_arr(nrows-1))

    end subroutine compute_fluxes


    ! Heating rate from flux divergence, accounting for liquid water heat capacity.
    subroutine compute_heating_rate_1d(nrows, ncols_max, nums, rad_box, dv, dF_dz, heating)
        integer(i4), intent(in) :: nrows, ncols_max
        integer(i4), intent(in) :: nums(nrows)
        real(dp), intent(in) :: rad_box(nrows, ncols_max), dv(nrows), dF_dz(nrows)
        real(dp), intent(out) :: heating(nrows)

        real(dp) :: w_d, vol_water
        integer :: i, j

        do i = 1, nrows
            vol_water = 0.0
            do j = 1, nums(i)
                vol_water = vol_water + &
                    pi_43 * rad_box(i,j)**3
            end do
            w_d = vol_water * rho_l / dv(i)
            heating(i) = -dF_dz(i) / (RHO_A * cp + w_d * c_l)
        end do

    end subroutine compute_heating_rate_1d


    ! Full 1D two-stream radiation solve.
    ! Populates module-level rad_F_net and rad_heating_rate.
    subroutine solve_1d(nrows, ncols_max, rad_box, nums, dx, dy, &
                        T_profile, T_bot, T_top, precomputed_kappa)
        integer(i4), intent(in) :: nrows, ncols_max
        real(dp), intent(in) :: rad_box(nrows, ncols_max)
        integer(i4), intent(in) :: nums(nrows)
        real(dp), intent(in) :: dx, dy
        real(dp), intent(in) :: T_profile(nrows), T_bot, T_top
        real(dp), intent(in), optional :: precomputed_kappa(:)

        real(dp) :: dz_arr(nrows), dv(nrows)
        real(dp) :: kappa_prof(nrows), kappa_d(nrows), tau(nrows)
        real(dp) :: I_plus(nrows), I_minus(nrows), dF_dz(nrows)

        call compute_dz_dv(dx, dy, dz_arr, dv)

        if (present(precomputed_kappa)) then
            kappa_prof = precomputed_kappa
        else
            call compute_kappa_prof(nrows, ncols_max, rad_box, nums, &
                                    dv, T_profile, kappa_prof)
        end if

        call compute_tau(nrows, kappa_prof, dz_arr, tau, kappa_d)
        call compute_fluxes(nrows, dz_arr, z, tau, kappa_d, T_profile, &
                            eps_top, eps_bot, T_bot, T_top, &
                            I_plus, I_minus, rad_F_net, dF_dz)
        call compute_heating_rate_1d(nrows, ncols_max, nums, rad_box, &
                                     dv, dF_dz, rad_heating_rate)

    end subroutine solve_1d


! =========================================================================
! 3D Monte Carlo Solver
! =========================================================================

    ! Linear interpolation from one 1D grid to another.
    subroutine interp1_linear(x_in, y_in, x_out, y_out, n_in, n_out)
        integer(i4), intent(in) :: n_in, n_out
        real(dp), intent(in) :: x_in(n_in), y_in(n_in), x_out(n_out)
        real(dp), intent(out) :: y_out(n_out)
        integer :: i, j
        real(dp) :: t_frac

        do j = 1, n_out
            if (x_out(j) <= x_in(1)) then
                y_out(j) = y_in(1); cycle
            end if
            if (x_out(j) >= x_in(n_in)) then
                y_out(j) = y_in(n_in); cycle
            end if
            do i = 1, n_in - 1
                if (x_out(j) >= x_in(i) .and. x_out(j) <= x_in(i+1)) then
                    t_frac = (x_out(j) - x_in(i)) / (x_in(i+1) - x_in(i))
                    y_out(j) = (1.0 - t_frac) * y_in(i) + t_frac * y_in(i+1)
                    exit
                end if
            end do
        end do

    end subroutine interp1_linear


    ! Wall photon emission probabilities for the 3D MC domain.
    ! Uses precomputed boundary areas (mc_A_bottom, mc_A_top, mc_A_sides).
    subroutine wall_photon_probabilities(T_bottom, T_top_mc, lT_side, &
                                         P_bottom, P_top_mc, P_sides, Q_total)
        real(dp), intent(in) :: T_bottom, T_top_mc, lT_side
        real(dp), intent(out) :: P_bottom, P_top_mc, P_sides, Q_total

        real(dp) :: Q_bottom, Q_top, Q_sides

        Q_bottom = SIGMA_SB * T_bottom**4 * mc_A_bottom
        Q_top = SIGMA_SB * T_top_mc**4 * mc_A_top
        Q_sides = SIGMA_SB * lT_side**4 * mc_A_sides

        Q_total = Q_bottom + Q_top + Q_sides
        P_bottom = Q_bottom / Q_total
        P_top_mc = Q_top / Q_total
        P_sides = Q_sides / Q_total

    end subroutine wall_photon_probabilities


    ! Sample emission point and outward normal on a wall.
    subroutine sample_point_and_normal(P_bottom, P_top_mc, P_sides, pt, norm_vec)
        real(dp), intent(in) :: P_bottom, P_top_mc, P_sides
        real(dp), intent(out) :: pt(3), norm_vec(3)
        real(dp) :: r, rx, ry, rz

        call random_number(r)

        if (r < P_bottom) then
            call random_number(rx); call random_number(ry)
            pt = (/ Lx_rad*rx, Ly_rad*ry, 0.0 /)
            norm_vec = (/ 0.0, 0.0, 1.0 /)
            return
        end if
        r = r - P_bottom

        if (r < P_top_mc) then
            call random_number(rx); call random_number(ry)
            pt = (/ Lx_rad*rx, Ly_rad*ry, H /)
            norm_vec = (/ 0.0, 0.0, -1.0 /)
            return
        end if
        r = (r - P_top_mc) / P_sides

        if (r < 0.25) then
            call random_number(ry); call random_number(rz)
            pt = (/ 0.0, Ly_rad*ry, H*rz /)
            norm_vec = (/ 1.0, 0.0, 0.0 /)
        else if (r < 0.50) then
            call random_number(ry); call random_number(rz)
            pt = (/ Lx_rad, Ly_rad*ry, H*rz /)
            norm_vec = (/ -1.0, 0.0, 0.0 /)
        else if (r < 0.75) then
            call random_number(rx); call random_number(rz)
            pt = (/ Lx_rad*rx, 0.0, H*rz /)
            norm_vec = (/ 0.0, 1.0, 0.0 /)
        else
            call random_number(rx); call random_number(rz)
            pt = (/ Lx_rad*rx, Ly_rad, H*rz /)
            norm_vec = (/ 0.0, -1.0, 0.0 /)
        end if

    end subroutine sample_point_and_normal


    ! Sample a Lambertian (cosine-weighted) direction from an axis-aligned normal.
    subroutine sample_lambertian_direction(norm_vec, d)
        real(dp), intent(in) :: norm_vec(3)
        real(dp), intent(out) :: d(3)
        real(dp) :: u, v, w, nx, ny, nz
        logical :: accepted

        nx = norm_vec(1); ny = norm_vec(2); nz = norm_vec(3)

        accepted = .false.
        do while (.not. accepted)
            call random_number(u); u = 2.0*u - 1.0
            call random_number(v); v = 2.0*v - 1.0
            accepted = (u*u + v*v <= 1.0)
        end do
        w = sqrt(1.0 - u*u - v*v)

        if (nz == 1.0) then
            d = (/ u, v, w /)
        else if (nz == -1.0) then
            d = (/ u, v, -w /)
        else if (nx == 1.0) then
            d = (/ w, u, v /)
        else if (nx == -1.0) then
            d = (/ -w, u, v /)
        else if (ny == 1.0) then
            d = (/ u, w, v /)
        else if (ny == -1.0) then
            d = (/ u, -w, v /)
        else
            d = (/ u, v, w /)
        end if

    end subroutine sample_lambertian_direction


    ! Compute distance from point pt along direction d to exit the box.
    subroutine compute_exit_distance(pt, d, t_exit)
        real(dp), intent(in) :: pt(3), d(3)
        real(dp), intent(out) :: t_exit
        real(dp) :: ts(3)
        real(dp) :: dx, dy, dz, px, py, pz
        integer :: n_ts

        px = pt(1); py = pt(2); pz = pt(3)
        dx = d(1); dy = d(2); dz = d(3)
        n_ts = 0

        if (dx > 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = (Lx_rad - px) / dx
        else if (dx < 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = -px / dx
        end if
        if (dy > 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = (Ly_rad - py) / dy
        else if (dy < 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = -py / dy
        end if
        if (dz > 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = (H - pz) / dz
        else if (dz < 0.0) then
            n_ts = n_ts + 1; ts(n_ts) = -pz / dz
        end if

        t_exit = minval(ts(1:n_ts))

    end subroutine compute_exit_distance


    ! Sort a small real array in ascending order (insertion sort).
    subroutine sort_real_array(a)
        real(dp), intent(inout) :: a(:)
        integer :: i, j
        real(dp) :: temp

        do i = 2, size(a)
            temp = a(i)
            j = i - 1
            do while (j >= 1 .and. a(j) > temp)
                a(j+1) = a(j)
                j = j - 1
            end do
            a(j+1) = temp
        end do

    end subroutine sort_real_array


    ! Monte Carlo photon path-length accumulation into vertical bins.
    subroutine mc_pathlength(z_edges, P_bottom, P_top_mc, P_sides, &
                             pathLen_sum, L_mean_MC)
        real(dp), intent(in) :: P_bottom, P_top_mc, P_sides
        real(dp), intent(in) :: z_edges(nBins+1)
        real(dp), intent(out) :: pathLen_sum(nBins)
        real(dp), intent(out) :: L_mean_MC

        real(dp) :: total_L, t_exit
        real(dp) :: p0(3), n_wall(3), d(3)
        real(dp) :: dz, z0, z_a, z_b, z_min, z_max
        real(dp), allocatable :: z_edges_in(:), t_list(:)
        real(dp) :: t1, t2, L_sub, z_mid, ze
        integer :: n_photon, j, jj, m, count_edges

        pathLen_sum = 0.0
        total_L = 0.0

        do n_photon = 1, nPhotons
            call sample_point_and_normal(P_bottom, P_top_mc, P_sides, p0, n_wall)
            call sample_lambertian_direction(n_wall, d)
            call compute_exit_distance(p0, d, t_exit)
            total_L = total_L + t_exit

            dz = d(3)

            ! Nearly horizontal rays
            if (abs(dz) < 1.0e-8) then
                z0 = p0(3)
                if (z0 >= 0.0 .and. z0 <= H .and. t_exit > 0.0) then
                    do j = 1, nBins
                        if (z0 >= z_edges(j) .and. z0 <= z_edges(j+1)) then
                            pathLen_sum(j) = pathLen_sum(j) + t_exit
                            exit
                        end if
                    end do
                end if
                cycle
            end if

            ! General case
            z0 = p0(3)
            z_a = z0
            z_b = z0 + dz * t_exit

            if ((z_a < 0.0 .and. z_b < 0.0) .or. &
                (z_a > H .and. z_b > H)) cycle

            z_min = max(min(z_a, z_b), 0.0)
            z_max = min(max(z_a, z_b), H)
            if (z_max <= z_min) cycle

            count_edges = 0
            do j = 1, nBins + 1
                if (z_edges(j) > z_min .and. z_edges(j) < z_max) &
                    count_edges = count_edges + 1
            end do

            allocate(z_edges_in(count_edges))
            m = 0
            do j = 1, nBins + 1
                if (z_edges(j) > z_min .and. z_edges(j) < z_max) then
                    m = m + 1
                    z_edges_in(m) = z_edges(j)
                end if
            end do

            allocate(t_list(count_edges + 2))
            t_list(1) = 0.0
            t_list(2) = t_exit
            do j = 1, count_edges
                ze = z_edges_in(j)
                t_list(j+2) = (ze - z0) / dz
            end do
            call sort_real_array(t_list)

            do jj = 1, size(t_list) - 1
                t1 = t_list(jj)
                t2 = t_list(jj+1)
                L_sub = t2 - t1
                if (L_sub <= 0.0) cycle

                z_mid = z0 + dz * 0.5 * (t1 + t2)
                if (z_mid < 0.0 .or. z_mid > H) cycle

                do j = 1, nBins
                    if (z_mid >= z_edges(j) .and. z_mid < z_edges(j+1)) then
                        pathLen_sum(j) = pathLen_sum(j) + L_sub
                        exit
                    end if
                end do
            end do

            deallocate(z_edges_in, t_list)
        end do

        L_mean_MC = total_L / nPhotons

    end subroutine mc_pathlength


    ! Compute radiative heating from MC path lengths and absorption.
    subroutine compute_mc_heating(z_edges, pathLen_sum, Q_total, &
                                  kappa_bins, T_bins, dTdt)
        real(dp), intent(in) :: z_edges(nBins+1)
        real(dp), intent(in) :: pathLen_sum(nBins)
        real(dp), intent(in) :: kappa_bins(nBins), T_bins(nBins)
        real(dp), intent(in) :: Q_total
        real(dp), intent(out) :: dTdt(nBins)

        real(dp) :: power_per_photon, dz_bin
        real(dp) :: Q_abs_bins(nBins), q_abs(nBins), q_emit(nBins), q_net(nBins)

        power_per_photon = Q_total / real(nPhotons, dp)
        Q_abs_bins = kappa_bins * power_per_photon * pathLen_sum
        dz_bin = z_edges(2) - z_edges(1)
        q_abs = Q_abs_bins / (mc_A_bottom * dz_bin)
        q_emit = 4.0 * kappa_bins * SIGMA_SB * T_bins**4
        q_net = q_abs - q_emit
        dTdt = q_net / (RHO_A * cp)

    end subroutine compute_mc_heating


    ! Full 3D Monte Carlo radiation solve.
    ! Populates module-level rad_heating_rate. rad_F_net is zeroed
    ! (not directly computed by the MC method).
    subroutine solve_3d(nrows, T_profile, T_bot, T_top, kappa_prof)
        integer(i4), intent(in) :: nrows
        real(dp), intent(in) :: T_profile(nrows), T_bot, T_top
        real(dp), intent(in) :: kappa_prof(nrows)

        real(dp) :: T_bins(nBins), kappa_bins(nBins)
        real(dp) :: P_bottom, P_top_mc, P_sides, Q_total
        real(dp) :: pathLen_sum(nBins), dTdt(nBins), L_mean_MC

        call interp1_linear(z, T_profile, mc_z_centers, T_bins, nrows, nBins)
        call interp1_linear(z, kappa_prof, mc_z_centers, kappa_bins, nrows, nBins)

        call wall_photon_probabilities(T_bot, T_top, T_side, &
                                       P_bottom, P_top_mc, P_sides, Q_total)
        call mc_pathlength(mc_z_edges, P_bottom, P_top_mc, P_sides, &
                           pathLen_sum, L_mean_MC)
        call compute_mc_heating(mc_z_edges, pathLen_sum, Q_total, &
                                kappa_bins, T_bins, dTdt)

        call interp1_linear(mc_z_centers, dTdt, z, rad_heating_rate, nBins, nrows)

        rad_F_net = 0.0

    end subroutine solve_3d




! =========================================================================
! Public interface routines
! =========================================================================

    ! Read RADIATION namelist, load Mie table, allocate output arrays.
    subroutine initialize_radiation()
        integer :: ierr, nml_unit, k
        character(256) :: resolved_path

        namelist /RADIATION/ radiation_method, mie_data_file, eps_top, eps_bot, &
            sky_temp, sky_cooling_flag, max_droplets_per_cell, rad_call_interval, &
            nPhotons, nBins, Lx_rad, Ly_rad, T_side

        open(newunit=nml_unit, file=namelist_path, action='read', status='old', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: cannot open namelist for RADIATION'
            stop 1
        end if
        read(nml=RADIATION, unit=nml_unit, iostat=ierr)
        close(nml_unit)

        if (mie_data_file == '') then
            write(0,*) 'Error: mie_data_file must be set when do_radiation = .true.'
            stop 1
        end if

        resolved_path = resolve_path(namelist_dir, mie_data_file)
        call load_mie_table(resolved_path)

        allocate(rad_F_net(N))
        allocate(rad_heating_rate(N))
        rad_F_net = 0.0
        rad_heating_rate = 0.0

        if (radiation_method == '3d') then
            allocate(mc_z_edges(nBins+1), mc_z_centers(nBins))
            do k = 1, nBins + 1
                mc_z_edges(k) = H * real(k-1, dp) / real(nBins, dp)
            end do
            do k = 1, nBins
                mc_z_centers(k) = 0.5 * (mc_z_edges(k) + mc_z_edges(k+1))
            end do
            mc_A_bottom = Lx_rad * Ly_rad
            mc_A_top = Lx_rad * Ly_rad
            mc_A_sides = 2.0 * (Lx_rad + Ly_rad) * H
        end if

        if (rad_call_interval > 0.0) then
            allocate(kappa_sum(N), T_sum(N))
            kappa_sum = 0.0
            T_sum = 0.0
            dt_accumulated = 0.0
            next_rad_time = rad_call_interval
        end if

    end subroutine initialize_radiation


    ! Main radiation driver — called from the time loop.
    !
    ! Bins droplets into grid cells, computes absorption, dispatches to
    ! the selected solver, and applies the heating rate to Tarr.
    ! When rad_call_interval > 0, accumulates kappa and T between calls
    ! and only fires the solver when enough time has elapsed.
    subroutine compute_radiation(Tarr, delta_t, current_time)
        real(dp), intent(inout) :: Tarr(:)
        real(dp), intent(in) :: delta_t, current_time

        real(dp) :: rad_box(N, max_droplets_per_cell)
        integer(i4) :: nums(N)
        real(dp) :: dx, T_bot, T_top
        real(dp) :: kappa_prof(N), dv(N), dz_arr(N)
        real(dp) :: kappa_avg(N), T_avg(N)
        integer :: i

        dx = sqrt(gridcell_volume / dz_length)

        call build_rad_box(rad_box, nums)

        if (rad_call_interval > 0.0) then
            call compute_dz_dv(dx, dx, dz_arr, dv)
            call compute_kappa_prof(N, max_droplets_per_cell, rad_box, nums, &
                                    dv, Tarr, kappa_prof)
            kappa_sum = kappa_sum + kappa_prof * delta_t
            T_sum = T_sum + Tarr * delta_t
            dt_accumulated = dt_accumulated + delta_t

            if (current_time >= next_rad_time) then
                kappa_avg = kappa_sum / dt_accumulated
                T_avg = T_sum / dt_accumulated

                call get_boundary_temps(T_avg, T_bot, T_top)

                if (radiation_method == '1d') then
                    call solve_1d(N, max_droplets_per_cell, rad_box, nums, &
                                  dx, dx, T_avg, T_bot, T_top, kappa_avg)
                else
                    call solve_3d(N, T_avg, T_bot, T_top, kappa_avg)
                end if

                do i = 1, N
                    Tarr(i) = Tarr(i) + rad_heating_rate(i) * dt_accumulated
                end do
                budget_radiation_delta_T = budget_radiation_delta_T + &
                    sum(rad_heating_rate) * dt_accumulated

                kappa_sum = 0.0
                T_sum = 0.0
                dt_accumulated = 0.0
                next_rad_time = next_rad_time + rad_call_interval
            end if
        else
            call get_boundary_temps(Tarr, T_bot, T_top)

            if (radiation_method == '1d') then
                call solve_1d(N, max_droplets_per_cell, rad_box, nums, &
                              dx, dx, Tarr, T_bot, T_top)
            else
                call compute_dz_dv(dx, dx, dz_arr, dv)
                call compute_kappa_prof(N, max_droplets_per_cell, rad_box, nums, &
                                        dv, Tarr, kappa_prof)
                call solve_3d(N, Tarr, T_bot, T_top, kappa_prof)
            end if

            do i = 1, N
                Tarr(i) = Tarr(i) + rad_heating_rate(i) * delta_t
            end do
            budget_radiation_delta_T = budget_radiation_delta_T + &
                sum(rad_heating_rate) * delta_t
        end if

    end subroutine compute_radiation


    ! Deallocate all radiation arrays.
    subroutine finalize_radiation()

        if (allocated(rad_F_net)) deallocate(rad_F_net)
        if (allocated(rad_heating_rate)) deallocate(rad_heating_rate)
        if (allocated(mie_table)) deallocate(mie_table)
        if (allocated(mie_wavelength)) deallocate(mie_wavelength)
        if (allocated(mie_delta_lambda)) deallocate(mie_delta_lambda)
        if (allocated(mie_radius)) deallocate(mie_radius)
        if (allocated(kappa_sum)) deallocate(kappa_sum)
        if (allocated(T_sum)) deallocate(T_sum)

    end subroutine finalize_radiation


! =========================================================================
! Internal helpers
! =========================================================================

    ! Set boundary temperatures from field endpoints, with sky cooling override.
    subroutine get_boundary_temps(Tarr, T_bot, T_top)
        real(dp), intent(in) :: Tarr(:)
        real(dp), intent(out) :: T_bot, T_top

        T_bot = Tarr(1)
        if (sky_cooling_flag) then
            T_top = sky_temp
        else
            T_top = Tarr(N)
        end if

    end subroutine get_boundary_temps


    ! Bin activated droplet radii into grid cells for radiation.
    subroutine build_rad_box(rad_box, nums)
        real(dp), intent(out) :: rad_box(:,:)
        integer(i4), intent(out) :: nums(:)
        integer :: i, gc, slot

        rad_box = 0.0
        nums = 0

        do i = 1, current_n_particles
            if (.not. particles(i)%activated) cycle
            gc = particles(i)%gridcell
            if (gc < 1 .or. gc > N) cycle
            if (nums(gc) >= max_droplets_per_cell) cycle
            nums(gc) = nums(gc) + 1
            slot = nums(gc)
            rad_box(gc, slot) = particles(i)%radius
        end do

    end subroutine build_rad_box


    ! Load Mie Q_abs table from text file (called once at init).
    subroutine load_mie_table(filepath)
        character(*), intent(in) :: filepath
        integer :: i, funit, ierr
        character(1) :: peek

        allocate(mie_table(NROWS_MIE, NCOLS_MIE))
        allocate(mie_wavelength(NCOLS_MIE - 1))
        allocate(mie_delta_lambda(NCOLS_MIE - 1))
        allocate(mie_radius(NROWS_MIE - 1))

        open(newunit=funit, file=trim(filepath), status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: cannot open mie_data_file: ', trim(filepath)
            stop 1
        end if
        ! Skip comment lines starting with #
        do
            read(funit, '(a1)', iostat=ierr) peek
            if (ierr /= 0) exit
            if (peek /= '#') then
                backspace(funit)
                exit
            end if
        end do
        do i = 1, NROWS_MIE
            read(funit, *) mie_table(i, :)
        end do
        close(funit)

        mie_wavelength = mie_table(1, 2:NCOLS_MIE)
        mie_radius = mie_table(2:NROWS_MIE, 1)

        mie_delta_lambda(1) = mie_wavelength(2) - mie_wavelength(1)
        do i = 2, NCOLS_MIE - 2
            mie_delta_lambda(i) = (mie_wavelength(i+1) - mie_wavelength(i-1)) / 2.0
        end do
        mie_delta_lambda(NCOLS_MIE - 1) = &
            mie_wavelength(NCOLS_MIE - 1) - mie_wavelength(NCOLS_MIE - 2)

    end subroutine load_mie_table


    ! Compute dz and dv arrays from the global z grid.
    subroutine compute_dz_dv(dx, dy, dz_arr, dv)
        real(dp), intent(in) :: dx, dy
        real(dp), intent(out) :: dz_arr(:), dv(:)
        integer :: i

        dz_arr(1) = z(2) - z(1)
        do i = 2, N
            dz_arr(i) = z(i) - z(i-1)
        end do
        dv = dx * dy * dz_arr

    end subroutine compute_dz_dv

end module radiation