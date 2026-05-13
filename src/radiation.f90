!> Longwave radiative heating for chamber-mode CODT.
!>
!> Two solvers available via radiation_method namelist parameter:
!>   '1d' — 1D two-stream (diffusivity-factor approximation)
!>   '3d' — 3D Monte Carlo photon path-length method
!>
!> Ported from Suryadev Singh's radiation module (Oct 2025).
module radiation
    use globals, only: dp, i4, N, H, z, T, gridcell_volume, dz_length, &
                       do_radiation, budget_radiation_delta_T, &
                       nc_verify, resolve_path, namelist_path, namelist_dir
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
    integer(i4) :: sky_cooling_flag = 0
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

    ! --- Mie table (read once at init) ---
    integer(i4), parameter :: NROWS_MIE = 501, NCOLS_MIE = 112
    real(dp), allocatable :: mie_table(:,:)
    real(dp), allocatable :: mie_wavelength(:), mie_delta_lambda(:), mie_radius(:)

    ! --- Interval-based calling accumulators ---
    real(dp), allocatable :: kappa_sum(:), T_sum(:)
    real(dp) :: dt_accumulated = 0.0
    real(dp) :: next_rad_time = 0.0

contains
    subroutine planck_lambda(wavelength_m, T, B_lambda)
    	implicit none
   	double precision, intent(in) :: wavelength_m, T
    	double precision, intent(out) :: B_lambda
    	double precision, parameter :: h = 6.62607015d-34
    	double precision, parameter :: c = 2.99792458d8
    	double precision, parameter :: k = 1.380649d-23
    	double precision :: exponent

    	exponent = h * c / (wavelength_m * k * T)
    	B_lambda = (2.d0 * h * c**2 / wavelength_m**5) / (dexp(exponent) - 1.d0)

    end subroutine planck_lambda 
    
    
    subroutine volume_calculation(z, dx, dy, dz, dv)
       implicit none
       double precision, intent(in) :: z(:)
       double precision, intent(in) :: dx, dy
       double precision, intent(out) :: dz(:), dv(:)
       integer :: i, n
       
       n = size(z)
       dz(1) = z(2) - z(1)
       do i = 2, n
          dz(i) = z(i) - z(i-1)
       end do
       dv = dx*dy*dz
   end subroutine volume_calculation
   
   
   
 !!!! reading mie file   
subroutine read_mie_data_fixed(nrows_mie, ncols_mie, mie_file, wavelength, delta_lambda, Radius)
    implicit none
    integer, intent(in) :: nrows_mie, ncols_mie
    double precision, intent(out) :: wavelength(ncols_mie-1), delta_lambda(ncols_mie-1), Radius(nrows_mie-1)

    double precision, intent(out) :: mie_file(nrows_mie, ncols_mie)
    integer :: i

    open(unit=10, file="../../../mie_qabs_wavelength_vs_radius.txt",&
     status="old", action="read")
    !! change with relative path.  
    !! if you want to update this dataset, see mie_scattering_data_file.ipynb in the folder and run accordingly.
    do i = 1, nrows_mie
        read(10, *) mie_file(i, :)
    end do
    close(10)

    wavelength = mie_file(1, 2:ncols_mie)

    delta_lambda = 0.0d0
    do i = 2, ncols_mie-2
        delta_lambda(i) = (wavelength(i+1) - wavelength(i-1)) / 2.0d0
    end do
    delta_lambda(1) = wavelength(2) - wavelength(1)
    delta_lambda(ncols_mie-1) = wavelength(ncols_mie-1) - wavelength(ncols_mie-2)

    Radius = mie_file(2:nrows_mie, 1)

end subroutine read_mie_data_fixed

subroutine compute_kappa_prof(nrows, ncols_max, nrows_mie, ncols_mie, rad_box, nums, dv, T_profile, &
                              Radius, wavelength, delta_lambda, mie_file, kappa_prof)
    implicit none
    integer, intent(in) :: nrows, ncols_max
    integer, intent(in) :: nrows_mie, ncols_mie
    double precision, intent(in) :: rad_box(nrows, ncols_max)
    integer, intent(in) :: nums(nrows)
    double precision, intent(in) :: dv(nrows)
    double precision, intent(in) :: T_profile(nrows)
    double precision, intent(in) :: Radius(nrows_mie-1)
    double precision, intent(in) :: wavelength(ncols_mie-1)
    double precision, intent(in) :: delta_lambda(ncols_mie-1)
    double precision, intent(in) :: mie_file(nrows_mie, ncols_mie)
    double precision, intent(out) :: kappa_prof(nrows)

    double precision :: C_abs(ncols_mie-1), kappa(ncols_mie-1), kappa_mean
    double precision :: B(ncols_mie-1), qabs(ncols_mie-1)
    integer :: i, j, k, idx
    integer :: n_wave, n_radius

    n_wave = ncols_mie - 1
    n_radius = nrows_mie - 1

    kappa_prof = 0.0d0

    do i = 1, nrows
        if (nums(i) > 0) then
            C_abs = 0.0d0
            do j = 1, nums(i)
                idx = 0
                do while (idx < n_radius .and. Radius(idx+1) <= rad_box(i,j))
                    idx = idx + 1
                end do
                qabs = mie_file(idx+1, 2:n_wave+1)
                C_abs = C_abs + 3.141592653589793d0 * rad_box(i,j)**2 * qabs
            end do

            
            kappa = C_abs / dv(i)

            do k = 1, n_wave
                call planck_lambda(wavelength(k), T_profile(i), B(k))
            end do

            ! Planck mean
            kappa_mean = sum(kappa * B * delta_lambda) / sum(B * delta_lambda)
            kappa_prof(i) = kappa_mean
        else
            kappa_prof(i) = 0.0d0
        end if
    end do

end subroutine compute_kappa_prof


!!! for 3-D model -----------------------------------------------------------------------------------------------------------------
! coarse interpolation
subroutine interp_column(n_in, z_in, var_in, n_out, z_out, var_out)
    implicit none
    integer, intent(in) :: n_in       ! number of input points
    double precision, intent(in) :: z_in(n_in), var_in(n_in)
    integer, intent(in) :: n_out      ! number of output points
    double precision, intent(in) :: z_out(n_out)
    double precision, intent(out) :: var_out(n_out,1)   ! column vector

    integer :: i, j
    real :: t

    j = 1
    do i = 1, n_out
        ! find the interval in input array
        do while (j < n_in .and. z_out(i) > z_in(j+1))
            j = j + 1
        end do

        ! linear interpolation
        if (j == n_in) then
            var_out(i,1) = var_in(n_in)
        else
            t = (z_out(i) - z_in(j)) / (z_in(j+1) - z_in(j))
            var_out(i,1) = (1.0 - t)*var_in(j) + t*var_in(j+1)
        end if
    end do

end subroutine interp_column


!! computes photon probability from each wall ............................................
subroutine wall_photon_probabilities(Lx, Ly, H_cloud, T_bottom, T_top, T_side, &
                                      P_bottom, P_top, P_sides, Q_total)
    implicit none
    double precision, intent(in) :: Lx, Ly, H_cloud
    double precision, intent(in) :: T_bottom, T_top, T_side
    double precision, intent(out) :: P_bottom, P_top, P_sides

    double precision :: A_bottom, A_top, A_side_x, A_side_y
    double precision :: A_sides_total, A_total
    double precision :: F_bottom, F_top, F_side
    double precision :: Q_bottom, Q_top, Q_sides
    double precision, intent(out) :: Q_total
    double precision, parameter :: sigmaSB = 5.670374e-8  !  [W/m^2/K^4]

    ! Wall areas
    A_bottom = Lx * Ly
    A_top    = Lx * Ly
    A_side_x = Ly * H_cloud         ! x = 0 or x = Lx
    A_side_y = Lx * H_cloud         ! y = 0 or y = Ly
    A_sides_total = 2.0d0 * A_side_x + 2.0d0 * A_side_y
    A_total = A_bottom + A_top + A_sides_total  ! not used further

    ! Wall fluxes (hemispheric)
    F_bottom = sigmaSB * T_bottom**4
    F_top    = sigmaSB * T_top**4
    F_side   = sigmaSB * T_side**4

    ! Emitted power from each wall
    Q_bottom = F_bottom * A_bottom
    Q_top    = F_top    * A_top
    Q_sides  = F_side   * A_sides_total

    Q_total = Q_bottom + Q_top + Q_sides

    ! Probabilities for picking emitting wall
    P_bottom = Q_bottom / Q_total
    P_top    = Q_top    / Q_total
    P_sides  = Q_sides  / Q_total   ! = 1 - P_bottom - P_top

end subroutine wall_photon_probabilities

!! sample point and normal calculation on wall
subroutine sample_point_and_normal(Lx, Ly, H_cloud, P_bottom, P_top, P_sides, p, n)
    implicit none
    ! Inputs
    double precision, intent(in) :: Lx, Ly, H_cloud
    double precision, intent(in) :: P_bottom, P_top, P_sides
    ! Outputs
    double precision, intent(out) :: p(3), n(3)

    ! Local variables
    double precision :: r, r_side
    double precision :: rx, ry, rz

    ! Sample uniform random number
    call random_number(r)

    ! --- Bottom wall ---
    if (r < P_bottom) then
        call random_number(rx)
        call random_number(ry)
        p = (/ Lx*rx, Ly*ry, 0.0d0 /)
        n = (/ 0.0d0, 0.0d0, 1.0d0 /)
        return
    end if
    r = r - P_bottom

    ! --- Top wall ---
    if (r < P_top) then
        call random_number(rx)
        call random_number(ry)
        p = (/ Lx*rx, Ly*ry, H_cloud /)
        n = (/ 0.0d0, 0.0d0, -1.0d0 /)
        return
    end if
    r = (r - P_top) / P_sides   ! normalize for sides

    ! --- Side walls ---
    if (r < 0.25d0) then
        ! x = 0, normal +x
        call random_number(ry)
        call random_number(rz)
        p = (/ 0.0d0, Ly*ry, H_cloud*rz /)
        n = (/ 1.0d0, 0.0d0, 0.0d0 /)
    else if (r < 0.50d0) then
        ! x = Lx, normal -x
        call random_number(ry)
        call random_number(rz)
        p = (/ Lx, Ly*ry, H_cloud*rz /)
        n = (/ -1.0d0, 0.0d0, 0.0d0 /)
    else if (r < 0.75d0) then
        ! y = 0, normal +y
        call random_number(rx)
        call random_number(rz)
        p = (/ Lx*rx, 0.0d0, H_cloud*rz /)
        n = (/ 0.0d0, 1.0d0, 0.0d0 /)
    else
        ! y = Ly, normal -y
        call random_number(rx)
        call random_number(rz)
        p = (/ Lx*rx, Ly, H_cloud*rz /)
        n = (/ 0.0d0, -1.0d0, 0.0d0 /)
    end if

end subroutine sample_point_and_normal

!! sampling lambertian direction-distribution.
subroutine sample_lambertian_direction(n, d)
    implicit none
    ! Input: unit normal vector along ±x, ±y, or ±z
    double precision, intent(in)  :: n(3)
    ! Output: sampled cosine-weighted direction
    double precision, intent(out) :: d(3)

    ! Local variables
    double precision :: u, v, w
    double precision :: nx, ny, nz
    double precision :: norm_d
    integer :: ok

    nx = n(1)
    ny = n(2)
    nz = n(3)

    ! --- Sample (u,v) inside unit disk ---
    ok = 0
    do while (ok == 0)
        call random_number(u)
        u = 2.0d0*u - 1.0d0
        call random_number(v)
        v = 2.0d0*v - 1.0d0
        if (u*u + v*v <= 1.0d0) ok = 1
    end do

    w = sqrt(max(0.0d0, 1.0d0 - u*u - v*v))

    ! --- Map to hemisphere around axis-aligned normal ---
    if (nx == 0.0d0 .and. ny == 0.0d0 .and. nz == 1.0d0) then
        ! +z
        d = (/ u, v, w /)
    else if (nx == 0.0d0 .and. ny == 0.0d0 .and. nz == -1.0d0) then
        ! -z
        d = (/ u, v, -w /)
    else if (nx == 1.0d0 .and. ny == 0.0d0 .and. nz == 0.0d0) then
        ! +x
        d = (/ w, u, v /)
    else if (nx == -1.0d0 .and. ny == 0.0d0 .and. nz == 0.0d0) then
        ! -x
        d = (/ -w, u, v /)
    else if (nx == 0.0d0 .and. ny == 1.0d0 .and. nz == 0.0d0) then
        ! +y
        d = (/ u, w, v /)
    else if (nx == 0.0d0 .and. ny == -1.0d0 .and. nz == 0.0d0) then
        ! -y
        d = (/ u, -w, v /)
    else
        ! fallback
        d = (/ u, v, w /)
    end if

    norm_d = sqrt(d(1)**2 + d(2)**2 + d(3)**2)
    d = d / norm_d

end subroutine sample_lambertian_direction

!! campute exit distance for a given ray. 
subroutine compute_exit_distance_sub(p, d, Lx, Ly, H_cloud, t_exit)
    implicit none
    ! Inputs
    double precision, intent(in) :: p(3)   ! starting point (x,y,z)
    double precision, intent(in) :: d(3)   ! direction vector (dx,dy,dz)
    double precision, intent(in) :: Lx, Ly, H_cloud
    ! Output
    double precision, intent(out) :: t_exit

    ! Local variables
    double precision :: ts(3)
    integer :: n, i
    double precision :: eps_t
    double precision :: dx, dy, dz, x, y, z
    double precision :: tmin

    x = p(1); y = p(2); z = p(3)
    dx = d(1); dy = d(2); dz = d(3)
    eps_t = 1.0d-9

    n = 0

    ! --- x planes ---
    if (dx > 0.0d0) then
        n = n + 1
        ts(n) = (Lx - x) / dx
    else if (dx < 0.0d0) then
        n = n + 1
        ts(n) = (0.0d0 - x) / dx
    end if

    ! --- y planes ---
    if (dy > 0.0d0) then
        n = n + 1
        ts(n) = (Ly - y) / dy
    else if (dy < 0.0d0) then
        n = n + 1
        ts(n) = (0.0d0 - y) / dy
    end if

    ! --- z planes ---
    if (dz > 0.0d0) then
        n = n + 1
        ts(n) = (H_cloud - z) / dz
    else if (dz < 0.0d0) then
        n = n + 1
        ts(n) = (0.0d0 - z) / dz
    end if

    ! Find the minimum positive t
    tmin = 1.0d300   ! large number
    do i = 1, n
        if (ts(i) > eps_t .and. ts(i) < tmin) tmin = ts(i)
    end do

    t_exit = tmin

end subroutine compute_exit_distance_sub

subroutine sort_real_array(a)
    implicit none
    double precision, intent(inout) :: a(:)
    integer :: i, j
    double precision :: temp

    do i = 1, size(a)-1
        do j = i+1, size(a)
            if (a(j) < a(i)) then
                temp = a(i)
                a(i) = a(j)
                a(j) = temp
            end if
        end do
    end do
end subroutine sort_real_array



subroutine monte_carlo_pathlength_sub(Lx, Ly, H_cloud, z_edges, P_bottom, P_top, P_sides, &
                                     pathLen_sum, L_mean_MC)
    implicit none
    ! Inputs
    double precision, intent(in) :: Lx, Ly, H_cloud
    double precision, intent(in) :: P_bottom, P_top, P_sides
    double precision, intent(in) :: z_edges(nBins+1)
    double precision, intent(out) :: pathLen_sum(nBins)
    double precision, intent(out) :: L_mean_MC
    
    integer :: n, j, jj
    double precision :: total_L, t_exit
    double precision :: p0(3), n_wall(3), d(3)
    double precision :: t_a, t_b, dz
    double precision :: z0, z_a, z_b
    double precision :: z_min, z_max
    double precision, allocatable :: z_edges_in(:)
    logical, allocatable :: mask(:)
    double precision, allocatable :: t_list(:)
    integer :: m, count_edges
    double precision :: t1, t2, L_sub, t_mid, z_mid
    double precision :: ze, t_e
    double precision :: L_mean_theory

    ! Initialize
    pathLen_sum = 0.0d0
    total_L = 0.0d0

    ! ------------------ Monte Carlo loop ------------------
    do n = 1, nPhotons

        ! 1. Pick emitting wall, point, normal
        call sample_point_and_normal(Lx, Ly, H_cloud, P_bottom, P_top, P_sides, p0, n_wall)

        ! 2. Sample Lambertian direction
        call sample_lambertian_direction(n_wall, d)

        ! 3. Distance to exit from the box
        call compute_exit_distance_sub(p0, d, Lx, Ly, H_cloud, t_exit)
        total_L = total_L + t_exit

        ! 4. Segment inside cloud
        t_a = 0.0d0
        t_b = t_exit
        dz = d(3)
        

! === Nearly horizontal rays ===
        if (abs(dz) < 1.0d-8) then
            z0 = p0(3)
            if (z0 >= 0.0d0 .and. z0 <= H_cloud) then
                L_sub = t_b - t_a
                if (L_sub > 0.0d0) then
                    ! locate bin
                    do j = 1, nBins
                        if (z0 >= z_edges(j) .and. z0 <= z_edges(j+1)) then
                            pathLen_sum(j) = pathLen_sum(j) + L_sub
                            exit
                        end if
                    end do
                end if
            end if
            cycle
        end if

        ! === General case ===
        z0 = p0(3)
        z_a = z0 + dz * t_a
        z_b = z0 + dz * t_b

        ! Clip
        if ((z_a < 0.0d0 .and. z_b < 0.0d0) .or. &
            (z_a > H_cloud .and. z_b > H_cloud)) cycle

        if (dz > 0.0d0) then
            z_min = max(min(z_a, z_b), 0.0d0)
            z_max = min(max(z_a, z_b), H_cloud)
        else
            z_min = max(min(z_a, z_b), 0.0d0)
            z_max = min(max(z_a, z_b), H_cloud)
        end if

        if (z_max <= z_min) cycle

        ! z-edges strictly between z_min and z_max
        allocate(mask(nBins+1))
        do j = 1, nBins+1
            mask(j) = (z_edges(j) > z_min) .and. (z_edges(j) < z_max)
        end do

        count_edges = count(mask)
        allocate(z_edges_in(count_edges))
        m = 0
        do j = 1, nBins+1
            if (mask(j)) then
                m = m + 1
                z_edges_in(m) = z_edges(j)
            end if
        end do
        deallocate(mask)

        ! t-list = [t_a, t_b] + converted edges
        allocate(t_list(count_edges+2))
        t_list(1) = t_a
        t_list(2) = t_b

        do j = 1, count_edges
            ze = z_edges_in(j)
            t_e = (ze - z0) / dz
            t_list(j+2) = t_e
        end do

        ! Sort t_list
        call sort_real_array(t_list)   ! you must provide this helper

        ! Accumulate segments
        do jj = 1, size(t_list)-1
            t1 = t_list(jj)
            t2 = t_list(jj+1)
            L_sub = t2 - t1
            if (L_sub <= 0.0d0) cycle

            t_mid = 0.5d0 * (t1 + t2)
            z_mid = z0 + dz * t_mid
            if (z_mid < 0.0d0 .or. z_mid > H_cloud) cycle

            ! bin index: searchsorted(z_edges, z_mid) - 1
            do j = 1, nBins
                if (z_mid >= z_edges(j) .and. z_mid < z_edges(j+1)) then
                    pathLen_sum(j) = pathLen_sum(j) + L_sub
                    exit
                end if
            end do

        end do

        deallocate(z_edges_in)
        deallocate(t_list)

    end do

    L_mean_MC = total_L / nPhotons
    !L_mean_theory = 4.0d0 * Lx * Ly * H_cloud / (2.0d0*Lx*Ly + 4.0d0*Lx*H_cloud)

end subroutine monte_carlo_pathlength_sub






!*************************************************************************************************************************
!! something seems wrong in below subroutine ............................................................................
subroutine compute_radiative_heating(Lx, Ly, z_edges, pathLen_sum, Q_total, &
                                     kappa_prof1, T_g, dTdt)
    implicit none
    ! Inputs
    double precision, intent(in) :: Lx, Ly
    double precision, intent(in) :: z_edges(nBins+1)
    double precision, intent(in) :: pathLen_sum(nBins)
    double precision, intent(in) :: kappa_prof1(nBins), T_g(nBins)
    double precision, intent(in) :: Q_total

    double precision, intent(out) :: dTdt(nBins)

    double precision      ::sigmaSB = 5.670374e-8  ! Stefan-Boltzmann constant [W/m^2/K^4]
    double precision      :: rho=1.2, cp=1007.0 
    double precision :: power_per_photon, dz_bin
    double precision :: Q_abs_bins(nBins), q_abs(nBins), q_emit(nBins), q_net(nBins)
    integer :: i

    power_per_photon = Q_total / dble(nPhotons)
    Q_abs_bins = kappa_prof1 * power_per_photon * pathLen_sum
    dz_bin = z_edges(2) - z_edges(1)
    q_abs = Q_abs_bins / (Lx * Ly * dz_bin)
    q_emit = 4.0d0 * kappa_prof1 * sigmaSB * T_g**4
    q_net = q_abs - q_emit
    dTdt = q_net / (rho * cp)
    !print*, dTdt*3600, "heating rate"

end subroutine compute_radiative_heating


subroutine interp1_linear(x, y, xin, yout, n, m)
    implicit none
    integer, intent(in) :: n, m
    double precision, intent(in)  :: x(n), y(n)     ! original grid
    double precision, intent(in)  :: xin(m)         ! new grid
    double precision, intent(out) :: yout(m)        ! interpolated result

    integer :: i, j

    do j = 1, m

        if (xin(j) <= x(1)) then
            yout(j) = y(1)
            cycle
        end if

        if (xin(j) >= x(n)) then
            yout(j) = y(n)
            cycle
        end if

        ! find interval x(i) <= xin(j) < x(i+1)
        do i = 1, n-1
            if (xin(j) >= x(i) .and. xin(j) <= x(i+1)) then
                yout(j) = y(i) + (y(i+1)-y(i)) *                 &
                        ( (xin(j)-x(i)) / (x(i+1)-x(i)) )
                exit
            end if
        end do

    end do
end subroutine interp1_linear


subroutine init_bins(H_cloud, z_edges, z_centers, nrows, z, T_profile, & 
                    kappa_prof, T_profile1, kappa_prof1)
    implicit none
    integer, intent(in)                    :: nrows
    double precision, intent(in)           :: H_cloud
    double precision, intent(in)           :: z(nrows), T_profile(nrows), kappa_prof(nrows)
    double precision, intent(out)        :: z_centers(nBins), z_edges(nBins+1)
    double precision, intent(out)        :: T_profile1(nBins, 1), kappa_prof1(nBins, 1)
    integer :: i
    
    do i = 1, nBins+1                            ! z bin edges from 0 to H_cloud
        z_edges(i) = H_cloud * real(i-1) / real(nBins)
    end do
    
    do i = 1, nBins                                   ! z bin centers
        z_centers(i) = 0.5 * (z_edges(i) + z_edges(i+1))
    end do
    
    
    call interp_column(nrows, z, T_profile, nBins, z_centers, T_profile1)
    call interp_column(nrows, z, kappa_prof, nBins, z_centers, kappa_prof1)
end subroutine init_bins


subroutine compute_kappa_local(nrows, max_droplets, nrows_mie, ncols_mie, rad_box, nums, &
                                     dx, dy, z, T_profile, kappa_prof)
    implicit none
    integer, intent(in) :: nrows, max_droplets
    integer, intent(in) :: nrows_mie, ncols_mie
    double precision, intent(in) :: rad_box(nrows,max_droplets)
    integer, intent(in) :: nums(nrows)
    double precision, intent(in)  :: dx, dy
    double precision, intent(in)  :: z(nrows), T_profile(nrows)
    double precision, intent(out) :: kappa_prof(nrows) 

    double precision :: dv(nrows), dz(nrows)
    double precision :: Radius(nrows_mie-1), wavelength(ncols_mie-1), delta_lambda(ncols_mie-1)
    double precision :: mie_file(nrows_mie,ncols_mie)

    !! below two function can be called during initialization. -- a cost reduction step.
    call volume_calculation(z, dx, dy, dz, dv)
    call read_mie_data_fixed(nrows_mie, ncols_mie, mie_file, wavelength, delta_lambda, Radius)

    call compute_kappa_prof(nrows, max_droplets, nrows_mie, ncols_mie, rad_box, nums, dv, T_profile, &
                            Radius, wavelength, delta_lambda, mie_file, kappa_prof)   
end subroutine compute_kappa_local



subroutine compute_MC_heating_profile(nrows, z, T_profile, T_bot, T_top, H_cloud, kappa_prof, &
                                      dTdt_intp)
    implicit none
    integer, intent(in) :: nrows
    double precision, intent(in) :: z(nrows), T_profile(nrows)
    double precision, intent(in) :: T_bot, T_top
    double precision, intent(in) :: kappa_prof(nrows)
    double precision, intent(in)               :: H_cloud
    double precision                           :: z_centers(nBins), z_edges(nBins+1)
    
    double precision                           :: T_profile1(nBins, 1), kappa_prof1(nBins, 1)
    double precision                           :: P_bottom, P_top, P_sides, Q_total
    double precision                           :: pathLen_sum(nBins), dTdt(nBins)
    double precision                           :: L_mean_MC
    double precision, intent(out)              :: dTdt_intp(nrows)

    call init_bins(H_cloud, z_edges, z_centers, nrows, z, T_profile, kappa_prof, T_profile1, kappa_prof1)
    call wall_photon_probabilities(Lx, Ly, H_cloud, T_bot, T_top, T_side, &
                                      P_bottom, P_top, P_sides, Q_total)
    call monte_carlo_pathlength_sub(Lx, Ly, H_cloud, z_edges, P_bottom, P_top, P_sides, &
                                     pathLen_sum, L_mean_MC)
                                     
    call compute_radiative_heating(Lx, Ly, z_edges, pathLen_sum, Q_total, kappa_prof1, T_profile1, dTdt)
    call interp1_linear(z_centers, dTdt, z, dTdt_intp, nBins, nrows)      
end subroutine compute_MC_heating_profile




! =========================================================================
! Public interface routines
! =========================================================================

    !> Read RADIATION namelist, load Mie table, allocate output arrays.
    subroutine initialize_radiation()
        integer :: ierr, nml_unit
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

        if (rad_call_interval > 0.0) then
            allocate(kappa_sum(N), T_sum(N))
            kappa_sum = 0.0
            T_sum = 0.0
            dt_accumulated = 0.0
            next_rad_time = rad_call_interval
        end if

    end subroutine initialize_radiation


    !> Main radiation driver — called from the time loop.
    !>
    !> Bins droplets into grid cells, computes absorption, dispatches to
    !> the selected solver, and applies the heating rate to Tarr.
    !> When rad_call_interval > 0, accumulates kappa and T between calls
    !> and only fires the solver when enough time has elapsed.
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


    !> Deallocate all radiation arrays.
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

    !> Set boundary temperatures from field endpoints, with sky cooling override.
    subroutine get_boundary_temps(Tarr, T_bot, T_top)
        real(dp), intent(in) :: Tarr(:)
        real(dp), intent(out) :: T_bot, T_top

        T_bot = Tarr(1)
        if (sky_cooling_flag > 0) then
            T_top = sky_temp
        else
            T_top = Tarr(N)
        end if

    end subroutine get_boundary_temps


    !> Bin activated droplet radii into grid cells for radiation.
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


    !> Load Mie Q_abs table from text file (called once at init).
    subroutine load_mie_table(filepath)
        character(*), intent(in) :: filepath
        integer :: i, funit, ierr

        allocate(mie_table(NROWS_MIE, NCOLS_MIE))
        allocate(mie_wavelength(NCOLS_MIE - 1))
        allocate(mie_delta_lambda(NCOLS_MIE - 1))
        allocate(mie_radius(NROWS_MIE - 1))

        open(newunit=funit, file=trim(filepath), status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: cannot open mie_data_file: ', trim(filepath)
            stop 1
        end if
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


    !> Compute dz and dv arrays from the global z grid.
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