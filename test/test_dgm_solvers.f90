program test_dgm_solvers
  use globals
  use microphysics, only: saturation_vapor_pressure
  use DGM, only: set_aerosol_properties, integrate_ODE, growth_jacobian, growth_rhs
  use ode_integrators, only: ros3_integrate, rkck45_integrate
  implicit none

  integer, parameter :: n_aero = 4
  integer, parameter :: n_ic   = 4
  integer, parameter :: nsteps_fine = 20
  real(dp), parameter :: dt_fine = 0.01
  real(dp), parameter :: dt_coarse = 0.1

  character(len=10) :: aero_labels(n_aero)
  real(dp) :: r_solutes(n_aero)

  character(len=10) :: ic_labels(n_ic)
  real(dp) :: ic_ss(n_ic)
  real(dp) :: ic_r_mult_all(n_ic, n_aero)
  real(dp) :: r_init

  real(dp) :: y_rk(3), y_ros(3), y0(3)
  real(dp) :: qv_sat, es_val, grid_rho, grid_mass, inv_grid_mass
  real(dp) :: tstart, tend, ss_out, ss_env
  real(dp) :: h_rk, h_ros
  integer(i4) :: istat, ia, ic, k

  real(dp), parameter :: perturb_r  = 1.0e-6
  real(dp), parameter :: ss_drop    = -0.70

  integer(i4), parameter :: species = 2  ! (NH4)2SO4
  real(dp), parameter :: rho_solute = 2165.0
  real(dp), parameter :: temp0 = 280.0
  real(dp), parameter :: h_init = 1.0e-3
  real(dp), parameter :: rtol_arr(3) = [1.0e-4, 1.0e-4, 1.0e-4]
  real(dp), parameter :: atol_arr(3) = [1.0e-10, 1.0e-10, 1.0e-10]

  real(dp) :: r_sol, sol_mass

  real(dp), parameter :: vol_scaling = 50.0

  ! Timing variables
  integer(i4), parameter :: n_bench = 10000
  real(dp) :: t_cpu_start, t_cpu_end
  real(dp) :: y_bench(3), h_bench
  integer(i4) :: ib

  ! Grid properties
  grid_rho = pres / (Rd * temp0)
  grid_mass = vol_scaling * domain_width * domain_width * (H / N) * grid_rho
  inv_grid_mass = 1.0 / grid_mass

  ! Aerosol sizes
  aero_labels(1) = "rs5nm";    r_solutes(1) = 5.0e-9
  aero_labels(2) = "rs50nm";   r_solutes(2) = 50.0e-9
  aero_labels(3) = "rs200nm";  r_solutes(3) = 200.0e-9
  aero_labels(4) = "rs1000nm"; r_solutes(4) = 1000.0e-9

  ic_labels(1) = "haze";     ic_ss(1) = -0.005
  ic_labels(2) = "critical"; ic_ss(2) = 0.001
  ic_labels(3) = "growing";  ic_ss(3) = 0.005
  ic_labels(4) = "large";    ic_ss(4) = 0.002

  ic_r_mult_all(1,:) = [2.0,   2.5,   3.0,   3.5]
  ic_r_mult_all(2,:) = [10.0,  12.0,  14.0,  16.0]
  ic_r_mult_all(3,:) = [50.0,  60.0,  70.0,  80.0]
  ic_r_mult_all(4,:) = [200.0, 220.0, 240.0, 260.0]

  ! ============================================================
  ! Part 1: Accuracy comparison (CSV output)
  ! ============================================================
  write(*,'(A)') "case,aerosol,ic,solver,step,time,radius,SS"

  do ia = 1, n_aero
    r_sol = r_solutes(ia)
    sol_mass = pi_43 * r_sol**3 * rho_solute

    do ic = 1, n_ic
      r_init = r_sol * ic_r_mult_all(ic, ia)
      es_val = saturation_vapor_pressure(temp0)
      qv_sat = 0.622 * es_val / (pres - es_val)
      y0(1) = r_init
      y0(2) = qv_sat * (1.0 + ic_ss(ic))
      y0(3) = temp0

      ss_out = compute_ss(y0(2), y0(3))
      call write_row(aero_labels(ia), ic_labels(ic), "RK45", 0, 0.0_dp, y0(1), ss_out)
      call write_row(aero_labels(ia), ic_labels(ic), "ROS3", 0, 0.0_dp, &
                     y0(1) * (1.0 + perturb_r), ss_out)

      y_rk  = y0
      y_ros = y0
      y_ros(1) = y0(1) * (1.0 + perturb_r)
      h_rk  = h_init
      h_ros = h_init

      ! --- Fine timesteps ---
      do k = 1, nsteps_fine
        tstart = (k - 1) * dt_fine
        tend   = k * dt_fine

        if (k < 4) then
          ss_env = ic_ss(ic)
        else if (k == 4) then
          ss_env = ss_drop
          es_val = saturation_vapor_pressure(y_rk(3))
          qv_sat = 0.622 * es_val / (pres - es_val)
          y_rk(2) = qv_sat * (1.0 + ss_drop)

          es_val = saturation_vapor_pressure(y_ros(3))
          qv_sat = 0.622 * es_val / (pres - es_val)
          y_ros(2) = qv_sat * (1.0 + ss_drop)
        end if

        ! RK45
        call set_aerosol_properties(species, sol_mass, r_sol, inv_grid_mass, ss_env)
        call integrate_ODE(y_rk, tstart, tend, h_rk, stat=istat)
        y_rk(1) = max(y_rk(1), r_sol * 1.01)
        ss_out = compute_ss(y_rk(2), y_rk(3))
        call write_row(aero_labels(ia), ic_labels(ic), "RK45", k, tend, y_rk(1), ss_out)

        ! ROS3
        call set_aerosol_properties(species, sol_mass, r_sol, inv_grid_mass, ss_env)
        call ros3_integrate(growth_rhs, growth_jacobian, 3, y_ros, tstart, tend, &
                              h_ros, rtol_arr, atol_arr, istat)
        y_ros(1) = max(y_ros(1), r_sol * 1.01)
        ss_out = compute_ss(y_ros(2), y_ros(3))
        call write_row(aero_labels(ia), ic_labels(ic), "ROS3", k, tend, y_ros(1), ss_out)
      end do



    end do
  end do

  ! ============================================================
  ! Part 2: Timing benchmark
  ! ============================================================
  write(0,'(/,A)') "=== TIMING BENCHMARK ==="
  write(0,'(A)') "Each solver called 10000 times for a single dt=0.01s integration."
  write(0,'(A)') ""

  ! Non-stiff case: 200nm aerosol, growing (SS=0.005), r=70*r_sol
  r_sol = 200.0e-9
  sol_mass = pi_43 * r_sol**3 * rho_solute
  es_val = saturation_vapor_pressure(temp0)
  qv_sat = 0.622 * es_val / (pres - es_val)
  y0(1) = r_sol * 70.0
  y0(2) = qv_sat * (1.0 + 0.005)
  y0(3) = temp0
  call set_aerosol_properties(species, sol_mass, r_sol, inv_grid_mass, 0.005_dp)

  write(0,'(A,ES10.3,A)') "Non-stiff case: r=", y0(1), " m, SS=0.5%"

  call cpu_time(t_cpu_start)
  do ib = 1, n_bench
    y_bench = y0
    h_bench = h_init
    call rkck45_integrate(growth_rhs, 3, y_bench, 0.0_dp, dt_fine, h_bench, &
                          rtol_arr, atol_arr, istat)
  end do
  call cpu_time(t_cpu_end)
  write(0,'(A,F8.4,A)') "  RK45:  ", t_cpu_end - t_cpu_start, " s"

  call cpu_time(t_cpu_start)
  do ib = 1, n_bench
    y_bench = y0
    h_bench = h_init
    call ros3_integrate(growth_rhs, growth_jacobian, 3, y_bench, 0.0_dp, dt_fine, &
                          h_bench, rtol_arr, atol_arr, istat)
  end do
  call cpu_time(t_cpu_end)
  write(0,'(A,F8.4,A)') "  ROS3:  ", t_cpu_end - t_cpu_start, " s"

  ! Stiff case: 5nm aerosol, haze (SS=-0.005), r=2*r_sol
  r_sol = 5.0e-9
  sol_mass = pi_43 * r_sol**3 * rho_solute
  y0(1) = r_sol * 2.0
  y0(2) = qv_sat * (1.0 - 0.005)
  y0(3) = temp0
  call set_aerosol_properties(species, sol_mass, r_sol, inv_grid_mass, -0.005_dp)

  write(0,'(/,A,ES10.3,A)') "Stiff case: r=", y0(1), " m, SS=-0.5%"

  call cpu_time(t_cpu_start)
  do ib = 1, n_bench
    y_bench = y0
    h_bench = h_init
    call rkck45_integrate(growth_rhs, 3, y_bench, 0.0_dp, dt_fine, h_bench, &
                          rtol_arr, atol_arr, istat)
    if (istat < 0) exit
  end do
  call cpu_time(t_cpu_end)
  if (istat < 0) then
    write(0,'(A,I0,A)') "  RK45:  FAILED after ", ib, " calls (max_steps exceeded)"
  else
    write(0,'(A,F8.4,A)') "  RK45:  ", t_cpu_end - t_cpu_start, " s"
  end if

  call cpu_time(t_cpu_start)
  do ib = 1, n_bench
    y_bench = y0
    h_bench = h_init
    call ros3_integrate(growth_rhs, growth_jacobian, 3, y_bench, 0.0_dp, dt_fine, &
                          h_bench, rtol_arr, atol_arr, istat)
  end do
  call cpu_time(t_cpu_end)
  write(0,'(A,F8.4,A)') "  ROS3:  ", t_cpu_end - t_cpu_start, " s"

contains

  subroutine write_row(aero, ic_tag, solver, step, time, radius, ss_val)
    character(len=*), intent(in) :: aero, ic_tag, solver
    integer, intent(in) :: step
    real(dp), intent(in) :: time, radius, ss_val

    write(*,'(A,"_",A,",",A,",",A,",",A,",",I0,",",ES14.6,",",ES14.6,",",ES14.6)') &
      trim(aero), trim(ic_tag), trim(aero), trim(ic_tag), trim(solver), step, time, radius, ss_val
  end subroutine write_row

  function compute_ss(qv, temp) result(ss)
    real(dp), intent(in) :: qv, temp
    real(dp) :: ss
    real(dp) :: es_loc, qv_sat_loc

    es_loc = saturation_vapor_pressure(temp)
    qv_sat_loc = 0.622 * es_loc / (pres - es_loc)
    ss = qv / qv_sat_loc - 1.0
  end function compute_ss

end program test_dgm_solvers
