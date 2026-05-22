program test_jac_magnitudes
  use globals
  use microphysics, only: saturation_vapor_pressure
  use DGM, only: set_aerosol_properties, growth_jacobian
  implicit none

  integer, parameter :: n_cases = 12
  real(dp) :: r_sol, sol_mass, es_val, qv_sat
  real(dp) :: y(3), jac(3,3)
  real(dp) :: inv_grid_mass, grid_rho, grid_mass
  integer :: i, j, k

  real(dp) :: radii(n_cases), ss_vals(n_cases), r_solutes(n_cases)
  character(len=40) :: labels(n_cases)

  grid_rho = pres / (Rd * 280.0)
  grid_mass = 50.0 * domain_width * domain_width * (H / N) * grid_rho
  inv_grid_mass = 1.0 / grid_mass

  labels(1)  = "5nm haze (r=10nm, S=-0.5%)";       r_solutes(1)  = 5e-9;    radii(1)  = 10e-9;   ss_vals(1)  = -0.005
  labels(2)  = "5nm critical (r=50nm, S=0.1%)";     r_solutes(2)  = 5e-9;    radii(2)  = 50e-9;   ss_vals(2)  = 0.001
  labels(3)  = "50nm haze (r=125nm, S=-0.5%)";      r_solutes(3)  = 50e-9;   radii(3)  = 125e-9;  ss_vals(3)  = -0.005
  labels(4)  = "50nm growing (r=3um, S=0.5%)";      r_solutes(4)  = 50e-9;   radii(4)  = 3e-6;    ss_vals(4)  = 0.005
  labels(5)  = "200nm haze (r=600nm, S=-0.5%)";     r_solutes(5)  = 200e-9;  radii(5)  = 600e-9;  ss_vals(5)  = -0.005
  labels(6)  = "200nm growing (r=14um, S=0.5%)";    r_solutes(6)  = 200e-9;  radii(6)  = 14e-6;   ss_vals(6)  = 0.005
  labels(7)  = "1000nm growing (r=80um, S=0.5%)";   r_solutes(7)  = 1000e-9; radii(7)  = 80e-6;   ss_vals(7)  = 0.005
  labels(8)  = "1000nm large (r=260um, S=0.2%)";    r_solutes(8)  = 1000e-9; radii(8)  = 260e-6;  ss_vals(8)  = 0.002
  ! Deeply subsaturated — entrainment-like
  labels(9)  = "5nm haze (r=10nm, S=-50%)";         r_solutes(9)  = 5e-9;    radii(9)  = 10e-9;   ss_vals(9)  = -0.50
  labels(10) = "50nm haze (r=125nm, S=-50%)";       r_solutes(10) = 50e-9;   radii(10) = 125e-9;  ss_vals(10) = -0.50
  labels(11) = "200nm growing (r=14um, S=-50%)";    r_solutes(11) = 200e-9;  radii(11) = 14e-6;   ss_vals(11) = -0.50
  labels(12) = "1000nm large (r=260um, S=-50%)";    r_solutes(12) = 1000e-9; radii(12) = 260e-6;  ss_vals(12) = -0.50

  write(*,'(A)') "Jacobian magnitudes: jac(i,j) for each case"
  write(*,'(A)') "================================================"

  do k = 1, n_cases
    r_sol = r_solutes(k)
    sol_mass = pi_43 * r_sol**3 * 2165.0
    es_val = saturation_vapor_pressure(280.0_dp)
    qv_sat = 0.622 * es_val / (pres - es_val)

    call set_aerosol_properties(2_i4, sol_mass, r_sol, inv_grid_mass, ss_vals(k))

    y(1) = radii(k)
    y(2) = qv_sat * (1.0 + ss_vals(k))
    y(3) = 280.0

    call growth_jacobian(0.0_dp, y, jac)

    write(*,'(/,A)') trim(labels(k))
    write(*,'(A)') "            dr           dqv          dT"
    write(*,'(A,3(ES12.4))') "  dr/dt: ", (jac(1,j), j=1,3)
    write(*,'(A,3(ES12.4))') "  dqv/dt:", (jac(2,j), j=1,3)
    write(*,'(A,3(ES12.4))') "  dT/dt: ", (jac(3,j), j=1,3)
  end do

end program
