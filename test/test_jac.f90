program test_jac
  use globals
  use DGM, only: set_aerosol_properties, growth_rhs, growth_jacobian
  use microphysics, only: saturation_vapor_pressure
  implicit none

  real(dp) :: y(3), jac(3,3), jac_fd(3,3)
  real(dp) :: f0(3), fp(3), yp(3)
  real(dp) :: es_val, qv_sat, delta
  integer :: j

  real(dp), parameter :: r_sol = 200.0e-9
  real(dp), parameter :: rho_solute = 2165.0
  real(dp), parameter :: sol_mass = 4.1887902047863905 * r_sol**3 * rho_solute
  real(dp), parameter :: temp0 = 280.0
  real(dp), parameter :: inv_gm = 1.0e-3

  es_val = saturation_vapor_pressure(temp0)
  qv_sat = 0.622 * es_val / (pres - es_val)

  y(1) = r_sol * 70.0  ! 14 um
  y(2) = qv_sat * 1.005
  y(3) = temp0

  call set_aerosol_properties(2, sol_mass, r_sol, inv_gm, 0.005_dp)

  call growth_rhs(0.0_dp, y, f0)
  call growth_jacobian(0.0_dp, y, jac)

  ! Finite difference Jacobian
  delta = 1.0e-8
  do j = 1, 3
    yp = y
    yp(j) = yp(j) + delta * max(abs(y(j)), 1.0e-15)
    call growth_rhs(0.0_dp, yp, fp)
    jac_fd(:,j) = (fp - f0) / (delta * max(abs(y(j)), 1.0e-15))
  end do

  write(*,'(A)') "RHS:"
  write(*,'(3ES16.6)') f0

  write(*,'(/,A)') "Analytical Jacobian:"
  do j = 1, 3
    write(*,'(3ES16.6)') jac(j,:)
  end do

  write(*,'(/,A)') "Finite-diff Jacobian:"
  do j = 1, 3
    write(*,'(3ES16.6)') jac_fd(j,:)
  end do

  write(*,'(/,A)') "Relative error (|ana-fd|/|fd|):"
  do j = 1, 3
    write(*,'(3ES16.6)') abs(jac(j,:) - jac_fd(j,:)) / max(abs(jac_fd(j,:)), 1.0e-30)
  end do

end program test_jac
