program test_ros_simple
  use globals, only: dp, i4
  use rosenbrock, only: ros3_integrate
  implicit none

  real(dp) :: y(3), h_last, rtol(3), atol(3)
  integer(i4) :: ierr

  rtol = [1.0e-4, 1.0e-4, 1.0e-4]
  atol = [1.0e-10, 1.0e-10, 1.0e-10]

  ! Test 1: dy/dt = [1, 0, 0], y(0) = [0, 0, 0], solve to t=1
  y = [0.0_dp, 0.0_dp, 0.0_dp]
  h_last = 0.1_dp
  call ros3_integrate(const_rhs, zero_jac, 3, y, 0.0_dp, 1.0_dp, h_last, rtol, atol, ierr)
  write(*,'(A,ES16.8,A,ES16.8,A,I0)') "dy/dt=1: y(1)=", y(1), " err=", abs(y(1) - 1.0), " ierr=", ierr

  ! Test 2: dy/dt = [y1, 0, 0], y(0) = [1, 0, 0], solve to t=1 (expect e)
  y = [1.0_dp, 0.0_dp, 0.0_dp]
  h_last = 0.1_dp
  call ros3_integrate(exp_rhs, exp_jac, 3, y, 0.0_dp, 1.0_dp, h_last, rtol, atol, ierr)
  write(*,'(A,ES16.8,A,ES16.8,A,I0)') "dy/dt=y: y(1)=", y(1), " err=", abs(y(1) - exp(1.0_dp)), " ierr=", ierr

  ! Test 3: Stiff test dy/dt = -1000*y, y(0) = 1, solve to t=0.01 (expect exp(-10))
  y = [1.0_dp, 0.0_dp, 0.0_dp]
  h_last = 0.001_dp
  call ros3_integrate(stiff_rhs, stiff_jac, 3, y, 0.0_dp, 0.01_dp, h_last, rtol, atol, ierr)
  write(*,'(A,ES16.8,A,ES16.8,A,I0)') "stiff:   y(1)=", y(1), " err=", abs(y(1) - exp(-10.0_dp)), " ierr=", ierr

contains

  subroutine const_rhs(t, yv, dydt)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: dydt(:)
    dydt = [1.0_dp, 0.0_dp, 0.0_dp]
  end subroutine

  subroutine zero_jac(t, yv, jac)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: jac(:,:)
    jac = 0.0_dp
  end subroutine

  subroutine exp_rhs(t, yv, dydt)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: dydt(:)
    dydt = [yv(1), 0.0_dp, 0.0_dp]
  end subroutine

  subroutine exp_jac(t, yv, jac)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: jac(:,:)
    jac = 0.0_dp
    jac(1,1) = 1.0_dp
  end subroutine

  subroutine stiff_rhs(t, yv, dydt)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: dydt(:)
    dydt = [-1000.0_dp * yv(1), 0.0_dp, 0.0_dp]
  end subroutine

  subroutine stiff_jac(t, yv, jac)
    real(dp), intent(in) :: t, yv(:)
    real(dp), intent(out) :: jac(:,:)
    jac = 0.0_dp
    jac(1,1) = -1000.0_dp
  end subroutine

end program test_ros_simple
