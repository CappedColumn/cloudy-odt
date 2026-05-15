module ode_integrators
  use globals, only: dp, i4
  implicit none

  private
  public :: ode_rhs, ode_integrate, rkck45_integrate, integrator_iface

  abstract interface
    subroutine ode_rhs(t, y, dydt)
      import dp
      real(dp), intent(in)  :: t
      real(dp), intent(in)  :: y(:)
      real(dp), intent(out) :: dydt(:)
    end subroutine ode_rhs
  end interface

  abstract interface
    subroutine integrator_iface(rhs, n, y, t_start, t_end, h, rtol, atol, ierr)
      import dp, i4, ode_rhs
      procedure(ode_rhs) :: rhs
      integer(i4), intent(in)    :: n
      real(dp),    intent(inout) :: y(n)
      real(dp),    intent(in)    :: t_start, t_end, h
      real(dp),    intent(in)    :: rtol(n), atol(n)
      integer(i4), intent(out)   :: ierr
    end subroutine integrator_iface
  end interface

  procedure(integrator_iface), pointer :: ode_integrate => rkck45_integrate

contains

! Cash-Karp embedded Runge-Kutta 4(5) adaptive integrator.
! Reference: Cash & Karp (1990), ACM Trans. Math. Software, 16, 201-222.
subroutine rkck45_integrate(rhs, n, y, t_start, t_end, h, rtol, atol, ierr)
  procedure(ode_rhs) :: rhs
  integer(i4), intent(in)    :: n
  real(dp),    intent(inout) :: y(n)
  real(dp),    intent(in)    :: t_start, t_end, h
  real(dp),    intent(in)    :: rtol(n), atol(n)
  integer(i4), intent(out)   :: ierr

  ! Butcher tableau — Cash-Karp coefficients
  real(dp), parameter :: a2 = 0.2,  a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875

  real(dp), parameter :: b21 = 0.2
  real(dp), parameter :: b31 = 3.0/40.0,     b32 = 9.0/40.0
  real(dp), parameter :: b41 = 0.3,           b42 = -0.9,           b43 = 1.2
  real(dp), parameter :: b51 = -11.0/54.0,    b52 = 2.5,            b53 = -70.0/27.0,    b54 = 35.0/27.0
  real(dp), parameter :: b61 = 1631.0/55296.0, b62 = 175.0/512.0, &
                         b63 = 575.0/13824.0,  b64 = 44275.0/110592.0, b65 = 253.0/4096.0

  ! 4th-order solution weights
  real(dp), parameter :: c1 = 37.0/378.0,  c3 = 250.0/621.0, &
                         c4 = 125.0/594.0, c6 = 512.0/1771.0

  ! Error estimate weights (4th - 5th order difference)
  real(dp), parameter :: dc1 = c1 - 2825.0/27648.0
  real(dp), parameter :: dc3 = c3 - 18575.0/48384.0
  real(dp), parameter :: dc4 = c4 - 13525.0/55296.0
  real(dp), parameter :: dc5 = -277.0/14336.0
  real(dp), parameter :: dc6 = c6 - 0.25

  ! Step control parameters
  real(dp), parameter :: safety    = 0.9
  real(dp), parameter :: grow_max  = 5.0
  real(dp), parameter :: shrink_min = 0.1
  integer,  parameter :: max_steps = 10000

  real(dp) :: k1(n), k2(n), k3(n), k4(n), k5(n), k6(n)
  real(dp) :: ytmp(n), yerr(n)
  real(dp) :: t, dt, dt_new, errmax, scale
  integer  :: i, step

  ierr = 0
  t  = t_start
  dt = min(h, t_end - t_start)

  do step = 1, max_steps
    if (t >= t_end) return
    if (t + dt > t_end) dt = t_end - t

    ! Stage 1
    call rhs(t, y, k1)
    ytmp = y + dt * b21 * k1

    ! Stage 2
    call rhs(t + a2 * dt, ytmp, k2)
    ytmp = y + dt * (b31 * k1 + b32 * k2)

    ! Stage 3
    call rhs(t + a3 * dt, ytmp, k3)
    ytmp = y + dt * (b41 * k1 + b42 * k2 + b43 * k3)

    ! Stage 4
    call rhs(t + a4 * dt, ytmp, k4)
    ytmp = y + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)

    ! Stage 5
    call rhs(t + a5 * dt, ytmp, k5)
    ytmp = y + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)

    ! Stage 6
    call rhs(t + a6 * dt, ytmp, k6)

    ! 4th-order solution and error estimate
    errmax = 0.0
    do i = 1, n
      ytmp(i) = y(i) + dt * (c1 * k1(i) + c3 * k3(i) + c4 * k4(i) + c6 * k6(i))
      yerr(i) = dt * (dc1 * k1(i) + dc3 * k3(i) + dc4 * k4(i) + dc5 * k5(i) + dc6 * k6(i))
      scale   = atol(i) + rtol(i) * max(abs(y(i)), abs(ytmp(i)))
      errmax  = max(errmax, abs(yerr(i)) / scale)
    end do

    if (errmax <= 1.0) then
      ! Accept step
      t = t + dt
      y = ytmp

      if (errmax > 1.0e-30) then
        dt_new = safety * dt * errmax**(-0.2)
        dt = min(dt_new, grow_max * dt)
      else
        dt = grow_max * dt
      end if
    else
      ! Reject step — reduce dt
      dt_new = safety * dt * errmax**(-0.25)
      dt = max(dt_new, shrink_min * dt)
    end if
  end do

  ierr = -1

end subroutine rkck45_integrate

end module ode_integrators
