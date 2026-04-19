module ode_integrators
    use globals, only: dp, i4
    implicit none

    abstract interface
        subroutine ode_rhs(t, y, dydt)
            import dp
            real(dp), intent(in) :: t
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dydt(:)
        end subroutine ode_rhs
    end interface

    procedure(integrator_iface), pointer :: ode_integrate => rkck45_integrate

    abstract interface
        subroutine integrator_iface(f, n, y, t1, t2, h, rtol, atol, ierr)
            import dp, i4, ode_rhs
            procedure(ode_rhs) :: f
            integer(i4), intent(in) :: n
            real(dp), intent(inout) :: y(n)
            real(dp), intent(in) :: t1, t2, h
            real(dp), intent(in) :: rtol(n), atol(n)
            integer(i4), intent(out) :: ierr
        end subroutine integrator_iface
    end interface

    private
    public :: ode_rhs, ode_integrate, rkck45_integrate, integrator_iface

contains

subroutine rkck45_integrate(f, n, y, t1, t2, h, rtol, atol, ierr)
    procedure(ode_rhs) :: f
    integer(i4), intent(in) :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(in) :: t1, t2, h
    real(dp), intent(in) :: rtol(n), atol(n)
    integer(i4), intent(out) :: ierr

    real(dp), parameter :: a2 = 0.2_dp
    real(dp), parameter :: a3 = 0.3_dp
    real(dp), parameter :: a4 = 0.6_dp
    real(dp), parameter :: a5 = 1.0_dp
    real(dp), parameter :: a6 = 0.875_dp

    real(dp), parameter :: b21 = 0.2_dp
    real(dp), parameter :: b31 = 3.0_dp/40.0_dp,  b32 = 9.0_dp/40.0_dp
    real(dp), parameter :: b41 = 0.3_dp,           b42 = -0.9_dp,          b43 = 1.2_dp
    real(dp), parameter :: b51 = -11.0_dp/54.0_dp, b52 = 2.5_dp,          b53 = -70.0_dp/27.0_dp, b54 = 35.0_dp/27.0_dp
    real(dp), parameter :: b61 = 1631.0_dp/55296.0_dp, b62 = 175.0_dp/512.0_dp, &
                           b63 = 575.0_dp/13824.0_dp,  b64 = 44275.0_dp/110592.0_dp, b65 = 253.0_dp/4096.0_dp

    real(dp), parameter :: c1 = 37.0_dp/378.0_dp, c3 = 250.0_dp/621.0_dp, &
                           c4 = 125.0_dp/594.0_dp, c6 = 512.0_dp/1771.0_dp

    real(dp), parameter :: dc1 = c1 - 2825.0_dp/27648.0_dp
    real(dp), parameter :: dc3 = c3 - 18575.0_dp/48384.0_dp
    real(dp), parameter :: dc4 = c4 - 13525.0_dp/55296.0_dp
    real(dp), parameter :: dc5 = -277.0_dp/14336.0_dp
    real(dp), parameter :: dc6 = c6 - 0.25_dp

    real(dp), parameter :: safety = 0.9_dp
    real(dp), parameter :: grow_max = 5.0_dp
    real(dp), parameter :: shrink_min = 0.1_dp
    integer, parameter  :: max_steps = 10000

    real(dp) :: k1(n), k2(n), k3(n), k4(n), k5(n), k6(n)
    real(dp) :: ytmp(n), yerr(n)
    real(dp) :: t, dt, dt_new, errmax, scale
    integer :: i, step

    ierr = 0
    t = t1
    dt = h

    if (dt > t2 - t1) dt = t2 - t1

    do step = 1, max_steps
        if (t >= t2) return

        if (t + dt > t2) dt = t2 - t

        call f(t, y, k1)
        do i = 1, n
            ytmp(i) = y(i) + dt*b21*k1(i)
        end do

        call f(t + a2*dt, ytmp, k2)
        do i = 1, n
            ytmp(i) = y(i) + dt*(b31*k1(i) + b32*k2(i))
        end do

        call f(t + a3*dt, ytmp, k3)
        do i = 1, n
            ytmp(i) = y(i) + dt*(b41*k1(i) + b42*k2(i) + b43*k3(i))
        end do

        call f(t + a4*dt, ytmp, k4)
        do i = 1, n
            ytmp(i) = y(i) + dt*(b51*k1(i) + b52*k2(i) + b53*k3(i) + b54*k4(i))
        end do

        call f(t + a5*dt, ytmp, k5)
        do i = 1, n
            ytmp(i) = y(i) + dt*(b61*k1(i) + b62*k2(i) + b63*k3(i) + b64*k4(i) + b65*k5(i))
        end do

        call f(t + a6*dt, ytmp, k6)

        errmax = 0.0_dp
        do i = 1, n
            ytmp(i) = y(i) + dt*(c1*k1(i) + c3*k3(i) + c4*k4(i) + c6*k6(i))
            yerr(i) = dt*(dc1*k1(i) + dc3*k3(i) + dc4*k4(i) + dc5*k5(i) + dc6*k6(i))
            scale = atol(i) + rtol(i) * max(abs(y(i)), abs(ytmp(i)))
            errmax = max(errmax, abs(yerr(i)) / scale)
        end do

        if (errmax <= 1.0_dp) then
            t = t + dt
            y = ytmp

            if (errmax > 1.0e-30_dp) then
                dt_new = safety * dt * errmax**(-0.2_dp)
                dt = min(dt_new, grow_max * dt)
            else
                dt = grow_max * dt
            end if
        else
            dt_new = safety * dt * errmax**(-0.25_dp)
            dt = max(dt_new, shrink_min * dt)
        end if
    end do

    ierr = -1

end subroutine rkck45_integrate

end module ode_integrators
