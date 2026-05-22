module ode_integrators
  use globals, only: dp, i4
  implicit none

  private
  public :: ode_rhs, ode_jacobian, rkck45_integrate, ros3_integrate

  abstract interface
    pure subroutine ode_rhs(t, y, dydt, ierr)
      import dp, i4
      real(dp),    intent(in)  :: t
      real(dp),    intent(in)  :: y(:)
      real(dp),    intent(out) :: dydt(:)
      integer(i4), intent(out) :: ierr
    end subroutine ode_rhs
  end interface

  abstract interface
    pure subroutine ode_jacobian(t, y, jac)
      import dp
      real(dp), intent(in)  :: t
      real(dp), intent(in)  :: y(:)
      real(dp), intent(out) :: jac(:,:)
    end subroutine ode_jacobian
  end interface

contains


! Cash-Karp embedded Runge-Kutta 4(5) adaptive integrator.
! Reference: Cash & Karp (1990), ACM Trans. Math. Software, 16, 201-222.
subroutine rkck45_integrate(rhs, n, y, t_start, t_end, h_last, rtol, atol, ierr)
  procedure(ode_rhs) :: rhs
  integer(i4), intent(in)    :: n
  real(dp),    intent(inout) :: y(n)
  real(dp),    intent(in)    :: t_start, t_end
  real(dp),    intent(inout) :: h_last
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
  real(dp) :: t, h, h_new, errmax, scale
  integer  :: i, step
  integer(i4) :: ierr_rhs

  ierr = 0
  t = t_start
  h = min(h_last, t_end - t_start)

  do step = 1, max_steps
    if (t >= t_end) then
      h_last = h
      return
    end if
    if (t + h > t_end) h = t_end - t

    ! Stage 1
    call rhs(t, y, k1, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    ytmp = y + h * b21 * k1

    ! Stage 2
    call rhs(t + a2 * h, ytmp, k2, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    ytmp = y + h * (b31 * k1 + b32 * k2)

    ! Stage 3
    call rhs(t + a3 * h, ytmp, k3, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    ytmp = y + h * (b41 * k1 + b42 * k2 + b43 * k3)

    ! Stage 4
    call rhs(t + a4 * h, ytmp, k4, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    ytmp = y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)

    ! Stage 5
    call rhs(t + a5 * h, ytmp, k5, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    ytmp = y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)

    ! Stage 6
    call rhs(t + a6 * h, ytmp, k6, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if

    ! 4th-order solution and error estimate
    errmax = 0.0
    do i = 1, n
      ytmp(i) = y(i) + h * (c1 * k1(i) + c3 * k3(i) + c4 * k4(i) + c6 * k6(i))
      yerr(i) = h * (dc1 * k1(i) + dc3 * k3(i) + dc4 * k4(i) + dc5 * k5(i) + dc6 * k6(i))
      scale   = atol(i) + rtol(i) * max(abs(y(i)), abs(ytmp(i)))
      errmax  = max(errmax, abs(yerr(i)) / scale)
    end do

    if (errmax <= 1.0) then
      t = t + h
      y = ytmp

      if (errmax > 1.0e-30) then
        h_new = safety * h * errmax**(-0.2)
        h = min(h_new, grow_max * h)
      else
        h = grow_max * h
      end if
    else
      h_new = safety * h * errmax**(-0.25)
      h = max(h_new, shrink_min * h)
    end if
  end do

  ierr = -1
  h_last = h

end subroutine rkck45_integrate


! 3-stage, L-stable, order 3(2) Rosenbrock method.
! Reference: ROS3 from Sandu et al. (1997), T-formulation (Hairer & Wanner).
subroutine ros3_integrate(rhs, jac, n, y, t_start, t_end, h_last, rtol, atol, ierr)
  procedure(ode_rhs)      :: rhs
  procedure(ode_jacobian)  :: jac
  integer(i4), intent(in)    :: n
  real(dp),    intent(inout) :: y(n)
  real(dp),    intent(in)    :: t_start, t_end
  real(dp),    intent(inout) :: h_last
  real(dp),    intent(in)    :: rtol(n), atol(n)
  integer(i4), intent(out)   :: ierr

  ! ROS3 coefficients (Sandu et al. 1997, KPP)
  real(dp), parameter :: gamma = 0.43586652150845899941601945119356

  real(dp), parameter :: c21 = -1.0156171083877702091517163265855
  real(dp), parameter :: c31 =  4.0759956452537699824805835358067
  real(dp), parameter :: c32 =  9.2076794298034434657925517058582

  real(dp), parameter :: m1 =  1.0
  real(dp), parameter :: m2 =  6.1697947043828245592553615689730
  real(dp), parameter :: m3 = -0.42772256543218573326238373806514

  real(dp), parameter :: e1 =  0.5
  real(dp), parameter :: e2 = -2.9079558716805469821718236208017
  real(dp), parameter :: e3 =  0.22354069897811569627360909276199

  real(dp), parameter :: safety    = 0.9
  real(dp), parameter :: grow_max  = 5.0
  real(dp), parameter :: shrink_min = 0.1
  integer,  parameter :: max_steps = 10000

  real(dp) :: J_mat(n,n), W(n,n)
  real(dp) :: f1(n), f2(n)
  real(dp) :: k1(n), k2(n), k3(n)
  real(dp) :: ytmp(n), ynew(n), yerr(n)
  real(dp) :: t, h, h_new, errmax, scale, inv_gamma_h
  integer  :: step, i
  integer(i4) :: ierr_rhs
  logical  :: new_jac

  ierr = 0
  t = t_start
  h = min(h_last, t_end - t_start)
  new_jac = .true.

  do step = 1, max_steps
    if (t >= t_end) then
      h_last = h
      return
    end if
    if (t + h > t_end) h = t_end - t

    if (new_jac) then
      call jac(t, y, J_mat)
      new_jac = .false.
    end if

    ! Form W = I/(h*gamma) - J
    inv_gamma_h = 1.0 / (gamma * h)
    W = -J_mat
    do i = 1, n
      W(i,i) = W(i,i) + inv_gamma_h
    end do

    ! Stage 1: W * k1 = f(y)
    call rhs(t, y, f1, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    k1 = f1
    call solve3(W, k1)

    ! Stage 2: W * k2 = f(y + k1) + c21/h * k1
    ytmp = y + k1
    call rhs(t + gamma * h, ytmp, f2, ierr_rhs)
    if (ierr_rhs /= 0) then; ierr = ierr_rhs; return; end if
    k2 = f2 + (c21 / h) * k1
    call solve3(W, k2)

    ! Stage 3: W * k3 = f2 + c31/h * k1 + c32/h * k2  (reuses f2)
    k3 = f2 + (c31 / h) * k1 + (c32 / h) * k2
    call solve3(W, k3)

    ! 3rd-order solution
    ynew = y + m1 * k1 + m2 * k2 + m3 * k3

    ! Error estimate (embedded 2nd-order)
    yerr = e1 * k1 + e2 * k2 + e3 * k3

    errmax = 0.0
    do i = 1, n
      scale = atol(i) + rtol(i) * max(abs(y(i)), abs(ynew(i)))
      errmax = max(errmax, abs(yerr(i)) / scale)
    end do

    if (errmax <= 1.0) then
      t = t + h
      y = ynew
      new_jac = .true.

      if (errmax > 1.0e-30) then
        h_new = safety * h * errmax**(-1.0 / 3.0)
        h = min(h_new, grow_max * h)
      else
        h = grow_max * h
      end if
    else
      h_new = safety * h * errmax**(-0.5)
      h = max(h_new, shrink_min * h)
    end if
  end do

  ierr = -1
  h_last = h

end subroutine ros3_integrate


! Direct 3x3 linear solve W*x = b with partial pivoting.
subroutine solve3(W, b)
  real(dp), intent(in)    :: W(3,3)
  real(dp), intent(inout) :: b(3)

  real(dp) :: A(3,3), tmp
  integer  :: piv, j, i

  A = W

  ! Column 1: pivot
  piv = 1
  if (abs(A(2,1)) > abs(A(piv,1))) piv = 2
  if (abs(A(3,1)) > abs(A(piv,1))) piv = 3
  if (piv /= 1) then
    do j = 1, 3; tmp = A(1,j); A(1,j) = A(piv,j); A(piv,j) = tmp; end do
    tmp = b(1); b(1) = b(piv); b(piv) = tmp
  end if
  do i = 2, 3
    A(i,1) = A(i,1) / A(1,1)
    do j = 2, 3
      A(i,j) = A(i,j) - A(i,1) * A(1,j)
    end do
    b(i) = b(i) - A(i,1) * b(1)
  end do

  ! Column 2: pivot rows 2-3
  if (abs(A(3,2)) > abs(A(2,2))) then
    do j = 2, 3; tmp = A(2,j); A(2,j) = A(3,j); A(3,j) = tmp; end do
    tmp = b(2); b(2) = b(3); b(3) = tmp
  end if
  A(3,2) = A(3,2) / A(2,2)
  A(3,3) = A(3,3) - A(3,2) * A(2,3)
  b(3)   = b(3) - A(3,2) * b(2)

  ! Back-substitution
  b(3) = b(3) / A(3,3)
  b(2) = (b(2) - A(2,3) * b(3)) / A(2,2)
  b(1) = (b(1) - A(1,2) * b(2) - A(1,3) * b(3)) / A(1,1)

end subroutine solve3

end module ode_integrators
