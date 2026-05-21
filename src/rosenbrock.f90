! 3-stage L-stable Rosenbrock solver for small stiff ODE systems.
! Coefficients: ROS3 from Sandu et al. (1997), as implemented in KPP.
! Formulation: T-type (Hairer & Wanner, Solving ODEs II, Section IV.7).
module rosenbrock
  use globals, only: dp, i4
  use ode_integrators, only: ode_rhs
  implicit none

  private
  public :: ros3_integrate, ode_jacobian

  abstract interface
    subroutine ode_jacobian(t, y, jac)
      import dp
      real(dp), intent(in)  :: t
      real(dp), intent(in)  :: y(:)
      real(dp), intent(out) :: jac(:,:)
    end subroutine ode_jacobian
  end interface

contains

! 3-stage, L-stable, order 3(2) Rosenbrock method (2 function evaluations per step).
!
! T-formulation stages:
!   W = I/(h*gamma) - J
!   W * k1 = f(y)
!   W * k2 = f(y + k1) + c21/h * k1
!   W * k3 = f2        + c31/h * k1 + c32/h * k2  (reuses f2, no new eval)
!   y_{n+1} = y + m1*k1 + m2*k2 + m3*k3
!
! h_last is inout for warm-start across calls.
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

  ! c_ij: coupling coefficients in stage equations
  real(dp), parameter :: c21 = -1.0156171083877702091517163265855
  real(dp), parameter :: c31 =  4.0759956452537699824805835358067
  real(dp), parameter :: c32 =  9.2076794298034434657925517058582

  ! m_i: solution weights
  real(dp), parameter :: m1 =  1.0
  real(dp), parameter :: m2 =  6.1697947043828245592553615689730
  real(dp), parameter :: m3 = -0.42772256543218573326238373806514

  ! e_i: error estimate weights
  real(dp), parameter :: e1 =  0.5
  real(dp), parameter :: e2 = -2.9079558716805469821718236208017
  real(dp), parameter :: e3 =  0.22354069897811569627360909276199

  ! Step control
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
    call rhs(t, y, f1)
    k1 = f1
    call solve3(W, k1)

    ! Stage 2: W * k2 = f(y + k1) + c21/h * k1
    ytmp = y + k1
    call rhs(t + gamma * h, ytmp, f2)
    k2 = f2 + (c21 / h) * k1
    call solve3(W, k2)

    ! Stage 3: W * k3 = f2 + c31/h * k1 + c32/h * k2  (reuses f2)
    k3 = f2 + (c31 / h) * k1 + (c32 / h) * k2
    call solve3(W, k3)

    ! 3rd-order solution
    ynew = y + m1 * k1 + m2 * k2 + m3 * k3

    ! Error estimate (embedded 2nd-order)
    yerr = e1 * k1 + e2 * k2 + e3 * k3

    ! Scaled error norm
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
! Overwrites b with the solution x.
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

end module rosenbrock
