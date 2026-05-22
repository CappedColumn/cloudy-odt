! Adaptive ODE integrators for small systems.
!
! Provides two solvers with shared step-control logic:
!   rkck45_integrate  — Cash-Karp RK4(5), explicit, for non-stiff systems
!   ros3_integrate    — ROS3 Rosenbrock 3(2), L-stable, for stiff systems
!
! Both integrate y'(t) = f(t,y) from t_start to t_end with local error control.
! Step size h_last is preserved across calls for warm-starting.
!
! Error control: at each step the local error is scaled as
!   err_i = |yerr_i| / (atol_i + rtol_i * max(|y_i|, |ynew_i|))
! The step is accepted when max(err_i) <= 1.
!
! On success ierr = 0; on max_steps exceeded ierr = -1;
! if the RHS returns nonzero, that value is forwarded as ierr.
module ode_integrators
  use globals, only: dp, i4
  implicit none

  private
  public :: ode_rhs, ode_rhs_jac, rkck45_integrate, ros3_integrate

  ! Step-size controller constants (shared by both solvers)
  real(dp), parameter :: SAFETY     = 0.9   ! safety factor on new step size
  real(dp), parameter :: GROW_MAX   = 5.0   ! max step size growth ratio
  real(dp), parameter :: SHRINK_MIN = 0.1   ! min step size shrink ratio
  integer,  parameter :: MAX_STEPS  = 10000

  abstract interface
    ! Right-hand side: dydt = f(t, y).  ierr /= 0 signals failure.
    pure subroutine ode_rhs(t, y, dydt, ierr)
      import dp, i4
      real(dp),    intent(in)  :: t
      real(dp),    intent(in)  :: y(:)
      real(dp),    intent(out) :: dydt(:)
      integer(i4), intent(out) :: ierr
    end subroutine ode_rhs

    ! Fused RHS + Jacobian: computes both dydt and jac at (t, y) in one call,
    ! sharing intermediate quantities (thermodynamic properties, etc.).
    pure subroutine ode_rhs_jac(t, y, dydt, jac, ierr)
      import dp, i4
      real(dp),    intent(in)  :: t
      real(dp),    intent(in)  :: y(:)
      real(dp),    intent(out) :: dydt(:)
      real(dp),    intent(out) :: jac(:,:)
      integer(i4), intent(out) :: ierr
    end subroutine ode_rhs_jac
  end interface

contains


! ==========================================================================
! Cash-Karp RK4(5) — explicit adaptive integrator
! Reference: Cash & Karp (1990), ACM Trans. Math. Software 16:201-222
!
! 6 RHS evaluations per step.  Embedded 4th/5th order pair for error
! estimation.  Step accepted at 4th-order accuracy.
!
! Arguments:
!   rhs      — RHS subroutine matching ode_rhs interface
!   nvar     — system size
!   y        — state vector, overwritten with solution at t_end
!   t_start  — integration start time
!   t_end    — integration end time
!   h_last   — step size hint on entry, last accepted step on exit
!   rtol     — relative error tolerance per component
!   atol     — absolute error tolerance per component
!   ierr     — 0 on success, -1 if max_steps exceeded, or RHS error code
! ==========================================================================
subroutine rkck45_integrate(rhs, nvar, y, t_start, t_end, h_last, rtol, atol, ierr)
  procedure(ode_rhs)          :: rhs
  integer(i4), intent(in)    :: nvar
  real(dp),    intent(inout) :: y(nvar)
  real(dp),    intent(in)    :: t_start, t_end
  real(dp),    intent(inout) :: h_last
  real(dp),    intent(in)    :: rtol(nvar), atol(nvar)
  integer(i4), intent(out)   :: ierr

  ! ---- Butcher tableau (Cash-Karp coefficients) ----
  ! Time fractions for stages 2-6
  real(dp), parameter :: a2 = 0.2,  a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875

  ! Stage coupling weights
  real(dp), parameter :: b21 = 0.2
  real(dp), parameter :: b31 = 3.0/40.0,       b32 = 9.0/40.0
  real(dp), parameter :: b41 = 0.3,             b42 = -0.9,             b43 = 1.2
  real(dp), parameter :: b51 = -11.0/54.0,      b52 = 2.5,              b53 = -70.0/27.0,      b54 = 35.0/27.0
  real(dp), parameter :: b61 = 1631.0/55296.0,  b62 = 175.0/512.0, &
                         b63 = 575.0/13824.0,   b64 = 44275.0/110592.0, b65 = 253.0/4096.0

  ! 4th-order solution weights
  real(dp), parameter :: c1 = 37.0/378.0,  c3 = 250.0/621.0, &
                         c4 = 125.0/594.0, c6 = 512.0/1771.0

  ! Error weights: difference between 4th and 5th order solutions
  real(dp), parameter :: dc1 = c1 - 2825.0/27648.0
  real(dp), parameter :: dc3 = c3 - 18575.0/48384.0
  real(dp), parameter :: dc4 = c4 - 13525.0/55296.0
  real(dp), parameter :: dc5 = -277.0/14336.0
  real(dp), parameter :: dc6 = c6 - 0.25

  ! ---- Work arrays ----
  real(dp) :: k1(nvar), k2(nvar), k3(nvar), k4(nvar), k5(nvar), k6(nvar)  ! stage derivatives
  real(dp) :: y_trial(nvar)  ! candidate solution
  real(dp) :: y_err(nvar)    ! local error estimate

  real(dp) :: t              ! current time
  real(dp) :: h              ! current step size
  real(dp) :: h_new          ! proposed next step size
  real(dp) :: err_max        ! worst-case scaled error across components
  real(dp) :: tol_scale      ! per-component error scaling
  integer  :: i, step
  integer(i4) :: rhs_err     ! RHS error flag

  ierr = 0
  t = t_start
  h = min(h_last, t_end - t_start)

  do step = 1, MAX_STEPS
    if (t >= t_end) then
      h_last = h
      return
    end if
    if (t + h > t_end) h = t_end - t

    ! Evaluate 6 stages of the Cash-Karp tableau
    call rhs(t, y, k1, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    y_trial = y + h * b21 * k1

    call rhs(t + a2 * h, y_trial, k2, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    y_trial = y + h * (b31 * k1 + b32 * k2)

    call rhs(t + a3 * h, y_trial, k3, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    y_trial = y + h * (b41 * k1 + b42 * k2 + b43 * k3)

    call rhs(t + a4 * h, y_trial, k4, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    y_trial = y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)

    call rhs(t + a5 * h, y_trial, k5, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    y_trial = y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)

    call rhs(t + a6 * h, y_trial, k6, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if

    ! Compute 4th-order solution and embedded error estimate
    err_max = 0.0
    do i = 1, nvar
      y_trial(i) = y(i) + h * (c1 * k1(i) + c3 * k3(i) + c4 * k4(i) + c6 * k6(i))
      y_err(i)   = h * (dc1 * k1(i) + dc3 * k3(i) + dc4 * k4(i) + dc5 * k5(i) + dc6 * k6(i))
      tol_scale  = atol(i) + rtol(i) * max(abs(y(i)), abs(y_trial(i)))
      err_max    = max(err_max, abs(y_err(i)) / tol_scale)
    end do

    ! Accept or reject, then adjust step size
    call adjust_step_size(h, err_max, 0.2_dp, 0.25_dp, h_new)
    if (err_max <= 1.0) then
      t = t + h
      y = y_trial
    end if
    h = h_new
  end do

  ierr = -1
  h_last = h

end subroutine rkck45_integrate


! ==========================================================================
! ROS3 Rosenbrock 3(2) — L-stable implicit adaptive integrator
! Reference: Sandu et al. (1997), T-formulation per Hairer & Wanner
!
! Per step: 1 fused RHS+Jacobian (Stage 1) + 1 RHS (Stage 2) + 3 solves.
! The fused call avoids redundant computation of shared thermodynamic
! intermediates between the RHS and Jacobian evaluations.
!
! T-formulation stages (W = I/(h*gamma) - J):
!   W * k1 = f(y)
!   W * k2 = f(y + k1) + (c21/h) * k1
!   W * k3 = f2        + (c31/h) * k1 + (c32/h) * k2   [reuses f2]
!   y_new  = y + m1*k1 + m2*k2 + m3*k3
!   y_err  =     e1*k1 + e2*k2 + e3*k3
!
! Arguments:
!   rhs      — RHS-only subroutine (ode_rhs) for Stage 2
!   rhs_jac  — fused RHS+Jacobian subroutine (ode_rhs_jac) for Stage 1
!   nvar     — system size
!   y        — state vector, overwritten with solution at t_end
!   t_start  — integration start time
!   t_end    — integration end time
!   h_last   — step size hint on entry, last accepted step on exit
!   rtol     — relative error tolerance per component
!   atol     — absolute error tolerance per component
!   ierr     — 0 on success, -1 if max_steps exceeded, or RHS error code
! ==========================================================================
subroutine ros3_integrate(rhs, rhs_jac, nvar, y, t_start, t_end, h_last, rtol, atol, ierr)
  procedure(ode_rhs)         :: rhs
  procedure(ode_rhs_jac)     :: rhs_jac
  integer(i4), intent(in)    :: nvar
  real(dp),    intent(inout) :: y(nvar)
  real(dp),    intent(in)    :: t_start, t_end
  real(dp),    intent(inout) :: h_last
  real(dp),    intent(in)    :: rtol(nvar), atol(nvar)
  integer(i4), intent(out)   :: ierr

  ! ---- ROS3 coefficients (Sandu et al. 1997, KPP) ----
  ! Implicit parameter: appears in W = I/(h*gamma) - J
  real(dp), parameter :: gamma = 0.43586652150845899941601945119356

  ! Stage coupling coefficients
  real(dp), parameter :: c21 = -1.0156171083877702091517163265855
  real(dp), parameter :: c31 =  4.0759956452537699824805835358067
  real(dp), parameter :: c32 =  9.2076794298034434657925517058582

  ! Solution weights (3rd-order)
  real(dp), parameter :: m1 =  1.0
  real(dp), parameter :: m2 =  6.1697947043828245592553615689730
  real(dp), parameter :: m3 = -0.42772256543218573326238373806514

  ! Error estimate weights (embedded 2nd-order)
  real(dp), parameter :: e1 =  0.5
  real(dp), parameter :: e2 = -2.9079558716805469821718236208017
  real(dp), parameter :: e3 =  0.22354069897811569627360909276199

  ! ---- Work arrays ----
  real(dp) :: J_mat(nvar,nvar)     ! Jacobian matrix J = df/dy
  real(dp) :: W_matrix(nvar,nvar)  ! iteration matrix W = I/(h*gamma) - J
  real(dp) :: f_current(nvar)      ! f(t, y) at current state
  real(dp) :: f_stage2(nvar)       ! f(t + gamma*h, y + k1)
  real(dp) :: k1(nvar), k2(nvar), k3(nvar)  ! stage vectors (solutions of W*k = rhs)
  real(dp) :: y_trial(nvar)        ! candidate solution at t + h
  real(dp) :: y_err(nvar)          ! local error estimate
  real(dp) :: y_stage(nvar)        ! intermediate state for RHS evaluation

  real(dp) :: t                    ! current time
  real(dp) :: h                    ! current step size
  real(dp) :: h_new                ! proposed next step size
  real(dp) :: err_max              ! worst-case scaled error across components
  real(dp) :: tol_scale            ! per-component error scaling
  real(dp) :: inv_gamma_h          ! 1/(gamma*h), diagonal term of W
  integer  :: step, i
  integer(i4) :: rhs_err           ! RHS error flag
  logical  :: need_jacobian        ! true when J must be recomputed

  ierr = 0
  t = t_start
  h = min(h_last, t_end - t_start)
  need_jacobian = .true.

  do step = 1, MAX_STEPS
    if (t >= t_end) then
      h_last = h
      return
    end if
    if (t + h > t_end) h = t_end - t

    ! Stage 1: evaluate f(y) and (if needed) J(y) at the current state.
    ! When the Jacobian is needed, the fused call computes both together,
    ! sharing thermodynamic intermediates.  On rejected steps the Jacobian
    ! is reused and only the RHS is recomputed.
    if (need_jacobian) then
      call rhs_jac(t, y, f_current, J_mat, rhs_err)
      need_jacobian = .false.
    else
      call rhs(t, y, f_current, rhs_err)
    end if
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if

    ! Build iteration matrix: W = I/(h*gamma) - J
    inv_gamma_h = 1.0 / (gamma * h)
    W_matrix = -J_mat
    do i = 1, nvar
      W_matrix(i,i) = W_matrix(i,i) + inv_gamma_h
    end do

    ! Stage 1 solve: W * k1 = f(y)
    k1 = f_current
    call solve3(W_matrix, k1)

    ! Stage 2: solve W * k2 = f(y + k1) + (c21/h) * k1
    y_stage = y + k1
    call rhs(t + gamma * h, y_stage, f_stage2, rhs_err)
    if (rhs_err /= 0) then; ierr = rhs_err; return; end if
    k2 = f_stage2 + (c21 / h) * k1
    call solve3(W_matrix, k2)

    ! Stage 3: solve W * k3 = f_stage2 + (c31/h)*k1 + (c32/h)*k2  [reuses f_stage2]
    k3 = f_stage2 + (c31 / h) * k1 + (c32 / h) * k2
    call solve3(W_matrix, k3)

    ! 3rd-order solution and embedded 2nd-order error estimate
    y_trial = y + m1 * k1 + m2 * k2 + m3 * k3
    y_err   = e1 * k1 + e2 * k2 + e3 * k3

    ! Compute worst-case scaled error
    err_max = 0.0
    do i = 1, nvar
      tol_scale = atol(i) + rtol(i) * max(abs(y(i)), abs(y_trial(i)))
      err_max   = max(err_max, abs(y_err(i)) / tol_scale)
    end do

    ! Accept or reject, then adjust step size
    call adjust_step_size(h, err_max, 1.0_dp/3.0_dp, 0.5_dp, h_new)
    if (err_max <= 1.0) then
      t = t + h
      y = y_trial
      need_jacobian = .true.
    end if
    h = h_new
  end do

  ierr = -1
  h_last = h

end subroutine ros3_integrate


! ==========================================================================
! Step size controller shared by both solvers.
!
! On accept (err_max <= 1): h_new = min(safety * h * err^(-grow_exp), grow_max * h)
! On reject (err_max >  1): h_new = max(safety * h * err^(-shrink_exp), shrink_min * h)
! ==========================================================================
pure subroutine adjust_step_size(h, err_max, grow_exp, shrink_exp, h_new)
  real(dp), intent(in)  :: h, err_max, grow_exp, shrink_exp
  real(dp), intent(out) :: h_new

  if (err_max <= 1.0) then
    if (err_max > 1.0e-30) then
      h_new = min(SAFETY * h * err_max**(-grow_exp), GROW_MAX * h)
    else
      h_new = GROW_MAX * h
    end if
  else
    h_new = max(SAFETY * h * err_max**(-shrink_exp), SHRINK_MIN * h)
  end if

end subroutine adjust_step_size


! ==========================================================================
! Direct 3x3 linear solve: A*x = b via LU with partial pivoting.
! Overwrites b with the solution x.  A is not modified.
! ==========================================================================
subroutine solve3(A_in, b)
  real(dp), intent(in)    :: A_in(3,3)
  real(dp), intent(inout) :: b(3)

  real(dp) :: A(3,3), tmp
  integer  :: piv, i, j

  A = A_in

  ! Forward elimination — column 1
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

  ! Forward elimination — column 2
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
