module DGM
  use globals
  use ode_integrators, only: ode_rhs, ode_integrate

  implicit none

  private
  public :: integrate_ODE, set_aerosol_properties, growth_jacobian, growth_rhs
  public :: kohler_equilibrium_radius, raoult_only_radius, bisection_equilibrium_radius
  public :: critical_radius
  public :: equilibrium_timescale, tau_ratio, kelvin_guard_threshold
  public :: growth_rate

  integer(i4), parameter :: nvar = 3

  ! Aerosol properties (set per-droplet by set_aerosol_properties)
  integer(i4) :: aerosol_type
  real(dp)    :: solute_mass
  real(dp)    :: solute_c7
  real(dp)    :: raoult_coeff
  real(dp)    :: flux_coeff
  real(dp)    :: r_floor
  real(dp)    :: ode_supersat

  ! Aerosol species parameters
  real(dp) :: molar_mass_solute
  real(dp) :: c7_solute
  real(dp) :: n_ions

  ! Species-specific c7 coefficients (solute density correction)
  real(dp), parameter :: c7_ammonium_sulfate = 0.4363021
  real(dp), parameter :: c7_sodium_chloride  = 0.5381062

  ! Surface tension of water
  real(dp), parameter :: sigma_w = 7.392730e-2

  ! Kinetic coefficients
  real(dp), parameter :: thermal_accom    = 1.0
  real(dp), parameter :: condensation_eff = 0.04

  ! Dry radius margin
  real(dp), parameter :: eps_r = 1.0e-2

  ! Equilibrium bypass threshold: bypass ODE when tau < dt * tau_ratio
  real(dp), parameter :: tau_ratio = 1.0

  ! Kelvin guard: skip Newton bypass for mild subsaturations (fraction, not %)
  real(dp), parameter :: kelvin_guard_threshold = -0.005

  ! ODE tolerances
  real(dp) :: ode_rtol(nvar) = [1.0e-4, 1.0e-4, 1.0e-4]
  real(dp) :: ode_atol(nvar) = [1.0e-10, 1.0e-10, 1.0e-10]


contains

subroutine set_aerosol_properties(species, mass, r_solute, grid_scale, supersat)
  integer(i4), intent(in) :: species
  real(dp), intent(in)    :: mass, r_solute, grid_scale, supersat

  aerosol_type  = species
  solute_mass   = mass
  r_floor       = r_solute * (1.0 + eps_r)
  ode_supersat  = supersat

  select case (species)
  case (1) ! NaCl
    molar_mass_solute = 58.4428e-3
    c7_solute         = c7_sodium_chloride
    n_ions            = 2.0
  case (2) ! (NH4)2SO4
    molar_mass_solute = 132.1395e-3
    c7_solute         = c7_ammonium_sulfate
    n_ions            = 3.0
  case (3) ! fumaric acid
    molar_mass_solute = 115.11e-3
    c7_solute         = c7_ammonium_sulfate
    n_ions            = 2.0
  case default
    error stop "set_aerosol_properties: unknown aerosol species"
  end select

  solute_c7    = solute_mass * c7_solute
  raoult_coeff = n_ions * (Mw / molar_mass_solute) * solute_mass
  flux_coeff   = pi_4 * grid_scale * rho_l

end subroutine set_aerosol_properties


subroutine integrate_ODE(ystart, t_start, t_end, h_init, stat, rtol, atol)
  real(dp), intent(inout) :: ystart(nvar)
  real(dp), intent(in)    :: t_start, t_end, h_init
  integer(i4), intent(out), optional :: stat
  real(dp), intent(in), optional :: rtol(nvar), atol(nvar)

  integer(i4) :: ierr
  real(dp) :: rtol_use(nvar), atol_use(nvar)

  if (present(rtol)) then
    rtol_use = rtol
  else
    rtol_use = ode_rtol
  end if
  if (present(atol)) then
    atol_use = atol
  else
    atol_use = ode_atol
  end if

  call ode_integrate(growth_rhs, nvar, ystart, t_start, t_end, h_init, &
                     rtol_use, atol_use, ierr)

  if (present(stat)) then
    stat = ierr
  else if (ierr < 0) then
    write(0,*) 'DGM ODE integration failed, ierr = ', ierr
    write(0,*) '  radius=', ystart(1), ' qv=', ystart(2), ' T=', ystart(3)
    error stop
  end if

end subroutine integrate_ODE


! RHS of the droplet growth ODE system.
! y(1) = radius, y(2) = water vapor mixing ratio, y(3) = temperature.
! Reference: Su et al. (1998), droplet growth with kinetic corrections.
subroutine growth_rhs(ltime, y, dydt)
  real(dp), intent(in)  :: ltime
  real(dp), intent(in)  :: y(:)
  real(dp), intent(out) :: dydt(:)

  real(dp) :: radius, qv, temp
  real(dp) :: es
  real(dp) :: cp_moist, Lv, K_therm, D_vapor
  real(dp) :: jump_thermal, jump_vapor
  real(dp) :: vent_thermal, vent_vapor
  real(dp) :: solution_density
  real(dp) :: kelvin, raoult, diffusion_denom

  radius = max(y(1), r_floor)
  qv     = y(2)
  temp   = y(3)

  if (qv < 0.0 .or. temp > 320.0 .or. temp < 193.0) then
    write(0, '(3(A,ES16.8))') "growth_rhs: T=", temp, " qv=", qv, " r=", radius
    error stop "growth_rhs: state out of bounds"
  end if

  cp_moist = cp * ((1.0 + cp_wv / cp * qv) / (1.0 + qv))
  Lv    = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
  K_therm  = 7.7e-5 * (temp - Tice) + 0.02399
  D_vapor  = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres

  es = esat(temp)

  ! Kinetic correction lengths (Fukuta & Walter 1970)
  jump_thermal = K_therm * sqrt(2.0 * pi * Ma * R_univ * temp) &
               / (thermal_accom * pres * (cv + R_univ / 2.0))
  jump_vapor   = sqrt(2.0 * pi * Mw / (Rv * temp)) * D_vapor / condensation_eff

  ! Ventilation coefficients
  vent_thermal = radius / (radius + jump_thermal)
  vent_vapor   = radius / (radius + jump_vapor)

  ! Solution density (accounts for dissolved solute)
  solution_density = (radius**3 * pi_43 * rho_l + solute_c7) / (radius**3 * pi_43)

  ! Köhler terms
  kelvin = 2.0 * sigma_w / (Rv * temp * solution_density * radius)
  raoult = raoult_coeff / (pi_43 * radius**3 * solution_density - solute_mass)

  ! Thermodynamic diffusion denominator
  diffusion_denom = solution_density &
    * (Rv * temp / (vent_vapor * D_vapor * es) &
     + Lv**2 / (vent_thermal * K_therm * Rv * temp**2))

  ! Radius tendency
  dydt(1) = (1.0 / radius) * (ode_supersat - kelvin + raoult) / diffusion_denom

  ! Vapor tendency (mass conservation with gridcell)
  dydt(2) = -flux_coeff * radius**2 * dydt(1)

  ! Limit evaporation to available vapor
  if (dydt(2) < 0.0 .and. abs(dydt(2)) > qv) then
    dydt(2) = -qv
    dydt(1) = -dydt(2) / (flux_coeff * radius**2)
  end if

  ! Temperature tendency (latent heating)
  dydt(3) = -Lv / cp_moist * dydt(2)

end subroutine growth_rhs


! Analytical Jacobian of growth_rhs: jac(i,j) = d(dydt_i)/d(y_j).
! Uses module-level aerosol state set by set_aerosol_properties.
! Ignores the vapor-limiter branch (non-smooth; Jacobian approximates
! the unconstrained system, which is correct whenever the limiter is inactive).
subroutine growth_jacobian(ltime, y, jac)
  real(dp), intent(in)  :: ltime
  real(dp), intent(in)  :: y(:)
  real(dp), intent(out) :: jac(:,:)

  real(dp) :: r, qv, temp
  real(dp) :: es, des_dT
  real(dp) :: cp_moist, Lv, K_therm, D_vapor
  real(dp) :: dLv_dT, dK_dT, dD_dT
  real(dp) :: jump_th, jump_vp, vent_th, vent_vp
  real(dp) :: djth_dr, djvp_dr, dvth_dr, dvvp_dr
  real(dp) :: djth_dT, djvp_dT, dvth_dT, dvvp_dT
  real(dp) :: rho_sol, drho_dr
  real(dp) :: kelvin, raoult_term, denom_water
  real(dp) :: dkelvin_dr, dkelvin_dT, draoult_dr
  real(dp) :: diff_denom, ddiff_dr, ddiff_dT
  real(dp) :: SS_eff, drdt, ddrdt_dr, ddrdt_dT
  real(dp) :: r3, r2, r4

  r    = max(y(1), r_floor)
  qv   = y(2)
  temp = y(3)

  r2 = r * r
  r3 = r2 * r
  r4 = r3 * r

  ! Thermodynamic quantities and their T-derivatives
  cp_moist = cp * ((1.0 + cp_wv / cp * qv) / (1.0 + qv))
  Lv       = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
  K_therm  = 7.7e-5 * (temp - Tice) + 0.02399
  D_vapor  = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres

  dLv_dT   = -0.00237e6
  dK_dT    = 7.7e-5
  dD_dT    = 1.57e-7 * 1.0e5 / pres

  es = esat(temp)
  des_dT = desat_dT(temp)

  ! Kinetic jump lengths and ventilation coefficients
  jump_th = K_therm * sqrt(2.0 * pi * Ma * R_univ * temp) &
          / (thermal_accom * pres * (cv + R_univ / 2.0))
  jump_vp = sqrt(2.0 * pi * Mw / (Rv * temp)) * D_vapor / condensation_eff

  vent_th = r / (r + jump_th)
  vent_vp = r / (r + jump_vp)

  ! d(jump_th)/dr = 0, d(jump_vp)/dr = 0
  dvth_dr = jump_th / (r + jump_th)**2
  dvvp_dr = jump_vp / (r + jump_vp)**2

  ! d(jump_th)/dT
  djth_dT = (dK_dT * sqrt(2.0 * pi * Ma * R_univ * temp) &
           + K_therm * pi * Ma * R_univ / sqrt(2.0 * pi * Ma * R_univ * temp)) &
           / (thermal_accom * pres * (cv + R_univ / 2.0))
  ! d(jump_vp)/dT
  djvp_dT = -0.5 * sqrt(2.0 * pi * Mw / (Rv * temp)) * D_vapor / (condensation_eff * temp) &
          + sqrt(2.0 * pi * Mw / (Rv * temp)) * dD_dT / condensation_eff

  dvth_dT = -r * djth_dT / (r + jump_th)**2
  dvvp_dT = -r * djvp_dT / (r + jump_vp)**2

  ! Solution density
  rho_sol  = (r3 * pi_43 * rho_l + solute_c7) / (r3 * pi_43)
  drho_dr  = -3.0 * solute_c7 / (pi_43 * r4)

  ! Köhler terms
  kelvin      = 2.0 * sigma_w / (Rv * temp * rho_sol * r)
  denom_water = pi_43 * r3 * rho_sol - solute_mass
  raoult_term = raoult_coeff / denom_water

  ! d(kelvin)/dr
  dkelvin_dr = -2.0 * sigma_w * (rho_sol + r * drho_dr) &
             / (Rv * temp * (rho_sol * r)**2)
  ! d(kelvin)/dT
  dkelvin_dT = -2.0 * sigma_w / (Rv * temp**2 * rho_sol * r)

  ! d(raoult)/dr: denom_water = pi_43*r3*rho_sol - m_s
  ! d(denom_water)/dr = 3*pi_43*r2*rho_sol + pi_43*r3*drho_dr
  draoult_dr = -raoult_coeff * (3.0 * pi_43 * r2 * rho_sol + pi_43 * r3 * drho_dr) &
             / denom_water**2

  ! Effective supersaturation driving growth
  SS_eff = ode_supersat - kelvin + raoult_term

  ! Diffusion denominator: D = rho_sol * (Rv*T/(fv*Dv*es) + Lv^2/(fth*K*Rv*T^2))
  diff_denom = rho_sol &
    * (Rv * temp / (vent_vp * D_vapor * es) &
     + Lv**2 / (vent_th * K_therm * Rv * temp**2))

  ! d(diff_denom)/dr
  ddiff_dr = drho_dr &
    * (Rv * temp / (vent_vp * D_vapor * es) &
     + Lv**2 / (vent_th * K_therm * Rv * temp**2)) &
    + rho_sol &
    * (-Rv * temp * dvvp_dr / (vent_vp**2 * D_vapor * es) &
     - Lv**2 * dvth_dr / (vent_th**2 * K_therm * Rv * temp**2))

  ! d(diff_denom)/dT
  ddiff_dT = rho_sol * ( &
    (Rv * vent_vp * D_vapor * es &
     - Rv * temp * (dvvp_dT * D_vapor * es + vent_vp * dD_dT * es + vent_vp * D_vapor * des_dT)) &
    / (vent_vp * D_vapor * es)**2 &
    + (2.0 * Lv * dLv_dT * vent_th * K_therm * Rv * temp**2 &
     - Lv**2 * (dvth_dT * K_therm * Rv * temp**2 &
              + vent_th * dK_dT * Rv * temp**2 &
              + vent_th * K_therm * Rv * 2.0 * temp)) &
    / (vent_th * K_therm * Rv * temp**2)**2 )

  ! dr/dt = (1/r) * SS_eff / diff_denom
  drdt = (1.0 / r) * SS_eff / diff_denom

  ! d(drdt)/dr = d/dr[(1/r) * SS_eff / D]
  !            = (-1/r²)(SS_eff/D) + (1/r)(dSS/dr)/D - (1/r)(SS_eff)(dD/dr)/D²
  ddrdt_dr = (-1.0 / r2) * SS_eff / diff_denom &
           + (1.0 / r) * (-dkelvin_dr + draoult_dr) / diff_denom &
           - (1.0 / r) * SS_eff * ddiff_dr / diff_denom**2

  ! d(drdt)/dT = (1/r) * [(-dkelvin_dT)/D - SS_eff * ddiff_dT / D²]
  ddrdt_dT = (1.0 / r) * (-dkelvin_dT / diff_denom &
           - SS_eff * ddiff_dT / diff_denom**2)

  ! --- Assemble Jacobian ---

  ! Row 1: d(drdt)/dy
  jac(1,1) = ddrdt_dr
  jac(1,2) = 0.0     ! qv does not enter drdt (SS is module-level constant)
  jac(1,3) = ddrdt_dT

  ! Row 2: d(dqv/dt)/dy where dqv/dt = -flux_coeff * r² * drdt
  ! d/dr: -flux_coeff * (2*r*drdt + r²*ddrdt_dr)
  jac(2,1) = -flux_coeff * (2.0 * r * drdt + r2 * ddrdt_dr)
  jac(2,2) = 0.0
  jac(2,3) = -flux_coeff * r2 * ddrdt_dT

  ! Row 3: d(dT/dt)/dy where dT/dt = -(Lv/cp_moist) * dqv/dt
  ! d/dr: -(Lv/cp_moist) * jac(2,1)
  jac(3,1) = -(Lv / cp_moist) * jac(2,1)
  ! d/dqv: -(Lv/cp_moist) * jac(2,2) + correction from d(cp_moist)/dqv
  ! cp_moist = cp*(1 + (cp_wv/cp)*qv)/(1+qv)
  ! d(cp_moist)/dqv = cp*(cp_wv/cp - 1)/(1+qv)² = (cp_wv - cp)/(1+qv)²
  ! dT/dt depends on cp_moist: d/dqv[-(Lv/cp_moist)*dqvdt]
  !   = (Lv * dqvdt / cp_moist²) * d(cp_moist)/dqv  (since dqvdt doesn't depend on qv)
  jac(3,2) = -(Lv / cp_moist**2) * flux_coeff * r2 * drdt &
           * (cp_wv - cp) / (1.0 + qv)**2
  jac(3,3) = -(Lv / cp_moist) * jac(2,3) &
           - (dLv_dT / cp_moist) * (-flux_coeff * r2 * drdt)

end subroutine growth_jacobian


! Derivative of saturation vapor pressure with respect to temperature [Pa/K].
! Flatau et al. (1992) polynomial derivative.
pure function desat_dT(temp) result(des)
  real(dp), intent(in) :: temp
  real(dp) :: des

  real(dp), parameter :: c(9) = [ &
    6.11239921d0, 0.443987641d0, 0.142986287d-1, &
    0.264847430d-3, 0.302950461d-5, 0.206739458d-7, &
    0.640689451d-10, -0.952447341d-13, -0.976195544d-15]
  real(dp) :: Tc
  integer  :: i

  Tc = max(-80.0, temp - Tice)
  des = 0.0
  do i = 9, 2, -1
    des = des * Tc + real(i - 1, dp) * c(i)
  end do
  des = des * 100.0

end function desat_dT


! Raoult-only equilibrium radius [m], ignoring the Kelvin (curvature) term.
! Cheap initial estimate used as the starting point for the full Newton solver.
function raoult_only_radius(supersat) result(r_eq)
  real(dp), intent(in) :: supersat  ! supersaturation [fraction]
  real(dp) :: r_eq

  r_eq = (solute_mass * (1.0 - c7_solute) + raoult_coeff / abs(supersat)) &
       / (pi_43 * rho_l)
  r_eq = r_eq**(1.0 / 3.0)

end function raoult_only_radius


! Critical radius [m] from Rogers & Yau Eqs. 6.7-6.8.
! Uses module-level aerosol state set by set_aerosol_properties.
function critical_radius(temp) result(r_crit)
  real(dp), intent(in) :: temp
  real(dp) :: r_crit

  real(dp) :: a, b

  a = a_RY / temp
  b = 4.3 * solute_mass * n_ions / molar_mass_solute

  r_crit = sqrt(3.0 * b / a) * m_per_cm

end function critical_radius


! Köhler equilibrium radius [m] on the stable (sub-critical) branch.
! Newton iteration on the full Köhler equation including solution density,
! Kelvin (curvature), and Raoult (solute) terms.
function kohler_equilibrium_radius(supersat, temp) result(r_eq)
  real(dp), intent(in) :: supersat  ! supersaturation [fraction]
  real(dp), intent(in) :: temp      ! temperature [K]
  real(dp) :: r_eq

  integer, parameter  :: max_iter = 20
  real(dp), parameter :: tol = 1.0e-30

  integer  :: iter
  real(dp) :: rho_sol, denom, kelvin, raoult_term
  real(dp) :: f, drho_dr, dk_dr, ddenom_dr, dra_dr, fprime

  r_eq = raoult_only_radius(supersat)

  do iter = 1, max_iter
    rho_sol = (r_eq**3 * pi_43 * rho_l + solute_c7) / (r_eq**3 * pi_43)
    denom   = pi_43 * r_eq**3 * rho_sol - solute_mass
    if (denom < tol) exit

    kelvin      = 2.0 * sigma_w / (Rv * temp * rho_sol * r_eq)
    raoult_term = raoult_coeff / denom
    f = supersat - kelvin + raoult_term

    ! Jacobian terms including solution density gradient
    drho_dr   = -3.0 * solute_c7 / (pi_43 * r_eq**4)
    dk_dr     = -2.0 * sigma_w * (rho_sol + r_eq * drho_dr) &
              / (Rv * temp * (rho_sol * r_eq)**2)
    ddenom_dr = 3.0 * pi_43 * r_eq**2 * rho_sol + pi_43 * r_eq**3 * drho_dr
    dra_dr    = -raoult_coeff * ddenom_dr / denom**2
    fprime    = -dk_dr + dra_dr

    if (abs(fprime) < tol) exit
    r_eq = r_eq - f / fprime
    if (r_eq < r_floor) r_eq = r_floor
  end do

end function kohler_equilibrium_radius


! Köhler equilibrium radius [m] via bisection on the stable branch.
! Guaranteed to converge where Newton may diverge (small aerosol, mild SS).
! Brackets: [r_floor, r_critical] where r_critical is the Köhler curve peak.
function bisection_equilibrium_radius(supersat, temp) result(r_eq)
  real(dp), intent(in) :: supersat  ! supersaturation [fraction]
  real(dp), intent(in) :: temp      ! temperature [K]
  real(dp) :: r_eq

  integer, parameter  :: max_iter = 60
  real(dp), parameter :: rtol = 1.0e-10

  real(dp) :: r_lo, r_hi, r_mid, f_lo, f_mid
  integer  :: iter

  ! Lower bound: dry radius (Raoult dominates, f > 0)
  r_lo = r_floor * (1.0 + eps_r)

  ! Upper bound: depends on sign of supersaturation
  if (supersat < 0.0) then
    r_hi = raoult_only_radius(supersat)
  else
    r_hi = critical_radius(temp)
  end if

  f_lo = kohler_residual(r_lo, supersat, temp)

  if (f_lo < 0.0) then
    r_eq = r_hi
    return
  end if

  ! No sign change: SS > SS_crit, droplet is activated
  if (kohler_residual(r_hi, supersat, temp) > 0.0) then
    r_eq = r_hi
    return
  end if

  do iter = 1, max_iter
    r_mid = 0.5 * (r_lo + r_hi)
    f_mid = kohler_residual(r_mid, supersat, temp)

    if (f_mid > 0.0) then
      r_lo = r_mid
    else
      r_hi = r_mid
    end if

    if ((r_hi - r_lo) < rtol * r_lo) exit
  end do

  r_eq = 0.5 * (r_lo + r_hi)

end function bisection_equilibrium_radius


! Köhler residual: f(r) = SS - kelvin(r) + raoult(r).
! Zero crossing on the stable branch gives equilibrium radius.
pure function kohler_residual(radius, supersat, temp) result(f)
  real(dp), intent(in) :: radius, supersat, temp
  real(dp) :: f

  real(dp), parameter :: tol = 1.0e-30
  real(dp) :: rho_sol, denom, kelvin, raoult_term

  rho_sol = (radius**3 * pi_43 * rho_l + solute_c7) / (radius**3 * pi_43)
  denom   = pi_43 * radius**3 * rho_sol - solute_mass
  if (denom < tol) then
    f = huge(1.0)
    return
  end if

  kelvin      = 2.0 * sigma_w / (Rv * temp * rho_sol * radius)
  raoult_term = raoult_coeff / denom
  f = supersat - kelvin + raoult_term

end function kohler_residual


! Linearized relaxation timescale [s] near Köhler equilibrium.
! Small tau means the droplet equilibrates much faster than the ODE timestep (stiff).
function equilibrium_timescale(radius, supersat, temp) result(tau)
  real(dp), intent(in) :: radius    ! droplet radius [m]
  real(dp), intent(in) :: supersat  ! supersaturation [fraction]
  real(dp), intent(in) :: temp      ! temperature [K]
  real(dp) :: tau

  real(dp), parameter :: tol = 1.0e-30

  real(dp) :: rho_sol, denom, kelvin_coeff, df_dr
  real(dp) :: es, Lv, K_therm, D_vapor, denom_thermo, eigenvalue

  rho_sol = (radius**3 * pi_43 * rho_l + solute_c7) / (radius**3 * pi_43)
  denom   = pi_43 * radius**3 * rho_sol - solute_mass
  if (denom < tol) then
    tau = 0.0
    return
  end if

  ! Köhler curve slope df/dr at current radius
  kelvin_coeff = 2.0 * sigma_w / (Rv * temp * rho_sol)
  df_dr = kelvin_coeff / radius**2 &
        - 3.0 * pi_43 * rho_sol * radius**2 * raoult_coeff / denom**2

  ! Thermodynamic diffusion resistance
  es      = esat(temp)
  Lv      = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
  K_therm = 7.7e-5 * (temp - Tice) + 0.02399
  D_vapor = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres

  denom_thermo = rho_sol * (Rv * temp / (D_vapor * es) &
               + Lv**2 / (K_therm * Rv * temp**2))

  eigenvalue = (1.0 / radius) * df_dr / denom_thermo

  if (abs(eigenvalue) < tol) then
    tau = huge(1.0)
  else
    tau = 1.0 / abs(eigenvalue)
  end if

end function equilibrium_timescale


! Radius growth rate dr/dt [m s-1] at the given state.
! Thin wrapper around the ODE RHS for diagnostic use.
! supersat and temp must match what was passed to set_aerosol_properties.
function growth_rate(radius, temp) result(drdt)
  real(dp), intent(in) :: radius, temp
  real(dp) :: drdt

  real(dp) :: y(3), dydt(3), es_loc, qv_sat, qv

  es_loc = esat(temp)
  qv_sat = 0.622 * es_loc / (pres - es_loc)
  qv     = qv_sat * (1.0 + ode_supersat)

  y(1) = radius
  y(2) = qv
  y(3) = temp
  call growth_rhs(0.0_dp, y, dydt)
  drdt = dydt(1)

end function growth_rate


! Saturation vapor pressure over liquid water [Pa].
! Flatau et al. (1992) polynomial, valid for Tc in [-80, 50] °C.
pure function esat(temp) result(es)
  real(dp), intent(in) :: temp
  real(dp) :: es

  real(dp), parameter :: c(9) = [ &
    6.11239921d0, 0.443987641d0, 0.142986287d-1, &
    0.264847430d-3, 0.302950461d-5, 0.206739458d-7, &
    0.640689451d-10, -0.952447341d-13, -0.976195544d-15]
  real(dp) :: Tc
  integer  :: i

  Tc = max(-80.0, temp - Tice)
  es = c(9)
  do i = 8, 1, -1
    es = es * Tc + c(i)
  end do
  es = es * 100.0

end function esat

end module DGM
