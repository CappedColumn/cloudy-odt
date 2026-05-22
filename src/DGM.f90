module DGM
  use globals
  use ode_integrators, only: ode_rhs
  use rosenbrock, only: ros3_integrate

  implicit none

  private
  public :: integrate_ODE, set_aerosol_properties, growth_jacobian, growth_rhs
  public :: critical_radius
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


subroutine integrate_ODE(ystart, t_start, t_end, h_last, stat, rtol, atol)
  real(dp), intent(inout) :: ystart(nvar)
  real(dp), intent(in)    :: t_start, t_end
  real(dp), intent(inout) :: h_last
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

  call ros3_integrate(growth_rhs, growth_jacobian, nvar, ystart, t_start, t_end, &
                      h_last, rtol_use, atol_use, ierr)

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
pure subroutine growth_rhs(ltime, y, dydt, ierr)
  real(dp),    intent(in)  :: ltime
  real(dp),    intent(in)  :: y(:)
  real(dp),    intent(out) :: dydt(:)
  integer(i4), intent(out) :: ierr

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
    ierr = 1
    dydt = 0.0
    return
  end if

  ierr = 0

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
pure subroutine growth_jacobian(ltime, y, jac)
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


! Radius growth rate dr/dt [m s-1] at the given state.
! Thin wrapper around the ODE RHS for diagnostic use.
! supersat and temp must match what was passed to set_aerosol_properties.
function growth_rate(radius, temp) result(drdt)
  real(dp), intent(in) :: radius, temp
  real(dp) :: drdt

  real(dp) :: y(3), dydt(3), es_loc, qv_sat, qv
  integer(i4) :: ierr_rhs

  es_loc = esat(temp)
  qv_sat = 0.622 * es_loc / (pres - es_loc)
  qv     = qv_sat * (1.0 + ode_supersat)

  y(1) = radius
  y(2) = qv
  y(3) = temp
  call growth_rhs(0.0_dp, y, dydt, ierr_rhs)
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
