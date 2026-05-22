! Droplet Growth Model (DGM)
!
! Solves the coupled ODE system for single-droplet condensational growth:
!   y(1) = droplet radius        [m]
!   y(2) = water vapor mixing ratio [kg/kg]
!   y(3) = temperature           [K]
!
! Physics: Su et al. (1998) growth equation with Fukuta & Walter (1970)
! kinetic corrections and full Köhler equilibrium (Kelvin + Raoult).
!
! The module maintains per-droplet aerosol state (solute mass, species
! properties, derived coefficients) via set_aerosol_properties, which must
! be called before each integration.
module DGM
  use globals
  use ode_integrators, only: ode_rhs, ros3_integrate

  implicit none

  private
  public :: integrate_ODE, set_aerosol_properties
  public :: growth_rhs, growth_jacobian
  public :: growth_rate, critical_radius

  integer(i4), parameter :: nvar = 3

  ! ---- Per-droplet aerosol state (set by set_aerosol_properties) ----
  integer(i4) :: aerosol_type       ! species index (1=NaCl, 2=(NH4)2SO4, 3=fumaric)
  real(dp)    :: solute_mass         ! dry solute mass [kg]
  real(dp)    :: solute_c7           ! solute_mass * c7 — density correction product
  real(dp)    :: raoult_coeff        ! Raoult prefactor: n_ions * (Mw/M_s) * m_s
  real(dp)    :: flux_coeff          ! vapor flux coupling: (4/3)*pi * (1/grid_mass) * rho_l
  real(dp)    :: r_floor             ! minimum radius: dry radius * (1 + eps_r)
  real(dp)    :: ode_supersat        ! ambient supersaturation [fraction, not %]

  ! ---- Species lookup tables ----
  real(dp) :: molar_mass_solute      ! molar mass of solute [kg/mol]
  real(dp) :: c7_solute              ! solution density correction factor
  real(dp) :: n_ions                 ! van't Hoff factor (effective ions per molecule)

  ! ---- Physical constants ----
  ! Solute density corrections (c7): ratio of apparent volume occupied by dissolved
  ! solute to its dry-particle volume.  Species-specific empirical values.
  real(dp), parameter :: c7_ammonium_sulfate = 0.4363021
  real(dp), parameter :: c7_sodium_chloride  = 0.5381062

  real(dp), parameter :: sigma_w = 7.392730e-2  ! surface tension of water [N/m]

  ! Kinetic accommodation coefficients (Fukuta & Walter 1970)
  real(dp), parameter :: thermal_accom    = 1.0   ! thermal accommodation coefficient
  real(dp), parameter :: condensation_eff = 0.04  ! condensation (mass) coefficient

  real(dp), parameter :: eps_r = 1.0e-2  ! dry radius margin: r_floor = r_dry * (1 + eps_r)

  ! ---- Default ODE tolerances ----
  real(dp) :: ode_rtol(nvar) = [1.0e-4, 1.0e-4, 1.0e-4]
  real(dp) :: ode_atol(nvar) = [1.0e-10, 1.0e-10, 1.0e-10]


contains


! ==========================================================================
! Set per-droplet aerosol properties before ODE integration.
!
! Must be called once per droplet before integrate_ODE.  Populates module
! state: solute mass, species-dependent constants, and derived coefficients
! used by growth_rhs and growth_jacobian.
!
! Arguments:
!   species    — aerosol type: 1=NaCl, 2=(NH4)2SO4, 3=fumaric acid
!   mass       — dry solute mass [kg]
!   r_solute   — dry solute radius [m]
!   grid_scale — 1 / (gridcell air mass) [kg^-1]
!   supersat   — ambient supersaturation [fraction, not %]
! ==========================================================================
subroutine set_aerosol_properties(species, mass, r_solute, grid_scale, supersat)
  integer(i4), intent(in) :: species
  real(dp), intent(in)    :: mass, r_solute, grid_scale, supersat

  aerosol_type = species
  solute_mass  = mass
  r_floor      = r_solute * (1.0 + eps_r)
  ode_supersat = supersat

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

  ! Precompute derived coefficients used by the RHS
  solute_c7    = solute_mass * c7_solute
  raoult_coeff = n_ions * (Mw / molar_mass_solute) * solute_mass
  flux_coeff   = pi_4 * grid_scale * rho_l

end subroutine set_aerosol_properties


! ==========================================================================
! Integrate the droplet growth ODE from t_start to t_end.
!
! Wraps ros3_integrate with default tolerances and error handling.
! Calls set_aerosol_properties first to configure module state.
!
! Arguments:
!   y        — state vector [radius, qv, T], overwritten with solution
!   t_start  — integration start time [s]
!   t_end    — integration end time [s]
!   h_last   — step size hint on entry, last accepted step on exit [s]
!   stat     — (optional) error code; if absent, failure causes error stop
!   rtol     — (optional) relative tolerance override per component
!   atol     — (optional) absolute tolerance override per component
! ==========================================================================
subroutine integrate_ODE(y, t_start, t_end, h_last, stat, rtol, atol)
  real(dp), intent(inout) :: y(nvar)
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

  call ros3_integrate(growth_rhs, growth_jacobian, nvar, y, t_start, t_end, &
                      h_last, rtol_use, atol_use, ierr)

  if (present(stat)) then
    stat = ierr
  else if (ierr < 0) then
    write(0,*) 'DGM ODE integration failed, ierr = ', ierr
    write(0,*) '  radius=', y(1), ' qv=', y(2), ' T=', y(3)
    error stop
  end if

end subroutine integrate_ODE


! ==========================================================================
! RHS of the droplet growth ODE: dydt = f(t, y).
!
! Computes radius, vapor, and temperature tendencies from the Su et al.
! (1998) growth equation with Fukuta & Walter (1970) kinetic corrections.
!
! State vector:
!   y(1) = radius [m],  y(2) = mixing ratio [kg/kg],  y(3) = temperature [K]
!
! The supersaturation driving growth comes from the module-level ode_supersat,
! not from y(2).  y(2) evolves via mass conservation with the gridcell.
! ==========================================================================
pure subroutine growth_rhs(t, y, dydt, ierr)
  real(dp),    intent(in)  :: t
  real(dp),    intent(in)  :: y(:)
  real(dp),    intent(out) :: dydt(:)
  integer(i4), intent(out) :: ierr

  real(dp) :: radius, qv, temp
  real(dp) :: es
  real(dp) :: cp_moist, latent_heat, thermal_cond, vapor_diff
  real(dp) :: jump_thermal, jump_vapor
  real(dp) :: vent_thermal, vent_vapor
  real(dp) :: solution_density
  real(dp) :: kelvin, raoult, diffusion_denom

  radius = max(y(1), r_floor)
  qv     = y(2)
  temp   = y(3)

  ! Guard against unphysical state
  if (qv < 0.0 .or. temp > 320.0 .or. temp < 193.0) then
    ierr = 1
    dydt = 0.0
    return
  end if

  ierr = 0

  ! Thermodynamic properties
  cp_moist    = cp * ((1.0 + cp_wv / cp * qv) / (1.0 + qv))
  latent_heat = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6        ! [J/kg]
  thermal_cond = 7.7e-5 * (temp - Tice) + 0.02399                 ! [W/(m·K)]
  vapor_diff  = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres  ! [m²/s]

  es = esat(temp)

  ! Kinetic correction: mean free path jump lengths (Fukuta & Walter 1970)
  jump_thermal = thermal_cond * sqrt(2.0 * pi * Ma * R_univ * temp) &
               / (thermal_accom * pres * (cv + R_univ / 2.0))
  jump_vapor   = sqrt(2.0 * pi * Mw / (Rv * temp)) * vapor_diff / condensation_eff

  ! Ventilation coefficients (transition regime correction)
  vent_thermal = radius / (radius + jump_thermal)
  vent_vapor   = radius / (radius + jump_vapor)

  ! Solution density: accounts for dissolved solute volume
  solution_density = (radius**3 * pi_43 * rho_l + solute_c7) / (radius**3 * pi_43)

  ! Köhler equilibrium terms
  kelvin = 2.0 * sigma_w / (Rv * temp * solution_density * radius)
  raoult = raoult_coeff / (pi_43 * radius**3 * solution_density - solute_mass)

  ! Thermodynamic diffusion resistance
  diffusion_denom = solution_density &
    * (Rv * temp / (vent_vapor * vapor_diff * es) &
     + latent_heat**2 / (vent_thermal * thermal_cond * Rv * temp**2))

  ! Radius tendency: dr/dt = (1/r) * (S - Kelvin + Raoult) / D
  dydt(1) = (1.0 / radius) * (ode_supersat - kelvin + raoult) / diffusion_denom

  ! Vapor tendency: mass conservation with gridcell
  dydt(2) = -flux_coeff * radius**2 * dydt(1)

  ! Limit evaporation to available vapor
  if (dydt(2) < 0.0 .and. abs(dydt(2)) > qv) then
    dydt(2) = -qv
    dydt(1) = -dydt(2) / (flux_coeff * radius**2)
  end if

  ! Temperature tendency: latent heating
  dydt(3) = -latent_heat / cp_moist * dydt(2)

end subroutine growth_rhs


! ==========================================================================
! Analytical Jacobian of growth_rhs: jac(i,j) = d(dydt_i)/d(y_j).
!
! Differentiates the unconstrained growth equations (ignores the vapor
! limiter branch, which is non-smooth; the Jacobian is correct whenever
! the limiter is inactive).
!
! The core quantity is dr/dt = (1/r) * SS_eff / D, where:
!   SS_eff = S_ambient - Kelvin(r,T) + Raoult(r)
!   D      = rho_sol * [Rv*T/(f_v*Dv*es) + Lv²/(f_th*K*Rv*T²)]
!
! qv does not appear in dr/dt (supersaturation is module-level), so
! d(dr/dt)/dqv = 0.  The only nonzero qv derivative is in dT/dt via
! the moist heat capacity cp_moist(qv).
!
! Uses module-level aerosol state set by set_aerosol_properties.
! ==========================================================================
pure subroutine growth_jacobian(t, y, jac)
  real(dp), intent(in)  :: t
  real(dp), intent(in)  :: y(:)
  real(dp), intent(out) :: jac(:,:)

  ! State
  real(dp) :: r, qv, temp
  real(dp) :: r2, r3, r4

  ! Base quantities (shared by forward and derivative calculations)
  real(dp) :: es, cp_moist, latent_heat, thermal_cond, vapor_diff
  real(dp) :: jump_th, jump_vp, vent_th, vent_vp
  real(dp) :: rho_sol, kelvin_term, raoult_term, denom_water
  real(dp) :: ss_eff, diff_denom, drdt

  ! --- Derivatives w.r.t. radius ---
  real(dp) :: d_rho_dr                        ! d(rho_sol)/dr
  real(dp) :: d_vth_dr, d_vvp_dr              ! d(ventilation)/dr
  real(dp) :: d_kelvin_dr, d_raoult_dr        ! d(Köhler terms)/dr
  real(dp) :: d_diff_dr                       ! d(diffusion denom)/dr
  real(dp) :: d_drdt_dr                       ! d(dr/dt)/dr

  ! --- Derivatives w.r.t. temperature ---
  real(dp) :: des_dt_val                      ! d(es)/dT
  real(dp) :: d_lv_dt, d_k_dt, d_dv_dt       ! d(Lv,K,Dv)/dT
  real(dp) :: d_jth_dt, d_jvp_dt              ! d(jump lengths)/dT
  real(dp) :: d_vth_dt, d_vvp_dt              ! d(ventilation)/dT
  real(dp) :: d_kelvin_dt                     ! d(Kelvin)/dT
  real(dp) :: d_diff_dt                       ! d(diffusion denom)/dT
  real(dp) :: d_drdt_dt                       ! d(dr/dt)/dT

  ! =====================================================================
  ! Unpack state
  ! =====================================================================
  r    = max(y(1), r_floor)
  qv   = y(2)
  temp = y(3)

  r2 = r * r
  r3 = r2 * r
  r4 = r3 * r

  ! =====================================================================
  ! Base quantities (needed by both forward evaluation and derivatives)
  ! =====================================================================

  ! Thermodynamic properties
  cp_moist     = cp * ((1.0 + cp_wv / cp * qv) / (1.0 + qv))
  latent_heat  = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
  thermal_cond = 7.7e-5 * (temp - Tice) + 0.02399
  vapor_diff   = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres
  es = esat(temp)

  ! Kinetic jump lengths (Fukuta & Walter 1970)
  jump_th = thermal_cond * sqrt(2.0 * pi * Ma * R_univ * temp) &
          / (thermal_accom * pres * (cv + R_univ / 2.0))
  jump_vp = sqrt(2.0 * pi * Mw / (Rv * temp)) * vapor_diff / condensation_eff

  ! Ventilation coefficients
  vent_th = r / (r + jump_th)
  vent_vp = r / (r + jump_vp)

  ! Solution density
  rho_sol = (r3 * pi_43 * rho_l + solute_c7) / (r3 * pi_43)

  ! Köhler terms
  kelvin_term = 2.0 * sigma_w / (Rv * temp * rho_sol * r)
  denom_water = pi_43 * r3 * rho_sol - solute_mass
  raoult_term = raoult_coeff / denom_water

  ! Effective supersaturation and diffusion denominator
  ss_eff = ode_supersat - kelvin_term + raoult_term
  diff_denom = rho_sol &
    * (Rv * temp / (vent_vp * vapor_diff * es) &
     + latent_heat**2 / (vent_th * thermal_cond * Rv * temp**2))

  ! Forward growth rate
  drdt = (1.0 / r) * ss_eff / diff_denom

  ! =====================================================================
  ! Derivatives w.r.t. radius (r)
  !
  ! Radius enters: rho_sol(r), ventilation(r), Kelvin(r), Raoult(r),
  ! diffusion_denom(r), and the 1/r prefactor.
  ! =====================================================================

  ! d(rho_sol)/dr
  d_rho_dr = -3.0 * solute_c7 / (pi_43 * r4)

  ! d(ventilation)/dr — jump lengths are independent of r
  d_vth_dr = jump_th / (r + jump_th)**2
  d_vvp_dr = jump_vp / (r + jump_vp)**2

  ! d(Kelvin)/dr
  d_kelvin_dr = -2.0 * sigma_w * (rho_sol + r * d_rho_dr) &
              / (Rv * temp * (rho_sol * r)**2)

  ! d(Raoult)/dr via chain rule on denom_water(r)
  d_raoult_dr = -raoult_coeff * (3.0 * pi_43 * r2 * rho_sol + pi_43 * r3 * d_rho_dr) &
              / denom_water**2

  ! d(diff_denom)/dr: contributions from rho_sol(r) and ventilation(r)
  d_diff_dr = d_rho_dr &
    * (Rv * temp / (vent_vp * vapor_diff * es) &
     + latent_heat**2 / (vent_th * thermal_cond * Rv * temp**2)) &
    + rho_sol &
    * (-Rv * temp * d_vvp_dr / (vent_vp**2 * vapor_diff * es) &
     - latent_heat**2 * d_vth_dr / (vent_th**2 * thermal_cond * Rv * temp**2))

  ! d(dr/dt)/dr = d/dr[(1/r) * SS_eff / D]
  !             = -(1/r²)(SS/D) + (1/r)(dSS/dr)/D - (1/r)(SS)(dD/dr)/D²
  d_drdt_dr = (-1.0 / r2) * ss_eff / diff_denom &
            + (1.0 / r) * (-d_kelvin_dr + d_raoult_dr) / diff_denom &
            - (1.0 / r) * ss_eff * d_diff_dr / diff_denom**2

  ! =====================================================================
  ! Derivatives w.r.t. temperature (T)
  !
  ! Temperature enters: es(T), Lv(T), K(T), Dv(T), jump lengths(T),
  ! ventilation(T via jumps), Kelvin(T), and diffusion_denom(T).
  ! Raoult has no T dependence.
  ! =====================================================================

  ! d(es)/dT
  des_dt_val = desat_dT(temp)

  ! d(thermodynamic properties)/dT — linear approximations
  d_lv_dt = -0.00237e6
  d_k_dt  = 7.7e-5
  d_dv_dt = 1.57e-7 * 1.0e5 / pres

  ! d(jump lengths)/dT
  d_jth_dt = (d_k_dt * sqrt(2.0 * pi * Ma * R_univ * temp) &
           + thermal_cond * pi * Ma * R_univ / sqrt(2.0 * pi * Ma * R_univ * temp)) &
           / (thermal_accom * pres * (cv + R_univ / 2.0))
  d_jvp_dt = -0.5 * sqrt(2.0 * pi * Mw / (Rv * temp)) * vapor_diff / (condensation_eff * temp) &
           + sqrt(2.0 * pi * Mw / (Rv * temp)) * d_dv_dt / condensation_eff

  ! d(ventilation)/dT via jump length dependence
  d_vth_dt = -r * d_jth_dt / (r + jump_th)**2
  d_vvp_dt = -r * d_jvp_dt / (r + jump_vp)**2

  ! d(Kelvin)/dT — only the 1/T factor
  d_kelvin_dt = -2.0 * sigma_w / (Rv * temp**2 * rho_sol * r)

  ! d(diff_denom)/dT: contributions from es(T), Lv(T), K(T), Dv(T), ventilation(T)
  d_diff_dt = rho_sol * ( &
    (Rv * vent_vp * vapor_diff * es &
     - Rv * temp * (d_vvp_dt * vapor_diff * es + vent_vp * d_dv_dt * es &
                  + vent_vp * vapor_diff * des_dt_val)) &
    / (vent_vp * vapor_diff * es)**2 &
    + (2.0 * latent_heat * d_lv_dt * vent_th * thermal_cond * Rv * temp**2 &
     - latent_heat**2 * (d_vth_dt * thermal_cond * Rv * temp**2 &
              + vent_th * d_k_dt * Rv * temp**2 &
              + vent_th * thermal_cond * Rv * 2.0 * temp)) &
    / (vent_th * thermal_cond * Rv * temp**2)**2 )

  ! d(dr/dt)/dT = (1/r) * [(-dKelvin/dT)/D - SS * (dD/dT)/D²]
  d_drdt_dt = (1.0 / r) * (-d_kelvin_dt / diff_denom &
            - ss_eff * d_diff_dt / diff_denom**2)

  ! =====================================================================
  ! Assemble Jacobian
  !
  ! The ODE system is:
  !   f1 = dr/dt  = (1/r) * SS_eff / D
  !   f2 = dqv/dt = -flux_coeff * r² * dr/dt
  !   f3 = dT/dt  = -(Lv/cp_moist) * dqv/dt
  !
  ! Column 1: d/dr    — all three equations depend on r
  ! Column 2: d/dqv   — only dT/dt depends on qv (via cp_moist)
  ! Column 3: d/dT    — all three equations depend on T
  ! =====================================================================

  ! Column 1: d(f)/dr
  jac(1,1) = d_drdt_dr
  jac(2,1) = -flux_coeff * (2.0 * r * drdt + r2 * d_drdt_dr)
  jac(3,1) = -(latent_heat / cp_moist) * jac(2,1)

  ! Column 2: d(f)/dqv
  ! dr/dt and dqv/dt do not depend on qv.
  ! dT/dt depends on qv through cp_moist:
  !   d(cp_moist)/dqv = (cp_wv - cp) / (1 + qv)²
  jac(1,2) = 0.0
  jac(2,2) = 0.0
  jac(3,2) = -(latent_heat / cp_moist**2) * flux_coeff * r2 * drdt &
           * (cp_wv - cp) / (1.0 + qv)**2

  ! Column 3: d(f)/dT
  jac(1,3) = d_drdt_dt
  jac(2,3) = -flux_coeff * r2 * d_drdt_dt
  jac(3,3) = -(latent_heat / cp_moist) * jac(2,3) &
           - (d_lv_dt / cp_moist) * (-flux_coeff * r2 * drdt)

end subroutine growth_jacobian


! ==========================================================================
! Derivative of saturation vapor pressure w.r.t. temperature [Pa/K].
! Flatau et al. (1992) polynomial derivative.
! ==========================================================================
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


! ==========================================================================
! Critical radius [m] from Rogers & Yau Eqs. 6.7-6.8.
! Uses module-level aerosol state set by set_aerosol_properties.
! ==========================================================================
function critical_radius(temp) result(r_crit)
  real(dp), intent(in) :: temp
  real(dp) :: r_crit

  real(dp) :: a, b

  a = a_RY / temp
  b = 4.3 * solute_mass * n_ions / molar_mass_solute

  r_crit = sqrt(3.0 * b / a) * m_per_cm

end function critical_radius


! ==========================================================================
! Radius growth rate dr/dt [m/s] at the given state.
! Thin wrapper around growth_rhs for diagnostic use.
! set_aerosol_properties must be called first.
! ==========================================================================
function growth_rate(radius, temp) result(drdt)
  real(dp), intent(in) :: radius, temp
  real(dp) :: drdt

  real(dp) :: y(nvar), dydt(nvar), es_loc, qv_sat, qv
  integer(i4) :: rhs_err

  es_loc = esat(temp)
  qv_sat = 0.622 * es_loc / (pres - es_loc)
  qv     = qv_sat * (1.0 + ode_supersat)

  y(1) = radius
  y(2) = qv
  y(3) = temp
  call growth_rhs(0.0_dp, y, dydt, rhs_err)
  drdt = dydt(1)

end function growth_rate


! ==========================================================================
! Saturation vapor pressure over liquid water [Pa].
! Flatau et al. (1992) polynomial, valid for Tc in [-80, 50] °C.
! ==========================================================================
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
