module DGM
  use globals
  use ode_integrators, only: ode_rhs, ode_integrate
  implicit none

  private
  public :: integrate_ODE, set_aerosol_properties
  public :: kohler_equilibrium_radius, equilibrium_timescale, tau_ratio

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
  real(dp), parameter :: tau_ratio = 0.01

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


subroutine integrate_ODE(ystart, t_start, t_end, h_init)
  real(dp), intent(inout) :: ystart(nvar)
  real(dp), intent(in)    :: t_start, t_end, h_init

  integer(i4) :: ierr

  call ode_integrate(growth_rhs, nvar, ystart, t_start, t_end, h_init, &
                     ode_rtol, ode_atol, ierr)

  if (ierr < 0) then
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


! Analytical Köhler equilibrium radius, neglecting the Kelvin term.
! Valid when S < 0 (subsaturated); Kelvin correction is ~3% for 50 nm aerosol.
function kohler_equilibrium_radius(supersat) result(r_eq)
  real(dp), intent(in) :: supersat
  real(dp) :: r_eq

  r_eq = (solute_mass * (1.0 - c7_solute) + raoult_coeff / abs(supersat)) &
       / (pi_43 * rho_l)
  r_eq = r_eq**(1.0 / 3.0)

end function kohler_equilibrium_radius


! Linearized relaxation timescale near Köhler equilibrium.
! Small tau means the droplet equilibrates much faster than the ODE timestep (stiff).
function equilibrium_timescale(radius, supersat, temp) result(tau)
  real(dp), intent(in) :: radius, supersat, temp
  real(dp) :: tau

  real(dp) :: D_water, df_dr, denom_thermo
  real(dp) :: kelvin_coeff, es, Lv, K_therm, D_vapor, eigenvalue

  kelvin_coeff = 2.0 * sigma_w / (Rv * temp * rho_l)

  D_water = pi_43 * radius**3 * rho_l - solute_mass * (1.0 - c7_solute)
  if (D_water < 1.0e-30) then
    tau = 0.0
    return
  end if

  ! df/dr at current radius (Köhler curve slope)
  df_dr = kelvin_coeff / radius**2 &
        - 3.0 * pi_43 * rho_l * radius**2 * raoult_coeff / D_water**2

  ! Thermodynamic resistance
  es = esat(temp)

  Lv   = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
  K_therm = 7.7e-5 * (temp - Tice) + 0.02399
  D_vapor = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres

  denom_thermo = rho_l * (Rv * temp / (D_vapor * es) &
               + Lv**2 / (K_therm * Rv * temp**2))

  eigenvalue = (1.0 / radius) * df_dr / denom_thermo

  if (abs(eigenvalue) < 1.0e-30) then
    tau = huge(1.0)
  else
    tau = 1.0 / abs(eigenvalue)
  end if

end function equilibrium_timescale

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
