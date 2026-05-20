module lsoda_integrator
  use globals
  use microphysics, only: saturation_vapor_pressure
  use odepack_mod, only: lsoda_class
  implicit none

  private
  public :: lsoda_solver

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

  integer, parameter :: nvar = 3

  type, extends(lsoda_class) :: lsoda_solver
    integer(i4) :: aerosol_type  = 0
    real(dp)    :: solute_mass   = 0.0
    real(dp)    :: solute_c7     = 0.0
    real(dp)    :: raoult_coeff  = 0.0
    real(dp)    :: flux_coeff    = 0.0
    real(dp)    :: r_floor       = 0.0
    real(dp)    :: ode_supersat  = 0.0
    real(dp)    :: molar_mass_solute = 0.0
    real(dp)    :: c7_solute     = 0.0
    real(dp)    :: n_ions        = 0.0
    logical     :: initialized   = .false.
    logical     :: warmed_up     = .false.
    real(dp)    :: global_time   = 0.0
  contains
    procedure :: set_properties => solver_set_properties
    procedure :: solve => solver_solve
  end type

contains


subroutine solver_set_properties(self, species, mass, r_solute, grid_scale, supersat)
  class(lsoda_solver), intent(inout) :: self
  integer(i4), intent(in) :: species
  real(dp), intent(in)    :: mass, r_solute, grid_scale, supersat

  self%aerosol_type = species
  self%solute_mass  = mass
  self%r_floor      = r_solute * (1.0 + eps_r)
  self%ode_supersat = supersat

  select case (species)
  case (1) ! NaCl
    self%molar_mass_solute = 58.4428e-3
    self%c7_solute         = c7_sodium_chloride
    self%n_ions            = 2.0
  case (2) ! (NH4)2SO4
    self%molar_mass_solute = 132.1395e-3
    self%c7_solute         = c7_ammonium_sulfate
    self%n_ions            = 3.0
  case (3) ! fumaric acid
    self%molar_mass_solute = 115.11e-3
    self%c7_solute         = c7_ammonium_sulfate
    self%n_ions            = 2.0
  case default
    error stop "solver_set_properties: unknown aerosol species"
  end select

  self%solute_c7    = self%solute_mass * self%c7_solute
  self%raoult_coeff = self%n_ions * (Mw / self%molar_mass_solute) * self%solute_mass
  self%flux_coeff   = pi_4 * grid_scale * rho_l

end subroutine solver_set_properties


subroutine solver_solve(self, ystart, t_start, t_end, stat)
  class(lsoda_solver), intent(inout) :: self
  real(dp), intent(inout) :: ystart(nvar)
  real(dp), intent(in)    :: t_start, t_end
  integer(i4), intent(out), optional :: stat

  integer :: istate
  real(dp) :: t, dt
  real(dp), parameter :: rtol = 1.0e-4
  real(dp), parameter :: atol(nvar) = [1.0e-10, 1.0e-10, 1.0e-10]

  if (.not. self%initialized) then
    call self%initialize(growth_rhs, nvar, istate=istate, iprint=0)
    if (istate < 0) then
      if (present(stat)) then
        stat = istate
        return
      else
        write(0,*) 'LSODA initialization failed: ', self%error_message
        error stop
      end if
    end if
    self%initialized = .true.
  end if

  dt = t_end - t_start
  t = self%global_time

  if (self%warmed_up) then
    istate = 3
  else
    istate = 1
    self%warmed_up = .true.
  end if

  call self%integrate(ystart, t, t + dt, rtol, atol, 1, istate)

  self%global_time = t

  if (istate < 0) then
    if (istate == -1) then
      istate = 1
      self%warmed_up = .false.
      call self%integrate(ystart, t, t + dt, rtol, atol, 1, istate)
      if (istate > 0) then
        self%warmed_up = .true.
        self%global_time = t
      end if
    end if
    if (istate < 0) then
      if (present(stat)) then
        stat = istate
        return
      else
        write(0,*) 'LSODA integration failed, istate = ', istate
        write(0,*) '  ', self%error_message
        write(0,*) '  radius=', ystart(1), ' qv=', ystart(2), ' T=', ystart(3)
        error stop
      end if
    end if
  end if

  if (present(stat)) stat = 0

end subroutine solver_solve


subroutine growth_rhs(self, neq, t, y, ydot, ierr)
  class(lsoda_class), intent(inout) :: self
  integer, intent(in)  :: neq
  real(dp), intent(in) :: t
  real(dp), intent(in) :: y(neq)
  real(dp), intent(out) :: ydot(neq)
  integer, intent(out) :: ierr

  real(dp) :: radius, qv, temp
  real(dp) :: es
  real(dp) :: cp_moist, Lv, K_therm, D_vapor
  real(dp) :: jump_thermal, jump_vapor
  real(dp) :: vent_thermal, vent_vapor
  real(dp) :: solution_density
  real(dp) :: kelvin, raoult, diffusion_denom

  select type (s => self)
  type is (lsoda_solver)

    radius = max(y(1), s%r_floor)
    qv     = y(2)
    temp   = y(3)

    if (qv < 0.0 .or. temp > 320.0 .or. temp < 193.0) then
      ierr = -1
      return
    end if

    cp_moist = cp * ((1.0 + cp_wv / cp * qv) / (1.0 + qv))
    Lv       = (2.501 - 0.00237 * (temp - Tice)) * 1.0e6
    K_therm  = 7.7e-5 * (temp - Tice) + 0.02399
    D_vapor  = (1.57e-7 * (temp - Tice) + 2.211e-5) * 1.0e5 / pres

    es = saturation_vapor_pressure(temp)

    jump_thermal = K_therm * sqrt(2.0 * pi * Ma * R_univ * temp) &
                 / (thermal_accom * pres * (cv + R_univ / 2.0))
    jump_vapor   = sqrt(2.0 * pi * Mw / (Rv * temp)) * D_vapor / condensation_eff

    vent_thermal = radius / (radius + jump_thermal)
    vent_vapor   = radius / (radius + jump_vapor)

    solution_density = (radius**3 * pi_43 * rho_l + s%solute_c7) / (radius**3 * pi_43)

    kelvin = 2.0 * sigma_w / (Rv * temp * solution_density * radius)
    raoult = s%raoult_coeff / (pi_43 * radius**3 * solution_density - s%solute_mass)

    diffusion_denom = solution_density &
      * (Rv * temp / (vent_vapor * D_vapor * es) &
       + Lv**2 / (vent_thermal * K_therm * Rv * temp**2))

    ydot(1) = (1.0 / radius) * (s%ode_supersat - kelvin + raoult) / diffusion_denom

    ydot(2) = -s%flux_coeff * radius**2 * ydot(1)

    if (ydot(2) < 0.0 .and. abs(ydot(2)) > qv) then
      ydot(2) = -qv
      ydot(1) = -ydot(2) / (s%flux_coeff * radius**2)
    end if

    ydot(3) = -Lv / cp_moist * ydot(2)

    ierr = 0

  class default
    ierr = -1
  end select

end subroutine growth_rhs

end module lsoda_integrator
