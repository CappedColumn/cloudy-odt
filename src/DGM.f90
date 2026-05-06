module DGM
    use globals

    use ode_integrators, only: ode_rhs, ode_integrate
    implicit none

    integer(i4), parameter :: nvar = 4

    integer(i4)            :: dmaxa
    real(dp)               :: grid_scale, solute_mass
    real(dp)               :: r_floor
    real(dp), parameter    :: eps_r = 1.0e-2_dp

    real(dp)               :: Ms_sp
    real(dp)               :: c7_sp
    real(dp)               :: nions_sp

    real(dp)               :: solute_c7
    real(dp)               :: raoult_coeff
    real(dp)               :: flux_coeff
    real(dp)               :: inv_grid_scale

    real(dp)               :: ode_supersat
    real(dp)               :: ode_press

    real(dp), parameter :: c7am   = 0.4363021
    real(dp), parameter :: c7nacl = 0.5381062

    real(dp) :: ode_rtol(nvar) = [1.0e-4_dp, 1.0e-4_dp, 1.0e-4_dp, 1.0e-4_dp]
    real(dp) :: ode_atol(nvar) = [1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]

    public :: integrate_ODE, set_aerosol_properties, ode_supersat, ode_press
    private

contains

subroutine set_aerosol_properties(dmax, m0_aerosol, r0_solute, gscale)
    integer(i4), intent(in) :: dmax
    real(dp), intent(in) :: m0_aerosol, r0_solute, gscale

    dmaxa = dmax
    solute_mass = m0_aerosol
    grid_scale = gscale
    r_floor = r0_solute * (1.0_dp + eps_r)

    select case (dmax)
    case (1)
        Ms_sp    = 58.4428e-3
        c7_sp    = c7nacl
        nions_sp = 2.0
    case (2)
        Ms_sp    = 132.1395e-3
        c7_sp    = c7am
        nions_sp = 3.0
    case (3)
        Ms_sp    = 115.11e-3
        c7_sp    = c7am
        nions_sp = 2.0
    case default
        write(*,*) "set_aerosol_properties: Error: aerosol type unknown"
        stop 1
    end select

    solute_c7      = solute_mass * c7_sp
    raoult_coeff   = nions_sp * (Mw / Ms_sp) * solute_mass
    flux_coeff     = pi_4 * grid_scale * rho_l
    inv_grid_scale = 1.0_dp / grid_scale

end subroutine set_aerosol_properties


subroutine integrate_ODE(ystart, x1, x2, h1)
    real(dp), intent(inout) :: ystart(nvar)
    real(dp), intent(in) :: x1, x2, h1

    integer(i4) :: ierr

    call ode_integrate(fcnkb, nvar, ystart, x1, x2, h1, ode_rtol, ode_atol, ierr)

    if (ierr < 0) then
        write(0,*) 'DGM ODE integration failed, ierr = ', ierr
        write(0,*) '  radius=', ystart(1), ' qv=', ystart(2), ' T=', ystart(3)
        stop 1
    end if

end subroutine integrate_ODE


subroutine fcnkb(ltime, y, dydt)
    real(dp), intent(in)  :: ltime
    real(dp), intent(in)  :: y(:)
    real(dp), intent(out) :: dydt(:)

    real(dp) :: ck, cr, denom
    real(dp) :: es, Tc_es
    real(dp) :: falpha, fbeta, rhol
    real(dp) :: radius, qv, temp, s, press
    real(dp) :: Lcond_temp, Ktemp, D, cpm
    real(dp) :: lalpha, lbeta
    real(dp), parameter :: alph = 1
    real(dp), parameter :: beta = 0.04
    real(dp), parameter :: sigma = 7.392730e-2
    ! Manual inlining of saturation_vapor_pressure() for speedup
    ! due to difficulty with cross-module inlining with some compilers
    double precision, parameter :: es_coeff(9) = [6.11239921d0, 0.443987641d0, 0.142986287d-1, &
        0.264847430d-3, 0.302950461d-5, 0.206739458d-7, 0.640689451d-10, &
        -0.952447341d-13, -0.976195544d-15]
    integer :: i_es

    radius = y(1)
    if (radius < r_floor) radius = r_floor
    qv     = y(2)
    temp   = y(3)
    s      = ode_supersat
    press  = ode_press

    cpm = cp*((1.0+cp_wv/cp*qv)/(1.0+qv))

    Lcond_temp = (2.501-0.00237*(temp-Tice))*1.e6
    Ktemp = 7.7e-5*(temp-Tice)+0.02399
    D     = 1.57e-7*(temp-Tice)+2.211e-5
    D     = D*1.e5/press

    Tc_es = max(-80.0_dp, temp - Tice)
    es = es_coeff(9)
    do i_es = 8, 1, -1
        es = es * Tc_es + es_coeff(i_es)
    end do
    es = es * 100.0_dp

    if (temp < 273 .or. qv < 0 .or. temp > 320) then
      write(*,*) ""
      write(*,'(3(A,E16.8))') "fcnb error: T: ", temp, " qv: ", qv,  " radius: " , radius
      write(*,*) "fcnb error: Ktemp: ", Ktemp
      write(*,*) "fcnb error: 2.0*pi*Md*R_uni : ", 2.0*pi*Ma*1e-3*R_univ
      flush(6)
      stop 1
    end if

    lalpha = Ktemp * SQRT(2.0*pi*Ma*R_univ*temp)/(alph*press*(cv+R_univ/2.0))
    lbeta  = SQRT(2.0*pi*Mw/(Rv*temp))*D/beta

    falpha = radius/(radius+lalpha)
    fbeta  = radius/(radius+lbeta)
    rhol   = (radius**3*pi_43*rho_l+solute_c7)/(radius**3*pi_43)

    ck    = (2*sigma/(Rv*temp*rhol*radius))
    cr    = raoult_coeff/(pi_43*radius**3*rhol-solute_mass)
    denom = rhol*(Rv*temp/(fbeta*D*es)+Lcond_temp**2/(falpha*Ktemp*Rv*temp**2))
    dydt(1) = 1.0/radius*(s-ck+cr)/denom

    dydt(2) = -flux_coeff*radius**2*dydt(1)

    if ( (dydt(2) < 0.0) .and. (abs(dydt(2)) > qv) ) then
      dydt(2) = - qv
      dydt(1) = - dydt(2)/(flux_coeff*radius**2)
    end if

    dydt(3) = -Lcond_temp/cpm*dydt(2)

    dydt(4) = -dydt(2)*inv_grid_scale

end subroutine fcnkb

end module DGM
