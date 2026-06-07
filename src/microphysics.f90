! Bulk thermodynamic relations shared by both simulation modes: saturation
! vapor pressure, saturation mixing ratio, virtual temperature, and the
! supersaturation field. These are the basic warm-cloud thermodynamics that the
! turbulence (ODT/LEM) and droplet-growth (DGM) modules build on.
module microphysics
    use globals
    implicit none

    private
    public :: saturation_vapor_pressure, saturation_mixing_ratio, virtual_temp, &
              update_supersat, calc_supersat

contains

    pure function saturation_vapor_pressure(lT) result(e_sat)
        ! Saturation vapor pressure over liquid water from an 8th-order polynomial
        ! in Celsius temperature (Horner evaluation). Valid roughly -80 to +50 C;
        ! the input is floored at -80 C to stay in the fitted range.
        ! Coefficients from Flatau, Walko & Cotton (1992), J. Appl. Meteor., 31, 1507.
        real(dp), intent(in) :: lT ! air temperature (K)
        real(dp) :: Tc             ! temperature in Celsius
        real(dp) :: e_sat          ! saturation vapor pressure (Pa, after *100)
        integer :: i
        ! Polynomial coefficients give e_sat in hPa (mb); converted to Pa below.
        double precision, parameter :: coeff(9) = [6.11239921,0.443987641,0.142986287e-1,0.264847430e-3, &
        0.302950461e-5, 0.206739458e-7,0.640689451e-10,-0.952447341e-13,-0.976195544e-15]

        Tc = max(-80., lT-Tice) ! (C), floored to the polynomial's valid range

        ! Horner's method: evaluate the polynomial from the highest-order term down
        e_sat = coeff(9)
        do i = 8,1,-1
            e_sat = e_sat * Tc + coeff(i)
        end do

        e_sat = e_sat * 100. ! hPa -> Pa

    end function saturation_vapor_pressure

    pure function saturation_mixing_ratio(lT, lpres) result(q_sat)
        ! Saturation mixing ratio over liquid water: q_sat = eps * e_sat/(p - e_sat),
        ! with eps = R_dry/R_vapor. The max() guards against p <= e_sat at very
        ! high temperature / low pressure (keeps the denominator >= e_sat).
        real(dp), intent(in) :: lT, lpres   ! temperature (K), pressure (Pa)
        double precision :: q_sat           ! saturation mixing ratio (kg/kg)
        double precision :: e_sat           ! saturation vapor pressure (Pa)

        e_sat = saturation_vapor_pressure(lT)
        q_sat = eps * e_sat / max(e_sat, (lpres-e_sat)) ! kg/kg

    end function saturation_mixing_ratio

    pure function virtual_temp(lT, lWV) result(Tvirt)
        ! Virtual temperature: the temperature dry air would need to match the
        ! density of this moist air. Uses the linearized form Tv = T*(1 + eps_tv*qv),
        ! where eps_tv = (1-eps)/eps; accurate for the small mixing ratios here and
        ! cheaper than the exact expression (commented below).
        real(dp), intent(in) :: lT, lWV   ! temperature (K), vapor mixing ratio (kg/kg)
        real(dp) :: Tvirt                 ! virtual temperature (K)

        !Tvirt = lT * ((lWV + eps) / (eps*(1 + lWv)))   ! exact form
        Tvirt = lT * (1 + (eps_tv * lWv))               ! linearized (fast path)

    end function virtual_temp

    pure function calc_supersat(temp, mr, p) result(supersat)
        ! Supersaturation as a percentage: 100*(qv/q_sat - 1). Positive => supersaturated.
        real(dp), intent(in) :: temp, mr, p   ! temperature (K), vapor mixing ratio (kg/kg), pressure (Pa)
        real(dp) :: sat_mr                    ! saturation mixing ratio (kg/kg)
        real(dp) :: supersat                  ! supersaturation (%)

        sat_mr = saturation_mixing_ratio(temp, p)
        supersat = 100*(mr/sat_mr - 1.)

    end function calc_supersat

    pure subroutine update_supersat(lT, lWV, lSS, lpres)
        ! Refreshes the whole supersaturation profile lSS from the current
        ! temperature and vapor fields (called after every diffusion/eddy update).
        real(dp), intent(in) :: lT(:), lWV(:)   ! temperature (K), vapor mixing ratio (kg/kg) profiles
        real(dp), intent(in) :: lpres           ! reference pressure (Pa)
        real(dp), intent(out) :: lSS(:)         ! supersaturation profile (%), overwritten

        integer(i4) :: k

        do concurrent (k = 1:N)
            lSS(k) = calc_supersat(lT(k), lWV(k), lpres)
        end do

    end subroutine update_supersat



end module microphysics