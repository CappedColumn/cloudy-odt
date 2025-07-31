module microphysics
    use globals
    implicit none

    private
    public :: saturation_vapor_pressure, saturation_mixing_ratio, virtual_temp, &
              update_dim_scalars, update_nondim_scalars, update_supersat, calc_supersat

contains

    pure function saturation_vapor_pressure(lT) result(e_sat)
        ! Saturation vapor pressure over water using an 8th-order polynomial approximation
        real(dp), intent(in) :: lT ! (K)
        real(dp) :: Tc, e_sat
        integer :: i
        double precision, parameter :: coeff(9) = [6.11239921,0.443987641,0.142986287e-1,0.264847430e-3, &
        0.302950461e-5, 0.206739458e-7,0.640689451e-10,-0.952447341e-13,-0.976195544e-15]

        Tc = max(-80., lT-Tice) ! (C)

        e_sat = coeff(9)
        do i = 8,1,-1
            e_sat = e_sat * Tc + coeff(i)
        end do

        e_sat = e_sat * 100. ! Convert to pascals

    end function saturation_vapor_pressure

    pure function saturation_mixing_ratio(lT, lpres) result(q_sat)
        ! Given a temperature and pressure, calculates the saturation mixing ratio over water.
        ! Returns a mixing ratio in kg/kg
        real(dp), intent(in) :: lT, lpres ! T (K), pres (pa)
        double precision :: q_sat, e_sat

        e_sat = saturation_vapor_pressure(lT)
        q_sat = eps * e_sat / max(e_sat, (lpres-e_sat)) ! kg/kg

    end function saturation_mixing_ratio

    pure function virtual_temp(lT, lWV) result(Tvirt)
        ! lT (K), lWV (kg/kg)
        ! Returns Tvirt (K)
        real(dp), intent(in) :: lT, lWV
        real(dp) :: Tvirt

        !Tvirt = lT * ((lWV + eps) / (eps*(1 + lWv)))
        Tvirt = lT * (1 + (eps_tv * lWv)) ! I wanna go fast

    end function virtual_temp

    pure function calc_supersat(temp, mr, p) result(supersat)
        ! Calculate the supersaturation (%)
        real(dp), intent(in) :: temp, mr, p
        real(dp) :: sat_mr, supersat

        sat_mr = saturation_mixing_ratio(temp, p)
        supersat = 100*(mr/sat_mr - 1.)

    end function calc_supersat

    ! FIX: Make Tv updated externally, then update dim fron ALL nondim
    pure subroutine update_dim_scalars(lW, lT, lWV, lTv, lWdim, lTdim, lWVdim, lTvdim)
        ! Takes the non-dimensional scalar arrays, and updates the dimensional
        ! arrays. Also calculates the virtual temperature
        real(dp), intent(in) :: lW(:), lT(:), lWV(:)
        real(dp), intent(out) ::  lWdim(:), lTdim(:), lWVdim(:), lTv(:), lTvdim(:)
        integer(i4) :: k

        ! Called after the non-dimensional T and WV fields are changed in some way
        ! Updates the dimensional values and calculates a new virtual temperature
        ! for the next eddy acceptance test

        ! Takes the non-dimensional arrays of T and WV and
        ! calculates the non-dim/dim virtual temperature

        ! Calculated dimensional scalar fields
        ! Note the subtraction because Tdiff and WVdiff are actually negative w.r.t height
        do concurrent (k = 1:N)
            lWdim(k) = lW(k) * w_dim_factor
            lTdim(k) = Tref - Tdiff * lT(k)
            lWVdim(k) = WVref - WVdiff * lWV(k)
        end do

        ! Calculate non-dimensional virtual temperature
        do concurrent (k = 1:N)
            lTvdim(k) = virtual_temp(lTdim(k), lWVdim(k))
            lTv(k) = (Tvref - lTvdim(k)) / Tdiff
        end do 

    end subroutine update_dim_scalars

    pure subroutine update_nondim_scalars(lTdim, lWVdim, lTvdim, lT, lWV, lTv)
        ! Updates the non-dimensional scalar fields based on the current value of 
        ! the dimension scalar fields
        real(dp), intent(in) :: lTdim(:), lWVdim(:), lTvdim(:)
        real(dp), intent(out) :: lT(:), lWV(:), lTv(:)
        integer(i4) :: k

        do concurrent (k = 1:N)
            lT(k) = -(lTdim(k) - Tref) / Tdiff
            lWV(k) = -(lWVdim(k) - WVref) / WVdiff
            lTv(k) = -(lTvdim(k) - Tvref) / Tvdiff
        end do

    end subroutine update_nondim_scalars

    pure subroutine update_supersat(lTdim, lWVdim, lSS, lpres)
        ! Calculates the supersaturation fields based on the current values of
        ! temperature and water vapor fields
        real(dp), intent(in) :: lTdim(:), lWVdim(:), lpres
        real(dp), intent(out) :: lSS(:)

        integer(i4) :: k

        do concurrent (k = 1:N)
            lSS(k) = calc_supersat(lTdim(k), lWVdim(k), lpres)
        end do

    end subroutine update_supersat



end module microphysics