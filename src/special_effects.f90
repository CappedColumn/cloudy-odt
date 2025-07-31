module special_effects
    use globals
    use microphysics, only: saturation_mixing_ratio, update_nondim_scalars

    ! Sidewall Parameters
    logical :: do_sidewalls = .false.
    real(dp) :: area_sw, area_bot, C_sw, T_sw, RH_sw, WV_sw, P_sw
    real(dp) :: Ra, Nuss, tau_sw
    real(dp) :: sw_iter, sw_nudging_time

    ! Stochastic Fallout parameters
    logical :: do_random_fallout = .false.
    real(dp) :: random_fallout_rate = 1.

    public :: initialize_special_effects, run_special_effects, do_random_fallout, &
    random_fallout_rate
    private

contains

    subroutine run_special_effects(Tarr, WVarr, delta_t)
        ! Interface subroutine to main.f90 to run special effects
        real(dp), intent(in) :: delta_t
        real(dp), intent(inout) :: Tarr(:), WVarr(:)

        ! Implement sidewall forcing and iteration
        if ( do_sidewalls ) then
            sw_iter = sw_iter + delta_t
            if ( sw_iter >= sw_nudging_time ) then
                call sidewall_fluxes(Tarr, WVarr, sw_iter)
                sw_iter = 0.
                call update_nondim_scalars(Tarr, WVarr, Tvdim, T, WV, Tv)
            end if
        end if

    end subroutine run_special_effects


    subroutine initialize_special_effects(filename)

        character(*), intent(in) :: filename
        integer :: nml_unit, ierr
        character(100) :: io_emsg, nml_line

        namelist /SPECIALEFFECTS/ do_sidewalls, area_sw, area_bot, C_sw, T_sw, RH_sw, P_sw, sw_nudging_time, &
        do_random_fallout, random_fallout_rate

        write(*,*) 'Initializing Special Effects...'
        open(newunit=nml_unit, file=trim(namelist_path), iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) io_emsg; stop
        end if
        read(nml=SPECIALEFFECTS, unit=nml_unit, iostat=ierr)
        ! Print value causing namelist read error
        if (ierr .ne. 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(*,'(a)') 'Invalid Namelist Parameter: '//trim(nml_line)
            stop
        end if
        close(nml_unit)

        ! Write namelist parameters
        open(newunit=nml_unit, file=trim(filename)//'_nml.txt', action='write', position='append')
        write(nml_unit, nml=SPECIALEFFECTS)
        close(nml_unit)

        ! Initialize enabled special effects
        if ( do_sidewalls ) call initialize_sidewalls()


    end subroutine initialize_special_effects


    subroutine initialize_sidewalls()

        real(dp) :: velocity_bot

        ! Set time tracker to zero
        sw_iter = 0.

        ! Convert sidewall temperature to Kelvin
        T_sw = T_sw + Tice

        ! Calculate rayleigh and nusselt numbers
        Ra = (g * Tdiff * H**3)/(Tref * nu * kT)

        Nuss = 0.124 * Ra**0.309

        ! Calculate the sidewall forcing timescale
        velocity_bot = P_sw * Nuss * kT / H
        tau_sw = (C_sw * velocity_bot) / ((H * area_bot) / area_sw)  ! inverse timescale for sidewall forcing

        ! Determine Water Vapor of sidewalls
        WV_sw = RH_sw * saturation_mixing_ratio(T_sw, pres)

    end subroutine initialize_sidewalls

    subroutine sidewall_fluxes(Tarr, WVarr, delta_t)
        ! Nudges scalar fields be the sidewall values
        real(dp), intent(inout) :: Tarr(:), WVarr(:)
        real(dp), intent(in) :: delta_t
        real(dp) :: time_ratio
        integer :: k

        time_ratio = delta_t * tau_sw

        do concurrent (k = 2:N-1)
            Tarr(k) = Tarr(k) + (T_sw - Tarr(k)) * time_ratio
            WVarr(k) = WVarr(k) + (WV_sw - WVarr(k)) * time_ratio
        end do

    end subroutine sidewall_fluxes

end module special_effects