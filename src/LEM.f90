module LEM
    use globals
    use microphysics, only: virtual_temp, update_supersat
    use droplets, only: particle, particles, current_n_particles, move_particles_by_cellmap
    implicit none

    private
    public :: initialize_LEM, lem_diffuse_step, lem_turbulence_step, lem_sync_after_physics
    public :: integral_length_scale, dissipation_rate
    public :: reynolds_number
    public :: actual_kolmogorov_scale, grid_eddy_scale, diffusivity_length_scale
    public :: smallest_eddy_scale, smallest_eddy_gridpoints, diffusivity_enhancement
    ! Exposed for the unit tests, which drive the derivation without a namelist.
    public :: derive_turbulence_scales, thermal_diffusivity, vapor_diffusivity

    ! LEM namelist parameters
    real(dp) :: integral_length_scale = 1 ! m
    real(dp) :: dissipation_rate = 0.01   ! m^2/s^3

    ! Derived LEM parameters (set in initialize_LEM)
    integer(i4) :: maps_per_event        ! Triplet maps per eddy event
    integer(i4) :: steps_between_events  ! Diffusion steps between eddy events
    integer(i4) :: iteration_count       ! Running iteration counter

    ! Turbulence scales and diffusivities (all set in initialize_LEM; see the
    ! naming table and derivation there before touching any of these).
    real(dp) :: integral_scale_diffusivity  ! 0.1*L_int^(4/3)*eps^(1/3) (m^2/s)
    real(dp) :: reynolds_number             ! (L_int / smallest_eddy_scale)^(4/3)
    real(dp) :: vapor_diffusivity           ! LEM diffusivity for water vapor (m^2/s)
    real(dp) :: thermal_diffusivity         ! LEM diffusivity for temperature (m^2/s)

    real(dp)    :: actual_kolmogorov_scale    ! (nu^3/eps)^(1/4), DIAGNOSTIC ONLY (m)
    real(dp)    :: grid_eddy_scale            ! min_eddy_gridpoints * dz_length (m)
    real(dp)    :: diffusivity_length_scale   ! where turbulent D == molecular (m)
    real(dp)    :: smallest_eddy_scale        ! smallest_eddy_gridpoints * dz (m) -- governs all
    integer(i4) :: smallest_eddy_gridpoints   ! multiple of 3, >= min_eddy_gridpoints
    real(dp)    :: diffusivity_enhancement    ! unitless, >= 1; scales kT and Dv alike

    ! The cell-label tracers used to compose the per-event triplet maps now live
    ! in globals alongside triplet_map, so ODT and LEM share one implementation.

contains

    subroutine initialize_LEM(domain_height)
        ! Sets up LEM for a parcel run: reads &TURBULENCE_LEM, derives the turbulent
        ! diffusivity and Reynolds number from the integral/Kolmogorov scales, then
        ! reconciles the molecular-diffusion stability step against the eddy-event
        ! rate (Krueger 1993) to fix dt, maps_per_event, and steps_between_events.
        real(dp), intent(in) :: domain_height
        real(dp) :: diffusion_timestep, convection_timestep
        real(dp) :: large_eddy_turnover_time, eddy_rate_per_length

        call read_lem_params()

        ! -------------------------------------------------------------------
        ! Turbulence scale naming -- FOUR length scales, three of them
        ! "Kolmogorov"-ish. Conflating them is the easy mistake here, so none is
        ! called plain `eta`.
        !
        !   actual_kolmogorov_scale   (nu^3/eps)^(1/4). The physical dissipation
        !                             scale. DIAGNOSTIC ONLY -- see below, it can
        !                             never govern.
        !   grid_eddy_scale           min_eddy_gridpoints * dz_length = 6*dz. The
        !                             smallest eddy the triplet map can represent
        !                             at all (below 6 cells each of the map's
        !                             three segments is one cell, so the map
        !                             carries no sub-eddy structure).
        !   diffusivity_length_scale  where the turbulent diffusivity falls to
        !                             the molecular one, i.e. the bottom of the
        !                             cascade as this closure sees it. Solving
        !                             0.1*eps^(1/3)*l^(4/3) = D gives
        !                             l_D = (D / (0.1*eps^(1/3)))^(3/4).
        !   smallest_eddy_scale       quantized max(grid, diffusivity). THE
        !                             GOVERNING SCALE: it sets the eddy sampler
        !                             bound, the Reynolds number, and the
        !                             diffusivity enhancement, so those three
        !                             cannot drift apart.
        !
        ! Why actual_kolmogorov_scale never wins: 0.1*eps^(1/3)*eta^(4/3) reduces
        ! to 0.1*nu identically (eps cancels), so
        ! diffusivity_length_scale/eta = (10*D/nu)^(3/4) ~ 7.7 for both kT and
        ! Dv. The physical eta is always ~8x below the scale at which this
        ! closure's turbulence actually hands off to molecular transport. It is
        ! reported for reference, not used.
        !
        ! THE TWO REGIMES, one formula. Let l_D = diffusivity_length_scale.
        !
        !   6*dz > l_D  the grid is coarser than the handoff, so eddies between
        !               l_D and 6*dz are real but unrepresentable. Their stirring
        !               is absorbed by scaling BOTH LEM diffusivities up by one
        !               common factor -- which preserves Pr and Sc exactly.
        !   6*dz < l_D  the grid is finer than the handoff, so eddies between
        !               6*dz and l_D sit below the diffusive cutoff and are
        !               meaningless. Raise the smallest eddy to l_D instead; the
        !               factor then comes out to 1 and the diffusivities stay
        !               molecular.
        !
        ! Both collapse into the branchless derivation below, because
        ! diffusivity_enhancement is 1 exactly when smallest_eddy_scale == l_D.
        !
        ! INVARIANT, the whole point of the scheme: diffusivity_enhancement >= 1,
        ! so thermal_diffusivity >= kT and vapor_diffusivity >= Dv always. LEM
        ! diffusion is never slower than molecular. `ceiling` (not `nint`) is
        ! what guarantees it -- rounding to nearest can land below l_D and drive
        ! the factor under 1.
        !
        ! SCOPE: these two diffusivities are private to LEM and reach only the
        ! diffusion stability step and diffuse_scalar_periodic. Droplet growth
        ! uses its own T- and p-dependent values (DGM.f90:204-205) and is
        ! deliberately untouched; the globals kT/Dv are not modified either,
        ! since Pr/Sc built from them drive chamber-mode ODT.
        !
        ! ORDERING IS LOAD-BEARING. Everything below -- the diffusivities,
        ! Reynolds number, both timesteps, dt, maps_per_event and
        ! steps_between_events -- derives from smallest_eddy_scale. Do not
        ! reorder.
        !
        ! Why the grid enters at all: in Krueger, eta IS the smallest rearranged
        ! segment, and EMPM enforces that by choosing its grid as
        ! ngrid = NINT(BL_inte/(eta/6))*(BL/BL_inte), i.e. 6*dx = eta exactly, so
        ! it can never reach either regime. CODT sets N independently of eta, so
        ! the scale has to be formed from the grid and reconciled against the
        ! diffusivities.
        ! -------------------------------------------------------------------

        call derive_turbulence_scales()

        ! -------------------------------------------------------------------
        ! Timestep naming -- four similar names, different quantities. Note in
        ! particular that `diffusion_step` and `diffusion_timestep` differ by one
        ! word and are NOT the same thing.
        !
        !   dt                  (global)    the model timestep; the main loop
        !                                   advances time = time + dt (CODT.f90:69)
        !   delta_time          (global)    time since the last physics update,
        !                                   time - last_time_updated (CODT.f90:73)
        !   diffusion_step      (global)    the backstop THRESHOLD compared
        !                                   against delta_time (CODT.f90:82);
        !                                   set to dt below
        !   diffusion_timestep  (LEM local) the diffusive STABILITY LIMIT
        !                                   0.2*dz^2/D, computed just below and
        !                                   used only to derive dt,
        !                                   maps_per_event, steps_between_events
        !
        ! diffusion_step == diffusion_timestep only in the
        ! diffusion_timestep >= convection_timestep branch, where dt is taken
        ! from diffusion_timestep. In the other branch dt comes from
        ! convection_timestep and the two differ.
        !
        ! The names change again at the call boundary: lem_turbulence_step takes
        ! dt and delta_time as dummies named ldt and ldelta_time.
        !
        ! See docs/known_issues.md, "Physics chain may run twice per iteration
        ! when an eddy is accepted" -- an open question about delta_time being
        ! computed once per iteration and not recomputed after the backstop.
        ! -------------------------------------------------------------------

        ! Diffusion stability timestep
        diffusion_timestep = 0.2 * dz_length**2 / max(vapor_diffusivity, thermal_diffusivity)

        ! Turbulent convection timestep (Krueger 1993)
        large_eddy_turnover_time = integral_length_scale**2 / integral_scale_diffusivity
        eddy_rate_per_length = (54./5. * reynolds_number**1.25) &
                               / (integral_length_scale * large_eddy_turnover_time)
        convection_timestep = 1. / (domain_height * eddy_rate_per_length)

        ! Determine maps_per_event and steps_between_events
        if (diffusion_timestep >= convection_timestep) then
            dt = diffusion_timestep
            maps_per_event = int(diffusion_timestep / convection_timestep)
            steps_between_events = 1
        else
            steps_between_events = int(convection_timestep / diffusion_timestep) + 1
            dt = convection_timestep / real(steps_between_events, dp)
            maps_per_event = 1
        end if

        diffusion_step = dt
        iteration_count = 0

        call begin_eddy_sequence()   ! allocates the shared tracers once

        write(*,*) '--- LEM Configuration ---'
        write(*,'(a)')          ' Turbulence length scales (m):'
        write(*,'(a,es12.4,a)') '   actual eta (nu^3/eps)^1/4  : ', actual_kolmogorov_scale, &
            '   (diagnostic only)'
        write(*,'(a,es12.4)')   '   grid eddy scale (6*dz)     : ', grid_eddy_scale
        write(*,'(a,es12.4)')   '   diffusivity length scale   : ', diffusivity_length_scale
        write(*,'(a,es12.4,a,i0,a,a)') '   smallest eddy scale        : ', smallest_eddy_scale, &
            '   (', smallest_eddy_gridpoints, ' cells)', &
            merge('   <- grid-limited       ', '   <- diffusivity-limited', &
                  grid_eddy_scale > diffusivity_length_scale)

        write(*,'(a)') ' Diffusivities (m^2/s):'
        write(*,'(a,es12.4)')   '   integral scale (DT)        : ', integral_scale_diffusivity
        write(*,'(a,f0.4,a)')   '   LEM enhancement factor     : ', diffusivity_enhancement, &
            '   (1 = molecular; DGM is never enhanced)'
        write(*,'(a,es12.4,a,es12.4)') '     heat  kT ', kT, '  ->', thermal_diffusivity
        write(*,'(a,es12.4,a,es12.4)') '     vapor Dv ', Dv, '  ->', vapor_diffusivity

        write(*,*) 'Re:                  ', reynolds_number
        write(*,*) 'dt (s):              ', dt
        write(*,*) 'maps_per_event:      ', maps_per_event
        write(*,*) 'steps_between_events:', steps_between_events
        write(*,*) '-------------------------'

    end subroutine initialize_LEM


    ! -----------------------------------------------
    ! Controller subroutines matching abstract interfaces
    ! -----------------------------------------------

    subroutine lem_diffuse_step(ldelta_time)
        real(dp), intent(in) :: ldelta_time
        real(dp) :: T_sum_before, WV_sum_before

        T_sum_before = sum(T)
        WV_sum_before = sum(WV)

        Nd = Nd + 1
        call diffuse_scalar_periodic(T, thermal_diffusivity, ldelta_time)
        call diffuse_scalar_periodic(WV, vapor_diffusivity, ldelta_time)
        call update_virtual_temperature()
        call update_supersat(T, WV, SS, pres)

        ! Should be ~0 for periodic BCs; kept for validation
        budget_diffusion_delta_T = budget_diffusion_delta_T + (sum(T) - T_sum_before)
        budget_diffusion_delta_WV = budget_diffusion_delta_WV + (sum(WV) - WV_sum_before)

    end subroutine lem_diffuse_step


    subroutine lem_turbulence_step(ldt, ltime, ldelta_time, &
                                   leddy_accepted, eddy_loc, eddy_len)
        ! Applies LEM eddy events on the schedule fixed at init. Every
        ! steps_between_events iterations, performs maps_per_event triplet maps:
        ! each samples an eddy size from the -5/3 inertial-subrange spectrum, snaps
        ! it to a multiple of 3 gridpoints, picks a random (periodic) location, and
        ! stirs T, WV, and particles. Unlike ODT there is no accept/reject test.
        real(dp), intent(inout) :: ldt
        real(dp), intent(in) :: ltime, ldelta_time
        logical, intent(out) :: leddy_accepted
        integer(i4), intent(out) :: eddy_loc, eddy_len

        integer(i4) :: j, c, gridpoints_raw, eddy_gridpoints, eddy_start
        real(dp) :: rand_size, rand_position, sampled_eddy_size

        iteration_count = iteration_count + 1
        leddy_accepted = .false.
        eddy_loc = 0
        eddy_len = 0

        if (mod(iteration_count, steps_between_events) /= 0) return

        ! Compose all maps_per_event maps into one net rearrangement, so droplets
        ! are displaced once per event rather than once per map.
        if (do_microphysics) call begin_eddy_sequence()

        do j = 1, maps_per_event
            ! Sample eddy size from the -5/3 inertial-subrange spectrum, over
            ! [smallest_eddy_scale, integral_length_scale].
            !
            ! Sampling below smallest_eddy_scale would double-count: those
            ! eddies would be floored back up here while the diffusivity
            ! enhancement set in initialize_LEM already stands in for them.
            call random_number(rand_size)
            sampled_eddy_size = (rand_size &
                * (integral_length_scale**(-5./3.) &
                   - smallest_eddy_scale**(-5./3.)) &
                + smallest_eddy_scale**(-5./3.)) ** (-3./5.)

            ! Convert to gridpoints, quantize to nearest multiple of 3
            gridpoints_raw = int(sampled_eddy_size / dz_length)
            eddy_gridpoints = nint(real(gridpoints_raw, dp) / 3.) * 3
            ! Floor at the DERIVED smallest_eddy_gridpoints, not the structural
            ! constant min_eddy_gridpoints. Flooring at 6 would let eddies land
            ! below the scale diffusivity_enhancement was computed for, so the
            ! enhancement and the sampled eddies would both be representing the
            ! same sub-scale stirring. smallest_eddy_gridpoints is already >= 6
            ! and a multiple of 3 (see initialize_LEM).
            eddy_gridpoints = max(smallest_eddy_gridpoints, min(eddy_gridpoints, N))

            ! Random starting position (periodic wrapping handles boundary)
            call random_number(rand_position)
            eddy_start = int(rand_position * N) + 1

            ! Apply periodic triplet map to dimensional scalars, and the identical
            ! permutation to the cell-label tracer so the maps compose.
            call triplet_map(eddy_gridpoints, eddy_start, T)
            call triplet_map(eddy_gridpoints, eddy_start, WV)
            if (do_microphysics) call accumulate_eddy(eddy_gridpoints, eddy_start)

            eddy_loc = eddy_start
            eddy_len = eddy_gridpoints
        end do

        if (do_microphysics) then
            call finalize_eddy_sequence()
            call move_particles_by_cellmap(particles, destination_cell)
        end if

        leddy_accepted = .true.

    end subroutine lem_turbulence_step


    subroutine lem_sync_after_physics()
        ! No-op: LEM operates entirely in dimensional space
    end subroutine lem_sync_after_physics


    ! -----------------------------------------------
    ! LEM-specific physics routines
    ! -----------------------------------------------

    subroutine diffuse_scalar_periodic(field, diffusivity, elapsed_time)
        ! One Crank-Nicolson diffusion step on a periodic field (parcel mode).
        ! Periodicity makes the implicit matrix *cyclic* tridiagonal (nonzero
        ! corners at (1,N) and (N,1)), which the plain Thomas algorithm cannot
        ! solve. The Sherman-Morrison formula writes the cyclic matrix as a base
        ! tridiagonal A' plus a rank-1 update u*v^T, solves two ordinary
        ! tridiagonal systems against A', and combines them, recovering O(N).
        ! Reference: Press et al., Numerical Recipes, "Cyclic Tridiagonal Systems".
        !
        !   field        - scalar to diffuse in place (e.g. T or WV)
        !   diffusivity  - diffusivity for this field (m^2/s). Molecular only
        !                  when the grid resolves the dissipation range; see
        !                  the scale test in initialize_LEM.
        !   elapsed_time - time step (s)
        real(dp), intent(inout) :: field(:)
        real(dp), intent(in) :: diffusivity, elapsed_time
        real(dp) :: De, gamma, correction
        real(dp) :: diag_val, off_diag
        real(dp) :: y_soln(N), q_soln(N)   ! solutions of the two A' systems
        real(dp) :: rhs(N)
        real(dp) :: lower(N), diag(N), upper(N)
        real(dp) :: q_rhs(N)
        integer(i4) :: k

        ! CN diffusion number; diag_val/off_diag are the implicit-matrix entries
        De = (elapsed_time * diffusivity) / (2.0 * dz_length**2)
        diag_val = 1.0 + 2.0 * De
        off_diag = -De

        ! Explicit (known) half of Crank-Nicolson, with neighbours wrapped periodically
        rhs(1) = (1.0 - 2.0*De)*field(1) + De*(field(2) + field(N))
        do k = 2, N-1
            rhs(k) = (1.0 - 2.0*De)*field(k) + De*(field(k+1) + field(k-1))
        end do
        rhs(N) = (1.0 - 2.0*De)*field(N) + De*(field(1) + field(N-1))

        ! gamma is the free Sherman-Morrison parameter; -diag_val is a stable choice.
        gamma = -diag_val

        ! A' = cyclic matrix with the two corner couplings removed via the rank-1
        ! update: subtract gamma from the (1,1) entry and off_diag^2/gamma from (N,N).
        diag(1) = diag_val - gamma
        do k = 2, N-1
            diag(k) = diag_val
        end do
        diag(N) = diag_val - off_diag * off_diag / gamma

        lower(:) = off_diag
        upper(:) = off_diag

        ! First solve: A' * y_soln = rhs (the physical RHS)
        call tridiagonal_periodic(lower, diag, upper, rhs, y_soln)

        ! Second solve: A' * q_soln = u, the rank-1 update vector
        ! u = [gamma, 0, ..., 0, off_diag] encodes the removed corner couplings.
        q_rhs(:) = 0.0
        q_rhs(1) = gamma
        q_rhs(N) = off_diag
        call tridiagonal_periodic(lower, diag, upper, q_rhs, q_soln)

        ! Sherman-Morrison correction factor = (v^T y) / (1 + v^T q), with
        ! v = [1, 0, ..., 0, off_diag/gamma] selecting the corner contributions.
        correction = (y_soln(1) + off_diag * y_soln(N) / gamma) &
                   / (1.0 + q_soln(1) + off_diag * q_soln(N) / gamma)

        ! Recombine the two solves into the cyclic-system solution
        do k = 1, N
            field(k) = y_soln(k) - correction * q_soln(k)
        end do

    end subroutine diffuse_scalar_periodic


    subroutine tridiagonal_periodic(l, d, u, rhs, x)
        ! Thomas algorithm for the base (non-cyclic) tridiagonal system A' x = rhs
        ! used by the Sherman-Morrison solve above. Forward elimination then back
        ! substitution; assumes A' is non-singular (diagonally dominant here).
        real(dp), intent(in) :: l(:), d(:), u(:)  ! sub-, main-, super-diagonals
        real(dp), intent(in) :: rhs(:)            ! right-hand side
        real(dp), intent(out) :: x(:)             ! solution
        real(dp) :: w(N), b                       ! super-diag/pivot ratios; running pivot
        integer(i4) :: k

        b = d(1)
        x(1) = rhs(1) / b
        do k = 2, N
            w(k) = u(k-1) / b
            b = d(k) - l(k) * w(k)
            x(k) = (rhs(k) - l(k) * x(k-1)) / b
        end do

        do k = N-1, 1, -1
            x(k) = x(k) - w(k+1) * x(k+1)
        end do

    end subroutine tridiagonal_periodic


    subroutine update_virtual_temperature()
        integer(i4) :: k
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
    end subroutine update_virtual_temperature


    subroutine read_lem_params()
        integer :: ierr, nml_unit
        character(256) :: nml_line, io_emsg

        namelist /TURBULENCE_LEM/ integral_length_scale, dissipation_rate

        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr /= 0) then
            write(error_unit,*) io_emsg; call exit(1)
        end if
        read(nml=TURBULENCE_LEM, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) call namelist_read_error(nml_unit, 'TURBULENCE_LEM')
        close(nml_unit)

    end subroutine read_lem_params


    subroutine derive_turbulence_scales()
        ! Derives every turbulence scale and both LEM diffusivities from
        ! dz_length, dissipation_rate and integral_length_scale. Separated from
        ! initialize_LEM so it can be exercised directly by the unit tests
        ! without a namelist or the eddy-tracer allocation; initialize_LEM calls
        ! it once, before anything that consumes its results.
        !
        ! See the naming table in initialize_LEM for what each scale means and
        ! why the ordering here is load-bearing.

        integral_scale_diffusivity = 0.1 * integral_length_scale**(4./3.) &
                                     * dissipation_rate**(1./3.)

        actual_kolmogorov_scale = (nu**3 / dissipation_rate) ** 0.25
        grid_eddy_scale         = min_eddy_gridpoints * dz_length

        ! max(kT, Dv) is the stricter reference: deriving from the larger
        ! molecular diffusivity guarantees BOTH fields end up at or above their
        ! own molecular value.
        diffusivity_length_scale = (max(kT, Dv) &
            / (0.1 * dissipation_rate**(1./3.))) ** 0.75

        ! Quantize up to a multiple of 3 cells (the triplet map's segment count)
        ! and never below the structural floor. `ceiling`, not `nint`: rounding
        ! to nearest can land below diffusivity_length_scale and drive the
        ! enhancement factor under 1.
        smallest_eddy_gridpoints = 3 * ceiling(max(grid_eddy_scale, diffusivity_length_scale) &
                                               / (3. * dz_length))
        smallest_eddy_gridpoints = max(min_eddy_gridpoints, smallest_eddy_gridpoints)
        smallest_eddy_scale      = smallest_eddy_gridpoints * dz_length

        ! >= 1 by construction (smallest_eddy_scale >= diffusivity_length_scale).
        ! Both fields take the SAME factor, so Pr and Sc are preserved exactly.
        diffusivity_enhancement = (smallest_eddy_scale / diffusivity_length_scale) ** (4./3.)
        thermal_diffusivity     = kT * diffusivity_enhancement
        vapor_diffusivity       = Dv * diffusivity_enhancement

        ! Reynolds number spans the integral scale down to the smallest eddy the
        ! model actually has, not the namelist input.
        reynolds_number = (integral_length_scale / smallest_eddy_scale) ** (4./3.)

        ! Fatal on configurations the eddy sampler cannot represent; must run
        ! before the timesteps in initialize_LEM consume reynolds_number.
        call validate_turbulence_scales()

    end subroutine derive_turbulence_scales


    subroutine validate_turbulence_scales()
        ! Guards the derived smallest_eddy_scale against configurations the eddy
        ! sampler cannot represent. Called from initialize_LEM immediately after
        ! smallest_eddy_scale is set and before anything derives from it.
        real(dp) :: scale_separation

        if (smallest_eddy_gridpoints > N) then
            write(error_unit,'(a)') 'Error: the domain cannot contain the smallest eddy.'
            write(error_unit,'(a,i0,a,i0,a)') &
                '  smallest eddy = ', smallest_eddy_gridpoints, ' cells, but N = ', N, '.'
            write(error_unit,'(a,es12.4,a)') &
                '  Required scale: ', smallest_eddy_scale, ' m (set by the larger of'
            write(error_unit,'(a)') &
                '  6*dz and the diffusivity length scale). Increase H, or increase N'
            write(error_unit,'(a)') &
                '  so 6*dz falls below it, or raise dissipation_rate.'
            call exit(1)
        end if

        ! The -5/3 sampler in lem_turbulence_step forms
        ! (L^(-5/3) - l_small^(-5/3)); if l_small >= L that bracket is <= 0 and
        ! raising it to -3/5 yields NaN silently.
        if (smallest_eddy_scale >= integral_length_scale) then
            write(error_unit,'(a)') 'Error: no inertial range -- the smallest eddy is not'
            write(error_unit,'(a)') '  smaller than the integral length scale.'
            write(error_unit,'(a,es12.4,a,es12.4,a)') &
                '  smallest_eddy_scale = ', smallest_eddy_scale, &
                ' m >= integral_length_scale = ', integral_length_scale, ' m.'
            write(error_unit,'(a)') &
                '  Increase integral_length_scale, or increase N to refine the grid.'
            call exit(1)
        end if

        scale_separation = integral_length_scale / smallest_eddy_scale
        if (scale_separation < 3.) then
            write(error_unit,'(a)') 'Warning: the inertial range is nearly absent.'
            write(error_unit,'(a,f0.2,a,f0.3,a)') &
                '  integral_length_scale / smallest_eddy_scale = ', scale_separation, &
                ' (Re = ', reynolds_number, ').'
            write(error_unit,'(a)') &
                '  Eddy statistics are drawn from a very narrow band; maps_per_event'
            write(error_unit,'(a)') &
                '  will be small and the turbulence poorly resolved. Consider a larger'
            write(error_unit,'(a)') &
                '  integral_length_scale or a finer grid.'
        end if

    end subroutine validate_turbulence_scales

end module LEM
