module ODT
    use globals
    use microphysics, only: update_dim_scalars, update_nondim_scalars, update_supersat
    use droplets, only: particles, move_particles_in_eddy
    use writeout, only: write_eddy
    implicit none

    ! Contains the subroutines and functions necessary for implementing the turbulent aspects
    ! of the model. This includes calculating the eddy probabilities, sizes, location, and
    ! the eddy rejection/acceptance method (the entire loop)

    private
    public :: initialize_ODT, diffusion, calc_eddy_length_cdf, eddy_acceptance_method
    public :: odt_diffuse_step, odt_turbulence_step, odt_sync_after_physics

    ! Non-dimensional time variables
    real(dp) :: dt_nd             ! Non-dimensional time step
    real(dp) :: time_conv_nd      ! Conversion factor: dimensional -> non-dim time (nu/H^2)

    ! Eddy acceptance/rejection counters
    integer(i4) :: Np = 0       ! Number of eddy probs. calculated
    integer(i4) :: Na = 0       ! Number of accepted eddies
    real(dp) :: Pa = 0.         ! Total acceptance probabilities

    ! Integrated eddy values (set in eddy_acceptance_prob, used in implement_eddy)
    real(dp) :: wK, TvK
    real(dp) :: pot_energy       ! Potential energy

contains

    subroutine initialize_ODT(H_domain)
        ! Initialize ODT module: set time conversion, timestep, and
        ! eddy probability parameters. Must be called after namelist
        ! read and thermodynamic boundary conditions are set.
        real(dp), intent(in) :: H_domain

        time_conv_nd = nu / (H_domain**2)
        dt_nd = 1.0 / (1.0 * N * N)

        diffusion_step = H_domain**2 / (nu * 1.0 * N * N)
        dt = diffusion_step

        buoy_nd = (8. * g * alpha * Tvdiff * C2 * H_domain**3) / (27. * nu * nu)
        LpD = 2 * Lprob
        Co = exp(-LpD / (1.*Lmin))
        Cm = exp(-LpD / (1.*Lmax))
        prob_coeff = (exp(-LpD/(1.*Lmax)) - exp(-LpD/(1.*Lmin))) * (N/(3.*LpD))

    end subroutine initialize_ODT

    subroutine diffusion(ldelta_time)
        ! Diffuses the non-dim scalar fields via Crank-Nicolson tridiagonal solver.
        ! Takes dimensional delta_time (seconds) and converts to non-dim internally.
        ! Caller is responsible for syncing dim/nondim fields afterward.
        real(dp), intent(in) :: ldelta_time
        real(dp) :: delta_time_nd

        delta_time_nd = ldelta_time * time_conv_nd

        call diffuse_scalar(W_nd, Ndnu, delta_time_nd)
        call diffuse_scalar(T_nd, Pr, delta_time_nd)
        call diffuse_scalar(WV_nd, Sc, delta_time_nd)

    end subroutine diffusion

    subroutine calc_eddy_length_cdf(prob_L)
        ! Calculates the assumed distribution of eddy lengths f(l).
        ! This distribution is sampled and tested against in the eddy
        ! acceptance method. See Appendix A of Kerstein 1999 for details.
        real(dp), intent(out) :: prob_L(:)

        ! integer(i4) :: L
        ! real(dp) :: norm_freq
        ! real(dp) :: length_pdf(N)

        ! ! Calculate Probability Distribution of eddy sizes
        ! length_pdf = 0.
        ! do concurrent (L = Lmin : Lmax)
        !     length_pdf(L) = exp(-LpD/(1.*L))*(exp(LpD/(L*(L+1.)))-1.)
        ! end do
        
        ! ! Create Normalized Cumulative Distribution function
        ! prob_L = 0.
        ! norm_freq = 1. / sum(length_pdf)
        ! do L = Lmin, Lmax
        !     prob_L(L) = prob_L(L-1) + (norm_freq*length_pdf(L))
        ! end do

        integer :: L
        real(dp) :: lC, lz
        lC = 0.d0
        do L = Lmin, Lmax
            lz = dexp(-LpD/(1.d0*L))*(dexp(LpD/(L*(L+1.d0)))-1.d0)
            lC = lC + lz
        enddo
        lC = 1.d0/lC
        do L = 1, Lmin-1
            prob_L(L) = 0.d0
        enddo
        do L = Lmin, Lmax
            lz = dexp(-LpD/(1.d0*L))*(dexp(LpD/(L*(L+1.d0)))-1.d0)
            prob_L(L) = prob_L(L-1) + (lC*lz)
        enddo
        do L = Lmax+1, size(prob_L)
            prob_L(L) = 0.d0
        enddo

        ! open(unit=999, file='prob.txt')
        ! do L = 1, size(prob_L)
        !     write(999,*) prob_L(L)
        ! end do
        ! close(999)
        

    end subroutine calc_eddy_length_cdf


    function sample_eddy_length(rand_num) result(L)
    ! This function might be really slow
    ! This function corresponds to equation 3.78 in section 3.3.5.1.3 of McDermott's Dissertation
        real(dp), intent(in) :: rand_num
        real(dp) :: initial_guess
        integer(i4) :: Pidx, L

        ! Equation (7) in appendix of ODTdoc1.pdf by Kerstein
        initial_guess = (-2. * Lmin) / log((Co*rand_num)+(Cm*(1. - rand_num)))
        Pidx = int(initial_guess) - 1

        do while (rand_num .gt. prob_eddy_length(Pidx))
            Pidx = Pidx + 1
        end do

        do while (rand_num .lt. prob_eddy_length(Pidx-1))
            Pidx = Pidx - 1
        end do

        L = 3 * Pidx

    end function sample_eddy_length


    pure function sample_eddy_location(Ng, L, rand_num) result(M)
        real(dp), intent(in) :: rand_num
        integer(i4), intent(in) :: Ng, L
        integer(i4) :: max_loc, rand_loc, M

        ! Given an eddy length (L), randomly finds a location
        ! for the eddy to occur in a valid manner

        ! Calculates upper most cell where eddy of length L will occur
        !max_loc = Ng - L - 1
        ! Locate a random, valid location for eddy to occur
        M = 1 + int(rand_num * (Ng - L))
        !M = 1 + min(max_loc, rand_loc)

    end function sample_eddy_location

    subroutine eddy_acceptance_prob(M, L, laccept_prob)
        integer(i4), intent(in) :: L, M
        real(dp), intent(out) :: laccept_prob
        real(dp) :: prob, Lnd, kin_energy

        real(dp) :: p

        ! Non-dimensional Eddy Size
        Lnd = (1.*L)/(1.*N)

        ! Integrate values across the eddy (fused single pass)
        call integrate_eddy_pair(L, M, W_nd, Tv_nd, wK, TvK)
        TvK = -TvK

        ! Calculate energies
        pot_energy = buoy_nd * TvK * Lnd
        kin_energy = wK*wK
        p = ((pot_energy + kin_energy)*Lnd*Lnd) - ZC2 ! REVIEW ! Different from Bodt
        !write(*,*) 'p: ', p
        ! Note: Bodt uses "disfac" and "ratefac" as parameters
        ! Bodt is VERY different
        if (p .gt. 0.) then
            prob = prob_coeff * (1.-Lnd) * sqrt(p) * exp(3.*LpD/(Lnd*N)) / (Lnd*Lnd)
        else
            prob = 0.
        end if

        laccept_prob = prob * dt_nd

        !write(9999, *) prob_coeff, Lnd, p, LpD, exp(3.*LpD/(Lnd*N))

    end subroutine eddy_acceptance_prob


    subroutine lower_dt(prob)
        ! The eddy acceptance/rejection method requires appropriate
        ! sampling of the eddy distribution. If the eddy acceptance
        ! prob. is too high, then we are essentially undersampling
        ! eddies of a certain size/location, and the implementation
        ! of eddies will not be physically accurate.

        real(dp), intent(inout) :: prob

        ! max_accept_prob is a parameter specified on initialization

        ! Test if the the acceptance prob is too high, If it is,
        ! then reduce the time frequency at which we sample to 
        ! induce a higher sampling rate
        if (prob .gt. 0) then
            if (prob .gt. max_accept_prob) then
                dt_nd = dt_nd * max_accept_prob / prob
                prob = max_accept_prob
            end if
            ! Used to calculate the average acceptance rate
            ! This will be used in raise_dt()
            Pa = Pa + prob
            Np = Np + 1
        end if

    end subroutine lower_dt

    subroutine raise_dt()
        real(dp) :: pmin
        
        pmin = 1.e-3
        pa = pa/Np
        if (pa .lt. (pmin/2.)) then
            dt_nd = dt_nd * 2.
        else
            dt_nd = dt_nd * pmin / pa
        end if
        pa = 0.
        Np = 0
        
    end subroutine raise_dt


    subroutine eddy_acceptance_method(ldt, M, L, eddy_flag)
        ! After sampling eddy size/length (L) and eddy location (M),
        ! computes an expected acceptance rate based on the current
        ! state of the velocity/vector fields.

        ! The acceptance probability is then tested, if too high an
        ! adjustment is made with by lowering the timestep.
        !
        ! Takes dimensional dt, converts to non-dim internally,
        ! and returns the (possibly adjusted) dimensional dt.

        real(dp), intent(inout) :: ldt
        integer(i4), intent(out) :: M, L
        logical, intent(out) :: eddy_flag
        real(dp) :: rand_num(3)
        real(dp) :: accept_prob

        ! Sync non-dim dt from dimensional input
        dt_nd = ldt * time_conv_nd

        ! Get 3 random numbers from a uniform distribution
        call random_number(rand_num)
        ! Use the random numbers to sample a eddy length and eddy location
        L = sample_eddy_length(rand_num(1))
        M = sample_eddy_location(N, L, rand_num(2))
        ! Calculated the probability of the randomly selected eddy
        call eddy_acceptance_prob(M, L, accept_prob)

        !if ( accept_prob > 0.0 ) write(*,*) accept_prob, L, M

        call lower_dt(accept_prob)

        ! Monte-Carlo test for this eddy
        if (rand_num(3) .lt. accept_prob) then
            ! Success
            eddy_flag = .true.
            call implement_eddy(L, M)

            Na = Na + 1

        else
            eddy_flag = .false.
        end if

        if (Np > 1e4) call raise_dt()

        ! Convert adjusted dt_nd back to dimensional
        ldt = dt_nd / time_conv_nd

    end subroutine eddy_acceptance_method


    pure function integrate_eddy(L, M, array) result(integral)
        ! Determines the total potential/kinetic energy over an eddy
        integer(i4), intent(in) :: L, M
        real(dp), intent(in) ::  array(:)
        
        integer(i4) :: k
        real(dp):: integral!, total, eddy_flip

        ! total = array(M)*L/2.
        ! do k = 1, L-1
        !     eddy_flip = L - (2*k)
        !     total = total + (array(M+k)*eddy_flip)
        ! end do
        ! total = total - (array(M+L)*L/2.)
        ! integral = 4. * total / (9.*L*L)

        integral = array(M)*L*.5
        do k=1, L-1
            integral = integral + (array(M+k) * (L-2*k))
        end do

        integral = integral - (array(M+L)*L*.5)
        integral = 4.*integral/(9.*L*L)

    end function integrate_eddy


    pure subroutine integrate_eddy_pair(L, M, array1, array2, integral1, integral2)
        integer(i4), intent(in) :: L, M
        real(dp), intent(in) :: array1(:), array2(:)
        real(dp), intent(out) :: integral1, integral2

        integer(i4) :: k, weight
        real(dp) :: half_L

        half_L = L * 0.5

        integral1 = array1(M) * half_L
        integral2 = array2(M) * half_L

        do k = 1, L-1
            weight = L - 2*k
            integral1 = integral1 + array1(M+k) * weight
            integral2 = integral2 + array2(M+k) * weight
        end do

        integral1 = integral1 - array1(M+L) * half_L
        integral2 = integral2 - array2(M+L) * half_L

        integral1 = 4. * integral1 / (9. * L * L)
        integral2 = 4. * integral2 / (9. * L * L)

    end subroutine integrate_eddy_pair


    subroutine implement_eddy(L, M)
        integer(i4), intent(in) :: L, M

        ! uK, vK, wK and PE are assigned in prob/eddy_acceptance_prob
        ! ci - amplitude constants
        ! qi - Pres. Scrambling kernel coefficient (for kinetic energy)
        real(dp) :: qw, cw

        ! Determine the w-scale parameters for energy conservation
        qw = sqrt(wK*wK + pot_energy)
        if (wK .gt. 0.) then
            cw = 6.75*(qw - wK)/(1.*L)
        else
            cw = -6.75*(qw + wK)/(1.*L)
        end if

        ! Create triplet copies
        call triplet_map(L, M, W_nd)
        call triplet_map(L, M, T_nd)
        call triplet_map(L, M, WV_nd)
        ! Scale for energy conservation
        call addK(L, M, W_nd, cw)
        
    end subroutine implement_eddy


    subroutine addK(L, M, ui, cui)
        ! Kernel adjustment for energy conservation
        integer(i4), intent(in) :: L, M
        real(dp), intent(in) :: cui
        real(dp), intent(inout) :: ui(:)
        integer(i4) :: Lseg, j, j1, j2, j3
        real(dp) :: y1, y2, y3

        Lseg = L/3
        do j = 1, Lseg
            y1 = -2.*j
            y2 = 4.*(j+Lseg) - 2.*L
            y3 = 2.*L - 2.*(j+Lseg+Lseg)
            j1 = M + j
            j2 = M + j + Lseg
            j3 = M + j + Lseg + Lseg
            ui(j1) = ui(j1) + cui*y1
            ui(j2) = ui(j2) + cui*y2
            ui(j3) = ui(j3) + cui*y3
        end do

    end subroutine addK

    subroutine diffuse_scalar(sclr, dim_num, ldelta_time_nd)
        real(dp), intent(in) :: dim_num
        real(dp), intent(in) :: ldelta_time_nd
        real(dp), intent(inout) :: sclr(:)
        real(dp) :: De, l(N+1), d(N+1), r(N+1), xsc(N+1)
        integer(i4) :: k

        De = (ldelta_time_nd*(N+1)*(N+1))/(2.*dim_num)

        l(1) = 0.
        d(1) = 1. + (2.*De)
        r(1) = -De
        do k = 2, N
            l(k) = -De
            d(k) = 1. + (2.*De)
            r(k) = -De
        end do
        l(N+1) = 0.
        d(N+1) = 1.
        r(N+1) = 0.

        xsc(1) = ((1.-(2.*De))*sclr(1)) + (De*sclr(2))
        do k = 2, N
            xsc(k) = ((1.-(2.*De))*sclr(k)) + (De*(sclr(k+1) + sclr(k-1)))
        end do
        xsc(N+1) = sclr(N+1)

        call tridiagonal(l, d, r, xsc, sclr)

    end subroutine diffuse_scalar


    subroutine tridiagonal(l, d, r, xui, ui)
        real(dp), intent(in) :: l(:), d(:), r(:), xui(:)
        real(dp), intent(inout) :: ui(:)
        real(dp) :: b, rdx(N)
        integer(i4) :: k

        b = d(1)
        ui(1) = xui(1)/b
        do k = 2, N+1
            rdx(k) = r(k-1)/b
            b = d(k) - (l(k)*rdx(k))
            if (b .eq. 0.) write(*,*) 'Tridiagonal: Failure'
            ui(k) = (xui(k) - (l(k)*ui(k-1)))/b
        end do

        do k = N, 1, -1
            ui(k) = ui(k) - (rdx(k+1)*ui(k+1))
        end do

    end subroutine tridiagonal


    ! -----------------------------------------------
    ! Controller subroutines matching abstract interfaces
    ! in globals.f90 — called via procedure pointers
    ! -----------------------------------------------

    subroutine odt_diffuse_step(ldelta_time)
        real(dp), intent(in) :: ldelta_time
        real(dp) :: T_sum_before, WV_sum_before

        T_sum_before = sum(T)
        WV_sum_before = sum(WV)

        Nd = Nd + 1
        call diffusion(ldelta_time)
        call update_dim_scalars(T_nd, WV_nd, Tv_nd, T, WV, Tv)
        call update_supersat(T, WV, SS, pres)

        budget_diffusion_delta_T = budget_diffusion_delta_T + (sum(T) - T_sum_before)
        budget_diffusion_delta_WV = budget_diffusion_delta_WV + (sum(WV) - WV_sum_before)

    end subroutine odt_diffuse_step


    subroutine odt_turbulence_step(ldt, ltime, ldelta_time, &
                                   leddy_accepted, eddy_loc, eddy_len)
        real(dp), intent(inout) :: ldt
        real(dp), intent(in) :: ltime, ldelta_time
        logical, intent(out) :: leddy_accepted
        integer(i4), intent(out) :: eddy_loc, eddy_len

        call eddy_acceptance_method(ldt, eddy_loc, eddy_len, leddy_accepted)

        if (.not. leddy_accepted) return

        if (write_eddies) call write_eddy(eddy_loc, eddy_len, ltime)

        if (do_microphysics) call move_particles_in_eddy(particles, eddy_loc, eddy_len)

    end subroutine odt_turbulence_step


    subroutine odt_sync_after_physics()
        call update_nondim_scalars(T, WV, Tv, T_nd, WV_nd, Tv_nd)
    end subroutine odt_sync_after_physics


end module ODT

! "Why won't you compile!", they cried.
!
! "I have sacrificed mind and bank account for the 
! holy 3-letter appendage to my name." 
!
! "Are you still not pleased?"
!
! "WHAT MORE CAN I GIVE!?"
!
! ...no response
!
! They sighed, and began to type
! 'GOTO' into their keyboard.