! Mode-agnostic blob entrainment mechanics.
!
! Implements the Krueger et al. (1997) blob method: at stochastic intervals,
! one or more contiguous regions ("blobs") of the 1-D domain are replaced with
! environmental air. Particles in the blob are detrained and fresh aerosols
! from the environmental population are entrained at Koehler equilibrium.
!
! This module owns the &ENTRAINMENT namelist and entrainment timing. It receives
! environmental state (T_env, qv_env, vel) from the caller, so the same mechanics
! can be driven by parcel mode (via parcel.f90) or chamber mode in the future.
!
! Reference: Krueger, S. K., Su, C.-W., & McMurtry, P. A. (1997),
!   J. Atmos. Sci., 54, 2697-2712.
module entrainment
    use globals
    use microphysics, only: virtual_temp, update_supersat
    use droplets, only: detrain_particles, entrain_particles, aerosol_concentration
    implicit none

    private
    public :: initialize_entrainment, apply_entrainment
    public :: ent_rate, n_blob, psigma, random_entrainment

    ! --- &ENTRAINMENT namelist variables ---
    real(dp) :: ent_rate = 2.0       ! fractional entrainment rate [1/m]
    integer(i4) :: n_blob = 1        ! number of blobs per entrainment event
    real(dp) :: psigma = 0.1         ! blob fraction of domain per blob [dimensionless]
    logical  :: random_entrainment = .true. ! Poisson-randomize entrainment timing

    ! --- Entrainment timing ---
    real(dp) :: t_next_entrain = 0.0 ! time of next entrainment event [s]

contains

    ! Read &ENTRAINMENT namelist, validate parameters, and schedule the first event.
    subroutine initialize_entrainment(vel)
        real(dp), intent(in) :: vel  ! current parcel/eddy velocity [m/s]
        integer :: nml_unit, ierr
        character(256) :: nml_line, io_emsg

        namelist /ENTRAINMENT/ ent_rate, n_blob, psigma, random_entrainment

        write(*,*) 'Reading ENTRAINMENT namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, &
             action='read', status='old')
        if (ierr /= 0) then
            write(0,*) io_emsg; stop 1
        end if
        read(nml=ENTRAINMENT, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) then
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(0,'(a)') 'Invalid ENTRAINMENT namelist parameter: '//trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        if (psigma <= 0.0 .or. psigma >= 1.0) then
            write(0,*) 'Error: psigma must be in (0, 1), got: ', psigma
            stop 1
        end if
        if (psigma * n_blob >= 1.0) then
            write(0,*) 'Error: psigma * n_blob must be < 1'
            stop 1
        end if
        if (ent_rate <= 0.0) then
            write(0,*) 'Error: ent_rate must be > 0, got: ', ent_rate
            stop 1
        end if

        t_next_entrain = compute_dt_entm(vel)

        write(*,*) 'entrainment:       ON'
        write(*,*) '  ent_rate:        ', ent_rate
        write(*,*) '  n_blob:          ', n_blob
        write(*,*) '  psigma:          ', psigma
        write(*,*) '  random_entrain:  ', random_entrainment
        write(*,*) '  first dt_entm:   ', t_next_entrain, ' s'

    end subroutine initialize_entrainment


    ! Execute one entrainment event if the scheduled time has been reached.
    !
    ! Sequence: (1) place blobs, (2) detrain particles in blob region,
    ! (3) replace scalars with environmental values, (4) update Tv and SS,
    ! (5) entrain fresh aerosols equilibrated to the new (entrained) air,
    ! (6) schedule the next event.
    !
    ! Detrainment precedes scalar replacement so budget counters capture the
    ! pre-entrainment particle state. Entrainment follows scalar replacement
    ! so new particles equilibrate to the entrained environmental air.
    subroutine apply_entrainment(T_env, qv_env, vel)
        real(dp), intent(in) :: T_env   ! environmental temperature [K]
        real(dp), intent(in) :: qv_env  ! environmental water vapor mixing ratio [kg/kg]
        real(dp), intent(in) :: vel     ! parcel/eddy velocity for timing [m/s]
        integer, parameter :: max_blobs = 10
        integer :: blob_start(max_blobs), blob_end(max_blobs), n_final
        integer :: i, k

        if (abs(vel) < 1.0e-30) return
        if (time < t_next_entrain) return

        call place_blobs(blob_start, blob_end, n_final)

        call detrain_particles(blob_start, blob_end, n_final)

        ! Replace scalar fields in blob regions with environmental air
        do i = 1, n_final
            do k = blob_start(i), blob_end(i)
                T(k) = T_env
                WV(k) = qv_env
            end do
        end do

        ! Recompute derived fields over the full domain
        do k = 1, N
            Tv(k) = virtual_temp(T(k), WV(k))
        end do
        call update_supersat(T, WV, SS, pres)

        call entrain_particles(blob_start, blob_end, n_final, aerosol_concentration)

        t_next_entrain = time + compute_dt_entm(vel)

    end subroutine apply_entrainment


    ! Compute the time interval until the next entrainment event [s].
    !
    ! Deterministic interval: dt = (n_blob / ent_rate) * (psigma / (1 - psigma)) / |vel|
    ! With random_entrainment, the interval is drawn from an exponential distribution
    ! (Poisson process) by multiplying by -ln(1 - U), U ~ Uniform(0,1).
    ! See Krueger et al. (1997), eq. (3).
    function compute_dt_entm(vel) result(dt_entm)
        real(dp), intent(in) :: vel  ! parcel/eddy velocity [m/s]
        real(dp) :: dt_entm, u

        dt_entm = (real(n_blob, dp) / ent_rate) * (psigma / (1.0 - psigma)) &
                  / abs(vel)

        if (random_entrainment) then
            call random_number(u)
            dt_entm = dt_entm * (-log(1.0 - u))
        end if
    end function compute_dt_entm


    ! Randomly place n_blob contiguous blobs on the periodic 1-D domain.
    !
    ! Each blob spans blob_size = psigma * N gridcells. Blobs that wrap past
    ! cell N are split into two contiguous segments (e.g., [s, N] and [1, remainder]),
    ! so n_final may exceed n_blob. The starts/ends arrays use gridcell indices [1, N].
    subroutine place_blobs(starts, ends, n_final)
        integer, intent(out) :: starts(:), ends(:), n_final
        integer :: xn, blob_size, s, i
        real(dp) :: u

        blob_size = int(psigma * N)
        xn = N - blob_size * n_blob  ! free cells available for random offset

        call random_number(u)
        s = int(u * xn) + 1

        n_final = 0
        do i = 1, n_blob
            if (i > 1) s = ends(n_final) + int(real(xn, dp) / n_blob) + 1

            ! Blob wraps around the periodic boundary — split into two segments
            if (s <= N .and. s + blob_size - 1 > N) then
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = N
                n_final = n_final + 1
                starts(n_final) = 1
                ends(n_final) = blob_size - (N - s + 1)
            else if (s > N) then
                n_final = n_final + 1
                starts(n_final) = mod(s - 1, N) + 1
                ends(n_final) = mod(s - 1 + blob_size - 1, N) + 1
                if (ends(n_final) < starts(n_final)) then
                    ends(n_final) = N
                    n_final = n_final + 1
                    starts(n_final) = 1
                    ends(n_final) = blob_size - (N - starts(n_final - 1) + 1)
                end if
            else
                n_final = n_final + 1
                starts(n_final) = s
                ends(n_final) = s + blob_size - 1
            end if
        end do
    end subroutine place_blobs

end module entrainment
