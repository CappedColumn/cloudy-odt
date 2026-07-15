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
    public :: initialize_entrainment, apply_entrainment, get_entrainment_params, &
              refresh_entrainment_schedule
    public :: ent_rate, n_blob, psigma, random_entrainment, time_varying_entrainment, &
              t_next_entrain

    ! --- &ENTRAINMENT namelist variables ---
    real(dp) :: ent_rate = 2.0       ! fractional entrainment rate [1/m]
    integer(i4) :: n_blob = 1        ! number of blobs per entrainment event
    real(dp) :: psigma = 0.1         ! blob fraction of domain per blob [dimensionless]
    logical  :: random_entrainment = .true. ! Poisson-randomize entrainment timing

    ! --- Entrainment timing ---
    real(dp) :: t_next_entrain = 0.0 ! time of next entrainment event [s]

    ! --- Optional time-varying schedule (parcel input v3) ---
    ! When time_varying_entrainment is set, ent_rate/n_blob/psigma above hold the
    ! value active at the current time, refreshed from these per-segment arrays by
    ! get_entrainment_params. The schedule shares the parcel velocity time axis.
    logical :: time_varying_entrainment = .false.
    integer(i4) :: n_ent_segments = 0
    integer(i4) :: active_segment = 1
    real(dp), allocatable :: sched_times(:)
    real(dp), allocatable :: sched_ent_rate(:)
    integer(i4), allocatable :: sched_n_blob(:)
    real(dp), allocatable :: sched_psigma(:)

    ! Maximum blobs per event; sizes the blob arrays in apply_entrainment and
    ! bounds n_blob during validation.
    integer(i4), parameter :: max_blobs = 10

contains

    ! Read &ENTRAINMENT namelist, validate parameters, and schedule the first event.
    !
    ! The schedule arguments are optional: when supplied (parcel input v3), the
    ! per-segment arrays make ent_rate/n_blob/psigma time-varying and override the
    ! namelist scalars (random_entrainment is always taken from the namelist). When
    ! absent, the namelist scalars are used unchanged for the whole run.
    subroutine initialize_entrainment(vel, sched_times_in, sched_ent_rate_in, &
                                      sched_n_blob_in, sched_psigma_in)
        real(dp), intent(in) :: vel  ! current parcel/eddy velocity [m/s]
        real(dp), intent(in), optional :: sched_times_in(:)     ! segment start times [s]
        real(dp), intent(in), optional :: sched_ent_rate_in(:)  ! ent_rate per segment [1/m]
        integer(i4), intent(in), optional :: sched_n_blob_in(:) ! n_blob per segment
        real(dp), intent(in), optional :: sched_psigma_in(:)    ! psigma per segment
        integer :: nml_unit, ierr, i
        character(256) :: nml_line, io_emsg

        namelist /ENTRAINMENT/ ent_rate, n_blob, psigma, random_entrainment

        write(*,*) 'Reading ENTRAINMENT namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, &
             action='read', status='old')
        if (ierr /= 0) then
            write(error_unit,*) io_emsg; call exit(1)
        end if
        read(nml=ENTRAINMENT, unit=nml_unit, iostat=ierr)
        if (ierr /= 0) call namelist_read_error(nml_unit, 'ENTRAINMENT')
        close(nml_unit)

        ! --- Optional time-varying schedule (parcel input v3) ---
        if (present(sched_times_in)) then
            time_varying_entrainment = .true.
            n_ent_segments = size(sched_times_in)
            allocate(sched_times(n_ent_segments), sched_ent_rate(n_ent_segments), &
                     sched_n_blob(n_ent_segments), sched_psigma(n_ent_segments))
            sched_times    = sched_times_in
            sched_ent_rate = sched_ent_rate_in
            sched_n_blob   = sched_n_blob_in
            sched_psigma   = sched_psigma_in
        end if

        ! --- Validate (every segment when time-varying, else the scalars) ---
        if (time_varying_entrainment) then
            do i = 1, n_ent_segments
                call validate_entrainment_params(sched_ent_rate(i), sched_n_blob(i), &
                                                 sched_psigma(i))
            end do
            ! Seed current values from the first (t=0) segment before scheduling.
            active_segment = 1
            ent_rate = sched_ent_rate(1)
            n_blob   = sched_n_blob(1)
            psigma   = sched_psigma(1)
        else
            call validate_entrainment_params(ent_rate, n_blob, psigma)
        end if

        t_next_entrain = compute_dt_entm(vel)

        write(*,*) 'entrainment:       ON'
        if (time_varying_entrainment) then
            write(*,*) '  schedule:        time-varying (', n_ent_segments, ' segments)'
        end if
        write(*,*) '  ent_rate:        ', ent_rate
        write(*,*) '  n_blob:          ', n_blob
        write(*,*) '  psigma:          ', psigma
        write(*,*) '  random_entrain:  ', random_entrainment
        write(*,*) '  first dt_entm:   ', t_next_entrain, ' s'

    end subroutine initialize_entrainment


    ! Abort with a message if any entrainment parameter is out of range. Shared by
    ! the scalar (namelist) and per-segment (v3 schedule) validation paths.
    subroutine validate_entrainment_params(er, nb, ps)
        real(dp), intent(in) :: er    ! entrainment rate [1/m]
        integer(i4), intent(in) :: nb ! number of blobs
        real(dp), intent(in) :: ps    ! blob fraction

        if (ps <= 0.0 .or. ps >= 1.0) then
            write(error_unit,*) 'Error: psigma must be in (0, 1), got: ', ps
            call exit(1)
        end if
        if (ps * nb >= 1.0) then
            write(error_unit,*) 'Error: psigma * n_blob must be < 1, got psigma, n_blob: ', ps, nb
            call exit(1)
        end if
        if (er <= 0.0) then
            write(error_unit,*) 'Error: ent_rate must be > 0, got: ', er
            call exit(1)
        end if
        if (nb < 1 .or. nb > max_blobs) then
            write(error_unit,*) 'Error: n_blob must be between 1 and ', max_blobs, ', got: ', nb
            call exit(1)
        end if

    end subroutine validate_entrainment_params


    ! Index of the schedule segment active at query_time (piecewise-constant;
    ! structural twin of parcel's get_velocity).
    pure function schedule_segment(query_time) result(idx)
        real(dp), intent(in) :: query_time   ! simulation time (s)
        integer :: idx, i

        idx = 1
        do i = 2, n_ent_segments
            if (query_time < sched_times(i)) exit
            idx = i
        end do
    end function schedule_segment


    ! Refresh ent_rate/n_blob/psigma to the segment active at query_time. No-op
    ! when the schedule is not time-varying, leaving the namelist scalars in place.
    subroutine get_entrainment_params(query_time)
        real(dp), intent(in) :: query_time   ! simulation time (s)
        integer :: i

        if (.not. time_varying_entrainment) return

        i = schedule_segment(query_time)
        ent_rate = sched_ent_rate(i)
        n_blob   = sched_n_blob(i)
        psigma   = sched_psigma(i)

    end subroutine get_entrainment_params


    ! Track the active schedule segment; on a segment change, adopt the new
    ! parameters and redraw the next-event time with them, so a rate increase
    ! takes effect immediately instead of waiting out an interval drawn under
    ! the old rate. Exact for the Poisson case (exponential waiting times are
    ! memoryless); the deterministic cadence simply restarts at the boundary.
    ! With a near-zero velocity the redraw is skipped (no meaningful interval
    ! exists); the next segment change redraws again.
    subroutine refresh_entrainment_schedule(query_time, vel)
        real(dp), intent(in) :: query_time   ! simulation time (s)
        real(dp), intent(in) :: vel          ! parcel/eddy velocity [m/s]
        integer :: seg

        if (.not. time_varying_entrainment) return

        seg = schedule_segment(query_time)
        if (seg == active_segment) return

        active_segment = seg
        call get_entrainment_params(query_time)
        if (abs(vel) >= 1.0e-30) then
            t_next_entrain = query_time + compute_dt_entm(vel)
        end if

    end subroutine refresh_entrainment_schedule


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
        ! Blobs that wrap the periodic boundary split in two, so up to 2*max_blobs
        ! contiguous segments can come back from place_blobs.
        integer :: blob_start(2*max_blobs), blob_end(2*max_blobs), n_final
        integer :: i, k

        ! Before the timing checks, so a segment change redraws t_next_entrain
        ! and the new rate takes effect immediately.
        call refresh_entrainment_schedule(time, vel)

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
