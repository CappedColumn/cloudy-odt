! Unit test for the time-varying entrainment schedule lookup
! (get_entrainment_params). Drives initialize_entrainment with a per-segment
! schedule, then checks that the active ent_rate/n_blob/psigma are the
! piecewise-constant values for a range of query times — including the segment
! breakpoints, where the lookup must switch exactly at t == segment start.
program test_entrainment_schedule
    use globals, only: dp, i4, namelist_path, time
    use entrainment, only: initialize_entrainment, get_entrainment_params, &
                           ent_rate, n_blob, psigma
    implicit none

    integer :: n_passed, n_failed, nml_unit
    real(dp) :: seg_times(3), seg_ent_rate(3), seg_psigma(3)
    integer(i4) :: seg_n_blob(3)
    real(dp) :: vel, q

    n_passed = 0
    n_failed = 0

    ! --- Minimal namelist so initialize_entrainment can read &ENTRAINMENT ---
    ! (random_entrainment off keeps compute_dt_entm deterministic; the scalar
    !  values are overridden by the schedule below.)
    namelist_path = 'test_entrainment_schedule.nml'
    open(newunit=nml_unit, file=namelist_path, status='replace', action='write')
    write(nml_unit, '(a)') '&ENTRAINMENT'
    write(nml_unit, '(a)') '  random_entrainment = .false.'
    write(nml_unit, '(a)') '/'
    close(nml_unit)

    ! --- Schedule: three segments on the shared time axis ---
    seg_times    = [0.0, 80.0, 160.0]
    seg_ent_rate = [2.0, 3.0, 1.5]
    seg_n_blob   = [1, 2, 1]
    seg_psigma   = [0.10, 0.05, 0.10]

    time = 0.0
    vel = 1.0
    call initialize_entrainment(vel, seg_times, seg_ent_rate, seg_n_blob, seg_psigma)

    ! After init, current values should be seeded from segment 1 (t=0).
    call check("init seeds segment 1", 1)

    ! Mid first segment.
    q = 40.0;    call get_entrainment_params(q); call check("t=40 -> segment 1", 1)
    ! Just before the first breakpoint stays in segment 1.
    q = 79.999;  call get_entrainment_params(q); call check("t=79.999 -> segment 1", 1)
    ! Exactly at the breakpoint switches to segment 2 (lookup uses query < start).
    q = 80.0;    call get_entrainment_params(q); call check("t=80 -> segment 2", 2)
    ! Mid second segment.
    q = 120.0;   call get_entrainment_params(q); call check("t=120 -> segment 2", 2)
    ! At and past the final breakpoint stays in the last segment.
    q = 160.0;   call get_entrainment_params(q); call check("t=160 -> segment 3", 3)
    q = 1.0e6;   call get_entrainment_params(q); call check("t=1e6 -> segment 3", 3)

    ! --- Cleanup ---
    open(newunit=nml_unit, file=namelist_path, status='old')
    close(nml_unit, status='delete')

    write(*,'(a,i0,a,i0,a)') 'test_entrainment_schedule: ', n_passed, ' passed, ', &
                             n_failed, ' failed'
    if (n_failed > 0) call exit(1)

contains

    ! Compare the current entrainment parameters against the expected segment.
    ! Both sides come from the same dp schedule arrays, so equality is exact.
    subroutine check(label, seg)
        character(*), intent(in) :: label
        integer, intent(in) :: seg
        logical :: ok

        ok = ent_rate == seg_ent_rate(seg) .and. &
             n_blob == seg_n_blob(seg) .and. &
             psigma == seg_psigma(seg)
        if (ok) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            write(*,'(a)') 'FAIL: '//label
            write(*,'(a,f10.5,a,i0,a,f10.5)') '  got   ent_rate=', ent_rate, &
                ' n_blob=', n_blob, ' psigma=', psigma
            write(*,'(a,f10.5,a,i0,a,f10.5)') '  want  ent_rate=', seg_ent_rate(seg), &
                ' n_blob=', seg_n_blob(seg), ' psigma=', seg_psigma(seg)
        end if
    end subroutine check

end program test_entrainment_schedule
