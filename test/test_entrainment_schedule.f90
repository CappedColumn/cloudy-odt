! Unit test for the vertically-varying entrainment schedule (parcel input v3).
!
! Drives initialize_entrainment with per-segment schedules on both vertical
! axes and checks:
!   - piecewise-constant lookup of ent_rate/n_blob/psigma vs height (ascending
!     coordinates) and vs pressure (descending coordinates), including exact
!     behavior at segment breakpoints;
!   - refresh_entrainment_schedule redraws t_next_entrain with the *new*
!     segment's parameters immediately on a segment change (the redraw is
!     anchored at the current simulation time), and does not redraw within a
!     segment. With random_entrainment off the redrawn interval is the
!     closed-form deterministic value, so it can be checked exactly.
program test_entrainment_schedule
    use globals, only: dp, i4, namelist_path, time, AXIS_HEIGHT, AXIS_PRESSURE
    use entrainment, only: initialize_entrainment, get_entrainment_params, &
                           refresh_entrainment_schedule, t_next_entrain, &
                           ent_rate, n_blob, psigma
    implicit none

    integer :: n_passed, n_failed, nml_unit
    real(dp) :: seg_coords(3), seg_ent_rate(3), seg_psigma(3)
    integer(i4) :: seg_n_blob(3)
    real(dp) :: vel, q, tn_before

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

    ! ========================================================================
    ! Height axis: three segments, coordinates in metres (ascending).
    ! Schedule ent_rate values are internal units (1/m), as passed by the
    ! parcel reader after its 1/km -> 1/m conversion.
    ! ========================================================================
    seg_coords   = [0.0, 800.0, 1600.0]
    seg_ent_rate = [2.0e-3, 3.0e-3, 1.5e-3]
    seg_n_blob   = [1, 2, 1]
    seg_psigma   = [0.10, 0.05, 0.10]

    time = 0.0
    vel = 1.0
    call initialize_entrainment(vel, AXIS_HEIGHT, seg_coords, seg_ent_rate, &
                                seg_n_blob, seg_psigma)

    ! After init, current values should be seeded from segment 1.
    call check("height: init seeds segment 1", 1)

    ! Mid first segment.
    q = 400.0;    call get_entrainment_params(q); call check("height: z=400 -> segment 1", 1)
    ! Just below the first breakpoint stays in segment 1.
    q = 799.99;   call get_entrainment_params(q); call check("height: z=799.99 -> segment 1", 1)
    ! Exactly at the breakpoint switches to segment 2 (lookup uses query < start).
    q = 800.0;    call get_entrainment_params(q); call check("height: z=800 -> segment 2", 2)
    ! Mid second segment.
    q = 1200.0;   call get_entrainment_params(q); call check("height: z=1200 -> segment 2", 2)
    ! At and past the final breakpoint stays in the last segment.
    q = 1600.0;   call get_entrainment_params(q); call check("height: z=1600 -> segment 3", 3)
    q = 1.0e6;    call get_entrainment_params(q); call check("height: z=1e6 -> segment 3", 3)

    ! --- Redraw on segment change -------------------------------------------
    ! Crossing a segment boundary must redraw t_next_entrain with the new
    ! parameters immediately, not wait out the interval drawn under the old
    ! ones. The redraw is anchored at the current simulation *time* even though
    ! the schedule coordinate is vertical. (The lookup calls above left the
    ! module params on segment 3, but the tracked active segment is still 1
    ! from initialization, so refreshing inside segment 2 is a genuine change.)
    time = 100.0
    q = 1000.0
    call refresh_entrainment_schedule(q, vel)
    call check("height: refresh at z=1000 -> segment 2", 2)
    call check_t_next("height: redraw uses segment 2 params", 2)

    ! No boundary crossed: a second refresh in the same segment must not redraw.
    tn_before = t_next_entrain
    time = 120.0
    q = 1100.0
    call refresh_entrainment_schedule(q, vel)
    call check_no_redraw("height: no redraw within a segment", tn_before)

    ! Crossing into segment 3 redraws again with its parameters.
    time = 250.0
    q = 2000.0
    call refresh_entrainment_schedule(q, vel)
    call check("height: refresh at z=2000 -> segment 3", 3)
    call check_t_next("height: redraw uses segment 3 params", 3)

    ! ========================================================================
    ! Pressure axis: descending coordinates in Pa. Re-initializes the module
    ! schedule (re-init is supported for exactly this purpose).
    ! ========================================================================
    seg_coords   = [90000.0, 80000.0, 70000.0]
    seg_ent_rate = [1.0e-3, 4.0e-3, 2.0e-3]
    seg_n_blob   = [1, 1, 2]
    seg_psigma   = [0.10, 0.20, 0.05]

    time = 0.0
    call initialize_entrainment(vel, AXIS_PRESSURE, seg_coords, seg_ent_rate, &
                                seg_n_blob, seg_psigma)

    call check("pressure: init seeds segment 1", 1)

    ! Above (higher pressure than) the first coordinate clamps to segment 1.
    q = 95000.0;  call get_entrainment_params(q); call check("pressure: p=95000 -> segment 1", 1)
    q = 85000.0;  call get_entrainment_params(q); call check("pressure: p=85000 -> segment 1", 1)
    ! Exactly at the breakpoint switches to segment 2 (lookup uses query > start).
    q = 80000.0;  call get_entrainment_params(q); call check("pressure: p=80000 -> segment 2", 2)
    q = 75000.0;  call get_entrainment_params(q); call check("pressure: p=75000 -> segment 2", 2)
    q = 70000.0;  call get_entrainment_params(q); call check("pressure: p=70000 -> segment 3", 3)
    q = 10000.0;  call get_entrainment_params(q); call check("pressure: p=10000 -> segment 3", 3)

    ! Redraw when the parcel's falling pressure crosses into segment 2.
    time = 500.0
    q = 78000.0
    call refresh_entrainment_schedule(q, vel)
    call check("pressure: refresh at p=78000 -> segment 2", 2)
    call check_t_next("pressure: redraw uses segment 2 params", 2)

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
            write(*,'(a,es12.4,a,i0,a,f10.5)') '  got   ent_rate=', ent_rate, &
                ' n_blob=', n_blob, ' psigma=', psigma
            write(*,'(a,es12.4,a,i0,a,f10.5)') '  want  ent_rate=', seg_ent_rate(seg), &
                ' n_blob=', seg_n_blob(seg), ' psigma=', seg_psigma(seg)
        end if
    end subroutine check


    ! Verify t_next_entrain equals the deterministic interval for segment seg's
    ! parameters, anchored at the current simulation time (random_entrainment
    ! is off).
    subroutine check_t_next(label, seg)
        character(*), intent(in) :: label
        integer, intent(in) :: seg
        real(dp) :: expected

        expected = time + (real(seg_n_blob(seg), dp) / seg_ent_rate(seg)) &
                   * (seg_psigma(seg) / (1.0 - seg_psigma(seg))) / abs(vel)
        if (abs(t_next_entrain - expected) < 1.0e-9) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            write(*,'(a)') 'FAIL: '//label
            write(*,'(a,es15.7,a,es15.7)') '  got t_next=', t_next_entrain, &
                ' want ', expected
        end if
    end subroutine check_t_next


    subroutine check_no_redraw(label, tn)
        character(*), intent(in) :: label
        real(dp), intent(in) :: tn

        if (t_next_entrain == tn) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            write(*,'(a)') 'FAIL: '//label
        end if
    end subroutine check_no_redraw

end program test_entrainment_schedule
