! Unit test for the per-leg entrainment schedule (parcel input v3).
!
! Drives initialize_entrainment with a per-leg schedule and checks:
!   - lookup: get_entrainment_params(leg) exposes exactly that leg's
!     ent_rate/n_blob/psigma;
!   - refresh_entrainment_schedule redraws t_next_entrain with the *new* leg's
!     parameters immediately on a leg change (anchored at the current
!     simulation time), and does not redraw while the leg is unchanged. With
!     random_entrainment off the redrawn interval is the closed-form
!     deterministic value, so it can be checked exactly.
program test_entrainment_schedule
    use globals, only: dp, i4, namelist_path, time
    use entrainment, only: initialize_entrainment, get_entrainment_params, &
                           refresh_entrainment_schedule, t_next_entrain, &
                           ent_rate, n_blob, psigma
    implicit none

    integer :: n_passed, n_failed, nml_unit
    real(dp) :: leg_ent_rate(3), leg_psigma(3)
    integer(i4) :: leg_n_blob(3)
    real(dp) :: vel, tn_before

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

    ! --- Three-leg schedule (ent_rate in internal 1/m, as passed by the
    !     parcel reader after its 1/km -> 1/m conversion) ---
    leg_ent_rate = [2.0e-3, 3.0e-3, 1.5e-3]
    leg_n_blob   = [1, 2, 1]
    leg_psigma   = [0.10, 0.05, 0.10]

    time = 0.0
    vel = 1.0
    call initialize_entrainment(vel, leg_ent_rate, leg_n_blob, leg_psigma)

    ! After init, current values should be seeded from leg 1.
    call check("init seeds leg 1", 1)

    ! Direct lookups follow the requested leg exactly.
    call get_entrainment_params(2_i4); call check("lookup leg 2", 2)
    call get_entrainment_params(3_i4); call check("lookup leg 3", 3)
    call get_entrainment_params(1_i4); call check("lookup leg 1", 1)

    ! --- Redraw on leg change --------------------------------------------
    ! Advancing to a new leg must redraw t_next_entrain with the new leg's
    ! parameters immediately, not wait out the interval drawn under the old
    ! ones. (The tracked active leg is still 1 from initialization.)
    time = 100.0
    call refresh_entrainment_schedule(2_i4, vel)
    call check("refresh to leg 2 adopts its params", 2)
    call check_t_next("redraw uses leg 2 params", 2)

    ! Same leg again: no redraw.
    tn_before = t_next_entrain
    time = 120.0
    call refresh_entrainment_schedule(2_i4, vel)
    if (t_next_entrain == tn_before) then
        n_passed = n_passed + 1
    else
        n_failed = n_failed + 1
        write(*,'(a)') 'FAIL: no redraw within a leg'
    end if

    ! Advancing to leg 3 redraws again with its parameters.
    time = 250.0
    call refresh_entrainment_schedule(3_i4, vel)
    call check("refresh to leg 3 adopts its params", 3)
    call check_t_next("redraw uses leg 3 params", 3)

    ! --- Cleanup ---
    open(newunit=nml_unit, file=namelist_path, status='old')
    close(nml_unit, status='delete')

    write(*,'(a,i0,a,i0,a)') 'test_entrainment_schedule: ', n_passed, ' passed, ', &
                             n_failed, ' failed'
    if (n_failed > 0) call exit(1)

contains

    ! Compare the current entrainment parameters against the expected leg.
    ! Both sides come from the same dp schedule arrays, so equality is exact.
    subroutine check(label, leg)
        character(*), intent(in) :: label
        integer, intent(in) :: leg
        logical :: ok

        ok = ent_rate == leg_ent_rate(leg) .and. &
             n_blob == leg_n_blob(leg) .and. &
             psigma == leg_psigma(leg)
        if (ok) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            write(*,'(a)') 'FAIL: '//label
            write(*,'(a,es12.4,a,i0,a,f10.5)') '  got   ent_rate=', ent_rate, &
                ' n_blob=', n_blob, ' psigma=', psigma
            write(*,'(a,es12.4,a,i0,a,f10.5)') '  want  ent_rate=', leg_ent_rate(leg), &
                ' n_blob=', leg_n_blob(leg), ' psigma=', leg_psigma(leg)
        end if
    end subroutine check


    ! Verify t_next_entrain equals the deterministic interval for leg `leg`'s
    ! parameters, anchored at the current simulation time (random_entrainment
    ! is off).
    !
    ! Deliberately independent of n_blob: psigma is the total fraction replaced
    ! per event, so event spacing depends on psigma alone and n_blob only
    ! subdivides that volume spatially. Leg 2 (n_blob = 2) is the case that
    ! distinguishes this from the pre-3.1.0 formula, which carried an n_blob
    ! factor here and so doubled leg 2's interval.
    subroutine check_t_next(label, leg)
        character(*), intent(in) :: label
        integer, intent(in) :: leg
        real(dp) :: expected

        expected = time + (leg_psigma(leg) / (1.0 - leg_psigma(leg))) &
                   / (leg_ent_rate(leg) * abs(vel))
        if (abs(t_next_entrain - expected) < 1.0e-9) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            write(*,'(a)') 'FAIL: '//label
            write(*,'(a,es15.7,a,es15.7)') '  got t_next=', t_next_entrain, &
                ' want ', expected
        end if
    end subroutine check_t_next

end program test_entrainment_schedule
