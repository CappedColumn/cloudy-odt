program test_collision_coalescence
    ! Standalone unit test for the event-driven collision-coalescence step.
    !
    ! Drives collision_coalescence_step directly on synthetic particle arrays,
    ! isolating pure collision-coalescence: no condensational growth (settling
    ! and microphysics are never invoked) and no fallout (parcel mode wraps
    ! periodically). This is the in-memory companion to the full-executable
    ! integration checks in test/scripts/test_cc_invariants.py.
    !
    ! The scenario is made deterministic (no dependence on the RNG sequence):
    !   - parcel mode           -> no fall events
    !   - tiny grid_area        -> geometric collision probability clamps to 1
    !   - 'unity' kernel        -> collection efficiency E == 1
    ! so every scheduled adjacent meeting (upper droplet falling faster than the
    ! lower one) results in a coalescence.
    use globals, only: dp, i4, H, N, simulation_mode, volume_scaling, pi_43, rho_l
    use particle_types, only: particle, calculate_terminal_velocity
    use collision_coalescence, only: collision_coalescence_step, &
        do_collisions, do_coalescence, write_collisions, &
        collisions_this_step, coalescences_this_step, fall_events_this_step
    use collection_efficiency, only: coalescence_kernel, set_kernel_selector
    implicit none

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    ! Global configuration shared by all subtests. grid_area inside the module is
    ! set once from domain_width**2 * volume_scaling on the first step call, so
    ! volume_scaling must be small enough to force p_coll -> 1 before any call.
    simulation_mode = 'parcel'      ! periodic: no fallout, droplets only merge
    N = 100                         ! gridcell count (unused by CC, set for sanity)
    volume_scaling = 1.0e-3_dp      ! grid_area = 1e-6 * 1e-3 = 1e-9 m^2 -> p_coll==1
    coalescence_kernel = 'unity'    ! collection efficiency E == 1 everywhere
    call set_kernel_selector()
    write_collisions = .false.

    call test_two_droplet_merge()
    call test_collision_without_coalescence()
    call test_separating_pair_does_not_collide()
    call test_pair_across_periodic_boundary()
    call test_two_droplet_laps()
    call test_dsd_box_conservation()

    print *
    if (n_failed > 0) then
        print '(a,i0,a,i0,a)', ' FAIL: ', n_failed, ' of ', n_passed + n_failed, ' checks failed'
        stop 1
    end if
    print '(a,i0,a)', ' PASS: all ', n_passed, ' checks passed'

contains

    ! -----------------------------------------------------------------------
    ! Subtest 1: a single deterministic two-droplet merge
    ! -----------------------------------------------------------------------
    subroutine test_two_droplet_merge()
        ! A large (fast) droplet placed just above a small (slow) one. The upper
        ! droplet catches the lower within the window and they coalesce. Verify
        ! the survivor conserves water volume, liquid mass and solute mass, and
        ! that the killed droplet is flagged out.
        type(particle) :: p(2)
        real(dp) :: v_pre, m_pre, s_pre, r_expected
        integer(i4) :: np

        H = 1.0_dp
        do_collisions = .true.
        do_coalescence = .true.

        ! Upper droplet (index 1, larger -> faster) and lower droplet (index 2).
        p(1) = make_droplet(1_i4, 30.0e-6_dp, 0.12_dp)
        p(2) = make_droplet(2_i4, 15.0e-6_dp, 0.10_dp)

        v_pre = p(1)%radius**3 + p(2)%radius**3
        m_pre = p(1)%water_liquid + p(2)%water_liquid
        s_pre = p(1)%solute_gross_mass + p(2)%solute_gross_mass
        r_expected = v_pre**(1.0_dp/3.0_dp)

        np = 2
        ! ltime is the window END, so a first window of 0.5 s ends at 0.5 s.
        call collision_coalescence_step(p, np, 0.5_dp, 0.5_dp)

        call check_int("two-drop: exactly one coalescence", coalescences_this_step, 1)
        call check_int("two-drop: one collision recorded", collisions_this_step, 1)
        call check_int("two-drop: no fall events (parcel)", fall_events_this_step, 0)

        ! handle_pair_event keeps the upper, larger droplet (index 1) and kills
        ! the lower, smaller one.
        call check_true("two-drop: survivor (1) not flagged coalesced", .not. p(1)%coalesced)
        call check_true("two-drop: killed (2) flagged coalesced", p(2)%coalesced)

        call check_close("two-drop: survivor radius = (r1^3+r2^3)^(1/3)", &
            p(1)%radius, r_expected, 1.0e-12_dp)
        call check_close("two-drop: water volume conserved (r^3)", &
            p(1)%radius**3, v_pre, 1.0e-12_dp)
        call check_close("two-drop: liquid water mass conserved", &
            p(1)%water_liquid, m_pre, 1.0e-12_dp)
        call check_close("two-drop: solute mass conserved", &
            p(1)%solute_gross_mass, s_pre, 1.0e-12_dp)
        call check_true("two-drop: DSD broadened (r grew)", p(1)%radius > 30.0e-6_dp)
    end subroutine test_two_droplet_merge


    ! -----------------------------------------------------------------------
    ! Subtest 2: collisions counted but do_coalescence off -> no merge
    ! -----------------------------------------------------------------------
    subroutine test_collision_without_coalescence()
        ! Same geometry, but coalescence disabled. The pair should register a
        ! collision and then pass (swap order) with both droplets unchanged.
        type(particle) :: p(2)
        real(dp) :: r1_pre, r2_pre
        integer(i4) :: np

        H = 1.0_dp
        do_collisions = .true.
        do_coalescence = .false.

        p(1) = make_droplet(1_i4, 30.0e-6_dp, 0.12_dp)
        p(2) = make_droplet(2_i4, 15.0e-6_dp, 0.10_dp)
        r1_pre = p(1)%radius
        r2_pre = p(2)%radius

        np = 2_i4
        call collision_coalescence_step(p, np, 0.5_dp, 0.5_dp)

        call check_true("no-coal: at least one collision counted", collisions_this_step >= 1)
        call check_int("no-coal: zero coalescences", coalescences_this_step, 0)
        call check_true("no-coal: droplet 1 survives", .not. p(1)%coalesced)
        call check_true("no-coal: droplet 2 survives", .not. p(2)%coalesced)
        call check_close("no-coal: radius 1 unchanged", p(1)%radius, r1_pre, 1.0e-15_dp)
        call check_close("no-coal: radius 2 unchanged", p(2)%radius, r2_pre, 1.0e-15_dp)
    end subroutine test_collision_without_coalescence


    ! -----------------------------------------------------------------------
    ! Subtest 3: a separating pair must not collide
    ! -----------------------------------------------------------------------
    subroutine test_separating_pair_does_not_collide()
        ! The two-droplet geometry reversed: the large (fast) droplet starts
        ! below the small (slow) one, so the gap between them grows. Nothing
        ! may collide within the window. Regression guard for the scheduler
        ! sign error fixed in d62c592, which queued exactly this pair.
        type(particle) :: p(2)
        integer(i4) :: np

        H = 1.0_dp
        do_collisions = .true.
        do_coalescence = .true.

        p(1) = make_droplet(1_i4, 30.0e-6_dp, 0.10_dp)
        p(2) = make_droplet(2_i4, 15.0e-6_dp, 0.12_dp)

        np = 2_i4
        call collision_coalescence_step(p, np, 0.5_dp, 0.5_dp)

        call check_int("separating: no collision", collisions_this_step, 0)
        call check_int("separating: no coalescence", coalescences_this_step, 0)
        call check_true("separating: droplet 1 survives", .not. p(1)%coalesced)
        call check_true("separating: droplet 2 survives", .not. p(2)%coalesced)
    end subroutine test_separating_pair_does_not_collide


    ! -----------------------------------------------------------------------
    ! Subtest 4: a pair that meets across the periodic boundary
    ! -----------------------------------------------------------------------
    subroutine test_pair_across_periodic_boundary()
        ! Parcel mode is periodic. A large (fast) droplet near the bottom leaves
        ! through the floor, re-enters at the top, and catches a small (slow)
        ! droplet near the top: the two are 0.03 m apart through the boundary
        ! and close at about 0.08 m/s, so they meet near 0.37 s. A neighbour
        ! list with two ends misses this pair; the ring finds it.
        type(particle) :: p(2)
        real(dp) :: v_pre
        integer(i4) :: np

        H = 1.0_dp
        do_collisions = .true.
        do_coalescence = .true.

        p(1) = make_droplet(1_i4, 30.0e-6_dp, 0.02_dp)
        p(2) = make_droplet(2_i4, 15.0e-6_dp, 0.99_dp)
        v_pre = p(1)%radius**3 + p(2)%radius**3

        np = 2_i4
        call collision_coalescence_step(p, np, 0.5_dp, 0.5_dp)

        call check_int("boundary: one collision across the boundary", collisions_this_step, 1)
        call check_int("boundary: one coalescence", coalescences_this_step, 1)
        call check_true("boundary: larger droplet (1) survives", .not. p(1)%coalesced)
        call check_true("boundary: smaller droplet (2) flagged coalesced", p(2)%coalesced)
        call check_close("boundary: water volume conserved (r^3)", p(1)%radius**3, v_pre, 1.0e-12_dp)
        call check_true("boundary: survivor position wrapped into [0, H)", &
            p(1)%position >= 0.0_dp .and. p(1)%position < H)
    end subroutine test_pair_across_periodic_boundary


    ! -----------------------------------------------------------------------
    ! Subtest 5: two droplets lapping each other
    ! -----------------------------------------------------------------------
    subroutine test_two_droplet_laps()
        ! Two droplets alone on the periodic ring, coalescence off, in one long
        ! window. The fast one passes the slow one, then laps the column and
        ! passes it again every H / (closing speed). Every pass is a collision
        ! here (p_coll == 1), so the count is known exactly.
        type(particle) :: p(2)
        real(dp), parameter :: window = 5.0_dp
        real(dp), parameter :: gap = 0.02_dp      ! fast droplet starts this far above the slow one
        real(dp) :: closing, t_first, lap_time
        integer :: expected_passes
        integer(i4) :: np

        H = 0.1_dp
        do_collisions = .true.
        do_coalescence = .false.

        p(1) = make_droplet(1_i4, 30.0e-6_dp, 0.05_dp + gap)
        p(2) = make_droplet(2_i4, 15.0e-6_dp, 0.05_dp)

        closing = abs(calculate_terminal_velocity(p(1))) - abs(calculate_terminal_velocity(p(2)))
        t_first = gap / closing
        lap_time = H / closing
        expected_passes = 1 + int((window - t_first) / lap_time)

        np = 2_i4
        call collision_coalescence_step(p, np, window, window)

        call check_true("laps: more than one pass expected in the window", expected_passes > 1)
        call check_int("laps: one collision per pass", collisions_this_step, expected_passes)
        call check_true("laps: both droplets survive", .not. p(1)%coalesced .and. .not. p(2)%coalesced)
        call check_true("laps: positions wrapped into [0, H)", &
            all(p%position >= 0.0_dp) .and. all(p%position < H))
    end subroutine test_two_droplet_laps


    ! -----------------------------------------------------------------------
    ! Subtest 6: many-droplet box -- DSD evolution under pure coalescence
    ! -----------------------------------------------------------------------
    subroutine test_dsd_box_conservation()
        ! A bidisperse population -- large "collector" droplets and small
        ! "collectee" droplets -- repeatedly collide and coalesce over many
        ! windows. Across the whole run, total water volume, liquid mass and
        ! solute mass must be conserved, the droplet count must be non-increasing,
        ! and the distribution must broaden (max radius grows). Coalesced droplets
        ! are compacted out between windows so they cannot re-collide.
        integer(i4), parameter :: n_large = 8     ! large collector droplets
        integer(i4), parameter :: n_small = 12    ! small collectee droplets
        integer(i4), parameter :: n0 = n_large + n_small
        integer(i4), parameter :: n_steps = 300
        real(dp), parameter :: dt = 0.05_dp
        real(dp), parameter :: r_large = 20.0e-6_dp
        real(dp), parameter :: r_small = 10.0e-6_dp
        type(particle) :: p(n0)
        real(dp) :: v0, m0, s0, vN, mN, sN, rmax0, rmaxN
        integer(i4) :: np, k, step

        H = 0.1_dp                  ! compact domain so adjacent pairs meet within dt
        do_collisions = .true.
        do_coalescence = .true.

        ! Two sizes only. Collectors (larger -> faster) fill the bottom of the
        ! domain, collectees (smaller -> slower) sit above, so the fast droplets
        ! sweep up the slow ones across periodic windows.
        do k = 1, n0
            if (k <= n_large) then
                p(k) = make_droplet(k, r_large, (real(k, dp) - 0.5_dp) * H / real(n0, dp))
            else
                p(k) = make_droplet(k, r_small, (real(k, dp) - 0.5_dp) * H / real(n0, dp))
            end if
        end do
        np = n0

        v0 = sum(p(1:np)%radius**3)
        m0 = sum(p(1:np)%water_liquid)
        s0 = sum(p(1:np)%solute_gross_mass)
        rmax0 = maxval(p(1:np)%radius)

        do step = 1, n_steps
            if (np <= 1) exit
            call collision_coalescence_step(p, np, dt, real(step, dp) * dt)
            call compact_survivors(p, np)
        end do

        vN = sum(p(1:np)%radius**3)
        mN = sum(p(1:np)%water_liquid)
        sN = sum(p(1:np)%solute_gross_mass)
        rmaxN = maxval(p(1:np)%radius)

        call check_true("box: at least one coalescence occurred", np < n0)
        call check_close("box: total water volume conserved", vN, v0, 1.0e-10_dp)
        call check_close("box: total liquid mass conserved", mN, m0, 1.0e-10_dp)
        call check_close("box: total solute mass conserved", sN, s0, 1.0e-10_dp)
        call check_true("box: DSD broadened (max radius grew)", rmaxN > rmax0)

        print '(a,i0,a,i0,a,es9.2,a,es9.2,a)', &
            '         (box: ', n0, ' -> ', np, ' droplets, r_max ', &
            rmax0, ' -> ', rmaxN, ' m)'
    end subroutine test_dsd_box_conservation


    ! -----------------------------------------------------------------------
    ! Helpers
    ! -----------------------------------------------------------------------
    function make_droplet(id, r, z) result(p)
        ! Minimal droplet with a finite terminal velocity (virt_temp > 0 avoids
        ! a divide-by-zero in the Stokes air-density term) and self-consistent
        ! liquid water and solute masses for conservation checks.
        integer(i4), intent(in) :: id
        real(dp), intent(in) :: r, z
        type(particle) :: p

        p%particle_id = id
        p%radius = r
        p%position = z
        p%virt_temp = 290.0_dp
        p%water_liquid = pi_43 * rho_l * r**3
        p%solute_gross_mass = 1.0e-18_dp
        p%solute_type%solute_density = 2160.0_dp   ! NaCl, kg/m^3
        p%coalesced = .false.
    end function make_droplet


    subroutine compact_survivors(p, np)
        ! Drop coalesced droplets (flagged by collision_coalescence_step) so the
        ! first np entries are the live population for the next window.
        type(particle), intent(inout) :: p(:)
        integer(i4), intent(inout) :: np
        integer(i4) :: src, dst

        dst = 0
        do src = 1, np
            if (.not. p(src)%coalesced) then
                dst = dst + 1
                if (dst /= src) p(dst) = p(src)
            end if
        end do
        np = dst
    end subroutine compact_survivors


    subroutine check_true(label, cond)
        character(*), intent(in) :: label
        logical, intent(in) :: cond

        if (cond) then
            n_passed = n_passed + 1
            print '(a,a)', '  PASS  ', label
        else
            n_failed = n_failed + 1
            print '(a,a)', '  FAIL  ', label
        end if
    end subroutine check_true


    subroutine check_int(label, got, expected)
        character(*), intent(in) :: label
        integer(i4), intent(in) :: got
        integer, intent(in) :: expected

        if (got == expected) then
            n_passed = n_passed + 1
            print '(a,a,a,i0,a)', '  PASS  ', label, '  (', got, ')'
        else
            n_failed = n_failed + 1
            print '(a,a,a,i0,a,i0,a)', '  FAIL  ', label, '  (got ', got, ', expected ', expected, ')'
        end if
    end subroutine check_int


    subroutine check_close(label, got, expected, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, rtol
        real(dp) :: relerr

        relerr = abs(got - expected) / max(abs(expected), 1.0e-30_dp)
        if (relerr <= rtol) then
            n_passed = n_passed + 1
            print '(a,a,a,es10.3,a)', '  PASS  ', label, '  (rel err ', relerr, ')'
        else
            n_failed = n_failed + 1
            print '(a,a,a,es12.5,a,es12.5,a,es10.3,a)', '  FAIL  ', label, &
                '  (got ', got, ', expected ', expected, ', rel err ', relerr, ')'
        end if
    end subroutine check_close

end program test_collision_coalescence
