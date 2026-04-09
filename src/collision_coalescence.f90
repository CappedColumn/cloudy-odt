module collision_coalescence
    ! Event-driven 1D collision-coalescence for Lagrangian cloud droplets.
    !
    ! Operates directly on type(particle) arrays. Within a time window dt,
    ! adjacent droplets (sorted by position) are checked for pair-passing events
    ! using a min-heap priority queue. When a pair meets:
    !   - With probability p = min(1, pi*(ri+rj)^2 / A), collision is accepted
    !   - On acceptance: merge the two particles (conserve water volume, liquid
    !     water mass, and solute mass)
    !   - Otherwise: swap adjacency order
    !
    ! Reference: coll-coal/collide_event_1d.f90, coll-coal/event_driven_1d_collision_high_level.pdf
    use globals, only: dp, i4, pi, pi_43, g, nu, rho_l, pres, Tref, Rd, N, H, domain_width, volume_scaling, &
                       simulation_mode
    use particle_types, only: particle, calculate_terminal_velocity
    implicit none
    private

    ! Namelist-controlled parameters
    logical, public :: do_collision_coalescence = .false.
    logical, public :: write_collisions = .false.
    real(dp), public :: wmax_collision = 10.0

    ! Per-step counters (reset each call to collision_coalescence_step)
    integer(i4), public :: collisions_this_step = 0
    integer(i4), public :: coalescences_this_step = 0
    integer(i4), public :: fall_events_this_step = 0

    ! Collision output file unit (unformatted stream binary)
    integer(i4) :: collision_unit

    public :: collision_coalescence_step, initialize_collision_file, seed_cc_rng

    ! Event types
    integer, parameter :: EV_PAIR = 1
    integer, parameter :: EV_FALL = 2

    type :: Event
        real(dp) :: t       ! event time
        integer  :: typ     ! EV_PAIR or EV_FALL
        integer  :: i, j    ! particle indices (j=0 for fall events)
    end type Event

    ! Independent xoshiro256** PRNG state (does not affect global random_number)
    integer(8) :: rng_s(4) = [123456789_8, 362436069_8, 521288629_8, 88675123_8]

contains

    subroutine collision_coalescence_step(lparticles, n_particles, ldt)
        ! Advance collision-coalescence for one time window ldt.
        ! Modifies particle radii and removes coalesced particles in-place.
        ! Positions are NOT written back — CODT settling is authoritative.
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(inout) :: n_particles
        real(dp), intent(in) :: ldt

        integer :: i, n, ncoll, head
        real(dp) :: A, tcur
        real(dp), allocatable :: zcur(:), w_fall(:), tstamp(:)
        logical, allocatable :: alive(:)
        integer, allocatable :: ord(:), prev(:), next(:)
        type(Event), allocatable :: heap(:)
        integer :: hsize
        type(Event) :: ev
        logical :: is_periodic

        n = n_particles
        collisions_this_step = 0
        coalescences_this_step = 0
        fall_events_this_step = 0

        ! Periodic boundary: parcel mode wraps via modulo(z, H); no fallout.
        is_periodic = (trim(simulation_mode) == 'parcel')

        if (n <= 1) then
            ! Nothing to collide — still advance settling so CC owns this path
            do i = 1, n
                call lparticles(i)%settling(ldt)
                if (is_periodic) lparticles(i)%position = modulo(lparticles(i)%position, H)
            end do
            return
        end if

        ! Cross-sectional area of the statistical volume
        A = domain_width**2 * volume_scaling
        if (A <= 0.0) return

        ! Allocate working arrays
        allocate(zcur(n), w_fall(n), tstamp(n), alive(n))
        allocate(ord(n), prev(n), next(n))
        allocate(heap(2*n))

        ! Initialize from particle data
        do i = 1, n
            zcur(i) = lparticles(i)%position
            w_fall(i) = abs(calculate_terminal_velocity(lparticles(i)))
            if (w_fall(i) > wmax_collision) w_fall(i) = wmax_collision
            tstamp(i) = 0.0
        end do

        alive = .true.
        tcur = 0.0
        ncoll = 0
        hsize = 0

        ! Build sorted ordering by position and linked list
        call build_order(zcur, n, ord)
        call build_links(ord, n, prev, next, head)

        ! Seed event queue with pair-passing events
        call push_all_pair_events(heap, hsize, zcur, tstamp, w_fall, alive, next, head, n, ldt, tcur)

        ! Seed fallout events (chamber mode only; parcel mode wraps via
        ! modulo at writeback and no particle exits the domain).
        if (.not. is_periodic) then
            do i = 1, n
                call push_fall_event(heap, hsize, i, zcur, tstamp, w_fall, alive, ldt, tcur)
            end do
        end if

        ! Main event loop
        do
            if (hsize <= 0) exit
            call heap_pop(heap, hsize, ev)
            if (ev%t > ldt) exit

            ! Validate event is still current
            if (ev%typ == EV_PAIR) then
                if (.not. is_pair_valid(ev%i, ev%j, alive, next, n)) cycle
            else if (ev%typ == EV_FALL) then
                if (ev%i < 1 .or. ev%i > n) cycle
                if (.not. alive(ev%i)) cycle
            end if

            tcur = ev%t

            if (ev%typ == EV_PAIR) then
                call handle_pair_event(ev%i, ev%j, tcur, lparticles, zcur, tstamp, w_fall, &
                                       alive, prev, next, head, n, A, ncoll, ldt, heap, hsize)
            else if (ev%typ == EV_FALL) then
                call handle_fall_event(ev%i, tcur, zcur, tstamp, w_fall, alive, prev, next, &
                                       head, n, ldt, heap, hsize)
            end if
        end do

        ! Position writeback: CC owns settling. Use update_working_pos to
        ! advance each alive particle's internal position (zcur) to t=ldt,
        ! accounting for piecewise velocity changes at merge events. Dead
        ! particles (from handle_fall_event) get z<0 so verify_particle_fallout
        ! removes them in the caller.
        do i = 1, n
            if (alive(i)) then
                call update_working_pos(i, ldt, zcur, tstamp, w_fall)
                if (is_periodic) then
                    lparticles(i)%position = modulo(zcur(i), H)
                else
                    lparticles(i)%position = zcur(i)
                end if
            else
                lparticles(i)%position = -1.0_dp
            end if
        end do

        deallocate(zcur, w_fall, tstamp, alive, ord, prev, next, heap)

    end subroutine collision_coalescence_step


    ! =========================================================================
    ! Event handlers
    ! =========================================================================

    subroutine handle_pair_event(i, j, tcur, lparticles, zcur, tstamp, w_fall, &
                                 alive, prev, next, head, n, A, ncoll, dt, heap, hsize)
        integer, intent(in) :: i, j, n
        real(dp), intent(in) :: tcur, A, dt
        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(inout) :: zcur(:), tstamp(:), w_fall(:)
        logical, intent(inout) :: alive(:)
        integer, intent(inout) :: prev(:), next(:), head, ncoll
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize

        real(dp) :: vrel, p_coll, u_coll, E_coal, r_keep, r_kill
        integer :: ip, jn, keep, kill

        call update_working_pos(i, tcur, zcur, tstamp, w_fall)
        call update_working_pos(j, tcur, zcur, tstamp, w_fall)

        vrel = w_fall(i) - w_fall(j)
        if (vrel <= 0.0) return

        ! Stage 1: Collision check (geometric kernel)
        p_coll = pi * (lparticles(i)%radius + lparticles(j)%radius)**2 / A
        if (p_coll > 1.0) p_coll = 1.0

        u_coll = cc_random()

        if (u_coll < p_coll) then
            ! Collision occurred
            collisions_this_step = collisions_this_step + 1
            lparticles(i)%n_collisions = lparticles(i)%n_collisions + 1
            lparticles(j)%n_collisions = lparticles(j)%n_collisions + 1

            ! Stage 2: Coalescence check (collection efficiency E)
            ! E = 1.0 for now (all collisions coalesce)
            E_coal = 0.0  ! TEMP: E=0 isolation test

            if (E_coal >= 1.0) then
                ! Coalescence: merge particles, conserving mass
                ncoll = ncoll + 1
                coalescences_this_step = coalescences_this_step + 1
                keep = i
                kill = j

                r_keep = lparticles(keep)%radius
                r_kill = lparticles(kill)%radius

                ! Record pre-coalescence state
                lparticles(keep)%radius_before_coalescence = r_keep
                lparticles(keep)%n_coalescences = lparticles(keep)%n_coalescences + 1

                ! Conserve water volume: r_new = (r1^3 + r2^3)^(1/3)
                lparticles(keep)%radius = (r_keep**3 + r_kill**3)**(1.0/3.0)

                ! Conserve liquid water mass
                lparticles(keep)%water_liquid = lparticles(keep)%water_liquid + lparticles(kill)%water_liquid

                ! Conserve solute mass and recompute solute radius
                lparticles(keep)%solute_gross_mass = lparticles(keep)%solute_gross_mass &
                                                    + lparticles(kill)%solute_gross_mass
                if (lparticles(keep)%solute_type%solute_density > 0.0) then
                    lparticles(keep)%solute_radius = &
                        (3.0 * lparticles(keep)%solute_gross_mass / &
                         (4.0 * pi * lparticles(keep)%solute_type%solute_density))**(1.0/3.0)
                end if

                ! Recalculate terminal velocity for merged droplet
                w_fall(keep) = abs(calculate_terminal_velocity(lparticles(keep)))
                if (w_fall(keep) > wmax_collision) w_fall(keep) = wmax_collision

                ! Update working position
                zcur(keep) = zcur(i)
                tstamp(keep) = tcur

                ! Write collision event log
                if (write_collisions) call write_collision( &
                    lparticles(keep)%particle_id, lparticles(kill)%particle_id, &
                    r_keep, r_kill, lparticles(keep)%radius, zcur(keep), tcur)

                ! Remove killed particle from linked list
                call unlink_particle(kill, alive, prev, next, head)
                alive(kill) = .false.

                ! Reschedule events for survivor and its neighbors
                ip = prev(keep)
                jn = next(keep)
                call push_pair_for_i(keep, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)
                if (ip > 0) call push_pair_for_i(ip, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)

                ! Reschedule fallout for survivor
                call push_fall_event(heap, hsize, keep, zcur, tstamp, w_fall, alive, dt, tcur)
                return
            end if
        end if

        ! No coalescence: swap adjacency order (they passed each other)
        call swap_adjacent(i, j, prev, next, head)

        ! Reschedule affected pairs
        ip = prev(j)  ! j is now before i after swap
        jn = next(i)
        if (ip > 0) call push_pair_for_i(ip, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)
        call push_pair_for_i(j, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)
        if (jn > 0) call push_pair_for_i(i, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)

    end subroutine handle_pair_event


    subroutine handle_fall_event(i, tcur, zcur, tstamp, w_fall, alive, prev, next, &
                                 head, n, dt, heap, hsize)
        integer, intent(in) :: i, n
        real(dp), intent(in) :: tcur, dt
        real(dp), intent(inout) :: zcur(:), tstamp(:), w_fall(:)
        logical, intent(inout) :: alive(:)
        integer, intent(inout) :: prev(:), next(:), head
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize

        real(dp) :: zi, tfall
        integer :: ip

        if (i < 1 .or. i > n) return
        if (.not. alive(i)) return

        call update_working_pos(i, tcur, zcur, tstamp, w_fall)
        zi = zcur(i)

        if (w_fall(i) <= 0.0) return

        ! Validate that fall time matches current time (stale check)
        tfall = tcur + zi / w_fall(i)
        if (abs(tfall - tcur) > 1.0e-8) return

        ! Remove particle
        ip = prev(i)
        call unlink_particle(i, alive, prev, next, head)
        alive(i) = .false.
        fall_events_this_step = fall_events_this_step + 1

        ! Reschedule neighbor pair
        if (ip > 0) call push_pair_for_i(ip, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)

    end subroutine handle_fall_event


    ! =========================================================================
    ! Event generation
    ! =========================================================================

    logical function is_pair_valid(i, j, alive, next, n)
        integer, intent(in) :: i, j, n
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:)
        is_pair_valid = .false.
        if (i < 1 .or. j < 1) return
        if (i > n .or. j > n) return
        if (.not. alive(i) .or. .not. alive(j)) return
        if (next(i) /= j) return
        is_pair_valid = .true.
    end function is_pair_valid


    subroutine push_all_pair_events(heap, hsize, zcur, tstamp, w_fall, alive, next, head, n, dt, tcur)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:), dt, tcur
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:), head, n
        integer :: i, count

        i = head
        count = 0
        do while (i /= 0 .and. count < n)
            call push_pair_for_i(i, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)
            i = next(i)
            count = count + 1
            if (i == 0) exit
        end do
    end subroutine push_all_pair_events


    subroutine push_pair_for_i(i, heap, hsize, zcur, tstamp, w_fall, alive, next, dt, tcur)
        integer, intent(in) :: i
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:), dt, tcur
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:)

        integer :: j
        real(dp) :: dz, vrel, te
        type(Event) :: ev

        if (i < 1 .or. i > size(alive)) return
        if (.not. alive(i)) return

        j = next(i)
        if (j <= 0) return
        if (.not. alive(j)) return

        call update_working_pos(i, tcur, zcur, tstamp, w_fall)
        call update_working_pos(j, tcur, zcur, tstamp, w_fall)

        ! Forward separation (j should be above i in sorted order)
        dz = zcur(j) - zcur(i)
        if (dz < 0.0) dz = dz + H

        vrel = w_fall(i) - w_fall(j)
        if (vrel <= 0.0) return

        te = tcur + dz / vrel
        if (te <= tcur) return
        if (te > dt) return

        ev%t = te
        ev%typ = EV_PAIR
        ev%i = i
        ev%j = j
        call heap_push(heap, hsize, ev)
    end subroutine push_pair_for_i


    subroutine push_fall_event(heap, hsize, i, zcur, tstamp, w_fall, alive, dt, tcur)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize
        integer, intent(in) :: i
        real(dp), intent(in) :: zcur(:), tstamp(:), w_fall(:), dt, tcur
        logical, intent(in) :: alive(:)
        type(Event) :: ev
        real(dp) :: zi, te

        if (i < 1 .or. i > size(alive)) return
        if (.not. alive(i)) return
        if (w_fall(i) <= 0.0) return

        ! Current position accounting for motion since last stamp
        zi = zcur(i) - w_fall(i) * (tcur - tstamp(i))
        if (zi <= 0.0) then
            te = tcur
        else
            te = tcur + zi / w_fall(i)
        end if
        if (te > dt) return

        ev%t = te
        ev%typ = EV_FALL
        ev%i = i
        ev%j = 0
        call heap_push(heap, hsize, ev)
    end subroutine push_fall_event


    ! =========================================================================
    ! Position updates
    ! =========================================================================

    subroutine update_working_pos(i, tcur, zcur, tstamp, w_fall)
        integer, intent(in) :: i
        real(dp), intent(in) :: tcur
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:)
        real(dp) :: dtloc

        dtloc = tcur - tstamp(i)
        if (dtloc /= 0.0) then
            zcur(i) = zcur(i) - w_fall(i) * dtloc
            tstamp(i) = tcur
        end if
    end subroutine update_working_pos


    ! =========================================================================
    ! Linked list operations
    ! =========================================================================

    subroutine build_links(ord, n, prev, next, head)
        integer, intent(in) :: ord(:), n
        integer, intent(out) :: prev(:), next(:), head
        integer :: k

        prev = 0
        next = 0

        do k = 1, n
            if (k > 1) prev(ord(k)) = ord(k-1)
            if (k < n) next(ord(k)) = ord(k+1)
        end do

        head = ord(1)
        ! Non-periodic: prev(head)=0, next(tail)=0 (already set)
    end subroutine build_links


    subroutine unlink_particle(kill, alive, prev, next, head)
        integer, intent(in) :: kill
        logical, intent(inout) :: alive(:)
        integer, intent(inout) :: prev(:), next(:), head
        integer :: ip, in_

        if (kill < 1 .or. kill > size(alive)) return
        if (.not. alive(kill)) return

        ip = prev(kill)
        in_ = next(kill)

        if (ip > 0) next(ip) = in_
        if (in_ > 0) prev(in_) = ip
        if (head == kill) head = in_

        prev(kill) = 0
        next(kill) = 0
    end subroutine unlink_particle


    subroutine swap_adjacent(i, j, prev, next, head)
        ! Swap adjacent nodes where next(i)=j: ... a-i-j-d ... -> ... a-j-i-d ...
        integer, intent(in) :: i, j
        integer, intent(inout) :: prev(:), next(:), head
        integer :: a, d

        a = prev(i)
        d = next(j)

        if (a /= 0) next(a) = j
        prev(j) = a

        next(j) = i
        prev(i) = j

        next(i) = d
        if (d /= 0) prev(d) = i

        if (head == i) head = j
    end subroutine swap_adjacent


    ! =========================================================================
    ! Sorting
    ! =========================================================================

    subroutine build_order(zcur, n, ord)
        real(dp), intent(in) :: zcur(:)
        integer, intent(in) :: n
        integer, intent(out) :: ord(:)
        integer :: i

        do i = 1, n
            ord(i) = i
        end do
        call sort_by_position(ord, zcur, 1, n)
    end subroutine build_order


    recursive subroutine sort_by_position(ord, key, lo, hi)
        integer, intent(inout) :: ord(:)
        real(dp), intent(in) :: key(:)
        integer, intent(in) :: lo, hi
        integer :: i, j, tmp
        real(dp) :: pval

        if (lo >= hi) return
        pval = key(ord((lo+hi)/2))
        i = lo
        j = hi
        do
            do while (key(ord(i)) < pval)
                i = i + 1
            end do
            do while (key(ord(j)) > pval)
                j = j - 1
            end do
            if (i <= j) then
                tmp = ord(i); ord(i) = ord(j); ord(j) = tmp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (lo < j) call sort_by_position(ord, key, lo, j)
        if (i < hi) call sort_by_position(ord, key, i, hi)
    end subroutine sort_by_position


    ! =========================================================================
    ! Min-heap (priority queue by event time)
    ! =========================================================================

    subroutine heap_push(heap, hsize, ev)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize
        type(Event), intent(in) :: ev
        integer :: k, parent_idx
        type(Event) :: tmp

        hsize = hsize + 1
        if (hsize > size(heap)) then
            ! Silently drop event if heap is full (should not happen in practice)
            hsize = hsize - 1
            return
        end if
        heap(hsize) = ev

        ! Sift up
        k = hsize
        do while (k > 1)
            parent_idx = k / 2
            if (heap(parent_idx)%t <= heap(k)%t) exit
            tmp = heap(parent_idx)
            heap(parent_idx) = heap(k)
            heap(k) = tmp
            k = parent_idx
        end do
    end subroutine heap_push


    subroutine heap_pop(heap, hsize, ev)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: hsize
        type(Event), intent(out) :: ev
        integer :: k, left, right, smallest
        type(Event) :: tmp

        if (hsize <= 0) then
            ev%t = huge(1.0_dp)
            ev%typ = 0
            ev%i = 0
            ev%j = 0
            return
        end if

        ev = heap(1)
        heap(1) = heap(hsize)
        hsize = hsize - 1

        ! Sift down
        k = 1
        do
            left = 2 * k
            right = 2 * k + 1
            smallest = k
            if (left <= hsize .and. heap(left)%t < heap(smallest)%t) smallest = left
            if (right <= hsize .and. heap(right)%t < heap(smallest)%t) smallest = right
            if (smallest == k) exit
            tmp = heap(k)
            heap(k) = heap(smallest)
            heap(smallest) = tmp
            k = smallest
        end do
    end subroutine heap_pop


    ! =========================================================================
    ! Collision event binary output
    ! =========================================================================

    subroutine initialize_collision_file(filename)
        character(*), intent(in) :: filename
        integer(i4) :: ierr

        open(newunit=collision_unit, file=trim(filename)//'_collisions.bin', &
             form='unformatted', access='stream', status='replace', iostat=ierr)
        if (ierr /= 0) then
            print *, "Error opening collision data file."
            stop 1
        end if

        ! Header: grid size, domain height, domain width, volume scaling
        write(collision_unit) N, H, domain_width, volume_scaling

    end subroutine initialize_collision_file


    subroutine write_collision(id_keep, id_kill, r_keep_m, r_kill_m, r_after_m, position_m, time_s)
        integer(i4), intent(in) :: id_keep, id_kill
        real(dp), intent(in) :: r_keep_m, r_kill_m, r_after_m, position_m, time_s

        write(collision_unit) id_keep, id_kill, r_keep_m, r_kill_m, r_after_m, position_m, time_s

    end subroutine write_collision


    ! =========================================================================
    ! Independent PRNG (xoshiro256** — does not affect global random_number)
    ! =========================================================================

    subroutine seed_cc_rng(seed_val)
        ! Seed the CC PRNG from a single integer. Uses splitmix64 to fill state.
        integer(8), intent(in) :: seed_val
        integer(8) :: z
        integer :: k

        z = seed_val
        do k = 1, 4
            z = z + 6364136223846793005_8
            z = ieor(z, ishft(z, -30)) * (-4658895280553007687_8)
            z = ieor(z, ishft(z, -27)) * (-7723592293110705685_8)
            z = ieor(z, ishft(z, -31))
            rng_s(k) = z
        end do
    end subroutine seed_cc_rng


    function cc_random() result(u)
        ! Return a uniform random number in [0, 1) from the CC-private PRNG.
        real(dp) :: u
        integer(8) :: res, t

        res = rng_s(2) * 5_8
        res = ior(ishft(res, 7), ishft(res, -57)) * 9_8

        t = ishft(rng_s(2), 17)

        rng_s(3) = ieor(rng_s(3), rng_s(1))
        rng_s(4) = ieor(rng_s(4), rng_s(2))
        rng_s(2) = ieor(rng_s(2), rng_s(3))
        rng_s(1) = ieor(rng_s(1), rng_s(4))

        rng_s(3) = ieor(rng_s(3), t)
        rng_s(4) = ior(ishft(rng_s(4), 45), ishft(rng_s(4), -19))

        ! Convert to [0, 1): use upper 53 bits for double precision
        u = real(ishft(res, -11), dp) * (1.0 / 9007199254740992.0)
    end function cc_random

end module collision_coalescence
