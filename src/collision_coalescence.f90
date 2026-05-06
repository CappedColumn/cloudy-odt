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
    use globals, only: dp, i4, i1, pi, pi_43, g, nu, rho_l, pres, Tref, Rd, N, H, domain_width, volume_scaling, &
                       simulation_mode
    use particle_types, only: particle, calculate_terminal_velocity
    use collection_efficiency, only: collection_efficiency_E
    implicit none
    private

    ! Namelist-controlled parameters
    logical, public :: do_collisions = .false.
    logical, public :: do_coalescence = .true.
    logical, public :: write_collisions = .false.
    real(dp), public :: wmax_collision = 10.0

    ! Per-step counters (reset each call to collision_coalescence_step)
    integer(i4), public :: collisions_this_step = 0
    integer(i4), public :: coalescences_this_step = 0
    integer(i4), public :: fall_events_this_step = 0

    ! Collision output file unit (unformatted stream binary)
    integer(i4) :: collision_unit

    public :: collision_coalescence_step, initialize_collision_file

    ! Event types
    integer, parameter :: EV_PAIR = 1
    integer, parameter :: EV_FALL = 2

    type :: Event
        real(dp) :: t
        integer  :: event_type
        integer  :: i, j
    end type Event

contains

    subroutine collision_coalescence_step(lparticles, n_particles, ldt)
        ! Advance collision-coalescence for one time window ldt.
        ! Modifies particle radii and removes coalesced particles in-place.
        ! Positions are NOT written back — CODT settling is authoritative.
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(inout) :: n_particles
        real(dp), intent(in) :: ldt

        integer :: i, n_active, n_coalescences, head
        real(dp) :: grid_area, current_time
        real(dp), allocatable :: zcur(:), w_fall(:), tstamp(:)
        logical, allocatable :: alive(:)
        integer, allocatable :: sort_order(:), prev(:), next(:)
        type(Event), allocatable :: heap(:)
        integer :: heap_size
        type(Event) :: current_event
        logical :: is_periodic

        n_active = n_particles
        collisions_this_step = 0
        coalescences_this_step = 0
        fall_events_this_step = 0

        ! Periodic boundary: parcel mode wraps via modulo(z, H); no fallout.
        is_periodic = (trim(simulation_mode) == 'parcel')

        if (n_active <= 1) then
            do i = 1, n_active
                call lparticles(i)%settling(ldt)
                if (is_periodic) lparticles(i)%position = modulo(lparticles(i)%position, H)
            end do
            return
        end if

        ! Cross-sectional area of the statistical volume
        grid_area = domain_width**2 * volume_scaling
        if (grid_area <= 0.0) return

        ! Allocate working arrays
        allocate(zcur(n_active), w_fall(n_active), tstamp(n_active), alive(n_active))
        allocate(sort_order(n_active), prev(n_active), next(n_active))
        allocate(heap(2*n_active))

        ! Initialize from particle data
        do i = 1, n_active
            zcur(i) = lparticles(i)%position
            w_fall(i) = abs(calculate_terminal_velocity(lparticles(i)))
            if (w_fall(i) > wmax_collision) w_fall(i) = wmax_collision
            tstamp(i) = 0.0
        end do

        alive = .true.
        current_time = 0.0
        n_coalescences = 0
        heap_size = 0

        ! Build sorted ordering by position and linked list
        call build_order(zcur, n_active, sort_order)
        call build_links(sort_order, n_active, prev, next, head)

        ! Seed event queue with pair-passing events
        call push_all_pair_events(heap, heap_size, zcur, tstamp, w_fall, alive, next, head, n_active, ldt, current_time)

        ! Seed fallout events (chamber mode only; parcel mode wraps via
        ! modulo at writeback and no particle exits the domain).
        if (.not. is_periodic) then
            do i = 1, n_active
                call push_fall_event(heap, heap_size, i, zcur, tstamp, w_fall, alive, ldt, current_time)
            end do
        end if

        ! Main event loop
        do
            if (heap_size <= 0) exit
            call heap_pop(heap, heap_size, current_event)
            if (current_event%t > ldt) exit

            ! Validate event is still current
            if (current_event%event_type == EV_PAIR) then
                if (.not. is_pair_valid(current_event%i, current_event%j, alive, next, n_active)) cycle
            else if (current_event%event_type == EV_FALL) then
                if (current_event%i < 1 .or. current_event%i > n_active) cycle
                if (.not. alive(current_event%i)) cycle
            end if

            current_time = current_event%t

            if (current_event%event_type == EV_PAIR) then
                call handle_pair_event(current_event%i, current_event%j, current_time, lparticles, zcur, tstamp, w_fall, &
                                       alive, prev, next, head, n_active, grid_area, n_coalescences, ldt, heap, heap_size)
            else if (current_event%event_type == EV_FALL) then
                call handle_fall_event(current_event%i, current_time, zcur, tstamp, w_fall, alive, prev, next, &
                                       head, n_active, ldt, heap, heap_size)
            end if
        end do

        ! Position writeback: CC owns settling. Use update_working_pos to
        ! advance each alive particle's internal position (zcur) to t=ldt,
        ! accounting for piecewise velocity changes at merge events. Dead
        ! particles (from handle_fall_event) get z<0 so verify_particle_fallout
        ! removes them in the caller.
        do i = 1, n_active
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

        deallocate(zcur, w_fall, tstamp, alive, sort_order, prev, next, heap)

    end subroutine collision_coalescence_step


    ! =========================================================================
    ! Event handlers
    ! =========================================================================

    subroutine handle_pair_event(i, j, current_time, lparticles, zcur, tstamp, w_fall, &
                                 alive, prev, next, head, n_active, grid_area, n_coalescences, dt, heap, heap_size)
        integer, intent(in) :: i, j, n_active
        real(dp), intent(in) :: current_time, grid_area, dt
        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(inout) :: zcur(:), tstamp(:), w_fall(:)
        logical, intent(inout) :: alive(:)
        integer, intent(inout) :: prev(:), next(:), head, n_coalescences
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size

        real(dp) :: vrel, p_coll, random_draw, coal_efficiency, r_keep, r_kill
        integer :: prev_neighbor, next_neighbor, keep, kill

        call update_working_pos(i, current_time, zcur, tstamp, w_fall)
        call update_working_pos(j, current_time, zcur, tstamp, w_fall)

        vrel = w_fall(i) - w_fall(j)
        if (vrel <= 0.0) return

        ! Stage 1: Collision check (geometric kernel)
        p_coll = pi * (lparticles(i)%radius + lparticles(j)%radius)**2 / grid_area
        if (p_coll > 1.0) p_coll = 1.0

        call random_number(random_draw)

        if (random_draw < p_coll) then
            ! Collision occurred
            collisions_this_step = collisions_this_step + 1
            lparticles(i)%n_collisions = lparticles(i)%n_collisions + 1
            lparticles(j)%n_collisions = lparticles(j)%n_collisions + 1

            if (.not. do_coalescence) then
                if (write_collisions) call write_collision( &
                    lparticles(i)%particle_id, lparticles(j)%particle_id, &
                    lparticles(i)%radius, lparticles(j)%radius, 0.0_dp, zcur(i), current_time, .false.)
            else
                coal_efficiency = collection_efficiency_E(lparticles(i)%radius, lparticles(j)%radius)

                call random_number(random_draw)
                if (random_draw <= coal_efficiency) then
                    ! Coalescence: merge particles, conserving mass
                    n_coalescences = n_coalescences + 1
                    coalescences_this_step = coalescences_this_step + 1
                    keep = i
                    kill = j

                    r_keep = lparticles(keep)%radius
                    r_kill = lparticles(kill)%radius

                    ! Record pre-coalescence state
                    lparticles(keep)%radius_before_coalescence = r_keep
                    lparticles(keep)%n_coalescences = lparticles(keep)%n_coalescences + 1

                    ! Conserve water volume: r_new = (r1^3 + r2^3)^(1/3)
                    lparticles(keep)%radius = (r_keep**3 + r_kill**3)**(1.0_dp/3.0_dp)

                    ! Conserve liquid water mass
                    lparticles(keep)%water_liquid = lparticles(keep)%water_liquid &
                                                   + lparticles(kill)%water_liquid

                    ! Conserve solute mass and recompute solute radius
                    lparticles(keep)%solute_gross_mass = lparticles(keep)%solute_gross_mass &
                                                        + lparticles(kill)%solute_gross_mass
                    if (lparticles(keep)%solute_type%solute_density > 0.0) then
                        lparticles(keep)%solute_radius = &
                            (3.0_dp * lparticles(keep)%solute_gross_mass / &
                             (4.0_dp * pi * lparticles(keep)%solute_type%solute_density))**(1.0_dp/3.0_dp)
                    end if

                    ! Recalculate terminal velocity for merged droplet
                    w_fall(keep) = abs(calculate_terminal_velocity(lparticles(keep)))
                    if (w_fall(keep) > wmax_collision) w_fall(keep) = wmax_collision

                    ! Update working position
                    zcur(keep) = zcur(i)
                    tstamp(keep) = current_time

                    if (write_collisions) call write_collision( &
                        lparticles(keep)%particle_id, lparticles(kill)%particle_id, &
                        r_keep, r_kill, lparticles(keep)%radius, zcur(keep), current_time, .true.)

                    ! Remove killed particle from linked list
                    call unlink_particle(kill, alive, prev, next, head)
                    alive(kill) = .false.
                    lparticles(kill)%coalesced = .true.

                    ! Reschedule events for survivor and its neighbors
                    prev_neighbor = prev(keep)
                    next_neighbor = next(keep)
                    call push_pair_for_i(keep, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)
                    if (prev_neighbor > 0) call push_pair_for_i(prev_neighbor, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)

                    ! Reschedule fallout for survivor
                    call push_fall_event(heap, heap_size, keep, zcur, tstamp, w_fall, alive, dt, current_time)
                    return
                else
                    if (write_collisions) call write_collision( &
                        lparticles(i)%particle_id, lparticles(j)%particle_id, &
                        lparticles(i)%radius, lparticles(j)%radius, 0.0_dp, zcur(i), current_time, .false.)
                end if
            end if
        end if

        ! No coalescence: swap adjacency order (they passed each other)
        call swap_adjacent(i, j, prev, next, head)

        ! Reschedule affected pairs
        prev_neighbor = prev(j)  ! j is now before i after swap
        next_neighbor = next(i)
        if (prev_neighbor > 0) call push_pair_for_i(prev_neighbor, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)
        call push_pair_for_i(j, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)
        if (next_neighbor > 0) call push_pair_for_i(i, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)

    end subroutine handle_pair_event


    subroutine handle_fall_event(i, current_time, zcur, tstamp, w_fall, alive, prev, next, &
                                 head, n_active, dt, heap, heap_size)
        integer, intent(in) :: i, n_active
        real(dp), intent(in) :: current_time, dt
        real(dp), intent(inout) :: zcur(:), tstamp(:), w_fall(:)
        logical, intent(inout) :: alive(:)
        integer, intent(inout) :: prev(:), next(:), head
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size

        real(dp) :: zi, tfall
        integer :: prev_neighbor

        if (i < 1 .or. i > n_active) return
        if (.not. alive(i)) return

        call update_working_pos(i, current_time, zcur, tstamp, w_fall)
        zi = zcur(i)

        if (w_fall(i) <= 0.0) return

        ! Validate that fall time matches current time (stale check)
        tfall = current_time + zi / w_fall(i)
        if (abs(tfall - current_time) > 1.0e-8) return

        ! Remove particle
        prev_neighbor = prev(i)
        call unlink_particle(i, alive, prev, next, head)
        alive(i) = .false.
        fall_events_this_step = fall_events_this_step + 1

        ! Reschedule neighbor pair
        if (prev_neighbor > 0) call push_pair_for_i(prev_neighbor, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)

    end subroutine handle_fall_event


    ! =========================================================================
    ! Event generation
    ! =========================================================================

    logical function is_pair_valid(i, j, alive, next, n_active)
        integer, intent(in) :: i, j, n_active
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:)
        is_pair_valid = .false.
        if (i < 1 .or. j < 1) return
        if (i > n_active .or. j > n_active) return
        if (.not. alive(i) .or. .not. alive(j)) return
        if (next(i) /= j) return
        is_pair_valid = .true.
    end function is_pair_valid


    subroutine push_all_pair_events(heap, heap_size, zcur, tstamp, w_fall, alive, next, head, n_active, dt, current_time)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:), dt, current_time
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:), head, n_active
        integer :: i, count

        i = head
        count = 0
        do while (i /= 0 .and. count < n_active)
            call push_pair_for_i(i, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)
            i = next(i)
            count = count + 1
            if (i == 0) exit
        end do
    end subroutine push_all_pair_events


    subroutine push_pair_for_i(i, heap, heap_size, zcur, tstamp, w_fall, alive, next, dt, current_time)
        integer, intent(in) :: i
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:), dt, current_time
        logical, intent(in) :: alive(:)
        integer, intent(in) :: next(:)

        integer :: j
        real(dp) :: dz, vrel, te
        type(Event) :: new_event

        if (i < 1 .or. i > size(alive)) return
        if (.not. alive(i)) return

        j = next(i)
        if (j <= 0) return
        if (.not. alive(j)) return

        call update_working_pos(i, current_time, zcur, tstamp, w_fall)
        call update_working_pos(j, current_time, zcur, tstamp, w_fall)

        ! Forward separation (j should be above i in sorted order)
        dz = zcur(j) - zcur(i)
        if (dz < 0.0) dz = dz + H

        vrel = w_fall(i) - w_fall(j)
        if (vrel <= 0.0) return

        te = current_time + dz / vrel
        if (te <= current_time) return
        if (te > dt) return

        new_event%t = te
        new_event%event_type = EV_PAIR
        new_event%i = i
        new_event%j = j
        call heap_push(heap, heap_size, new_event)
    end subroutine push_pair_for_i


    subroutine push_fall_event(heap, heap_size, i, zcur, tstamp, w_fall, alive, dt, current_time)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size
        integer, intent(in) :: i
        real(dp), intent(in) :: zcur(:), tstamp(:), w_fall(:), dt, current_time
        logical, intent(in) :: alive(:)
        type(Event) :: new_event
        real(dp) :: zi, te

        if (i < 1 .or. i > size(alive)) return
        if (.not. alive(i)) return
        if (w_fall(i) <= 0.0) return

        ! Current position accounting for motion since last stamp
        zi = zcur(i) - w_fall(i) * (current_time - tstamp(i))
        if (zi <= 0.0) then
            te = current_time
        else
            te = current_time + zi / w_fall(i)
        end if
        if (te > dt) return

        new_event%t = te
        new_event%event_type = EV_FALL
        new_event%i = i
        new_event%j = 0
        call heap_push(heap, heap_size, new_event)
    end subroutine push_fall_event


    ! =========================================================================
    ! Position updates
    ! =========================================================================

    subroutine update_working_pos(i, current_time, zcur, tstamp, w_fall)
        integer, intent(in) :: i
        real(dp), intent(in) :: current_time
        real(dp), intent(inout) :: zcur(:), tstamp(:)
        real(dp), intent(in) :: w_fall(:)
        real(dp) :: dtloc

        dtloc = current_time - tstamp(i)
        if (dtloc /= 0.0) then
            zcur(i) = zcur(i) - w_fall(i) * dtloc
            tstamp(i) = current_time
        end if
    end subroutine update_working_pos


    ! =========================================================================
    ! Linked list operations
    ! =========================================================================

    subroutine build_links(sort_order, n_active, prev, next, head)
        integer, intent(in) :: sort_order(:), n_active
        integer, intent(out) :: prev(:), next(:), head
        integer :: k

        prev = 0
        next = 0

        do k = 1, n_active
            if (k > 1) prev(sort_order(k)) = sort_order(k-1)
            if (k < n_active) next(sort_order(k)) = sort_order(k+1)
        end do

        head = sort_order(1)
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

    subroutine build_order(zcur, n_active, sort_order)
        real(dp), intent(in) :: zcur(:)
        integer, intent(in) :: n_active
        integer, intent(out) :: sort_order(:)
        integer :: i

        do i = 1, n_active
            sort_order(i) = i
        end do
        call sort_by_position(sort_order, zcur, 1, n_active)
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

    subroutine heap_push(heap, heap_size, new_event)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size
        type(Event), intent(in) :: new_event
        integer :: k, parent_idx
        type(Event) :: tmp

        heap_size = heap_size + 1
        if (heap_size > size(heap)) then
            heap_size = heap_size - 1
            return
        end if
        heap(heap_size) = new_event

        ! Sift up
        k = heap_size
        do while (k > 1)
            parent_idx = k / 2
            if (heap(parent_idx)%t <= heap(k)%t) exit
            tmp = heap(parent_idx)
            heap(parent_idx) = heap(k)
            heap(k) = tmp
            k = parent_idx
        end do
    end subroutine heap_push


    subroutine heap_pop(heap, heap_size, popped_event)
        type(Event), intent(inout) :: heap(:)
        integer, intent(inout) :: heap_size
        type(Event), intent(out) :: popped_event
        integer :: k, left, right, smallest
        type(Event) :: tmp

        if (heap_size <= 0) then
            popped_event%t = huge(1.0_dp)
            popped_event%event_type = 0
            popped_event%i = 0
            popped_event%j = 0
            return
        end if

        popped_event = heap(1)
        heap(1) = heap(heap_size)
        heap_size = heap_size - 1

        ! Sift down
        k = 1
        do
            left = 2 * k
            right = 2 * k + 1
            smallest = k
            if (left <= heap_size .and. heap(left)%t < heap(smallest)%t) smallest = left
            if (right <= heap_size .and. heap(right)%t < heap(smallest)%t) smallest = right
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


    subroutine write_collision(id_i, id_j, r_i_m, r_j_m, r_after_m, position_m, time_s, coalesced)
        integer(i4), intent(in) :: id_i, id_j
        real(dp), intent(in) :: r_i_m, r_j_m, r_after_m, position_m, time_s
        logical, intent(in) :: coalesced
        integer(i1) :: flag

        flag = merge(1_i1, 0_i1, coalesced)
        write(collision_unit) id_i, id_j, r_i_m, r_j_m, r_after_m, position_m, time_s, flag

    end subroutine write_collision

end module collision_coalescence
