! Lagrangian droplet management: the bridge between the Eulerian scalar fields
! and the individual particles. Owns the particle array and handles injection
! (chamber) / pre-loading (parcel), gravitational settling and fallout,
! eddy displacement, droplet growth (dispatching to DGM per particle), aerosol
! initialization from NetCDF, entrainment/detrainment, and the droplet-size-
! distribution diagnostics written to output. The growth ODE itself lives in DGM;
! the collision-coalescence events live in collision_coalescence.
module droplets
    use netcdf
    use globals
    use particle_types
    use collision_coalescence, only: do_collisions, do_coalescence, &
                                     collision_coalescence_step, wmax_collision, &
                                     write_collisions, collisions_this_step, coalescences_this_step
    use collection_efficiency, only: coalescence_kernel, set_kernel_selector
    use DGM, only: integrate_ODE, set_aerosol_properties
    use special_effects, only: do_random_fallout, random_fallout_rate
    use microphysics
    implicit none

    private
    public :: particles, particle, current_n_particles, total_n_particles, total_n_fellout, n_injected
    public :: initialize_microphysics, update_droplets
    public :: move_particles_by_cellmap
    public :: calculate_droplet_statistics, bin_droplet_radii
    public :: particle_bin_edges, size_distribution, n_DSD_bins, n_aer_category
    public :: dsd_varid, aerDSD_varids
    public :: write_trajectories, trajectory_start, trajectory_end, trajectory_timer
    public :: initial_wet_radius, init_drop_each_gridpoint, expected_Ndrops_per_gridpoint
    public :: aerosol_concentration
    public :: do_collisions, do_coalescence, wmax_collision, write_collisions, coalescence_kernel
    public :: aerosol_file
    public :: detrain_particles, entrain_particles
    public :: do_seeding, seed_hydration, seed_growth_time, n_seeded
    public :: read_aerosol_netcdf, sample_radius, aerosols, aerosol_bin_type, &
              aerosol_radii, aerosol_partition, aerosol_bin_freq, aerosol_size_edges, &
              injection_times, injection_rates
    public :: n_seed_bins, seed_radii, seed_partition, seed_bin_type, seed_bin_freq, &
              seed_event_coord, seed_event_conc

    ! Counters to track particles, used for statistics and array indexing
    integer(i4) :: current_n_particles = 0
    integer(i4) :: total_n_particles = 0
    integer(i4) :: total_n_fellout = 0

    ! Collision-coalescence accumulators (accumulate between writes, reset on stats output)
    integer(i4) :: collisions_since_write = 0
    integer(i4) :: coalescences_since_write = 0
    real(dp), allocatable :: particle_bin_edges(:), particle_bins(:)
    integer(i4), allocatable :: size_distribution(:,:) !(No. DSDs, rbins)
    integer(i4) :: n_DSD_bins, n_aer_category

    ! Cached NetCDF variable IDs for DSD variables
    integer :: dsd_varid
    integer, allocatable :: aerDSD_varids(:)

    ! Arrays to hold all particles and aerosol types
    type(particle), allocatable :: particles(:) ! allocated in initialize_microphysics()
    type(aerosol), allocatable :: aerosols(:) ! allocated in read_aerosol_netcdf()
    integer(i4) :: particle_array_expansion = 1000

    ! Aerosol injection variables
    real(dp), allocatable :: injection_times(:), injection_rates(:) ! allocated in read_aerosol_netcdf()
    real(dp), allocatable :: aerosol_size_edges(:), aerosol_bin_freq(:,:) ! allocated in read_aerosol_netcdf()
    real(dp), allocatable :: aerosol_radii(:) ! allocated in read_aerosol_netcdf()
    integer(i4), allocatable :: aerosol_partition(:)
    integer(i4), allocatable :: aerosol_bin_type(:) ! per-bin aerosol type; all 1 unless the file has bin_type

    ! Seed aerosol: a size distribution of its own, sampled and injected separately
    ! from the background rather than sharing its CDF, so the file states the two
    ! populations independently and neither implies the other. Present only when the
    ! aerosol file carries the seed group (see docs/data_formats.md).
    real(dp), allocatable :: seed_size_edges(:), seed_radii(:)
    real(dp), allocatable :: seed_bin_freq(:,:) ! (seed bin, event), CDF -> 1
    integer(i4), allocatable :: seed_partition(:), seed_bin_type(:)
    integer(i4) :: n_seed_bins = 0

    ! Seeding events. Each event names a point on the schedule axis and a
    ! concentration to release when the run reaches it: one burst, fired once.
    real(dp), allocatable :: seed_event_coord(:)  ! axis value that triggers the event
    real(dp), allocatable :: seed_event_conc(:)   ! concentration released [cm-3]
    logical, allocatable :: seed_event_fired(:)
    integer(i4) :: seed_event_idx = 1  ! event being injected; selects its CDF row
    real(dp) :: last_injection_time
    real(dp) :: injection_dt
    integer(i4) :: n_injected, inj_time_idx
    logical :: update_inj_rate
    real(dp) :: initial_wet_radius
    logical :: init_drop_each_gridpoint = .false.
    real(dp) :: expected_Ndrops_per_gridpoint = 1
    real(dp) :: aerosol_concentration = 0.0
    character(256) :: aerosol_file = ''

    ! Aerosol seeding. One driver serves both modes; only the schedule axis differs.
    ! Chamber seeds on time [s], parcel on the vertical coordinate set by &PARCEL
    ! vertical_axis (m or Pa). Either way an event fires once, when the run first
    ! reaches its level. See docs/data_formats.md.
    logical :: do_seeding = .false.
    character(16) :: seed_hydration = 'equilibrium' ! 'equilibrium' | 'double_growth' | 'dry'
    real(dp) :: seed_growth_time = 5.0 ! growth time a seed is allowed at injection, s
    real(dp), parameter :: seed_RH_cap = 0.99 ! keeps the Kohler solve on the stable branch
    integer(i4) :: n_seeded = 0
    real(dp) :: prev_seed_coord         ! axis value at the previous step, to detect crossings
    logical :: seed_coord_set = .false.

    ! Particle I/O Handling
    logical :: write_trajectories = .false.
    real(dp) :: trajectory_start = 0., trajectory_end = 0.
    real(dp) :: trajectory_timer = 1.

contains

    subroutine update_droplets(ltime, ldt, coord)
        ! Interface subroutine to main.f90. Does aerosol injection, settling,
        ! droplet-environment property update, and droplet growth.
        ! Caller is responsible for syncing nondim fields afterward.
        !
        ! coord is the seeding schedule axis. Parcel runs pass the parcel's vertical
        ! coordinate, which this module cannot read for itself: parcel depends on
        ! entrainment, which depends on this module, so a use here would be circular.
        ! Chamber runs omit it and seed on time.
        real(dp), intent(in) :: ltime, ldt
        real(dp), intent(in), optional :: coord
        real(dp) :: current_coord
        integer :: i

        if (simulation_mode == 'chamber') call injection_controller(time, particles)

        if (do_seeding) then
            current_coord = ltime
            if (present(coord)) current_coord = coord
            call seeding_controller(particles, current_coord)
        end if

        if (do_collisions) then
            ! CC owns settling across the ldt window (writes back final
            ! particle positions). Fallout removal and gridcell updates
            ! match the tail of move_particles_by_gravity.
            call collision_coalescence_step(particles, current_n_particles, ldt)
            collisions_since_write = collisions_since_write + collisions_this_step
            coalescences_since_write = coalescences_since_write + coalescences_this_step
            call verify_particle_fallout(particles, current_n_particles)
            do i = 1, current_n_particles
                call particles(i)%update_gridcell()
            end do
        else
            call move_particles_by_gravity(particles, ldt)
        end if

        call update_all_particles(particles, T, WV, Tv, SS)
        call droplet_growth_model(particles, ltime, ldt)

    end subroutine update_droplets

    subroutine seeding_controller(lparticles, coord)
        ! Fires any seeding event the run has just reached. Both modes share this
        ! driver; only the axis differs. Chamber runs pass time, so an event means
        ! "seed at t = 300 s"; parcel runs pass the vertical coordinate, so it means
        ! "seed at z = 600 m". Each event releases its concentration in a single
        ! burst and then never fires again, even if the parcel returns to the level.
        !
        ! An event triggers when the step brackets its coordinate, which needs no
        ! notion of which way the run is moving along the axis: a parcel descending
        ! onto a level crosses it exactly as one rising onto it does.
        !
        ! Input:
        ! lparticles - array of particle types
        ! coord - current value on the schedule axis (time, height, or pressure)
        type(particle), allocatable, intent(inout) :: lparticles(:)
        real(dp), intent(in) :: coord
        integer(i4) :: i, j, inject_n

        ! The first call only establishes where the run starts on the axis: with no
        ! previous value there is no interval yet, and nothing can have been crossed.
        if (.not. seed_coord_set) then
            prev_seed_coord = coord
            seed_coord_set = .true.
            return
        end if

        do i = 1, size(seed_event_coord)
            if (seed_event_fired(i)) cycle
            if ((prev_seed_coord - seed_event_coord(i)) * (coord - seed_event_coord(i)) > 0.0) cycle

            ! seed_event_idx selects this event's row of the seed CDF, so each event
            ! may release a different size distribution.
            seed_event_idx = i
            inject_n = nint(seed_event_conc(i) * 1.0e6 * domain_volume)

            do j = 1, inject_n
                call inject_particle(lparticles, T, WV, Tv, SS, seed=.true.)
                n_seeded = n_seeded + 1
            end do

            seed_event_fired(i) = .true.
            write(*,'(a,i0,a,es10.3,a,i0,a)') ' Seeding event ', i, ' fired at coord ', &
                  coord, ': ', inject_n, ' particles injected'
        end do

        prev_seed_coord = coord

    end subroutine seeding_controller

    subroutine injection_controller(ltime, lparticles)
        ! Controller subroutine to manage particle injection
        ! into the simulation. Particles are injected at a
        ! specified rate and are initialized with properties
        ! from the gridcell in which they are injected.
        ! Note, will under-inject by 1 particle. I'll take that.
        !
        ! Input:
        ! ltime - current model time, dimensional
        ! lparticles - array of particle types
        !
        ! Output:
        ! particles - array of particle types with new particles injected
        real(dp), intent(in) :: ltime
        type(particle), allocatable, intent(inout) :: lparticles(:)
        real(dp) :: time_since_last_injection
        integer(i4) :: i, inject_n

        ! Determine time since last injection
        time_since_last_injection = ltime - last_injection_time

        ! Inject particles if enough time has passed
        if ( time_since_last_injection >= injection_dt ) then

            ! Determine number of particles to inject
            inject_n = int(time_since_last_injection / injection_dt)
            do i = 1, inject_n
                call inject_particle(lparticles, T, WV, Tv, SS)
                n_injected = n_injected + 1
            end do

            ! For any remaining time, count it for the next injection iteration
            ! this will roll over any 'partial' particle injections excluded from above do loop
            last_injection_time = ltime - mod(time_since_last_injection, injection_dt)

        end if

        ! Determine if injection_rate needs to be updated
        if ( update_inj_rate ) then
            if ( ltime >= injection_times(inj_time_idx+1) ) then

                ! Update injection rates
                inj_time_idx = inj_time_idx + 1
                injection_dt = calculate_injection_rate(injection_rates(inj_time_idx), domain_volume)
                write(*,*) 'Injection rate updated: ', injection_dt

                ! Stop updating if at last injection rate index
                if ( inj_time_idx == size(injection_times) ) then
                    update_inj_rate = .false.
                end if

            end if
        end if

    end subroutine injection_controller

    subroutine inject_particle(lparticles_array, Temp, Vapor, VirtTemp, Supersat, seed)
        ! Injects a new particle into the simulation at a random location. Particle inherets
        ! the properties of the gridcell which it originates. The caller chooses only
        ! which population to draw from; the material within it comes from the file.
        type(particle), allocatable, intent(inout) :: lparticles_array(:)
        real(dp), intent(in) :: Temp(:), Vapor(:), VirtTemp(:), Supersat(:)
        logical, intent(in), optional :: seed
        type(particle) :: injected_particle
        type(particle), allocatable :: temp_array(:)
        real(dp) :: random_position
        integer(i4) :: grid_idx

        ! Increment the particle count
        current_n_particles = current_n_particles + 1
        total_n_particles = total_n_particles + 1

        ! Allocate space for the new particle, if needed
        ! Uses 'particle_array_expansion' to add x more elements in array
        if (size(lparticles_array) < current_n_particles) then
            ! Grow the array in place, preserving existing particles.
            ! move_alloc transfers the allocation (single copy, no leftover temp).
            allocate(temp_array(current_n_particles + particle_array_expansion))
            temp_array(:size(lparticles_array)) = lparticles_array
            call move_alloc(temp_array, lparticles_array)
            write(*,*) "WARNING: Allocated new particle array"
        end if

        ! Generate a random position for the new particle   
        ! and determine its gridcell
        call random_number(random_position)
        grid_idx = int(random_position * N) + 1

        ! Ensure the gridcell index is within the valid range
        if (grid_idx > N) grid_idx = N
        if (grid_idx < 1) grid_idx = 1

        ! Initialize the new particle with the properties of the gridcell
        random_position = grid_idx * (H / N)

        ! Initialize the new particle using properties from the gridcell
        call particle_initialize(injected_particle, total_n_particles, random_position, grid_idx, &
                    Temp(grid_idx), Vapor(grid_idx), VirtTemp(grid_idx), Supersat(grid_idx), seed)

        call injected_particle%update_gridcell()

        ! Accumulate budget
        budget_inject_solute_mass = budget_inject_solute_mass + injected_particle%solute_gross_mass
        budget_inject_liquid_mass = budget_inject_liquid_mass + injected_particle%water_liquid
        budget_n_injected = budget_n_injected + 1

        ! Index newly injected particle object into array of current particles
        lparticles_array(current_n_particles) = injected_particle

    end subroutine inject_particle


    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINES TO MOVE PARTICLES - SETTLING AND EDDIES !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine move_particles_by_gravity(lparticles, ldt)
        ! Moves particles based on their settling velocity and current model
        ! time increment. Particle fallout is then determined and the gridcell
        ! indices are then updated for remaining particles.
        !
        ! Input:
        ! lparticles - array of particle types
        ! ldt - model time increment, dimensional
        !  
        ! Output:
        ! lparticles - array of particle types with each particle position updated

        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(in) :: ldt
        integer :: i
        
        ! Move each particle based on settling velocity
        do i = 1, current_n_particles
            call lparticles(i)%settling(ldt)
        end do

        if (simulation_mode == 'parcel') then
            do i = 1, current_n_particles
                lparticles(i)%position = modulo(lparticles(i)%position, H)
            end do
        else
            call verify_particle_fallout(lparticles, current_n_particles)
        end if

        ! For remaining particles, update gridcell index and properties
        do i = 1, current_n_particles
            call lparticles(i)%update_gridcell()
        end do

    end subroutine move_particles_by_gravity

    subroutine verify_particle_fallout(lparticle_array, n_particles)
        ! Loops through each particle in the particle array and deteremines
        ! if the particle has fallen through the bottom boundary due to settling
        ! (i.e. this is called from the end of 'move_particles_by_gravity')
        ! 
        ! Considers the special effect 'random_fallout' and reindexes the number
        ! of active particles which remain in the domain
        type(particle), intent(inout) :: lparticle_array(:)
        integer(i4), intent(inout) :: n_particles
        integer :: i, n_removed
        real(dp) :: r, new_position

        ! Update fallout flag for particles which settled out of the domain
        ! Coalesced particles (killed by CC) are already marked and skip
        ! the random_fallout pathway since their mass is in the survivor.
        do i = 1, n_particles
            if ( lparticle_array(i)%coalesced ) cycle
            if ( lparticle_array(i)%position < 0.0_dp ) then
                ! I contemplate this immortal outcome
                ! A timeless trace of future times to come
                ! I know I'm more than the span of my life
                ! I know I'm more than my own local strife
                if ( do_random_fallout ) then
                    call random_fallout(lparticle_array(i), random_fallout_rate, H)
                else
                    lparticle_array(i)%fellout = .true.
                    lparticle_array(i)%position = 0.
                end if
            end if
        end do

        ! Remove dead particles (fellout or coalesced) and compact array.
        ! Budget counters are incremented separately for each removal type.
        n_removed = 0
        do i = 1, n_particles
            if ( lparticle_array(i)%coalesced ) then
                n_removed = n_removed + 1
                budget_n_coalesced = budget_n_coalesced + 1
            else if ( lparticle_array(i)%fellout ) then
                total_n_fellout = total_n_fellout + 1
                n_removed = n_removed + 1
                budget_fallout_liquid_mass = budget_fallout_liquid_mass + lparticle_array(i)%water_liquid
                budget_fallout_solute_mass = budget_fallout_solute_mass + lparticle_array(i)%solute_gross_mass
                budget_n_fellout = budget_n_fellout + 1
            else
                if (n_removed > 0) then
                    lparticle_array(i - n_removed) = lparticle_array(i)
                end if
            end if
        end do

        ! Update number of particles in domain
        n_particles = n_particles - n_removed

    end subroutine verify_particle_fallout

    subroutine random_fallout(lparticle, fallout_rate, height)
        ! Special effects-related function. When particle falls through bottom
        ! boundary, determines whether it is removed or falls back through top
        ! of domain
        type(particle), intent(inout) :: lparticle
        real(dp), intent(in) :: fallout_rate, height
        real(dp) :: r
        
        ! Stochastic test to see if droplet fell out of domain
        call random_number(r)
        if ( r < random_fallout_rate ) then
            lparticle%fellout = .true.
            lparticle%position = 0.
        else
            lparticle%position = height + lparticle%position
        end if

    end subroutine random_fallout


    ! Remove all particles whose gridcell falls within any blob segment.
    !
    ! Uses the same shift-down compaction pattern as verify_particle_fallout,
    ! but tracks detrainment budget counters separately. Does NOT call
    ! verify_particle_fallout to avoid double-counting as fallout.
    subroutine detrain_particles(blob_starts, blob_ends, n_blobs)
        integer(i4), intent(in) :: blob_starts(:), blob_ends(:), n_blobs
        integer :: i, j, n_removed
        logical :: in_blob

        n_removed = 0
        do i = 1, current_n_particles
            in_blob = .false.
            do j = 1, n_blobs
                if (particles(i)%gridcell >= blob_starts(j) .and. &
                    particles(i)%gridcell <= blob_ends(j)) then
                    in_blob = .true.
                    exit
                end if
            end do

            if (in_blob) then
                n_removed = n_removed + 1
                budget_detrain_liquid_mass = budget_detrain_liquid_mass + particles(i)%water_liquid
                budget_detrain_solute_mass = budget_detrain_solute_mass + particles(i)%solute_gross_mass
                budget_n_detrained = budget_n_detrained + 1
            else
                if (n_removed > 0) then
                    particles(i - n_removed) = particles(i)
                end if
            end if
        end do

        current_n_particles = current_n_particles - n_removed

    end subroutine detrain_particles


    ! Inject a single particle at a random gridcell within the blob region.
    !
    ! Similar to inject_particle, but restricts placement to blob segments.
    ! A flat random index over total blob cells is mapped to the actual gridcell
    ! by walking through the blob segments sequentially.
    subroutine inject_particle_in_region(lparticles_array, Temp, Vapor, VirtTemp, Supersat, &
                                         blob_starts, blob_ends, n_blobs, n_blob_cells)
        type(particle), allocatable, intent(inout) :: lparticles_array(:)
        real(dp), intent(in) :: Temp(:), Vapor(:), VirtTemp(:), Supersat(:)
        integer(i4), intent(in) :: blob_starts(:), blob_ends(:), n_blobs, n_blob_cells
        type(particle) :: injected_particle
        type(particle), allocatable :: temp_array(:)
        real(dp) :: rval
        integer(i4) :: grid_idx, cell_count, j

        current_n_particles = current_n_particles + 1
        total_n_particles = total_n_particles + 1

        if (size(lparticles_array) < current_n_particles) then
            allocate(temp_array(size(lparticles_array)))
            temp_array = lparticles_array
            deallocate(lparticles_array)
            allocate(lparticles_array(current_n_particles + particle_array_expansion))
            lparticles_array(:size(temp_array)) = temp_array
            deallocate(temp_array)
        end if

        ! Draw a flat random index over total blob cells, then map to gridcell
        call random_number(rval)
        cell_count = int(rval * n_blob_cells) + 1
        if (cell_count > n_blob_cells) cell_count = n_blob_cells

        grid_idx = 0
        do j = 1, n_blobs
            if (cell_count <= blob_ends(j) - blob_starts(j) + 1) then
                grid_idx = blob_starts(j) + cell_count - 1
                exit
            end if
            cell_count = cell_count - (blob_ends(j) - blob_starts(j) + 1)
        end do

        call particle_initialize(injected_particle, total_n_particles, &
                    grid_idx * (H / N), grid_idx, &
                    Temp(grid_idx), Vapor(grid_idx), VirtTemp(grid_idx), Supersat(grid_idx))
        call injected_particle%update_gridcell()

        budget_entrain_solute_mass = budget_entrain_solute_mass + injected_particle%solute_gross_mass
        budget_entrain_liquid_mass = budget_entrain_liquid_mass + injected_particle%water_liquid
        budget_n_entrained = budget_n_entrained + 1

        lparticles_array(current_n_particles) = injected_particle

    end subroutine inject_particle_in_region


    ! Inject fresh aerosols into the blob region at the given number concentration.
    !
    ! The number of particles to inject is computed from the blob volume fraction:
    !   n_inject = nint(concentration [cm^-3] * 1e6 * (n_blob_cells / N) * domain_volume)
    ! Concentration is passed as an argument (not read from the module variable) so
    ! callers can supply height-dependent values in the future.
    !
    ! The particle array is pre-expanded once before the injection loop to avoid
    ! repeated reallocations. After injection, only the newly added particles are
    ! equilibrated to Koehler equilibrium via equilibrate_particles.
    subroutine entrain_particles(blob_starts, blob_ends, n_blobs, concentration)
        integer(i4), intent(in) :: blob_starts(:), blob_ends(:), n_blobs
        real(dp), intent(in) :: concentration  ! environmental aerosol concentration [cm^-3]
        integer(i4) :: n_blob_cells, n_inject, old_count, i, j
        type(particle), allocatable :: temp_array(:)
        integer(i4) :: new_capacity

        n_blob_cells = 0
        do j = 1, n_blobs
            n_blob_cells = n_blob_cells + (blob_ends(j) - blob_starts(j) + 1)
        end do

        n_inject = nint(concentration * 1.0e6 * (real(n_blob_cells, dp) / N) * domain_volume)
        if (n_inject <= 0) return

        ! Pre-expand particle array to avoid per-particle reallocation
        new_capacity = current_n_particles + n_inject
        if (size(particles) < new_capacity) then
            allocate(temp_array(size(particles)))
            temp_array = particles
            deallocate(particles)
            allocate(particles(new_capacity + particle_array_expansion))
            particles(:size(temp_array)) = temp_array
            deallocate(temp_array)
        end if

        old_count = current_n_particles

        do i = 1, n_inject
            call inject_particle_in_region(particles, T, WV, Tv, SS, &
                                           blob_starts, blob_ends, n_blobs, n_blob_cells)
        end do

        ! Equilibrate only the newly entrained particles to Koehler equilibrium
        call equilibrate_particles(particles(old_count+1:), n_inject)

    end subroutine entrain_particles


    subroutine move_particles_by_cellmap(lparticles, dest_cell)
        ! Displace every particle once, using the net cell rearrangement produced
        ! by a whole sequence of composed triplet maps. dest_cell(c) is the grid
        ! cell that the fluid originally in cell c has ended up in after the
        ! sequence, built by begin_eddy_sequence / accumulate_eddy /
        ! finalize_eddy_sequence in globals. Both ODT (one eddy per turbulence
        ! step) and LEM (several) reach this routine through that same API, so
        ! droplets follow the fluid through the composition
        ! c -> map1(c) -> map2(map1(c)) -> ... and are touched exactly once per
        ! turbulence step.
        !
        ! The displacement is a whole number of cells, so each particle keeps its
        ! sub-cell offset. Because the offset is bounded (0 <= offset < dz_length)
        ! and dest_cell is always in [1, N], the new position lands inside
        ! [0, H) by construction, including for eddies that wrap the periodic
        ! boundary -- no modulo is needed. gridcell is set directly rather than
        ! re-derived, so position and gridcell stay consistent at all times.
        !
        ! (Argument is named dest_cell, not destination_cell, to avoid masking
        ! the module variable of that name in globals.)
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(in) :: dest_cell(:)
        integer(i4) :: i, start_cell, end_cell

        do i = 1, current_n_particles
            start_cell = lparticles(i)%gridcell
            end_cell = dest_cell(start_cell)
            if (end_cell /= start_cell) then
                lparticles(i)%position = z(end_cell) &
                    + (lparticles(i)%position - z(start_cell))
                lparticles(i)%gridcell = end_cell
            end if
        end do

    end subroutine move_particles_by_cellmap

    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SUBROUTINES TO INITIALIZE PARTICLES AND AEROSOLS !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine particle_initialize(this, ln_particles, pos, grid_idx, temp, vapor, virt_temp, supersat, seed)
        ! Draws a dry aerosol from a size distribution at the current schedule row
        ! and hydrates it. The seed flag selects which population to draw from; the
        ! material then follows from the sampled bin, so a file with several
        ! aerosol types populates the domain with the mixture its CDF describes.
        type(particle), intent(out) :: this
        integer(i4), intent(in) :: ln_particles
        real(dp), intent(in) :: pos
        integer(i4), intent(in) :: grid_idx
        real(dp), intent(in) :: temp
        real(dp), intent(in) :: vapor
        real(dp), intent(in) :: virt_temp
        real(dp), intent(in) :: supersat
        logical, intent(in), optional :: seed  ! draw seed material rather than background
        logical :: is_seed
        integer(i4) :: type_idx  ! composition row the sampled bin is made of

        ! Background is the default: callers that never seed (chamber injection,
        ! entrainment, the parcel preload) simply omit the flag.
        is_seed = .false.
        if (present(seed)) is_seed = seed

        ! Each of us, a cell of awareness
        ! Imperfect and incomplete
        ! Genetic blends, with uncertain ends
        ! On a fortune hunt that's far too fleet
        this%particle_id = ln_particles
        this%position = pos
        this%gridcell = grid_idx
        this%temperature = temp
        this%water_vapor = vapor
        this%virt_temp = virt_temp
        this%supersaturation = supersat

        ! Currently in domain
        this%fellout = .false.

        ! Determine solute properties (sampled from aerosol input) and initial radius.
        ! Each population is drawn from its own distribution at its own schedule row.
        if (is_seed) then
            call sample_radius(seed_bin_freq(:,seed_event_idx), seed_radii, seed_partition, &
                               seed_bin_type, this%solute_radius, this%aerosol_category, type_idx)
        else
            call sample_radius(aerosol_bin_freq(:,inj_time_idx), aerosol_radii, aerosol_partition, &
                               aerosol_bin_type, this%solute_radius, this%aerosol_category, type_idx)
        end if
        this%solute_type = aerosols(type_idx)
        this%solute_gross_mass = ( pi_43 * this%solute_type%solute_density ) * this%solute_radius**3
        this%radius = initial_wet_radius * this%solute_radius

        ! Determine liquid water content of the particle
        call this%calculate_water_content()

        ! Determine the critical radius and supersaturation for this particle
        call this%critical_kohler()

        ! Seed material is hydrated on its own terms; background keeps
        ! initial_wet_radius. Testing the material rather than the local is_seed says
        ! what the rule actually is. The two agree in any case: build_aerosol_table
        ! rejects a composition row that is both background and seed.
        if (this%solute_type%is_seed) call hydrate_seed_particle(this)

    end subroutine particle_initialize

    subroutine hydrate_seed_particle(this)
        ! Sets the wet radius of a freshly injected seed particle according to
        ! seed_hydration. Background particles never reach here: they keep the
        ! initial_wet_radius multiple of their dry radius.
        type(particle), intent(inout) :: this
        real(dp) :: RH

        select case (trim(seed_hydration))
        case ('dry')
            ! Injected as a bare solute core, leaving the growth model to wet it
            this%radius = this%solute_radius

        case ('double_growth')
            ! Fixed wet radius at twice the dry radius, independent of humidity
            this%radius = 2.0 * this%solute_radius

        case ('equilibrium')
            ! Haze equilibrium with the humidity of the cell it lands in, bounded by
            ! what the droplet could actually grow to in seed_growth_time. The RH cap
            ! keeps the solve on the stable branch even in a supersaturated cell,
            ! where a large seed would otherwise have no equilibrium radius at all.
            !
            ! The growth bound matters for GCCN: their equilibrium radius is tens of
            ! microns, and injecting them there would condense water they need
            ! minutes to collect. Anything small enough to equilibrate quickly
            ! reaches its equilibrium radius and the bound does not bind.
            !
            ! TODO: the growth bound is a stopgap. seed_growth_time is a free
            ! parameter with no physical derivation, and a bounded seed lands at a
            ! radius that is neither an equilibrium nor a state the flow produced,
            ! so a GCCN's initial size still depends on a tuning knob. Worth
            ! revisiting: inject GCCN dry and let the DGM supply the transient, or
            ! derive the bound from the growth timescale at the injection RH.
            RH = 1.0 + this%supersaturation / 100.0
            this%radius = min(this%kohler_equilibrium_radius(min(RH, seed_RH_cap)), &
                              grown_radius(this, seed_growth_time))
        end select

        call this%calculate_water_content()

    end subroutine hydrate_seed_particle

    function grown_radius(this, t_grow) result(radius)
        ! Radius a particle would reach by growing at the humidity of its own
        ! gridcell for t_grow seconds, starting from the initial wet radius. This is
        ! a hypothetical: it asks how far the droplet could get, and the caller wants
        ! only the answer, so the growth runs on a throwaway copy.
        !
        ! single_droplet_growth accumulates into the global condensation and
        ! temperature budgets, and this growth never happens -- the droplet is about
        ! to be born at the resulting radius, and the water it holds is accounted for
        ! by the injection budget. So the budgets are restored afterwards, or every
        ! seed release would book condensation the run never performed.
        type(particle), intent(in) :: this
        real(dp), intent(in) :: t_grow
        real(dp) :: radius
        type(particle) :: trial
        real(dp) :: t_elapsed, condensation_before, delta_T_before
        real(dp), parameter :: growth_dt = 0.01

        trial = this
        condensation_before = budget_condensation
        delta_T_before = budget_dgm_delta_T

        ! Growth depends on the interval, not on where the run sits in time, so the
        ! integration starts from zero like the other off-clock growth in this module.
        t_elapsed = 0.0
        do while (t_elapsed < t_grow)
            call single_droplet_growth(trial, 0.0_dp, growth_dt)
            call update_particle(trial)
            t_elapsed = t_elapsed + growth_dt
        end do

        budget_condensation = condensation_before
        budget_dgm_delta_T = delta_T_before

        radius = trial%radius

    end function grown_radius

    subroutine sample_radius(bin_freq, radii, partition, bin_types, radius, aer_partition, type_idx)
        ! Selects an aerosol size from a CDF in the aerosol input file. The caller
        ! passes one group's arrays (background or seed), which are indexed
        ! alike, so the two populations are drawn by the same code from their own
        ! distributions. The sampled bin carries its output category and its
        ! material type along with the radius.
        real(dp), intent(in) :: bin_freq(:), radii(:)
        integer(i4), intent(in) :: partition(:), bin_types(:)
        real(dp), intent(out) :: radius
        integer(i4), intent(out) :: aer_partition
        integer(i4), intent(out) :: type_idx
        real(dp) :: random_num
        integer :: idx

        ! Generate a random number between 0 and 1
        call random_number(random_num)

        idx = 1
        do while ( random_num .gt. bin_freq(idx) )
            idx = idx + 1
        end do

        radius = radii(idx) * m_per_nm
        aer_partition = partition(idx)
        type_idx = bin_types(idx)

    end subroutine sample_radius

    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINES TO UPDATE PARTICLE PROPERTIES !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINES TO INTERFACE WITH DROPLET GROWTH MODEL !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine droplet_growth_model(lparticles, ltime, ldt)
        ! Contoller subroutine to initiate the droplet growth model for each particle
        ! in an array of particles. Droplet properties are updated after call to the DGM,
        ! then scalar fields are updated for the gridcell in which the particle resides,
        ! before moving on to next dparticle.
        !
        ! Input:
        ! lparticles - array of particle types
        ! ltime - current model time, dimensional
        ! ldt - model time increment, dimensional
        !
        ! Output:
        ! lparticles - array of particle types with updated properties
        ! scalar fields - updated sequentially for each particle based on DGM
        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(in) :: ltime, ldt
        integer :: i

        do i = 1, current_n_particles
            call update_particle(lparticles(i))
            call single_droplet_growth(lparticles(i), ltime, ldt)
            !call lparticles(i)%verify_radius()
            call update_scalar_fields_DGM(lparticles(i), T, WV, Tv, SS)
        end do

    end subroutine droplet_growth_model


    subroutine single_droplet_growth(droplet, ltime, ldt)
        ! Interface to droplet growth model. Packages an individual particle's
        ! properties and calls the ODE integrator in DGM.f90.
        type(particle), intent(inout) :: droplet
        real(dp), intent(in) :: ltime, ldt
        real(dp) :: y_arr(3), y_before(3), grid_mass, inverse_grid_mass, grid_rho
        real(dp) :: wl_before, T_before

        ! Determine mass of air in gridcell
        grid_rho = pres / (Rd * droplet%virt_temp)
        grid_mass = gridcell_volume * grid_rho
        inverse_grid_mass = 1.0 / grid_mass

        ! Save pre-growth state for budget tracking
        wl_before = droplet%water_liquid
        T_before = droplet%temperature

        ! set odeint parameters for different aerosol mass of each droplet
        call set_aerosol_properties(1, droplet%solute_gross_mass, droplet%solute_radius, &
                                    inverse_grid_mass, droplet%supersaturation/100)

        ! Package droplet properties into an array for ode solver
        y_arr(1) = droplet%radius
        y_arr(2) = droplet%water_vapor
        y_arr(3) = droplet%temperature

        y_before = y_arr

        call integrate_ODE(y_arr, ltime, ltime + ldt, droplet%h_last)

        ! Unpack droplet properties from array
        droplet%radius = y_arr(1)
        droplet%water_vapor = y_arr(2)
        droplet%temperature = y_arr(3)
        droplet%water_liquid = wl_before - (y_arr(2) - y_before(2)) * grid_mass
        droplet%supersaturation = calc_supersat(droplet%temperature, droplet%water_vapor, pres)
        droplet%virt_temp = virtual_temp(droplet%temperature, droplet%water_vapor)

        ! Accumulate condensation/evaporation budget
        budget_condensation = budget_condensation + (droplet%water_liquid - wl_before)
        budget_dgm_delta_T = budget_dgm_delta_T + (droplet%temperature - T_before)

    end subroutine single_droplet_growth



    subroutine update_scalar_fields_DGM(droplet, lT, lWV, lTv, lSS)
        ! Once droplet growth model is complete/solved, scalar fields need to be updated
        ! to values determine from DGM
        type(particle), intent(in) :: droplet
        real(dp), intent(inout) :: lT(:), lWV(:), lTv(:), lSS(:)

        integer(i4) :: idx

        ! Locate gridcell to update based on droplet position
        idx = droplet%gridcell

        ! Update scalar fields for the gridcell
        lT(idx) = droplet%temperature
        lWV(idx) = droplet%water_vapor
        lSS(idx) = droplet%supersaturation
        lTv(idx) = droplet%virt_temp

    end subroutine update_scalar_fields_DGM

    subroutine update_all_particles(lparticles, lT, lWV, lTv, lSS)
        ! Controller subroutine to update particle properties based on the local properties
        ! of the gridcell. Called after any positional updates to the particles.
        !
        ! Input:
        ! lparticles - array of particle types
        ! lT - temperature at each gridcell, C
        ! lWV - water vapor content at each gridcell, kg/kg
        ! lTv - virtual temperature at each gridcell, C
        ! lSS - supersaturation at each gridcell, %
        !
        ! Output:
        ! lparticles - array of particle types with updated properties
        type(particle), intent(inout) :: lparticles(:)
        real(dp), intent(in) :: lT(:), lWV(:), lTv(:), lSS(:)
        integer :: i, idx

        do i = 1, current_n_particles
            idx = lparticles(i)%gridcell
            call lparticles(i)%update_scalars(lT(idx), lWV(idx), lTv(idx), lSS(idx))
            call lparticles(i)%verify_activation
        end do

    end subroutine update_all_particles

    subroutine update_particle(lparticle)
        ! Adjust particle properties to its local grid box
        class(particle), intent(inout) :: lparticle
        integer :: idx

        idx = lparticle%gridcell
        call lparticle%update_scalars(T(idx), WV(idx), Tv(idx), SS(idx))
        call lparticle%verify_activation

    end subroutine update_particle

    subroutine print_particle_properties(droplet)

        type(particle), intent(in) :: droplet

        write(*,*) "Particle Properties"
        write(*,*) "Particle ID: ", droplet%particle_id
        write(*,*) "Gridcell: ", droplet%gridcell
        write(*,*) "Position: ", droplet%position
        write(*,*) "Temperature: ", droplet%temperature
        write(*,*) "Water Vapor: ", droplet%water_vapor
        write(*,*) "Virtual Temp: ", droplet%virt_temp
        write(*,*) "Supersaturation: ", droplet%supersaturation
        write(*,*) "Radius: ", droplet%radius
    
    end subroutine print_particle_properties


    !!! INITIALIZATION !!!

    subroutine initialize_microphysics()
        ! Initialization of the MICROPHYSICS namelist and related parameters
        integer     :: ierr, nml_unit, i
        character(256) :: nml_line, io_emsg

        namelist /MICROPHYSICS/ init_drop_each_gridpoint, expected_Ndrops_per_gridpoint, aerosol_file, &
        write_trajectories, trajectory_start, trajectory_end, trajectory_timer, initial_wet_radius, &
        do_collisions, do_coalescence, wmax_collision, write_collisions, coalescence_kernel, &
        aerosol_concentration, do_seeding, seed_hydration, seed_growth_time

        ! Read in microphysical namelist parameters
        write(*,*) 'Reading MICROPHYSICS namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(error_unit,*) io_emsg; call exit(1)
        end if
        read(nml=MICROPHYSICS, unit=nml_unit, iostat=ierr)
        if (ierr .ne. 0) call namelist_read_error(nml_unit, 'MICROPHYSICS')
        close(nml_unit)

        call set_kernel_selector()

        ! Collision-coalescence consistency checks
        if (do_coalescence .and. .not. do_collisions) then
            write(error_unit,*) 'ERROR: do_coalescence requires do_collisions.'
            call exit(1)
        end if
        if (write_collisions .and. .not. do_collisions) then
            write(error_unit,*) 'ERROR: write_collisions requires do_collisions.'
            call exit(1)
        end if
        if (do_collisions .and. .not. do_coalescence) then
            write(*,*) 'NOTE: Collisions enabled without coalescence (collisions-only mode).'
        end if

        ! Namelist is copied to output directory in initialize_params
        

        ! Verify that the initial wet radius will be greater than dry radius
        if ( initial_wet_radius <= 1.) then
            write(*,*) 'WARNING: initial_wet_radius must be > 1.0 (got ', initial_wet_radius, ')'
            write(*,*) '         Auto-correcting to 1.1x dry radius.'
            initial_wet_radius = 1.1
        end if

        ! Resolve input paths relative to namelist directory
        aerosol_file = resolve_path(namelist_dir, trim(aerosol_file))

        ! Set up aerosol type, injection forcings, and DSD bin edges
        call read_aerosol_netcdf(trim(aerosol_file))

        call validate_seeding_params()

        ! Set up DSD arrays
        ! Seed bins carry their own categories, and a per-category DSD is allocated
        ! for each, so the count has to span both populations.
        n_aer_category = maxval(aerosol_partition)
        if (n_seed_bins > 0) n_aer_category = max(n_aer_category, maxval(seed_partition))
        if (n_aer_category > 1) then
            allocate(size_distribution(1 + n_aer_category, n_DSD_bins))
        else
            allocate(size_distribution(1, n_DSD_bins))
        end if
        size_distribution = 0

        ! Calculate mid-point radii of DSD
        allocate(particle_bins(n_DSD_bins))
        do i = 1, n_DSD_bins
            particle_bins(i) = 0.5*(particle_bin_edges(i) + particle_bin_edges(i+1))
        end do

        ! setup variable in netCDF
        call netcdf_add_DSD(ncid, particle_bins)

        ! Make multiple DSD variables for each aerosol category if applicable
        if ( n_aer_category > 1 ) then
            call netcdf_add_aerDSD(ncid, n_aer_category)
        end if

        if (simulation_mode == 'parcel') then
            call initialize_parcel_aerosol()
        else
            call initialize_chamber_aerosol()
        end if

    end subroutine initialize_microphysics


    subroutine validate_seeding_params()
        ! Checks the seeding namelist settings against each other and against the
        ! aerosol file. Collects every problem before exiting so one run surfaces
        ! them all.
        logical :: has_error
        integer(i4) :: i

        has_error = .false.

        select case (trim(seed_hydration))
        case ('equilibrium', 'double_growth', 'dry')
            ! valid
        case default
            write(error_unit,*) 'ERROR: seed_hydration must be equilibrium, double_growth, or dry. Got: ', &
                                trim(seed_hydration)
            has_error = .true.
        end select

        if (seed_growth_time <= 0.0) then
            write(error_unit,*) 'ERROR: seed_growth_time must be > 0. Got: ', seed_growth_time
            has_error = .true.
        end if

        ! do_seeding is the absolute controller. Enabling it requires a seed group
        ! to act on; without one there is nothing to seed, so that is fatal. With
        ! seeding off the group is never read (n_seed_bins == 0 here regardless of
        ! what the file contains), so a dormant seed group in the file is simply
        ! ignored rather than rejected.
        if (do_seeding .and. n_seed_bins == 0) then
            write(error_unit,*) 'ERROR: do_seeding is set but the aerosol file has no seed group.', &
                                ' See docs/data_formats.md.'
            has_error = .true.
        end if

        if (do_seeding .and. n_seed_bins > 0) then
            ! Events fire independently on crossing, so they need no ordering. Two
            ! events at the same coordinate would both fire on the same step, which
            ! is a confusing way to write one larger event.
            do i = 2, size(seed_event_coord)
                if (any(seed_event_coord(:i-1) == seed_event_coord(i))) then
                    write(error_unit,*) 'ERROR: duplicate seed_coord at event ', i, &
                                        '. Combine them into a single event instead.'
                    has_error = .true.
                    exit
                end if
            end do

            if (any(seed_event_conc < 0.0)) then
                write(error_unit,*) 'ERROR: seed_concentration must be >= 0.'
                has_error = .true.
            end if
        end if

        if (has_error) call exit(1)

    end subroutine validate_seeding_params


    subroutine initialize_chamber_aerosol()
        integer(i4) :: i

        call initialize_injection(injection_rates)
        allocate(particles(int(expected_Ndrops_per_gridpoint*N)))
        if (init_drop_each_gridpoint) then
            do i = 1, N
                call inject_particle(particles, T, WV, Tv, SS)
            end do
        end if

    end subroutine initialize_chamber_aerosol


    subroutine initialize_parcel_aerosol()
        integer(i4) :: n_total, i

        inj_time_idx = 1
        n_total = nint(aerosol_concentration * 1.0e6 * domain_volume)

        write(*,'(a,i0,a,f0.1,a)') ' Expected particles: ', n_total, &
              ' (', aerosol_concentration, ' cm-3)'

        allocate(particles(max(n_total, 1)))

        do i = 1, n_total
            call inject_particle(particles, T, WV, Tv, SS)
        end do

        call equilibrate_particles(particles, n_total)

    end subroutine initialize_parcel_aerosol


    subroutine equilibrate_particles(lparticles, n_particles)
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(in) :: n_particles
        real(dp) :: eq_dt, r_before, max_dr
        integer(i4) :: i, iter
        integer(i4), parameter :: max_iter = 100
        real(dp), parameter :: eq_tol = 1.0e-12

        eq_dt = 0.01

        do iter = 1, max_iter
            max_dr = 0.0
            do i = 1, n_particles
                r_before = lparticles(i)%radius
                call single_droplet_growth(lparticles(i), 0.0_dp, eq_dt)
                call update_particle(lparticles(i))
                max_dr = max(max_dr, abs(lparticles(i)%radius - r_before))
            end do

            if (max_dr < eq_tol) exit
        end do

        write(*,'(a,i0,a,es9.2)') ' Equilibrated particles in ', iter, &
              ' iterations, max dr = ', max_dr

    end subroutine equilibrate_particles


    subroutine netcdf_add_DSD(lncid, r_bins)
        ! Adds a DSD variable to the netcdf file
        integer, intent(in) :: lncid
        real(dp), intent(in) :: r_bins(:)

        integer :: t_dimid, r_dimid, r_varid, nbins, dimids(2)
        integer :: re_dimid, re_varid

        nbins = size(r_bins)

        call nc_verify( nf90_inq_dimid(lncid, "time", t_dimid))

        ! Open netcdf in definition mode
        call nc_verify( nf90_redef(lncid), "nf90_redef: DSD" )

        ! Create radius dimension and variable (bin centers)
        call nc_verify( nf90_def_dim(lncid, "radius", nbins, r_dimid), "nf90_def_dim: radius")
        call nc_verify( nf90_def_var(lncid, "radius", NF90_FLOAT, r_dimid, r_varid), "nf90_def_var: radius")
        call nc_verify( nf90_put_att(lncid, r_varid, "long_name", "Droplet Bin Centers"), "nf90_put_att: radius, name")
        call nc_verify( nf90_put_att(lncid, r_varid, "units", "microns"), "nf90_put_att: radius, units")

        ! Create radius_edges variable (bin edges)
        call nc_verify( nf90_def_dim(lncid, "radius_edges", nbins + 1, re_dimid), "nf90_def_dim: radius_edges")
        call nc_verify( nf90_def_var(lncid, "radius_edges", NF90_FLOAT, re_dimid, re_varid), "nf90_def_var: radius_edges")
        call nc_verify( nf90_put_att(lncid, re_varid, "units", "microns"), "nf90_put_att: radius_edges, units")
        call nc_verify( nf90_put_att(lncid, re_varid, "long_name", "Droplet Bin Edges"), "nf90_put_att: radius_edges, name")

        ! Create Droplet Size Distribution variable
        dimids = (/ r_dimid, t_dimid /)
        call nc_verify( nf90_def_var(lncid, "DSD", NF90_INT, dimids, dsd_varid, &
                        deflate_level=1, shuffle=.true.), "nf90_def_var: DSD")
        call nc_verify( nf90_put_att(lncid, dsd_varid, "long_name", "Droplet Size Distribution"), "nf90_put_att: DSD, name")
        call nc_verify( nf90_put_att(lncid, dsd_varid, "units", "#"), "nf90_put_att: DSD, units")

        call nc_verify( nf90_enddef(lncid), "nf90_enddef: DSD")

        ! Populate radius variables
        call nc_verify( nf90_put_var(lncid, r_varid, r_bins), "nf90_put_var: radius")
        call nc_verify( nf90_put_var(lncid, re_varid, particle_bin_edges), "nf90_put_var: radius_edges")

    end subroutine netcdf_add_DSD

    subroutine netcdf_add_aerDSD(lncid, n_DSDs)

        integer, intent(in) :: lncid, n_DSDs
        integer :: i
        character(100) :: name, strint

        integer :: t_dimid, r_dimid, nbins, dimids(2)

        ! Get dimension IDs
        call nc_verify( nf90_inq_dimid(lncid, "time", t_dimid))
        call nc_verify( nf90_inq_dimid(lncid, "radius", r_dimid))

        ! Open netcdf in definition mode, and create a DSD for each aerosol partition
        dimids = (/ r_dimid, t_dimid /)
        call nc_verify( nf90_redef(lncid), "nf90_redef: DSD_aerr" )
        allocate(aerDSD_varids(n_DSDs))
        do i = 1, n_DSDs
            write(strint,*) i
            name = "DSD_" // adjustl(strint)
            call nc_verify( nf90_def_var(lncid, trim(name), NF90_INT, dimids, aerDSD_varids(i), &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: DSD_aer" )
            name = "Droplet Size Distribution - " // adjustl(strint)
            call nc_verify( nf90_put_att(lncid, aerDSD_varids(i), "long_name", trim(name)), "nf90_put_att: DSD_aer, name")
            call nc_verify( nf90_put_att(lncid, aerDSD_varids(i), "units", "#"), "nf90_put_att: DSD_aer, units")

        end do
        call nc_verify( nf90_enddef(lncid), "nf90_enddef: DSD_aer")

    end subroutine netcdf_add_aerDSD

    subroutine read_aerosol_netcdf(filepath)
        ! Reads the background aerosol from a NetCDF file (CODT_aerosol_input_v1
        ! schema): the composition table, the per-bin size distribution and labels,
        ! and the injection schedule. Delegates the optional seed group to
        ! read_seed_group and the composition table to build_aerosol_table.
        !
        ! The schema has grown backward-compatibly under the v1 string: bin_type and
        ! the whole seed group are optional, and a file written before they existed
        ! reads as a single background material, which is what it is.
        character(*), intent(in) :: filepath
        integer :: i, aer_ncid, varid, dimid, status
        integer :: n_types, n_bins, n_edges, n_times, n_dsd_edges
        character(64) :: conventions
        character(20) :: aer_name
        integer, allocatable :: aer_n_ions(:)
        real(dp), allocatable :: aer_molar_mass(:), aer_density(:)

        write(*,*) 'Reading aerosol data from: ', trim(filepath)

        ! Open and validate schema version
        call nc_verify(nf90_open(trim(filepath), NF90_NOWRITE, aer_ncid), &
                       'opening aerosol file')
        call nc_verify(nf90_get_att(aer_ncid, NF90_GLOBAL, 'conventions', conventions), &
                       'reading conventions attribute')
        if (trim(conventions) /= 'CODT_aerosol_input_v1') then
            write(error_unit,*) 'Error: expected CODT_aerosol_input_v1, got: ', trim(conventions)
            call exit(1)
        end if

        ! Read dimensions
        call nc_verify(nf90_inq_dimid(aer_ncid, 'aerosol_type', dimid), 'finding aerosol_type dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_types), 'reading aerosol_type dim')
        call nc_verify(nf90_inq_dimid(aer_ncid, 'bin', dimid), 'finding bin dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_bins), 'reading bin dim')
        call nc_verify(nf90_inq_dimid(aer_ncid, 'edge', dimid), 'finding edge dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_edges), 'reading edge dim')
        call nc_verify(nf90_inq_dimid(aer_ncid, 'time', dimid), 'finding time dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_times), 'reading time dim')

        if (n_edges /= n_bins + 1) then
            write(error_unit,*) 'Error: edge dimension must equal bin + 1'
            call exit(1)
        end if

        ! Allocate and read size distribution arrays
        allocate(aerosol_size_edges(n_edges))
        call nc_verify(nf90_inq_varid(aer_ncid, 'edge_radii', varid), 'finding edge_radii')
        call nc_verify(nf90_get_var(aer_ncid, varid, aerosol_size_edges), 'reading edge_radii')

        allocate(injection_times(n_times))
        call nc_verify(nf90_inq_varid(aer_ncid, 'injection_time', varid), 'finding injection_time')
        call nc_verify(nf90_get_var(aer_ncid, varid, injection_times), 'reading injection_time')

        allocate(injection_rates(n_times))
        call nc_verify(nf90_inq_varid(aer_ncid, 'injection_rate', varid), 'finding injection_rate')
        call nc_verify(nf90_get_var(aer_ncid, varid, injection_rates), 'reading injection_rate')

        allocate(aerosol_partition(n_bins))
        call nc_verify(nf90_inq_varid(aer_ncid, 'category', varid), 'finding category')
        call nc_verify(nf90_get_var(aer_ncid, varid, aerosol_partition), 'reading category')

        ! bin_type maps each bin to a row of the composition table. It is optional:
        ! a single-material file omits it and every bin is type 1, exactly as before.
        allocate(aerosol_bin_type(n_bins))
        status = nf90_inq_varid(aer_ncid, 'bin_type', varid)
        if (status == NF90_NOERR) then
            call nc_verify(nf90_get_var(aer_ncid, varid, aerosol_bin_type), 'reading bin_type')
            if (minval(aerosol_bin_type) < 1 .or. maxval(aerosol_bin_type) > n_types) then
                write(error_unit,*) 'Error: bin_type values must lie in 1..', n_types
                call exit(1)
            end if
        else
            aerosol_bin_type = 1
        end if

        allocate(aerosol_bin_freq(n_bins, n_times))
        call nc_verify(nf90_inq_varid(aer_ncid, 'cumulative_frequency', varid), &
                       'finding cumulative_frequency')
        call nc_verify(nf90_get_var(aer_ncid, varid, aerosol_bin_freq), &
                       'reading cumulative_frequency')

        ! Read the composition table: one entry per aerosol type
        allocate(aer_n_ions(n_types), aer_molar_mass(n_types), aer_density(n_types))
        call nc_verify(nf90_inq_varid(aer_ncid, 'n_ions', varid), 'finding n_ions')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_n_ions), 'reading n_ions')

        call nc_verify(nf90_inq_varid(aer_ncid, 'molar_mass', varid), 'finding molar_mass')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_molar_mass), 'reading molar_mass')

        call nc_verify(nf90_inq_varid(aer_ncid, 'solute_density', varid), 'finding solute_density')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_density), 'reading solute_density')

        ! Read aerosol name from global attribute
        call nc_verify(nf90_get_att(aer_ncid, NF90_GLOBAL, 'aerosol_name', aer_name), &
                       'reading aerosol_name')

        ! Read DSD bin edges (microns)
        call nc_verify(nf90_inq_dimid(aer_ncid, 'dsd_edge', dimid), 'finding dsd_edge dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_dsd_edges), 'reading dsd_edge dim')
        n_DSD_bins = n_dsd_edges - 1
        allocate(particle_bin_edges(n_dsd_edges))
        call nc_verify(nf90_inq_varid(aer_ncid, 'dsd_bin_edges', varid), 'finding dsd_bin_edges')
        call nc_verify(nf90_get_var(aer_ncid, varid, particle_bin_edges), 'reading dsd_bin_edges')

        ! do_seeding is the sole gate on the seed group: when off, the group is
        ! not read at all (n_seed_bins stays 0), so a seed group present in the
        ! file is ignored and never reaches build_aerosol_table, n_aer_category,
        ! or any output. When on, read_seed_group requires the group to be
        ! present (its absence is caught by validate_seeding_params).
        if (do_seeding) then
            call read_seed_group(aer_ncid, n_types)
        else if (nf90_inq_dimid(aer_ncid, 'seed_bin', dimid) == NF90_NOERR) then
            ! Seed data present but seeding disabled: ignored, but warn so the
            ! user is not surprised that the file's seed population does nothing.
            write(*,*) 'Warning: seeding input detected but do_seeding is false; ' // &
                       'the seed group is ignored.'
        end if

        call nc_verify(nf90_close(aer_ncid), 'closing aerosol file')

        ! Calculate midpoint radii (nanometers)
        allocate(aerosol_radii(n_bins))
        do i = 1, n_bins
            aerosol_radii(i) = (aerosol_size_edges(i) + aerosol_size_edges(i+1)) / 2.0_dp
        end do

        call build_aerosol_table(aer_name, aer_n_ions, aer_molar_mass, aer_density)

    end subroutine read_aerosol_netcdf

    subroutine read_seed_group(aer_ncid, n_types)
        ! Reads the optional seed aerosol group: its own bins, size distribution,
        ! and event schedule, kept separate from the background so the file states
        ! each population explicitly and neither implies the other.
        !
        ! The group is present as a whole or not at all. Absent, n_seed_bins stays
        ! zero and nothing downstream seeds. seed_coord is time [s] in chamber mode
        ! and the vertical coordinate (m or Pa) in parcel mode.
        integer, intent(in) :: aer_ncid, n_types
        integer :: varid, dimid, status, i
        integer :: n_seed_edges, n_seed_events

        status = nf90_inq_dimid(aer_ncid, 'seed_bin', dimid)
        if (status /= NF90_NOERR) then
            n_seed_bins = 0
            return
        end if

        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_seed_bins), 'reading seed_bin dim')
        call nc_verify(nf90_inq_dimid(aer_ncid, 'seed_edge', dimid), 'finding seed_edge dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_seed_edges), 'reading seed_edge dim')
        call nc_verify(nf90_inq_dimid(aer_ncid, 'seed_event', dimid), 'finding seed_event dim')
        call nc_verify(nf90_inquire_dimension(aer_ncid, dimid, len=n_seed_events), 'reading seed_event dim')

        if (n_seed_edges /= n_seed_bins + 1) then
            write(error_unit,*) 'Error: seed_edge dimension must equal seed_bin + 1'
            call exit(1)
        end if

        allocate(seed_size_edges(n_seed_edges))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_edge_radii', varid), 'finding seed_edge_radii')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_size_edges), 'reading seed_edge_radii')

        allocate(seed_partition(n_seed_bins))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_category', varid), 'finding seed_category')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_partition), 'reading seed_category')

        ! Required for the seed group: this mapping is what marks a composition row
        ! as seed material, so there is no sensible default for it.
        allocate(seed_bin_type(n_seed_bins))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_bin_type', varid), 'finding seed_bin_type')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_bin_type), 'reading seed_bin_type')
        if (minval(seed_bin_type) < 1 .or. maxval(seed_bin_type) > n_types) then
            write(error_unit,*) 'Error: seed_bin_type values must lie in 1..', n_types
            call exit(1)
        end if

        allocate(seed_bin_freq(n_seed_bins, n_seed_events))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_frequency', varid), 'finding seed_frequency')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_bin_freq), 'reading seed_frequency')

        allocate(seed_event_coord(n_seed_events))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_coord', varid), 'finding seed_coord')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_event_coord), 'reading seed_coord')

        allocate(seed_event_conc(n_seed_events))
        call nc_verify(nf90_inq_varid(aer_ncid, 'seed_concentration', varid), 'finding seed_concentration')
        call nc_verify(nf90_get_var(aer_ncid, varid, seed_event_conc), 'reading seed_concentration')

        allocate(seed_event_fired(n_seed_events))
        seed_event_fired = .false.

        ! Midpoint radii (nanometers), as for the background bins
        allocate(seed_radii(n_seed_bins))
        do i = 1, n_seed_bins
            seed_radii(i) = (seed_size_edges(i) + seed_size_edges(i+1)) / 2.0_dp
        end do

    end subroutine read_seed_group

    subroutine build_aerosol_table(base_name, n_ions, molar_mass, density)
        ! Fills the aerosol composition table from the per-type arrays read out of
        ! the aerosol input file, then checks it against the bin -> type mappings.
        ! The global aerosol_name attribute names type 1; further types are named
        ! from their index, since the v1 schema carries no per-type names.
        !
        ! A type is seed material iff the seed group points at it: the file cannot
        ! state that separately and so cannot contradict itself. Background types
        ! are whatever bin_type references, which is what lets a file carry several
        ! background materials alongside the seed.
        character(*), intent(in) :: base_name
        integer, intent(in) :: n_ions(:)
        real(dp), intent(in) :: molar_mass(:), density(:)
        character(20) :: type_name
        character(8) :: strint
        logical :: is_background, is_seed
        integer :: i, n_aerosol_types

        n_aerosol_types = size(n_ions)
        allocate(aerosols(n_aerosol_types))

        do i = 1, n_aerosol_types
            if (i == 1) then
                type_name = trim(base_name)
            else
                write(strint,'(i0)') i
                type_name = trim(base_name) // '_type' // trim(strint)
            end if
            call aerosols(i)%aerosol_initialize(trim(type_name), i, n_ions(i), &
                                                molar_mass(i), density(i), is_seed_type(i))
        end do

        do i = 1, n_aerosol_types
            is_background = any(aerosol_bin_type == i)
            is_seed = is_seed_type(i)

            ! An unreachable type describes material that can never be sampled, and
            ! a type in both groups has no answer to whether seed_hydration applies
            ! to it. Both are writer bugs worth catching at read time. The
            ! unreachable check is skipped when do_seeding is off: the seed group
            ! is then unread, so its composition row is legitimately unreferenced
            ! (ignored, not a bug), which is what lets one file serve both modes.
            if (do_seeding .and. .not. is_background .and. .not. is_seed) then
                write(error_unit,*) 'Error: aerosol type ', i, ' is referenced by neither bin_type nor seed_bin_type'
                call exit(1)
            end if
            if (is_background .and. is_seed) then
                write(error_unit,*) 'Error: aerosol type ', i, ' is referenced as both background and seed.', &
                                    ' Duplicate the composition row so each group has its own.'
                call exit(1)
            end if
        end do

    end subroutine build_aerosol_table

    logical function is_seed_type(itype)
        ! A composition row holds seed material iff a seed bin maps to it.
        integer, intent(in) :: itype

        is_seed_type = .false.
        if (n_seed_bins > 0) is_seed_type = any(seed_bin_type == itype)

    end function is_seed_type

    subroutine initialize_injection(inj_rate)
        ! 
        real(dp), intent(in) :: inj_rate(:)

        inj_time_idx = 1
        ! Get frequency of aerosol injection in seconds
        injection_dt = calculate_injection_rate(inj_rate(inj_time_idx), domain_volume)
        
        ! Note if the injection rate will change through the simulation,
        ! if so, this flag will check for changing Inj. Rate during runtime
        if ( size(injection_times) > 1 ) then
            update_inj_rate = .true.
        else
            update_inj_rate = .false.
        end if

        ! Initialize injection variables
        last_injection_time = 0.0
        n_injected = 0

    end subroutine initialize_injection

    pure function calculate_injection_rate(inj_rate, dom_vol) result(injection_time)
        ! Determines the injection rate, i.e. after how many seconds 1 aerosol should be injected
        real(dp), intent(in) :: inj_rate, dom_vol
        real(dp) :: injection_time

        if ( inj_rate == 0. ) then
            injection_time = 9999999.
        else
            injection_time = 1.0 / (dom_vol * inj_rate) ! seconds per number
        end if

    end function calculate_injection_rate

    subroutine calculate_droplet_statistics(droplets, stats)
        ! Summary droplet statistics for the output time series. Fills stats(:) by
        ! index (the order must match the writeout reader):
        !   1 = total particle count
        !   2 = activated count        3 = unactivated count
        !   4 = mean radius (um)       5 = liquid water content (g/m3)
        !   6 = collisions since last write   7 = coalescences since last write
        ! Resets the collision/coalescence accumulators after reading them.
        type(particle), intent(in) :: droplets(:)
        real(dp), intent(out) :: stats(:)
        integer(i4) :: i, Nact
        real(dp) :: r_sum, lwc_sum

        ! Nall
        stats(1) = current_n_particles

        ! Nactivated
        Nact = 0 ! Nactivated
        r_sum = 0. ! Integral radius
        lwc_sum = 0.
        do i = 1, current_n_particles
            r_sum = r_sum + droplets(i)%radius
            lwc_sum = lwc_sum + droplets(i)%water_liquid
            if ( droplets(i)%activated ) then
                Nact = Nact + 1
            end if
        end do
        stats(2) = Nact
        stats(3) = current_n_particles - Nact ! Unactivated N
        stats(4) = (r_sum / current_n_particles) * um_per_m ! r_bar
        stats(5) = lwc_sum * g_per_kg / domain_volume ! g/m3

        ! Collision-coalescence counts since last write (0 when CC is off)
        stats(6) = collisions_since_write
        stats(7) = coalescences_since_write
        collisions_since_write = 0
        coalescences_since_write = 0

    end subroutine calculate_droplet_statistics

    subroutine bin_droplet_radii(droplets, bin_edges, histogram)
        ! Sets up array of particle radii for binning algorithm
        ! no need to set histogram size as already determined in bin_data and read_binning_data
        type(particle), intent(in) :: droplets(:)
        real(dp), intent(in) :: bin_edges(:)
        integer(i4), intent(out) :: histogram(:,:)
        real(dp) :: radii(current_n_particles)
        real(dp), allocatable :: cat_radii(:)
        logical :: include_cat(current_n_particles)
        integer(i4) :: i, j, n_particles


        ! First calculate the total DSD for all drops
        ! Create array to feed binning function
        do concurrent (i = 1:current_n_particles)
            radii(i) = droplets(i)%radius * um_per_m ! microns and diameter bins
        end do

        histogram(1,:) = bin_data(bin_edges, radii)

        ! Then cycle through each aerosol category and create a DSD for that
        if ( n_aer_category > 1 ) then
            do j = 1, n_aer_category
                n_particles = 0
                include_cat = .false.

                ! Determine which droplets belong to which category
                do concurrent (i = 1:current_n_particles)
                    if ( droplets(i)%aerosol_category == j ) then
                        include_cat(i) = .true.
                    end if
                end do
                n_particles = count(include_cat)

                ! Create array of those droplet radii and bin them
                allocate(cat_radii(n_particles))
                cat_radii(:) = pack(radii, include_cat)
                histogram(j+1,:) = bin_data(bin_edges, cat_radii)
                deallocate(cat_radii)

            end do
        end if

    end subroutine bin_droplet_radii



end module droplets