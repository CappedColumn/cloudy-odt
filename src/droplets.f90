module droplets
    use netcdf
    use globals
    use DGM, only: integrate_ODE, set_aerosol_properties
    use special_effects, only: do_random_fallout, random_fallout_rate
    use microphysics
    implicit none

    type :: aerosol
        ! Aerosol user-defined type. Aerosol is used by the 
        ! particle type and attributes are used in calculations
        character(20) :: aerosol_name = 'NA'
        integer(i4) :: aerosol_id = 0
        integer(i4) :: n_ions = 0
        real(dp) :: solute_molar_mass = 0.0
        real(dp) :: solute_density = 0.0 ! kg/m3
    contains
        ! All type-bound procedures begin with "aerosol_"
        ! Renamed for ease of reading when calling procedures
        procedure :: aerosol_initialize ! not renamed due to "initialize" procedure for particle below
    end type aerosol

    type, extends(aerosol) :: particle
        ! Droplet particle user-defined type, contains the unique attributes
        ! for each particle in the simulation
        integer(i4) :: particle_id

        ! Particle properties shared with gridcell
        real(dp) :: position = 0.0          ! z coordinate, m
        integer(i4) :: gridcell = 0      ! Gridcell index, i
        real(dp) :: temperature = 0.0       ! Temperature of the particle, K
        real(dp) :: water_vapor = 0.0      ! Water vapor content, kg/kg
        real(dp) :: virt_temp = 0.0         ! Virtual temperature, K
        real(dp) :: supersaturation = 0.0   ! Supersaturation level, %

        ! particle properties unique to each particle
        real(dp) :: critical_radius = 0.0
        real(dp) :: critical_supersaturation = 0.0
        real(dp) :: water_liquid = 0.0      ! Absolute liquid water, kg
        real(dp) :: radius = 0.0            ! Radius of the particle, m

        ! Particle properties pertaining to dissolved aerosol
        type(aerosol) :: solute_type = aerosol()    ! Type of aerosol particle
        real(dp) :: solute_gross_mass = 0.0      ! Gross mass of the solute in the particle (kg)
        real(dp) :: solute_radius = 0.0        ! Radius of the dry-solute in the particle (m)
        integer(i4) :: aerosol_category = 1

        logical :: activated = .false.
        logical :: fellout = .false.

    contains
        ! All type-bound procedures begin with "particle_"
        ! Renamed for ease of reading when calling procedures
        procedure :: initialize => particle_initialize
        procedure :: update_scalars => particle_update_scalars
        procedure :: critical_kohler => particle_critical_kohler
        procedure :: settling => particle_settling
        procedure :: update_gridcell => particle_update_gridcell
        procedure :: calculate_water_content => particle_calculate_water_content
        procedure :: verify_activation => particle_verify_activation
    end type particle

    ! Counters to track particles, used for statistics and array indexing
    integer(i4) :: current_n_particles = 0
    integer(i4) :: total_n_particles = 0
    integer(i4) :: total_n_fellout = 0
    real(dp), allocatable :: particle_bin_edges(:), particle_bins(:)
    integer(i4), allocatable :: size_distribution(:,:) !(No. DSDs, rbins)
    integer(i4) :: n_DSD_bins, n_aer_category

    ! Arrays to hold all particles and aerosol types
    type(particle), allocatable :: particles(:) ! allocated in initialize_microphysics()
    type(aerosol), allocatable :: aerosols(:) ! allocated in read_aerosol_netcdf()
    integer(i4) :: particle_array_expansion = 1000

    ! Aerosol injection variables
    real(dp), allocatable :: injection_times(:), injection_rates(:) ! allocated in read_aerosol_netcdf()
    real(dp), allocatable :: aerosol_size_edges(:), aerosol_bin_freq(:,:) ! allocated in read_aerosol_netcdf()
    real(dp), allocatable :: aerosol_radii(:) ! allocated in read_aerosol_netcdf()
    integer(i4), allocatable :: aerosol_partition(:)
    real(dp) :: last_injection_time
    real(dp) :: injection_dt
    integer(i4) :: n_injected, inj_time_idx
    logical :: update_inj_rate
    real(dp) :: initial_wet_radius
    logical :: init_drop_each_gridpoint = .true.
    real(dp) :: expected_Ndrops_per_gridpoint = 1

    ! DGM-Controlling variables
    real(dp), parameter :: RK5_min_timestep = 0.01

    ! Particle I/O Handling
    logical :: write_trajectories
    real(dp) :: trajectory_start = 0., trajectory_end = 0.
    real(dp) :: trajectory_timer = 1.

        ! private
    !public :: particles, aerosol, particle, inject_particle, initialize_aerosol_type, &
    !total_n_fellout, droplet_growth_model, initialize_injection_rate, n_injected

    ! Used in Main - write_trajectories

contains

    subroutine update_droplets(ltime, ldt)
        ! Interface subroutine to main.f90. Does aerosol injection, settling,
        ! droplet-environment property update, and droplet growth.
        ! Caller is responsible for syncing nondim fields afterward.
        real(dp), intent(in) :: ltime, ldt

        call injection_controller(time, particles)
        call move_particles_by_gravity(particles, ldt)
        call update_all_particles(particles, T, WV, Tv, SS)
        call droplet_growth_model(particles, ltime, ldt)

    end subroutine update_droplets

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
                call inject_particle(lparticles, T, WV, Tv, SS, aerosols(1))
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

    subroutine inject_particle(lparticles_array, Temp, Vapor, VirtTemp, Supersat, aerosol_type)
        ! Injects a new particle into the simulation at a random location. Particle inherets
        ! the properties of the gridcell which it originates
        type(particle), allocatable, intent(inout) :: lparticles_array(:)
        real(dp), intent(in) :: Temp(:), Vapor(:), VirtTemp(:), Supersat(:)
        type(aerosol), intent(in) :: aerosol_type
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
            ! Set up temporary array to hold current particles
            allocate(temp_array(size(lparticles_array)))
            temp_array = lparticles_array
            deallocate(lparticles_array)
            ! Array Expansion and reassignment
            allocate(lparticles_array(current_n_particles + particle_array_expansion))
            lparticles_array(:size(temp_array)) = temp_array
            deallocate(temp_array)
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
        call injected_particle%initialize(aerosol_type, total_n_particles, random_position, grid_idx, Temp(grid_idx), &
                    Vapor(grid_idx), VirtTemp(grid_idx), Supersat(grid_idx))

        call injected_particle%update_gridcell()

        ! Accumulate budget
        budget_inject_solute_mass = budget_inject_solute_mass + injected_particle%solute_gross_mass
        budget_inject_liquid_mass = budget_inject_liquid_mass + injected_particle%water_liquid
        budget_n_injected = budget_n_injected + 1

        ! Index newly injected particle object into array of current particles
        lparticles_array(current_n_particles) = injected_particle

    end subroutine inject_particle


    pure subroutine particle_critical_kohler(this)
        ! Calculate the critical radius and supersaturation for the particle
        ! based on Equations 6.7-6.8 in Rogers & Yau
        class(particle), intent(inout) :: this
        real(dp) :: radius, supersaturation
        real(dp) :: a, b ! m

        ! Determine numerical approximations
        a = a_RY / (this%temperature)
        b = 4.3 * this%solute_gross_mass * this%solute_type%n_ions / this%solute_type%solute_molar_mass

        ! Calculate the critical radius and supersaturation
        radius = sqrt(3 * b / a) * m_per_cm ! meters
        supersaturation = sqrt(4 * a**3 / (27 * b)) * 100 ! %

        ! Assign to particle properties
        this%critical_radius = radius
        this%critical_supersaturation = supersaturation

    end subroutine particle_critical_kohler

    pure subroutine particle_verify_activation(this)
        ! Changes the activated flag for a particle, depending on its
        ! current size and critical radius.
        class(particle), intent(inout) :: this

        if ( this%radius >= this%critical_radius ) then
            this%activated = .true.
        else
            this%activated = .false.
        end if

    end subroutine particle_verify_activation

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
        do concurrent (i = 1:current_n_particles)
            call lparticles(i)%settling(ldt)
        end do

        ! Verify if particle fellout of domain
        ! Updates number of particles currently in domain
        call verify_particle_fallout(lparticles, current_n_particles)

        ! For remaining particles, update gridcell index and properties
        do concurrent (i = 1:current_n_particles)
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
        integer :: i, n_just_fellout
        real(dp) :: r, new_position

        ! Update fallout flag for particles which settled out of the domain
        do concurrent (i = 1:n_particles)
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

        ! Count the number of particles that fellout
        ! and reassign array positions for particles still in domain
        n_just_fellout = 0 ! How many particles fellout since last check
        do i = 1, n_particles
            ! If particle fellout, increment counters
            if ( lparticle_array(i)%fellout ) then
                total_n_fellout = total_n_fellout + 1
                n_just_fellout = n_just_fellout + 1
                budget_fallout_liquid_mass = budget_fallout_liquid_mass + lparticle_array(i)%water_liquid
                budget_fallout_solute_mass = budget_fallout_solute_mass + lparticle_array(i)%solute_gross_mass
                budget_n_fellout = budget_n_fellout + 1
            else
                ! Reassign array position of particles still in domain to 
                ! fill in gaps of particles that fell out, saves space 
                ! and removes particles that fellout from further computation
                if (n_just_fellout > 0) then
                    lparticle_array(i - n_just_fellout) = lparticle_array(i)
                end if
            end if
        end do

        ! Update number of particles in domain
        n_particles = n_particles - n_just_fellout

    end subroutine verify_particle_fallout

    pure subroutine random_fallout(lparticle, fallout_rate, height)
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

    pure subroutine particle_settling(this, ldt)
        ! Determines particle fall speed and updates particle position
        ! with explicit scheme
        class(particle), intent(inout) :: this
        real(dp), intent(in) :: ldt
        real(dp) :: terminal_velocity

        ! Determine the terminal velocity of the particle
        terminal_velocity = calculate_terminal_velocity(this)

        ! Update the particle position
        this%position = this%position + (terminal_velocity * ldt)

    end subroutine particle_settling

    pure function calculate_terminal_velocity(this) result(terminal_velocity)
        ! Calculate the stokes terminal fall velocity of a particle using
        ! virtual effects
        use globals, only: g, nu_stokes, pres, Rd, rho_l
        class(particle), intent(in) :: this
        real(dp) :: terminal_velocity
        real(dp) :: coeff, rho_air

        rho_air = pres / (this%virt_temp * Rd)
        coeff = 2. * g * rho_l / (9.0 * nu * rho_air)
        terminal_velocity = -coeff * this%radius**2


    end function calculate_terminal_velocity

    subroutine move_particles_in_eddy(lparticles, M, L)
        ! Move particles within an eddy based on the triplet map.
        ! The triplet map moves gridcell z1 -> z2. A particle at pz1 in that
        ! cell carries its offset: pz2 = z2 + (pz1 - z1), then wraps periodic.
        type(particle), intent(inout) :: lparticles(:)
        integer(i4), intent(in) :: M, L
        real(dp) :: mapped_z(N)
        integer :: i, gc

        mapped_z = z
        call triplet_map(L, M, mapped_z)

        do i = 1, current_n_particles
            gc = lparticles(i)%gridcell
            if (mapped_z(gc) /= z(gc)) then
                lparticles(i)%position = modulo( &
                    mapped_z(gc) + (lparticles(i)%position - z(gc)), H)
            end if
        end do

    end subroutine move_particles_in_eddy

    pure subroutine particle_update_gridcell(this)
        ! Given the position of a particle, update the particle's gridcell
        ! index to reflect that location
        !
        ! Input:
        ! this - particle type
        !
        ! Output:
        ! this - particle type with updated gridcell index
        use globals, only: dz_length
        class(particle), intent(inout) :: this

        this%gridcell = int(this%position / dz_length) + 1

        ! Particle position can be moved to 1.0,
        ! this ensures the gridcell index is not out of bounds
        if ( this%gridcell > N ) this%gridcell = N

    end subroutine particle_update_gridcell

    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! SUBROUTINES TO INITIALIZE PARTICLES AND AEROSOLS !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine aerosol_initialize(this, name, id, n_ions, molar_mass, density)
        ! Type-bound procedure for aerosol. Initializes aerosol properties with
        ! passed values.
        class(aerosol), intent(out) :: this
        character(len=*), intent(in) :: name
        integer(i4), intent(in) :: id
        integer(i4), intent(in) :: n_ions
        real(dp), intent(in) :: molar_mass
        real(dp), intent(in) :: density

        this%aerosol_name = name
        this%aerosol_id = id
        this%n_ions = n_ions
        this%solute_molar_mass = molar_mass
        this%solute_density = density

    end subroutine aerosol_initialize


    subroutine particle_initialize(this, solute, ln_particles, pos, grid_idx, temp, vapor, virt_temp, supersat)

        class(particle), intent(out) :: this
        integer(i4), intent(in) :: ln_particles
        real(dp), intent(in) :: pos
        integer(i4), intent(in) :: grid_idx
        real(dp), intent(in) :: temp
        real(dp), intent(in) :: vapor
        real(dp), intent(in) :: virt_temp
        real(dp), intent(in) :: supersat
        class(aerosol), intent(in) :: solute

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

        ! Assign the selected aerosol type to the particle
        this%solute_type = solute

        ! Currently in domain
        this%fellout = .false.

        ! Determine solute properties (sampled from aerosol input) and initial radius
        call sample_radius(aerosol_bin_freq(:,inj_time_idx), aerosol_radii, this%solute_radius, this%aerosol_category)
        this%solute_gross_mass = ( pi_43 * this%solute_type%solute_density ) * this%solute_radius**3
        this%radius = initial_wet_radius * this%solute_radius

        ! Determine liquid water content of the particle
        call this%calculate_water_content()

        ! Determine the critical radius and supersaturation for this particle
        call this%critical_kohler()

    end subroutine particle_initialize

    subroutine sample_radius(bin_freq, radii, radius, aer_partition)
        ! Selects an aerosol size from CDF distribution in aerosol input file
        real(dp), intent(in) :: bin_freq(:), radii(:)
        real(dp), intent(out) :: radius
        integer(i4), intent(out) :: aer_partition
        real(dp) :: random_num
        integer :: idx

        ! Generate a random number between 0 and 1
        call random_number(random_num)
        
        idx = 1
        do while ( random_num .gt. bin_freq(idx) )
            idx = idx + 1
        end do

        radius = radii(idx) * m_per_nm
        aer_partition = aerosol_partition(idx)

    end subroutine sample_radius

    !-----------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SUBROUTINES TO UPDATE PARTICLE PROPERTIES !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pure subroutine particle_update_scalars(this, Temp, Vapor, VirtTemp, Supersat)
        ! Update the particle properties based on the gridcell properties
        ! at the particle's location.
        !
        ! Input:
        ! this - particle type
        ! Temp - temperature at particle location, K
        ! Vapor - water vapor content at particle location, kg/kg
        ! VirtTemp - virtual temperature at particle location, K
        ! Supersat - supersaturation at particle location, %
        !
        ! Output:
        ! this - particle type with updated properties
        class(particle), intent(inout) :: this
        real(dp), intent(in) :: Temp, Vapor, VirtTemp, Supersat

        this%temperature = Temp
        this%water_vapor = Vapor
        this%virt_temp = VirtTemp
        this%supersaturation = Supersat

        ! Determine the critical radius and supersaturation for this particle
        call this%critical_kohler()

        ! Determine liquid water content of the particle
        call this%calculate_water_content()

    end subroutine particle_update_scalars

    pure subroutine particle_calculate_water_content(this)
        ! Determine particle liquid water content
        class(particle), intent(inout) :: this

        this%water_liquid = pi_43 * (this%radius**3 - this%solute_radius**3) * rho_l

    end subroutine particle_calculate_water_content
    
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
        use DGM, only: integrate_ODE
        type(particle), intent(inout) :: droplet
        real(dp), intent(in) :: ltime, ldt
        real(dp) :: time_start, time_stop, time_iterate
        real(dp) :: y_arr(8), y_before(8), inverse_grid_mass, grid_rho
        logical :: substep_flag

        ! Determine mass of air in gridcell
        grid_rho = pres / (Rd * droplet%virt_temp)
        inverse_grid_mass = 1.0 / (gridcell_volume*grid_rho)

        ! set odeint parameters for different aerosol mass of each droplet
        call set_aerosol_properties(1, droplet%solute_gross_mass, inverse_grid_mass)

        ! Package droplet properties into an array for ode solver
        y_arr(1) = droplet%radius
        y_arr(2) = droplet%water_vapor
        y_arr(3) = droplet%temperature
        y_arr(4) = droplet%supersaturation/100
        y_arr(5) = pres
        y_arr(6) = 0.0 ! vertical velocity
        y_arr(7) = 0.0 ! height
        y_arr(8) = droplet%water_liquid
        y_before = y_arr

        ! In an attempt to please the continuum of time
        ! within a world of discrete tendencies,
        ! we may have to iterate the iterator.
        ! Eddies occur so frequently that this rarely happens.

        ! Test to ensure timestep is not too large for Runge-Kutta solver
        ! must be <1 sec, Mani - chose 0.01 sec to be safe
        substep_flag = .false.
        if ( dt >= RK5_min_timestep ) then
            substep_flag = .true.
        end if

        if ( substep_flag ) then 
            ! Set time bounds for iterative calls to ode solver
            time_stop = ltime + ldt
            time_start = ltime
            time_iterate = ltime + RK5_min_timestep
            ! Iterate through time in RK5_min_timestep steps
            do while ( time_iterate < time_stop )
                ! Call with RK5_min_timestep dt
                call integrate_ODE(y_arr, time_start, time_iterate, RK5_min_timestep)
                time_start = time_iterate
                time_iterate = time_iterate + RK5_min_timestep
            end do
            ! Call ode solver for remainder of time
            if ( time_stop > time_start .and. time_stop < time_iterate ) then
                call integrate_ODE(y_arr, time_start, time_stop, time_stop - time_start)
            end if

        else ! Call me maybe
            call integrate_ODE(y_arr, ltime, ltime + ldt, ldt)
        end if

        ! Unpack droplet properties from array
        droplet%radius = y_arr(1)
        droplet%water_vapor = y_arr(2)
        droplet%temperature = y_arr(3)
        droplet%water_liquid = y_arr(8)
        droplet%supersaturation = calc_supersat(droplet%temperature, droplet%water_vapor, pres)
        droplet%virt_temp = virtual_temp(droplet%temperature, droplet%water_vapor)

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

        do concurrent (i = 1:current_n_particles)
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
        
        character(256):: aerosol_file, bin_data_file ! Not allocatable since namelist-specified variable

        namelist /MICROPHYSICS/ init_drop_each_gridpoint, expected_Ndrops_per_gridpoint, aerosol_file, &
        bin_data_file, write_trajectories, trajectory_start, trajectory_end, trajectory_timer, initial_wet_radius

        ! Read in microphysical namelist parameters
        write(*,*) 'Reading MICROPHYSICS namelist values...'
        open(newunit=nml_unit, file=namelist_path, iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) io_emsg; stop 1
        end if
        read(nml=MICROPHYSICS, unit=nml_unit, iostat=ierr)
        ! Print value causing namelist read error
        if (ierr .ne. 0) then
            write(*,*) ierr
            backspace(nml_unit)
            read(nml_unit,'(a)') nml_line
            write(*,'(a)') 'Invalid Namelist Parameter: ', trim(nml_line)
            stop 1
        end if
        close(nml_unit)

        ! Namelist is copied to output directory in initialize_params
        

        ! Verify that the initial wet radius will be greater than dry radius
        if ( initial_wet_radius <= 1.) then
            write(*,*) 'WARNING: initial_wet_radius must be > 1.0 (got ', initial_wet_radius, ')'
            write(*,*) '         Auto-correcting to 1.1x dry radius.'
            initial_wet_radius = 1.1
        end if

        ! Resolve input paths relative to namelist directory
        aerosol_file = resolve_path(namelist_dir, trim(aerosol_file))
        bin_data_file = resolve_path(namelist_dir, trim(bin_data_file))

        ! Copy aerosol input to output directory
        i = scan(trim(aerosol_file), '/', back=.true.)
        call copy_file(trim(aerosol_file), trim(sim_output_dir)//trim(aerosol_file(i+1:)))

        ! Set up aerosol type and injection forcings
        call read_aerosol_netcdf(trim(aerosol_file))

        ! Bring in binning data
        n_aer_category = maxval(aerosol_partition)
        call read_binning_data(trim(bin_data_file), n_aer_category, particle_bin_edges, size_distribution)
        
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

        ! Set up starting injection rate
        call initialize_injection(injection_rates)

        ! Size particle array based on the expected number of droplets
        allocate(particles(int(expected_Ndrops_per_gridpoint*N)))

        ! Initialize particles in each gridpoint (approximately)
        if ( init_drop_each_gridpoint ) then
            do i = 1, N
                call inject_particle(particles, T, WV, Tv, SS, aerosols(1))
            end do
        end if

    end subroutine initialize_microphysics

    subroutine netcdf_add_DSD(lncid, r_bins)
        ! Adds a DSD variable to the netcdf file
        integer, intent(in) :: lncid
        real(dp), intent(in) :: r_bins(:)

        integer :: t_dimid, r_dimid, r_varid, dsd_varid, nbins, dimids(2)
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

        integer :: t_dimid, r_dimid, dsd_varid, nbins, dimids(2)

        ! Get dimension IDs
        call nc_verify( nf90_inq_dimid(lncid, "time", t_dimid))
        call nc_verify( nf90_inq_dimid(lncid, "radius", r_dimid))

        ! Open netcdf in definition mode, and create a DSD for each aerosol partition
        dimids = (/ r_dimid, t_dimid /)
        call nc_verify( nf90_redef(lncid), "nf90_redef: DSD_aerr" )
        do i = 1, n_DSDs
            write(strint,*) i
            name = "DSD_" // adjustl(strint)
            call nc_verify( nf90_def_var(lncid, trim(name), NF90_INT, dimids, dsd_varid, &
                            deflate_level=1, shuffle=.true.), "nf90_def_var: DSD_aer" )
            name = "Droplet Size Distribution - " // adjustl(strint)
            call nc_verify( nf90_put_att(lncid, dsd_varid, "long_name", trim(name)), "nf90_put_att: DSD_aer, name")
            call nc_verify( nf90_put_att(lncid, dsd_varid, "units", "#"), "nf90_put_att: DSD_aer, units")

        end do
        call nc_verify( nf90_enddef(lncid), "nf90_enddef: DSD_aer")

    end subroutine netcdf_add_aerDSD

    subroutine read_binning_data(location, n_cat, bin_edges, DSD)
        ! Acquires the bin edges for particle/droplet size distribution calculations
        character(*), intent(in) :: location
        integer, intent(in), value :: n_cat
        real(dp), allocatable, intent(out) :: bin_edges(:)
        integer, allocatable, intent(out) :: DSD(:,:)
        integer :: i, ierr, file_unit
        character(100) :: io_emsg
        logical :: file_exists

        integer :: n_bin_edges

        ! Open the bin data file
        write(*,*) 'Reading bin data from: ', location
        inquire(file=location, exist=file_exists)
        if ( .not. file_exists ) write(*,*) 'Bin data file does not exist: ', location
        open(newunit=file_unit, file=trim(location), iostat=ierr, iomsg=io_emsg, action='read', status='old')
        if (ierr .ne. 0) then
            write(*,*) 'Error opening injection data file: ', io_emsg
            stop 1
        end if

        read(file_unit, *) ! N Bin-Edges
        read(file_unit, *) n_bin_edges
        n_DSD_bins = n_bin_edges - 1
        if ( n_cat > 1 ) n_cat = n_cat + 1 ! Account for total and categorical DSDs
        allocate(DSD(n_cat, n_DSD_bins))
        DSD = 0
        
        ! Allocate and read in bin edge values
        allocate(bin_edges(n_bin_edges))
        read(file_unit, *) ! Bin-Edges
        do i = 1, n_bin_edges
            read(file_unit, *) bin_edges(i)
        end do

        close(file_unit)

    end subroutine read_binning_data

    subroutine read_aerosol_netcdf(filepath)
        ! Reads aerosol properties, injection schedule, and size distribution
        ! from a NetCDF file (CODT_aerosol_input_v1 schema).
        character(*), intent(in) :: filepath
        integer :: i, aer_ncid, varid, dimid
        integer :: n_types, n_bins, n_edges, n_times
        character(64) :: conventions
        character(20) :: aer_name
        integer :: aer_n_ions
        real(dp) :: aer_molar_mass, aer_density

        write(*,*) 'Reading aerosol data from: ', trim(filepath)

        ! Open and validate schema version
        call nc_verify(nf90_open(trim(filepath), NF90_NOWRITE, aer_ncid), &
                       'opening aerosol file')
        call nc_verify(nf90_get_att(aer_ncid, NF90_GLOBAL, 'conventions', conventions), &
                       'reading conventions attribute')
        if (trim(conventions) /= 'CODT_aerosol_input_v1') then
            write(*,*) 'Error: expected CODT_aerosol_input_v1, got: ', trim(conventions)
            stop 1
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
            write(*,*) 'Error: edge dimension must equal bin + 1'
            stop 1
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

        allocate(aerosol_bin_freq(n_bins, n_times))
        call nc_verify(nf90_inq_varid(aer_ncid, 'cumulative_frequency', varid), &
                       'finding cumulative_frequency')
        call nc_verify(nf90_get_var(aer_ncid, varid, aerosol_bin_freq), &
                       'reading cumulative_frequency')

        ! Read aerosol properties for type 1
        call nc_verify(nf90_inq_varid(aer_ncid, 'n_ions', varid), 'finding n_ions')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_n_ions, start=[1]), 'reading n_ions')

        call nc_verify(nf90_inq_varid(aer_ncid, 'molar_mass', varid), 'finding molar_mass')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_molar_mass, start=[1]), 'reading molar_mass')

        call nc_verify(nf90_inq_varid(aer_ncid, 'solute_density', varid), 'finding solute_density')
        call nc_verify(nf90_get_var(aer_ncid, varid, aer_density, start=[1]), 'reading solute_density')

        ! Read aerosol name from global attribute
        call nc_verify(nf90_get_att(aer_ncid, NF90_GLOBAL, 'aerosol_name', aer_name), &
                       'reading aerosol_name')

        call nc_verify(nf90_close(aer_ncid), 'closing aerosol file')

        ! Calculate midpoint radii (nanometers)
        allocate(aerosol_radii(n_bins))
        do i = 1, n_bins
            aerosol_radii(i) = (aerosol_size_edges(i) + aerosol_size_edges(i+1)) / 2.0_dp
        end do

        ! Initialize aerosol type
        allocate(aerosols(1))
        call aerosols(1)%aerosol_initialize(trim(aer_name), 1, aer_n_ions, &
                                            aer_molar_mass, aer_density)

    end subroutine read_aerosol_netcdf

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
        ! Calculate a few droplet statistics to be included in the netcdf file
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