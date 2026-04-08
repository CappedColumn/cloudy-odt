module particle_types
    use globals
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

        ! Collision-coalescence history
        integer(i4) :: n_collisions = 0              ! times this particle has collided
        integer(i4) :: n_coalescences = 0            ! times this particle was the keeper in a merge
        real(dp)    :: radius_before_coalescence = 0.0  ! radius before most recent coalescence

    contains
        ! All type-bound procedures begin with "particle_"
        ! Renamed for ease of reading when calling procedures
        procedure :: update_scalars => particle_update_scalars
        procedure :: critical_kohler => particle_critical_kohler
        procedure :: settling => particle_settling
        procedure :: update_gridcell => particle_update_gridcell
        procedure :: calculate_water_content => particle_calculate_water_content
        procedure :: verify_activation => particle_verify_activation
    end type particle

contains

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
        class(particle), intent(in) :: this
        real(dp) :: terminal_velocity
        real(dp) :: coeff, rho_air

        rho_air = pres / (this%virt_temp * Rd)
        coeff = 2. * g * rho_l / (9.0 * nu * rho_air)
        terminal_velocity = -coeff * this%radius**2

    end function calculate_terminal_velocity

    pure subroutine particle_update_gridcell(this)
        ! Given the position of a particle, update the particle's gridcell
        ! index to reflect that location
        class(particle), intent(inout) :: this

        this%gridcell = int(this%position / dz_length) + 1

        ! Particle position can be moved to 1.0,
        ! this ensures the gridcell index is not out of bounds
        if ( this%gridcell > N ) this%gridcell = N

    end subroutine particle_update_gridcell

    pure subroutine particle_update_scalars(this, Temp, Vapor, VirtTemp, Supersat)
        ! Update the particle properties based on the gridcell properties
        ! at the particle's location.
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

end module particle_types
