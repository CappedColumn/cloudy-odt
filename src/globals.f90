module globals
    use netcdf, only: nf90_noerr, nf90_strerror
    implicit none
    public

    ! Sets the default parameters of the model.
    ! Some parameters can be changed by the namelist
    ! during the initialize_simulation() subroutine.
    ! Also establishes variables for the simulation.

    ! There are a lot of global variables, primarily for
    ! convienence. We must trust our future selves to 
    ! act responsibly. 

    ! -----------------------------------------------
    ! ---- ESTABLISHING DATA TYPES FOR SIMULATION ---
    ! -----------------------------------------------

    !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
    integer, parameter :: sp = selected_real_kind(6, 37)
    !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
    integer, parameter :: dp = selected_real_kind(15, 307)
 
    !> Char length for integers, range -2⁷ to 2⁷-1; 8 bits
    integer, parameter :: i1 = selected_int_kind(2)
    !> Short length for integers, range -2¹⁵ to 2¹⁵-1; 16 bits
    integer, parameter :: i2 = selected_int_kind(4)
    !> Length of default integers, range -2³¹ to 2³¹-1; 32 bits
    integer, parameter :: i4 = selected_int_kind(9)
    
    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ------ Constants -------------------------------
    ! Change if you are not on planet earth
    ! or are working with strange chemicals.
    ! -----------------------------------------------

    real(dp), parameter :: nu = 1.488e-5           ! Kinematic Viscosity
    real(dp), parameter :: kT = 1.96e-5           ! Heat Diffusivity (Dt is timestep...)
    real(dp), parameter :: Dv = 2.2705e-5          ! Mass Diffusivity of water vapor
    real(dp), parameter :: Ndnu = 1.                ! For nondim diffusion of velocity
    real(dp), parameter :: Pr = nu/kT               ! Prandtl Number
    real(dp), parameter :: Sc = nu/Dv              ! Schmidt Number

    real(dp), parameter :: Ma = 28.96535e-3             ! Molecular mass of dry air - kg/mol
    real(dp), parameter :: Mw = 18.01528e-3             ! Molecular mass of water - kg/mol
    real(dp), parameter :: eps = Mw/Ma            ! For density calculations
    real(dp), parameter :: eps_tv = 0.608           ! For virtual temperature calculations
    real(dp), parameter :: Tice = 273.15            ! Freezing temperature - K

    real(dp), parameter :: g = 9.81               ! Gravitational Constant m/s2
    real(dp), parameter :: cp = 1005. !1004.             ! Specific heat of air J/kg/K - const. pressure
    real(dp), parameter :: cv = 718. !717.5             ! Specific heat of air J/kg/K - const. volume
    real(dp), parameter :: cp_wv = 1875.          ! Specific heat of water vapor J/kg/K  - const. pressure
    real(dp), parameter :: c_l = 4190.            ! Specific heat of water J/kg/K
    real(dp), parameter :: Lcond = 2.5e6 !2.5104e6       ! Latent heat of condensation J/kg
    real(dp), parameter :: Rv = 461.5             ! Individual gas constant of water vapor J/kg/K
    real(dp), parameter :: Rd = 287.0             ! Individual gas constant of dry air J/kg/K
    real(dp), parameter :: R_univ = 8.1344598     ! Universal gas constant
    real(dp), parameter :: rho_l = 1000.0         ! Density of water kg/m3

    real(dp), parameter :: pi = 3.1415926535897931      ! duh
    real(dp), parameter :: pi_4 = 12.566370614359172    ! 4*pi
    real(dp), parameter :: pi_43 = 4.1887902047863905   ! 4/3*pi

    real(dp), parameter :: alpha = 3.5e-3           ! Thermal Expansion Coefficient
    real(dp), parameter :: C2 = 1.5e3                ! Turbulent strength in ODT (dimensionless, squared): see Eq. 2.9 - Wunsch and Kerstein 2005
    real(dp), parameter :: C = sqrt(C2)
    real(dp) :: w_dim_factor
    real(dp), parameter :: ZC2 = 1.0e5                 ! ODT viscous cut-off parameter (dimensionless): see same equation/paper

    ! Aerosol Calculations/Constants
    real(dp), parameter :: a_RY = 3.3e-5  ! Rogers & Yau Eq. 6.7 alpha parameter
    real(dp), parameter :: nu_stokes = 1./(9. * nu) ! For fall-speed calculations


    ! Conversions
    real(dp), parameter :: m_per_cm = 1e-2
    real(dp), parameter :: m_per_mm = 1e-3
    real(dp), parameter :: m_per_mu = 1e-6
    real(dp), parameter :: um_per_m = 1e6
    real(dp), parameter :: m_per_nm = 1e-9
    real(dp), parameter :: nm_per_m = 1e9
    real(dp), parameter :: g_per_kg = 1e3
    real(dp), parameter :: kg_per_g = 1e-3

    ! Will likely change method of getting directory locations,
    ! but this will do for now
    character(*), parameter :: namelist_path = 'input/params.nml' ! Location of namelist file
    character(100) :: output_directory  ! Location of output direction

    ! -----------------------------------------------
    ! -----------------------------------------------


    ! -------------- Domain Parameters --------------
    ! These are default values, to be changed by 
    ! params.nml namelist in the initialization module
    ! -----------------------------------------------

    integer(i4) :: N = 6000         ! Number of Grid Cells
    integer(i4) :: Lmin = 6        ! 1/3 of Smallest Eddy Size (gridpoints)
    integer(i4) :: Lprob= 18        ! 1/3 of Most Probable Eddy Size (gridpoints)
    integer(i4) :: eddy_location = 1 ! Eddy index Location (M)
    integer(i4) :: eddy_length = 1   ! Eddy index Length (L)
    logical :: eddy_accepted = .false. ! Eddy acceptance flag
    logical :: write_eddies = .false.

    ! Note default values will yeild a domain of 1 cm^3, area_frac=2 gives 2 cm^3...
    real(dp), parameter :: domain_width = 0.001    ! Implied Domain Width (m)
    real(dp) :: volume_scaling = 1      ! Cross-Sectional Area scaling - controls volume
    real(dp) :: domain_volume           ! Volume of domain (m^3)
    real(dp) :: gridcell_volume         ! Volume of each grid cell (m^3)

    real(dp) :: Tdiff = 10.        ! Top-Bottom Temperature Difference (Celsius)
    real(dp) :: Tref = 15.          ! Bottom temperature used for reference
    real(dp) :: Ttop
    real(dp) :: WVref, WVtop, WVdiff
    real(dp) :: Tvref, Tvtop, Tvdiff
    real(dp) :: pres = 1.00e5       ! Pressure (Pa)
    real(dp) :: H = 1.            ! Domain Height (meters)
    real(dp) :: dz_length    ! Length of each grid cell (meters)

    real(dp) :: max_accept_prob = 0.1      ! Upper constraint for stability in eddy accpt. method

    logical :: same_random = .false.    ! Will use random numbers seeded from same state if true 

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! -------------- Microphysics Parameters --------------
    ! These are default values, to be changed by 
    ! params.nml namelist in the initialization module
    ! -----------------------------------------------------

    logical :: do_turbulence = .true.
    logical :: do_microphysics = .true.
    real(dp) :: injection_rate

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! -------------- Special Effects --------------
    ! Called from primary namelists, then initialized
    ! and run from the special effects namelist
    ! -----------------------------------------------------

    logical :: do_special_effects = .false.

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ------------------ Iterators ------------------
    ! Variables that are used to progress simulation
    ! and test for conditions
    ! -----------------------------------------------

    real(dp) :: tmax = 30            ! Maximum simulation time (seconds)
    real(dp) :: tmax_nd     ! Maxmimum simulation time (non-dimensional, params.nml)
    real(dp) :: time    ! Dimensional time
    real(dp) :: time_nd        ! Non-dimensional time
    real(dp) :: last_time       ! Time lagging by dt
    real(dp) :: dt                 ! Dimensional time step
    real(dp) :: dt_nd             ! Non-Dimensional time step
    real(dp) :: delta_time_nd
    real(dp) :: delta_time
    real(dp) :: diffusion_step             ! Difussive Time Step
    integer(i4) :: Nt = 0       ! Number of timesteps
    integer(i4) :: Nd = 0       ! Number of Diffusion Calls

    integer(i4) :: Np = 0       ! Number of eddy probs. calculated (acc/rej method)
    integer(i4) :: Na = 0       ! Number of accepted eddies
    real(dp) :: Pa = 0.         ! Total acceptance probabilitiess

    ! -----------------------------------------------
    ! -----------------------------------------------


    ! ----------- Calculated Parameters -------------
    ! Will be calculated in initialization module once
    ! namelist parameters are fully updated
    ! -----------------------------------------------

    integer(i4) :: Lmax      ! Largest Eddy Size (1/3 of domain, in gridpoints)
    integer(i4) :: LpD          ! Twice the most probable length
    real(dp) :: time_conv_nd ! Convert between dimensional and nd-time
    real(dp) :: buoy_nd           ! Dimensionless Buoyancy
    real(dp) :: prob_coeff          ! Used in calculation of eddy acceptance probability
    real(dp) :: Co, Cm              ! Used for initial eddy sample guess
    real(dp) :: t_diff       ! Diffusion Time Step

    ! Eddy acceptance/rejection related
    real(dp), allocatable :: prob_eddy_length(:)      ! Probability of eddy sizes
    real(dp) :: accept_prob                 ! Eddy Acceptance Prob. for eddy w/ (L, M)
    real(dp) :: pot_energy                  ! Potential Energy
    real(dp) :: tot_accept_prob             ! Tracker of the aggregate acceptance probabilities
    integer(i4) :: num_accept_prob          ! Tracker of the number of probabilities calculated
    
    ! -----------------------------------------------
    ! -----------------------------------------------
    
    ! ------------------- ARRAYS --------------------

    ! Velocity arrays
    real(dp), allocatable :: W(:), Wdim(:)   ! Velocity components

    ! Positional arrays
    real(dp), allocatable :: z(:)

    ! Scalar arrays
    ! Temperature, Water Vapor (dim and non-dim) and Virt. Temp, Supersaturation
    real(dp), allocatable :: T(:), WV(:), Tv(:), Tdim(:), WVdim(:), Tvdim(:)
    real(dp), allocatable :: SS(:)

    ! Statistics
    real(dp) :: statistics(5) ! N, Na, Nu, r_bar, LWC

    ! Integrated Eddy Values
    real(dp) :: wK, TvK

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! Read/Write Arrays

    integer :: ncid, ncid_particles
    ! writout iterators
    real(dp) :: write_time_iter = 0. ! iterator for write out
    real(dp) :: write_timer


contains

    subroutine nc_verify(status, error_msg)
        integer, intent(in) :: status
        character(*), intent(in), optional :: error_msg

        ! All netcdf function calls return a status code, rather
        ! than assigning the status to a variable name, wrapping
        ! this subroutine around the netCDF function call will verify
        ! its execution, and return useful error messages

        if (status /= nf90_noerr) then
            write(*,'(a)') "Error in netCDF procedure..."
            write(*,'((a), (a), (i4))') error_msg, ' :: ', status
            write(*,'(a)') trim(nf90_strerror(status))
            stop
        end if

    end subroutine nc_verify

    function bin_data(bin_edges, data) result(histogram)
        ! Given the bin edges and data values, returns and histogram
        ! with the number of data points found in each bin
        real(dp), intent(in) :: bin_edges(:), data(:)
        integer(i4) :: i, j, n_edges, n_data
        integer(i4), allocatable :: data_bin(:), histogram(:)

        n_edges = size(bin_edges)
        n_data = size(data)
        
        allocate(histogram(n_edges-1))
        histogram = 0

        ! Assign the bin for which each data point will be assigned
        do i = 1, n_data
            do j = 1, n_edges-1
                if ( (data(i)>=bin_edges(j)) .and. (data(i)<bin_edges(j+1)) ) then
                    histogram(j) = histogram(j) + 1
                    exit
                end if
            end do
        end do

        ! FASTER??
        ! allocate(data_bin(n_data))
        ! do concurrent (i=1:n_data, j=1:n_edges-1)
        !     if ( (data(i) >= bin_edges(j)) .and. (data(i) < bin_edges(j+1)) ) then
        !         data_bin(i) = j
        !     end if
        ! end do
        ! ! Calculate frequency for each bin
        ! do i = 1, n_data
        !     histogram(data_bin(i)) = histogram(data_bin(i)) + 1
        ! end do

    end function




end module globals