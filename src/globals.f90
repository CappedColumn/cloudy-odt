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

    character(256) :: namelist_path    ! Path to namelist file (set from command line)
    character(256) :: namelist_dir     ! Parent directory of namelist file (for resolving relative paths)
    character(256) :: output_directory ! Base output directory from namelist (absolute path)
    character(100) :: simulation_name  ! Simulation name from namelist
    character(256) :: sim_output_dir   ! {output_directory}/{simulation_name}/ — where all output files live
    character(256) :: file_prefix      ! {sim_output_dir}{simulation_name} — base path for output files (.nc, .nml, etc.)

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
    logical :: overwrite = .false.      ! Allow overwriting existing output files

    ! Simulation mode: 'chamber' (ODT, fixed BCs) or 'parcel' (LEM, periodic BCs)
    character(7) :: simulation_mode = 'chamber'

    ! LEM parameters (only used when simulation_mode = 'parcel')
    real(dp) :: integral_length_scale = 0.01_dp       ! Largest eddy size (m)
    real(dp) :: kolmogorov_length_scale = 0.001_dp    ! Smallest eddy size (m)
    real(dp) :: dissipation_rate = 0.01_dp            ! TKE dissipation rate (m^2/s^3)

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
    real(dp) :: time    ! Dimensional time
    real(dp) :: last_time_updated  ! Time of last physics update (diffusion/eddy event)
    real(dp) :: dt                 ! Dimensional time step
    real(dp) :: delta_time         ! Dimensional time since last diffusion
    real(dp) :: diffusion_step     ! Diffusive time step (dimensional, seconds)
    integer(i4) :: Nt = 0       ! Number of timesteps
    integer(i4) :: Nd = 0       ! Number of Diffusion Calls

    ! -----------------------------------------------
    ! -----------------------------------------------


    ! ----------- Calculated Parameters -------------
    ! Will be calculated in initialization module once
    ! namelist parameters are fully updated
    ! -----------------------------------------------

    integer(i4) :: Lmax      ! Largest Eddy Size (1/3 of domain, in gridpoints)
    integer(i4) :: LpD          ! Twice the most probable length
    real(dp) :: buoy_nd           ! Dimensionless Buoyancy
    real(dp) :: prob_coeff          ! Used in calculation of eddy acceptance probability
    real(dp) :: Co, Cm              ! Used for initial eddy sample guess

    ! Eddy acceptance/rejection related
    real(dp), allocatable :: prob_eddy_length(:)      ! Probability of eddy sizes
    
    ! -----------------------------------------------
    ! -----------------------------------------------
    
    ! ------------------- ARRAYS --------------------

    ! Velocity arrays
    real(dp), allocatable :: W_nd(:)   ! ODT velocity (nondim, energy conservation)

    ! Positional arrays
    real(dp), allocatable :: z(:)

    ! Scalar arrays
    ! Temperature, Water Vapor (dim and non-dim) and Virt. Temp, Supersaturation
    real(dp), allocatable :: T_nd(:), WV_nd(:), Tv_nd(:), T(:), WV(:), Tv(:)
    real(dp), allocatable :: SS(:)

    ! Statistics
    real(dp) :: statistics(7) ! N, Na, Nu, r_bar, LWC, N_coll, N_coal

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! Read/Write Arrays

    integer :: ncid
    ! writout iterators
    real(dp) :: write_timer

    ! ----------- Budget Accumulators -----------------
    ! Accumulated over each write interval, then reset.
    ! -------------------------------------------------
    integer(i4), parameter :: n_budgets = 12
    real(dp) :: budget_inject_solute_mass = 0.0
    real(dp) :: budget_inject_liquid_mass = 0.0
    real(dp) :: budget_fallout_liquid_mass = 0.0
    real(dp) :: budget_fallout_solute_mass = 0.0
    real(dp) :: budget_condensation = 0.0
    real(dp) :: budget_dgm_delta_T = 0.0
    real(dp) :: budget_diffusion_delta_T = 0.0
    real(dp) :: budget_diffusion_delta_WV = 0.0
    real(dp) :: budget_sidewall_delta_T = 0.0
    real(dp) :: budget_sidewall_delta_WV = 0.0
    integer(i4) :: budget_n_injected = 0
    integer(i4) :: budget_n_fellout = 0

    ! ----------- Turbulence Dispatch ----------------
    ! Abstract interfaces for mode-agnostic turbulence calls.
    ! Pointers are set once in initialize_simulation().
    ! ------------------------------------------------

    abstract interface
        subroutine diffuse_iface(ldelta_time)
            import :: dp
            real(dp), intent(in) :: ldelta_time
        end subroutine

        subroutine turbulence_iface(ldt, ltime, ldelta_time, &
                                    leddy_accepted, eddy_loc, eddy_len)
            import :: dp, i4
            real(dp), intent(inout) :: ldt
            real(dp), intent(in) :: ltime, ldelta_time
            logical, intent(out) :: leddy_accepted
            integer(i4), intent(out) :: eddy_loc, eddy_len
        end subroutine

        subroutine sync_iface()
        end subroutine
    end interface

    procedure(diffuse_iface), pointer :: diffuse_step => null()
    procedure(turbulence_iface), pointer :: turbulence_step => null()
    procedure(sync_iface), pointer :: sync_after_physics => null()

    ! -----------------------------------------------
    ! -----------------------------------------------

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
            stop 1
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

    function parent_directory(filepath) result(dir)
        ! Returns the parent directory of a file path.
        ! e.g. "/home/user/input/params.nml" -> "/home/user/input/"
        character(*), intent(in) :: filepath
        character(256) :: dir
        integer :: idx

        idx = scan(trim(filepath), '/', back=.true.)
        if (idx > 0) then
            dir = filepath(1:idx)
        else
            dir = './'
        end if
    end function parent_directory

    function resolve_path(basedir, filepath) result(full_path)
        ! If filepath is absolute, return it as-is.
        ! If relative, prepend basedir.
        character(*), intent(in) :: basedir, filepath
        character(256) :: full_path

        if (filepath(1:1) == '/') then
            full_path = filepath
        else
            full_path = trim(basedir) // trim(filepath)
        end if
    end function resolve_path

    subroutine copy_file(source, destination)
        ! Copies a file byte-for-byte from source to destination using stream I/O
        character(*), intent(in) :: source, destination
        integer :: in_unit, out_unit, ierr
        character :: byte

        open(newunit=in_unit, file=trim(source), status='old', access='stream', &
             form='unformatted', action='read', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: could not open source file: ', trim(source)
            stop 1
        end if

        open(newunit=out_unit, file=trim(destination), status='replace', access='stream', &
             form='unformatted', action='write', iostat=ierr)
        if (ierr /= 0) then
            write(0,*) 'Error: could not open destination file: ', trim(destination)
            close(in_unit)
            stop 1
        end if

        do
            read(in_unit, iostat=ierr) byte
            if (ierr /= 0) exit
            write(out_unit) byte
        end do

        close(in_unit)
        close(out_unit)
    end subroutine copy_file


    subroutine triplet_map(eddy_length, eddy_start, field)
        ! Applies the triplet map rearrangement to field.
        ! Uses mod indexing so wrapping eddies on periodic domains are handled
        ! automatically. For non-periodic domains the mod is a no-op.
        integer(i4), intent(in) :: eddy_length, eddy_start
        real(dp), intent(inout) :: field(:)

        real(dp) :: mapped_values(eddy_length)
        integer(i4) :: j, source_index, dest_index, segment_length

        segment_length = eddy_length / 3

        ! Segment 1: every 3rd element, forward
        do j = 1, segment_length
            source_index = mod(eddy_start + 3*(j-1) - 1, N) + 1
            mapped_values(j) = field(source_index)
        end do

        ! Segment 2: every 3rd element, reversed (block inversion)
        do j = 1, segment_length
            source_index = mod(eddy_start + eddy_length - 3*j, N) + 1
            mapped_values(j + segment_length) = field(source_index)
        end do

        ! Segment 3: every 3rd element, forward offset by 2
        do j = 1, segment_length
            source_index = mod(eddy_start + 3*j - 2, N) + 1
            mapped_values(j + 2*segment_length) = field(source_index)
        end do

        ! Write rearranged values back
        do j = 1, eddy_length
            dest_index = mod(eddy_start + j - 2, N) + 1
            field(dest_index) = mapped_values(j)
        end do

    end subroutine triplet_map


    subroutine reset_budgets()
        budget_inject_solute_mass = 0.0
        budget_inject_liquid_mass = 0.0
        budget_fallout_liquid_mass = 0.0
        budget_fallout_solute_mass = 0.0
        budget_condensation = 0.0
        budget_dgm_delta_T = 0.0
        budget_diffusion_delta_T = 0.0
        budget_diffusion_delta_WV = 0.0
        budget_sidewall_delta_T = 0.0
        budget_sidewall_delta_WV = 0.0
        budget_n_injected = 0
        budget_n_fellout = 0
    end subroutine reset_budgets

end module globals