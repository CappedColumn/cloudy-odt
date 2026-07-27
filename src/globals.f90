! Central definitions shared across the whole model: precision kinds, physical
! constants, the global field arrays and scalar state, the gridcell triplet map,
! the abstract interfaces for the turbulence/diffusion/sync procedure pointers,
! and small utilities (nc_verify, resolve_path, namelist_read_error). Almost
! every other module uses this one; keep it dependency-light.
module globals
    use iso_fortran_env, only: error_unit
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

    ! ------ Vertical-axis identifiers ---------------
    ! Coordinate the parcel trajectory-leg targets are specified on
    ! (&PARCEL vertical_axis).
    integer(i4), parameter :: AXIS_HEIGHT = 1   ! metres
    integer(i4), parameter :: AXIS_PRESSURE = 2 ! Pa

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
    real(dp), parameter :: Pa_per_mb = 100.0
    real(dp), parameter :: m_per_km = 1e3

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

    integer(i4) :: N = 2000         ! Number of Grid Cells
    logical :: write_eddies = .false.

    ! Note default values will yeild a domain of 1 cm^3, area_frac=2 gives 2 cm^3...
    real(dp), parameter :: domain_width = 0.001    ! Implied Domain Width (m)
    real(dp) :: volume_scaling = 10      ! Cross-Sectional Area scaling - controls volume
    real(dp) :: domain_volume           ! Volume of domain (m^3)
    real(dp) :: gridcell_volume         ! Volume of each grid cell (m^3)

    real(dp) :: Tref = 20.          ! Bottom temperature used for reference
    real(dp) :: pres = 1.00e5       ! Pressure (Pa)
    real(dp) :: H = 1.            ! Domain Height (meters)
    real(dp) :: dz_length    ! Length of each grid cell (meters)

    logical :: same_random = .false.    ! Will use random numbers seeded from same state if true
    logical :: overwrite = .false.      ! Allow overwriting existing output files

    ! Simulation mode: 'chamber' (ODT, fixed BCs) or 'parcel' (LEM, periodic BCs)
    character(7) :: simulation_mode = 'chamber'

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
    logical :: do_radiation = .false.
    logical :: do_entrainment = .false.

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ------------------ Iterators ------------------
    ! Variables that are used to progress simulation
    ! and test for conditions
    ! -----------------------------------------------

    real(dp) :: tmax = 100            ! Maximum simulation time (seconds)
    real(dp) :: time    ! Dimensional time
    real(dp) :: last_time_updated  ! Time of last physics update (diffusion/eddy event)
    real(dp) :: dt                 ! Dimensional time step
    real(dp) :: delta_time         ! Dimensional time since last diffusion
    real(dp) :: diffusion_step     ! Diffusive time step (dimensional, seconds)
    integer(i4) :: Nt = 0       ! Number of timesteps
    integer(i4) :: Nd = 0       ! Number of Diffusion Calls

    ! -----------------------------------------------
    ! -----------------------------------------------

    ! ------------------- ARRAYS --------------------

    ! Positional arrays
    real(dp), allocatable :: z(:)

    ! Scalar arrays
    ! Temperature, Water Vapor, Virtual Temperature, Supersaturation
    real(dp), allocatable :: T(:), WV(:), Tv(:)
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

    ! Field budgets (always active)
    integer(i4), parameter :: n_field_budgets = 4
    real(dp) :: budget_diffusion_delta_T = 0.0
    real(dp) :: budget_diffusion_delta_WV = 0.0
    real(dp) :: budget_sidewall_delta_T = 0.0
    real(dp) :: budget_sidewall_delta_WV = 0.0

    ! Microphysics budgets (only when do_microphysics = .true.)
    integer(i4), parameter :: n_micro_budgets = 9
    real(dp) :: budget_inject_solute_mass = 0.0
    real(dp) :: budget_inject_liquid_mass = 0.0
    real(dp) :: budget_fallout_liquid_mass = 0.0
    real(dp) :: budget_fallout_solute_mass = 0.0
    real(dp) :: budget_condensation = 0.0
    real(dp) :: budget_dgm_delta_T = 0.0
    integer(i4) :: budget_n_injected = 0
    integer(i4) :: budget_n_fellout = 0
    integer(i4) :: budget_n_coalesced = 0

    ! Entrainment budgets — accumulated per output interval, reset in reset_budgets.
    ! Only active when do_entrainment = .true.
    integer(i4), parameter :: n_entrain_budgets = 6
    real(dp) :: budget_detrain_liquid_mass = 0.0    ! liquid water removed by detrainment [kg]
    real(dp) :: budget_detrain_solute_mass = 0.0    ! solute mass removed by detrainment [kg]
    real(dp) :: budget_entrain_liquid_mass = 0.0    ! liquid water added by entrainment [kg]
    real(dp) :: budget_entrain_solute_mass = 0.0    ! solute mass added by entrainment [kg]
    integer(i4) :: budget_n_detrained = 0           ! particles removed by detrainment
    integer(i4) :: budget_n_entrained = 0           ! particles added by entrainment

    ! Radiation budget (only when do_radiation = .true.)
    real(dp) :: budget_radiation_delta_T = 0.0

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

    ! Generic triplet map: the real variant rearranges scalar fields (T, WV, ...);
    ! the integer variant rearranges a cell-label tracer using the identical
    ! permutation, so droplet transport can compose maps by advecting labels.
    interface triplet_map
        module procedure triplet_map_real, triplet_map_int
    end interface triplet_map

    ! Cell-label tracer backing the eddy-sequence API below. Shared by ODT (one
    ! eddy per turbulence step) and LEM (several), so both move droplets through
    ! exactly the same code.
    !   origin_cell(j)      = the cell whose fluid now sits at j   (gather view)
    !   destination_cell(c) = where the fluid originally in c ended up (forward)
    integer(i4), allocatable :: origin_cell(:)
    integer(i4), allocatable :: destination_cell(:)

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
            write(error_unit,'(a)') 'Error in netCDF procedure...'
            write(error_unit,'((a), (a), (i4))') error_msg, ' :: ', status
            write(error_unit,'(a)') trim(nf90_strerror(status))
            call exit(1)
        end if

    end subroutine nc_verify


    subroutine namelist_read_error(nml_unit, group_name)
        integer, intent(in) :: nml_unit
        character(*), intent(in) :: group_name

        character(256) :: bad_line, scan_line
        character(64) :: upper_group
        logical :: group_found
        integer :: i, ic, io_stat

        ! Grab the offending line before we rewind
        backspace(nml_unit)
        read(nml_unit,'(a)', iostat=io_stat) bad_line

        ! Check whether the group exists in the file
        upper_group = group_name
        do i = 1, len_trim(upper_group)
            ic = iachar(upper_group(i:i))
            if (ic >= iachar('a') .and. ic <= iachar('z')) &
                upper_group(i:i) = achar(ic - 32)
        end do

        group_found = .false.
        rewind(nml_unit)
        do
            read(nml_unit, '(a)', end=10) scan_line
            scan_line = adjustl(scan_line)
            if (scan_line(1:1) == '&') then
                do i = 2, len_trim(scan_line)
                    ic = iachar(scan_line(i:i))
                    if (ic >= iachar('a') .and. ic <= iachar('z')) &
                        scan_line(i:i) = achar(ic - 32)
                end do
                if (trim(scan_line(2:)) == trim(upper_group)) then
                    group_found = .true.
                    exit
                end if
            end if
        end do
        10 continue
        close(nml_unit)

        if (.not. group_found) then
            write(error_unit,'(a,a,a)') &
                'Error: namelist group &', trim(group_name), ' not found in namelist file.'
            write(error_unit,'(a)') 'Check that the group exists and is spelled correctly.'
        else
            write(error_unit,'(a,a,a)') &
                'Invalid parameter in &', trim(group_name), ':'
            write(error_unit,'(a,a)') '  ', trim(adjustl(bad_line))
        end if
        call exit(1)

    end subroutine namelist_read_error


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


    subroutine triplet_map_real(eddy_length, eddy_start, field)
        ! Applies the triplet map rearrangement to a real scalar field.
        ! Uses mod indexing so wrapping eddies on periodic domains are handled
        ! automatically. For non-periodic domains the mod is a no-op.
        ! Invoke via the generic name triplet_map.
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

    end subroutine triplet_map_real

    subroutine triplet_map_int(eddy_length, eddy_start, field)
        ! Integer counterpart of triplet_map_real, applying the identical
        ! permutation to a cell-label tracer. Keep the index arithmetic in
        ! lockstep with triplet_map_real: the whole point of the tracer is that
        ! it undergoes exactly the same rearrangement as the scalar fields.
        ! Invoke via the generic name triplet_map.
        integer(i4), intent(in) :: eddy_length, eddy_start
        integer(i4), intent(inout) :: field(:)

        integer(i4) :: mapped_values(eddy_length)
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

    end subroutine triplet_map_int


    ! -----------------------------------------------
    ! Eddy sequence: composing one or more triplet maps into the single net cell
    ! rearrangement used to displace droplets.
    !
    ! Usage, identical for ODT (one eddy) and LEM (several):
    !
    !     call begin_eddy_sequence()
    !     ... for each accepted eddy, beside the scalar triplet_map calls:
    !         call accumulate_eddy(eddy_length, eddy_start)
    !     call finalize_eddy_sequence()
    !     call move_particles_by_cellmap(particles, destination_cell)
    ! -----------------------------------------------

    subroutine begin_eddy_sequence()
        ! Start a new sequence with the tracer at the identity.
        integer(i4) :: c

        ! Reallocate if N has changed since the last sequence. In a run N is
        ! fixed, so this costs one size comparison per event; the unit tests do
        ! vary N between cases.
        if (allocated(origin_cell)) then
            if (size(origin_cell) /= N) deallocate(origin_cell)
        end if
        if (allocated(destination_cell)) then
            if (size(destination_cell) /= N) deallocate(destination_cell)
        end if
        if (.not. allocated(origin_cell)) allocate(origin_cell(N))
        if (.not. allocated(destination_cell)) allocate(destination_cell(N))

        do c = 1, N
            origin_cell(c) = c
        end do

    end subroutine begin_eddy_sequence


    subroutine accumulate_eddy(eddy_length, eddy_start)
        ! Fold one eddy into the sequence. Call with the same arguments as the
        ! triplet_map calls on the scalar fields, so the tracer and the scalars
        ! stay in lockstep.
        !
        ! Composition order needs no special handling: applying the maps to the
        ! tracer in the order the eddies occurred yields the correct composed
        ! gather automatically. Written out explicitly the nesting runs in
        ! reverse (eddies 1,2,3 compose as s1(s2(s3(z)))), which is easy to get
        ! backwards -- so do not replace this with hand-rolled index arithmetic.
        integer(i4), intent(in) :: eddy_length, eddy_start

        call triplet_map(eddy_length, eddy_start, origin_cell)

    end subroutine accumulate_eddy


    subroutine finalize_eddy_sequence()
        ! Invert the composed tracer into the forward map used to move droplets.
        !
        ! triplet_map is a gather: it fills each cell with the value pulled from
        ! its source, so after a sequence origin_cell(j) says what arrived at j.
        ! That receiver's view is exactly what a scalar field needs -- every cell
        ! is filled with whatever landed in it.
        !
        ! A droplet asks the opposite question: not "what arrived here?" but
        ! "where did my fluid go?". That is the sender's view, and it is the
        ! inverse permutation, which is what this builds.
        !
        ! The two coincide when the permutation is its own inverse -- true for
        ! eddy lengths 3 and 6, false for 9 and above.
        !
        ! origin_cell is a permutation of 1..N, so destination_cell is too: the
        ! droplet population is neither lost nor duplicated.
        integer(i4) :: j

        do j = 1, N
            destination_cell(origin_cell(j)) = j
        end do

    end subroutine finalize_eddy_sequence


    subroutine reset_budgets()
        budget_diffusion_delta_T = 0.0
        budget_diffusion_delta_WV = 0.0
        budget_sidewall_delta_T = 0.0
        budget_sidewall_delta_WV = 0.0
        if (do_microphysics) then
            budget_inject_solute_mass = 0.0
            budget_inject_liquid_mass = 0.0
            budget_fallout_liquid_mass = 0.0
            budget_fallout_solute_mass = 0.0
            budget_condensation = 0.0
            budget_dgm_delta_T = 0.0
            budget_n_injected = 0
            budget_n_fellout = 0
            budget_n_coalesced = 0
        end if
        if (do_radiation) then
            budget_radiation_delta_T = 0.0
        end if
        if (do_entrainment) then
            budget_detrain_liquid_mass = 0.0
            budget_detrain_solute_mass = 0.0
            budget_entrain_liquid_mass = 0.0
            budget_entrain_solute_mass = 0.0
            budget_n_detrained = 0
            budget_n_entrained = 0
        end if
    end subroutine reset_budgets

end module globals