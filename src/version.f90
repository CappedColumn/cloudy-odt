module version
    implicit none

    character(*), parameter :: code_version = '0.6.0'
    character(*), parameter :: git_commit   = '894cafc'

contains

    subroutine print_usage()
        write(*,'(a)') 'CODT - Cloudy One-Dimensional Turbulence - v', trim(code_version), ' (', trim(git_commit), ')'
        write(*,'(a)') ''
        write(*,'(a)') 'Usage: codt <namelist_path>'
        write(*,'(a)') '       codt --help | --version'
        write(*,'(a)') ''
        write(*,'(a)') 'Run a simulation configured by the Fortran namelist at <namelist_path>.'
        write(*,'(a)') 'The path must include a directory (e.g. ./params.nml).'
        write(*,'(a)') ''
        write(*,'(a)') 'Namelist groups:'
        write(*,'(a)') '  &PARAMETERS       Domain, timing, output, and physics switches'
        write(*,'(a)') '  &TURBULENCE_ODT   ODT turbulence parameters (chamber mode)'
        write(*,'(a)') '  &TURBULENCE_LEM   LEM turbulence parameters (parcel mode)'
        write(*,'(a)') '  &MICROPHYSICS     Aerosol injection, collision-coalescence, trajectories'
        write(*,'(a)') '  &PARCEL           Parcel ascent and initial conditions (parcel mode)'
        write(*,'(a)') '  &ENTRAINMENT      Entrainment/detrainment parameters (parcel mode)'
        write(*,'(a)') '  &SPECIALEFFECTS   Sidewall nudging and stochastic fallout'
        write(*,'(a)') '  &RADIATION        Radiative transfer parameters'
        write(*,'(a)') ''

    end subroutine print_usage

end module version
