# Cloudy One-Dimensional Turbulence

A numerical model for simulating convection cloud chambers, particularly the [Michigan Tech Pi-Chamber](https://doi.org/10.1175/BAMS-D-15-00203.1), a NSF Community Instruments and Facility (CIF) resource. CODT is designed for studying warm-cloud microphysics in the presence of turbulence. One spatial dimension allows for kolmogorov-scale resolution of turbulent fluctuations, coupled to a droplet growth model for individual aerosols species. Since CODT replicates the physics (and microphysics) of the Pi-Chamber reasonably well, while requiring basic compute resources (*i.e.* 1 CPU), it is best suited for scientists looking to test/understand multiple experimental cases which may not be possible to perform in the laboratory. 


## Features

- Simulates turbulent Rayleigh-Benard Convection at the kolmogorov scale
- Specifiable aerosol species and concentrations
- Lagrangian microphysics (not-so-super droplets)
- Written in [Modern Fortran](https://fortran-lang.org) with netCDF compatability


## Installation and Dependencies

CODT uses [netCDF](https://docs.unidata.ucar.edu/netcdf-fortran/current/) for certain data output, but otherwise has no dependencies. Post-processing jupyter notebooks are provided in the `utilities` folder, and require [numpy](https://numpy.org) and [netCDF4](https://unidata.github.io/netcdf4-python/) to generate ragged-array netCDF files of droplet trajectories.

For ease of build, it is highly recommended to use the [Fortran Package Manager](https://fpm.fortran-lang.org). The source files in `app` and `src` are highly modular, and FPM tracks build dependencies between each source file. If you choose to use FPM, the provided `fpm_env` bash script can be used to set the location of your fortran compilers and netcdf build. Simply run the script to establish the fpm build environment prior to running the `fpm build` command.

## Running CODT

Once built, simulations are setup entirely from `input/params.nml` (yes, a namelist, the original YAML/TOML/XML/CSON/JSON configuration file). This includes model parameters/physics, write locations and simulation naming conventions. CODT can then be executed from the command line by running the executable with no arguments, `./path/to/executable/CODT > stdout/file/SIM_NAME.log`.

## Resources and Documentation

This implementation of One-Dimensional Turbulence (ODT) is based on [Wunsch & Kerstein (2005)](https://doi.org/10.1017/S0022112004003258), and was further developed for moist conditions by [Chandrakar et al. (2020)](https://doi.org/10.1017/jfm.2019.895). The original code has been rewritten to modern fortran, and the cloud droplet growth model of [Su et al. (1998)](https://doi.org/10.1016/S0169-8095(98)00039-8) was added as the microphysics scheme. Full physical documentation, including comparison between CODT and DNS, will be available in a forthcoming paper.

## Contributing

CODT was developed with other researcher's interests in mind, and there is the potential to add new physics such as ice microphysics, aqueous-phase chemistry, new aerosol species, etc. If you would like to contribute, please use GitHub pull requests.

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)