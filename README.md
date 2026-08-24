# Cloudy One-Dimensional Turbulence

[![Fortran](https://img.shields.io/badge/Fortran-734f96.svg?logo=fortran&logoColor=white)](https://fortran-lang.org)
[![netCDF](https://img.shields.io/badge/netCDF-output-1f72b8.svg)](https://www.unidata.ucar.edu/software/netcdf/)

A numerical model for warm-cloud microphysics in turbulent flows. CODT was built to simulate convection cloud chambers, particularly the [Michigan Tech Pi-Chamber](https://doi.org/10.1175/BAMS-D-15-00203.1), a NSF Community Instruments and Facility (CIF) resource, and also supports adiabatic parcel ascent. One spatial dimension allows for kolmogorov-scale resolution of turbulent fluctuations, coupled to a droplet growth model for individual aerosol species. Since CODT replicates the physics (and microphysics) of the Pi-Chamber reasonably well, while requiring basic compute resources (*i.e.* 1 CPU), it is best suited for scientists looking to test/understand multiple experimental cases which may not be possible to perform in the laboratory.


## Features

- **Chamber mode** — turbulent Rayleigh-Benard convection at the kolmogorov scale via One-Dimensional Turbulence (ODT)
- **Parcel mode** — adiabatic ascent with Linear Eddy Model (LEM) turbulence and optional entrainment
- Lagrangian microphysics (not-so-super droplets) with explicit droplet growth
- Collision-coalescence
- Optional radiative effects
- Specifiable aerosol species and concentrations
- Written in [Modern Fortran](https://fortran-lang.org) with netCDF compatability


## Installation and Dependencies

CODT uses [netCDF](https://docs.unidata.ucar.edu/netcdf-fortran/current/) for data output, but otherwise has no dependencies. Post-processing is handled by the companion **CODT_tools** Python package, which reads CODT's netCDF output (profiles and time series, droplet trajectories) and the collision/eddy binary streams.

For ease of build, it is highly recommended to use the [Fortran Package Manager](https://fpm.fortran-lang.org). The source files in `app` and `src` are highly modular, and FPM tracks build dependencies between each source file.

Build flags live in the fpm manifest (`fpm.toml`) as **profiles**, selected with `--profile`:

| Profile | Purpose | gfortran flags |
|---------|---------|----------------|
| `release`  | optimized production build | `-O2 -finline-functions` |
| `debug`    | no optimization, runtime checks | `-g -O0 -fcheck=all -fbacktrace` |
| `profiled` | release plus gprof instrumentation | `-O2 … -pg` |

Two files hold the machine-specific settings; copy each template and edit it for your site. Both copies are gitignored, so your local paths and module names never get committed:

- **`fpm.toml`** (from `fpm.toml.template`) — put your netCDF include/lib paths in the `netcdf-local` feature. All optimization/debug/profiling flags are already generic and need no editing.
- **`build.sh`** (from `build.sh.template`) — maps each `fpm.toml` profile to the compiler module it needs, then builds. Edit the `case` block near the top.

```bash
cp fpm.toml.template fpm.toml        # then edit the netcdf-local paths
cp build.sh.template build.sh        # then edit the profile -> module case block
chmod +x build.sh

./build.sh                           # release + gfortran (the default), version-stamped
./build.sh debug                     # debug build (bounds/backtrace)
./build.sh release nvfortran         # pick profile and compiler positionally
```

`build.sh` takes `[profile] [compiler]`, loads the compiler module that profile needs, stamps the version into `src/version.f90`, and runs `fpm build`. For anything beyond those two knobs, call fpm directly — but pass `--profile`:

```bash
fpm build --profile release --compiler gfortran --verbose
```

> A bare `fpm build` (no `--profile`) is an **unoptimized** debug build, and also omits the `netcdf-local` feature, so it fails to find `netcdf.mod`. Always pass `--profile`, or just use `./build.sh`.

Every profile in `fpm.toml` needs a matching entry in `build.sh`, and vice versa: gfortran `.mod` files are not portable across gcc major versions, so a profile must be built under the same compiler its netCDF was.

`build.sh` derives the version from git tags: a released build (a clean checkout of a tagged commit) reports a bare version such as `v1.0.0`, while development and fork builds self-identify with a commit hash and a `-dirty` suffix. These values are written into the global attributes of CODT's netCDF output for provenance.

## Running CODT

Simulations are setup entirely from a Fortran namelist (yes, a namelist, the original YAML/TOML/XML/CSON/JSON configuration file). A template lives at `input/params.nml`; it sets model parameters/physics, write locations, and simulation naming conventions. CODT takes the namelist path as its single argument:

```bash
fpm run --profile release --compiler gfortran -- input/params.nml   # via fpm
CODT input/params.nml            # installed executable
CODT --help                      # usage and namelist groups
CODT --version                   # version and git commit
```

CODT redirects its own stdout to `{simulation_name}.log` in the output directory, so no shell redirection is needed. Running the executable with no argument prints usage and exits.

## Resources and Documentation

This implementation of One-Dimensional Turbulence (ODT) is based on [Wunsch & Kerstein (2005)](https://doi.org/10.1017/S0022112004003258), and was further developed for moist conditions by [Chandrakar et al. (2020)](https://doi.org/10.1017/jfm.2019.895). The original code has been rewritten to modern fortran, and the cloud droplet growth model of [Su et al. (1998)](https://doi.org/10.1016/S0169-8095(98)00039-8) was added as the microphysics scheme. The parcel mode with entrainment is based off of the Explicit Mixing Parcel Model (EMPM; [Krueger et al. 1997](https://journals.ametsoc.org/view/journals/atsc/54/23/1520-0469_1997_054_2697_meafmi_2.0.co_2.xml), [Su et al. (1998)](https://doi.org/10.1016/S0169-8095(98)00039-8)). Full physical documentation, including comparison between CODT and DNS, will be available in a forthcoming paper.

## Contributing

CODT was developed with other researcher's interests in mind, and there is the potential to add new physics such as ice microphysics, aqueous-phase chemistry, new aerosol species, etc. If you would like to contribute, please fork the repository and open a GitHub pull request. Release versions are tagged by the maintainer; development builds are automatically marked as such (see *Installation and Dependencies*).

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
