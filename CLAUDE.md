# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

CODT uses the Fortran Package Manager (fpm). Before building, source the environment script to set compiler and NetCDF paths:

```bash
source fpm_env          # sets FPM_FC, FPM_FFLAGS, FPM_LDFLAGS
fpm build               # compiles to build/ directory
```

The `fpm_env` script is configured for CHPC (University of Utah) modules. Edit the `compiler` variable and `base_netcdf_path` for other systems. Supported compilers: `nvfortran`, `gfortran`, `ifort`.

## Running

```bash
./build/*/app/CODT > output.log
```

All simulation parameters are set in `input/params.nml` (three namelists: `PARAMETERS`, `MICROPHYSICS`, `SPECIALEFFECTS`). The executable takes no arguments.

## Architecture

**Dual representation of fields:** Scalar fields exist in both non-dimensional (for ODT numerics) and dimensional (for physics) forms. Arrays `T`, `WV`, `Tv`, `W` are non-dimensional; `Tdim`, `WVdim`, `Tvdim`, `Wdim` are dimensional. Always update both after modifying either (see `update_dim_scalars` in `microphysics.f90`).

**Time loop structure** (`app/main.f90`): Each iteration advances `time_nd` by `dt_nd`. Two event types trigger physics:
1. **Diffusion step** — fires when `delta_time_nd >= diffusion_step`. Runs Crank-Nicolson tridiagonal diffusion on all scalar fields, then updates particles.
2. **Eddy step** — Monte-Carlo accept/reject of a triplet-map eddy. On acceptance: diffusion + triplet map + particle rearrangement.

Both paths call `update_droplets()` when microphysics is enabled.

**Particle system** (`src/droplets.f90`): `aerosol` base type holds solute properties; `particle` extends it with position, radius, thermodynamic state. Particles are stored in an allocatable array and track their grid cell index. Injection reads from `input/injection_data.txt`; binning uses `input/bin_data.txt`.

**Droplet Growth Model** (`src/DGM.f90`): Cash-Karp RK45 adaptive ODE integrator (Numerical Recipes heritage). The state vector `ystart(1:8)` carries: radius, qv, temperature, supersaturation, pressure, vertical velocity, height, liquid water. `set_aerosol_properties()` must be called before `integrate_ODE()` to configure solute type (`dmaxa`: 1=NaCl, 2=(NH4)2SO4, 3=ammonium bisulfate).

**Module dependency chain:** `globals` → `microphysics` → `ODT` → `droplets` → `DGM`. Most state lives in `globals` as module-level variables. `special_effects` is optional and controlled by namelist flags.

**Output:** Buffered NetCDF writes (`src/writeout.f90`) for profiles and droplet size distributions. Optional particle trajectory output (`src/write_particle.f90`) controlled by `write_trajectories` namelist flag.

## Fortran Conventions

- Precision kinds defined in `globals.f90`: `dp` (double), `sp` (single), `i4` (32-bit int), `i1`/`i2` for smaller ints.
- Physical constants are `parameter` values in `globals.f90` — do not duplicate them elsewhere.
- NetCDF calls should be wrapped with `nc_verify()` from `globals.f90` for error handling.
- The codebase uses implicit none everywhere; maintain this.


## Interface Contract with codt_tools

This section defines how the Python `codt_tools` package invokes the compiled
executable. **Both projects must honor this contract.** The corresponding spec
lives in `codt_tools/CLAUDE.md` section 2.6.

### Invocation

```
codt <NAMELIST_PATH>
```

- The executable reads the namelist from the provided path.
- All relative paths inside the namelist (`inj_data_path`, `bin_data_path`)
  are resolved **relative to the namelist's parent directory**.
- `output_directory` must be an **absolute path**.
- If no argument is provided, print an error message and stop with a
  non-zero exit code. Do **not** silently fall back to `./params.nml`.

### Fortran Implementation Requirements

1. Check `command_argument_count()` on startup
2. If an argument is present, use `get_command_argument(1, namelist_path)`
3. Open and read the namelist from `namelist_path`
4. Resolve `inj_data_path` and `bin_data_path` relative to the namelist's
   parent directory (not CWD)
5. If no argument is provided, print a usage message to stderr and
   `stop 1`

### Path Resolution Example

Given invocation: `codt /scratch/runs/sim01/input/params.nml`

And namelist contents:
```
inj_data_path = "injection_data.txt"
bin_data_path = "bin_data.txt"
output_directory = "/scratch/output"
```

The executable resolves:
- `inj_data_path` → `/scratch/runs/sim01/input/injection_data.txt`
- `bin_data_path` → `/scratch/runs/sim01/input/bin_data.txt`
- `output_directory` → `/scratch/output` (used as-is)

### Status

**Implemented.** The executable requires a namelist path as its first
argument. Relative input paths resolve from the namelist's parent directory.
`output_directory` must be absolute. Bare filenames (no `/`) are rejected.

---

## Planned Modifications

### ~~1. Command-line namelist path (interface contract above)~~ Done

### ~~2. Add namelist parameters as netCDF global attributes~~ Done

### ~~3. Store bin edges in netCDF~~ Done

Bin edges stored as `radius_edges` (201 values) alongside existing
`radius` bin centers (200 values).

### ~~4. Fix "celcius" typo~~ Done

### ~~5. Convert eddy output to binary stream~~ Done

Eddy output converted from formatted text (`_eddies.txt`) to unformatted
binary stream (`_eddies.bin`). Stores raw grid indices (M, L) and
dimensional time for direct replay, with a header containing N, H, C2,
ZC2, Tdiff, Tref.

### 6. Add predetermined eddies mode

Add a simulation mode where eddies are read from an existing `_eddies.bin`
file instead of generated via Monte Carlo acceptance. This enables running
microphysics on pre-computed turbulence fields. Requires:
- A new namelist flag (e.g., `use_predetermined_eddies`) and a path
  parameter (e.g., `eddy_file`) in the `PARAMETERS` namelist
- On startup, read the eddy file header and validate that `N` matches
- In the time loop, replace `eddy_acceptance_method` with sequential
  reads from the eddy file — call `implement_eddy(L, M)` when simulation
  time reaches the next eddy's timestamp
- Skip `lower_dt`/`raise_dt` since there is no acceptance sampling
- `write_eddies` and `use_predetermined_eddies` should be mutually exclusive

### 7. Add eddy data reader to codt_tools

Add an `eddy_io.py` module to the Python `codt_tools` package to read
`_eddies.bin` files. Should return the header as a dict and eddy records
as a structured numpy array using `numpy.frombuffer`. Integrate into
`CODTSimulation` as a `load_eddies()` method.

---

## Output File Reference

| File | Format | Description |
|------|--------|-------------|
| `{name}.nc` | netCDF4 | Main output: profiles + time series |
| `{name}_particles.nc` | netCDF4 | Particle-level data (if `write_trajectories=.true.`) |
| `{name}_eddies.bin` | Binary (stream) | Accepted eddy events (if `write_eddies=.true.`) |
| `{name}.nml` | ASCII | Copy of namelist |

All output is written to `{output_directory}/{simulation_name}/`.

---

## Cross-Project Sync

When committing changes to CODT, update `~/dev/CODT_tools/CLAUDE.md` at the end of the session to reflect any changes that affect the codt_tools interface. This includes changes to:
- Output file formats, names, or structure (sections 2.2–2.4 in codt_tools)
- Namelist parameters or their behavior
- Executable invocation or path resolution (section 2.5 in codt_tools)
- NetCDF variable names, dimensions, attributes, or units
- Fortran modification status (section 5 in codt_tools)

Update the relevant sections in `~/dev/CODT_tools/CLAUDE.md` to match what was actually implemented. Keep both CLAUDE.md files consistent.