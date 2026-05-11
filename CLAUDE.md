# CLAUDE.md

## Habits

BE CONSCIOUS OF TOKEN USAGE!

## Build & Run

```bash
source fpm_env          # sets compiler/NetCDF paths (once per shell)
fpm build               # compiles to build/
fpm test                # unit tests (paths, file copy, aerosol reader) — no full simulations
fpm run -- /path/to/params.nml   # run simulation
```

## Reference Test (Reftest)

Reftest files live outside the repo — ask the user for directory locations. Use the `run-reftest` and `archive-reftest` skills for execution and promotion.

`identical` = safe. `DIFFERS` = investigate. Uses `same_random = .true.` for determinism.

## Architecture

**Simulation modes** (`simulation_mode` namelist):
- `'chamber'` (default) — ODT turbulence, Dirichlet BCs, nondimensional scalars, adaptive dt.
- `'parcel'` — LEM turbulence, periodic BCs, dimensional scalars, fixed dt. Adiabatic ascent driven by `parcel_file` (piecewise-constant velocity). Particles pre-loaded from aerosol distribution at init (not injected over time), equilibrated to Köhler. Positions wrap periodically (no fallout).

### Chamber vs. Parcel Mode

| Aspect | Chamber | Parcel |
|--------|---------|--------|
| **Turbulence** | ODT (`&TURBULENCE_ODT`) | LEM (`&TURBULENCE_LEM`) |
| **Boundary conditions** | Dirichlet (fixed T at top/bottom) | Periodic |
| **Scalar representation** | Dual: nondimensional + dimensional | Dimensional only |
| **Time step** | Adaptive (eddy-driven) | Fixed |
| **Forcing** | Temperature gradient (`Tdiff`) | Adiabatic ascent (`parcel_file`) |
| **Particle initialization** | Injected over time (`injection_time`, `injection_rate` from aerosol NC) | Pre-loaded at init (`aerosol_concentration`) + Köhler equilibration |
| **Particle fallout** | Gravitational removal (`verify_particle_fallout`) | Periodic wrapping (`modulo(position, H)`) |
| **`Tref`** | Bottom boundary temperature (K, converted from °C) | Uniform initial temperature (K, converted from °C) |
| **`pres`** | Constant reference pressure | Initial pressure, evolves hydrostatically |
| **`initial_RH`** | Ignored | Initial relative humidity (0–1), sets WV field |
| **`aerosol_concentration`** | Ignored | Number concentration (cm⁻³) |
| **`parcel_file`** | Not used | NetCDF with velocity segments for ascent rate |
| **`Tdiff`** | Top-bottom ΔT driving convection | Not used |
| **Special effects** | Sidewalls, stochastic fallout | Sidewalls use Re placeholder (needs work) |

**Turbulence dispatch:** Abstract interfaces in `globals.f90` (`diffuse_iface`, `turbulence_iface`, `sync_iface`). Procedure pointers set at init, called from `main.f90`.

**Fields (chamber):** Dual representation — nondimensional (`T_nd`, `WV_nd`, `Tv_nd`, `W_nd`) for ODT numerics; dimensional (`T`, `WV`, `Tv`, `SS`) for physics. Parcel mode uses dimensional only.

**Time loop** (`app/main.f90`): Each step advances `time` by `dt`:
1. Output (before physics)
2. Diffusion → if fields updated: update_droplets → special_effects → sync
3. Turbulence → if eddy accepted: update_droplets → special_effects → sync

**When modifying the time loop, preserve the diffusion→droplets→sync and turbulence→droplets→sync ordering.**

**Key modules:**
- `globals.f90` — constants, arrays, triplet map, abstract interfaces, utilities (`nc_verify`, `resolve_path`)
- `ODT.f90` — eddy accept/reject, nondim Crank-Nicolson (Dirichlet), triplet map + addK kernel
- `LEM.f90` — periodic Crank-Nicolson (Sherman-Morrison), -5/3 eddy sampling, periodic triplet map
- `microphysics.f90` — thermodynamic functions, dim/nondim conversion
- `particle_types.f90` — `aerosol`→`particle` type hierarchy, Köhler theory, terminal velocity
- `droplets.f90` — Lagrangian tracking, DSD binning, sequential/batched DGM dispatch
- `collision_coalescence.f90` — event-driven 1D collision-coalescence (min-heap, linked list)
- `collection_efficiency.f90` — coalescence kernels: Hall (1980), Long (1974), unity
- `ode_integrators.f90` — Cash-Karp RK45 adaptive ODE integrator
- `DGM.f90` — droplet growth model RHS and ODE driver
- `special_effects.f90` — sidewall nudging, stochastic fallout
- `writeout.f90` — buffered NetCDF output
- `write_particle.f90` — particle trajectory NetCDF output
- `parcel.f90` — reads piecewise-constant velocity from NetCDF, applies adiabatic forcing (dT, dp)
- `initialize.f90` — namelist I/O, domain setup, pointer assignment

**Dependency chain:** `globals` → `microphysics` → `particle_types` → `droplets` → `DGM`. `ode_integrators` → `DGM`. `collection_efficiency` → `collision_coalescence` → `droplets`. ODT/LEM use `globals`, `microphysics`, `droplets`, `writeout`.

## Fortran Conventions

- Precision kinds in `globals.f90`: `dp` (double), `sp` (single), `i4` (32-bit int), `i1`/`i2` for smaller.
- Physical constants are `parameter` values in `globals.f90` — do not duplicate.
- Wrap NetCDF calls with `nc_verify()`.
- `implicit none` everywhere.

## Interface Contract with codt_tools

Executable invocation: `codt <NAMELIST_PATH>`. Relative path `aerosol_file` resolves from namelist's parent directory. `output_directory` must be absolute. No argument → error + `stop 1`. Corresponding spec in `~/dev/CODT_tools/CLAUDE.md`.

## Known Fragilities

**CC ↔ fallout interface:** `collision_coalescence_step` signals fallout by setting `position = -1.0`, then `verify_particle_fallout` detects this and does bookkeeping (`%fellout`, array compaction, `total_n_fellout`). Fragile because: (1) nothing between CC writeback and `verify_particle_fallout` may read `%position` or `%fellout`, (2) `do_random_fallout` can silently recycle CC-removed particles, (3) fellout counting is split across two modules. Future fix: explicit event interface instead of sentinel values.

**`copy_file` self-clobber:** When input dir == output dir, `copy_file` truncates the file to 0 bytes. Workaround: keep inputs in a subdirectory.

**DGM at low RH:** The RK45 solver hangs when droplets encounter strongly negative supersaturation (e.g. RH<50%). The 1/r amplification near `r_floor` creates stiffness that shrinks dt to near-zero. Parcel runs should start at RH≥0.9 until this is fixed.

**Sidewall Ra in parcel mode:** `initialize_special_effects` receives a Rayleigh number from the caller. In chamber mode this is the true Ra = gΔTH³/(T_ref·ν·κ). In parcel mode the LEM Reynolds number is substituted as a placeholder — this needs a proper formulation.

## Planned Modifications

- **Predetermined eddies mode:** Read eddies from `_eddies.bin` instead of Monte Carlo. New namelist flags `use_predetermined_eddies` + `eddy_file`. Mutually exclusive with `write_eddies`.
- **Eddy data reader (codt_tools):** `eddy_io.py` module to read `_eddies.bin`. Integrate as `CODTSimulation.load_eddies()`.

## Output Files

All output to `{output_directory}/{simulation_name}/`:

| File | Format | Description |
|------|--------|-------------|
| `{name}.nc` | netCDF4 | Profiles + time series |
| `{name}_particles.nc` | netCDF4 | Particle data (if `write_trajectories=.true.`) |
| `{name}_collisions.bin` | Binary stream | Collision/coalescence events (if `write_collisions=.true.`) |
| `{name}_eddies.bin` | Binary stream | Eddy events (if `write_eddies=.true.`) |
| `{name}.nml` | ASCII | Namelist copy |

## Cross-Project Sync

When committing CODT changes, update `~/dev/CODT_tools/CLAUDE.md` if any of these change: output formats/structure, namelist parameters, invocation/path resolution, NetCDF variables/dimensions/attributes.
