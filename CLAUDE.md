# CLAUDE.md

## Habits

BE CONSCIOUS OF TOKEN USAGE!

Before every commit, invoke the `codt-versioning` skill to determine whether a version bump is needed and at what level.

## Build & Run

```bash
source fpm_env          # sets compiler/NetCDF paths (once per shell)
./build.sh              # injects version+git hash into writeout.f90, then runs fpm build
fpm build               # plain build (no version injection — use build.sh for releases)
fpm test                # unit tests (paths, aerosol reader) — no full simulations
fpm run -- /path/to/params.nml   # run simulation
```

## Versioning

Semantic versioning in `fpm.toml` (MAJOR.MINOR.PATCH):
- **PATCH** — bug fixes, perf, internal refactors. Same output for same input.
- **MINOR** — new features, namelist params, output variables. Old inputs/tools still work.
- **MAJOR** — breaking changes. Conventions string bump, incompatible namelist or output format.

Bump the version in `fpm.toml` when merging to main. During development, the git commit hash (`-dirty` suffix for uncommitted changes) in output files provides traceability.

`build.sh` injects `code_version` and `git_commit` into `writeout.f90` at build time, written as global attributes in both profile and particle NetCDF files. The injected values stay in `writeout.f90` after the build (not restored to placeholders), so `fpm run` uses the same binary without recompiling. The `conventions` string (`CODT_output_v1`, etc.) is the format contract — bump it only when the file schema changes in a way that would break an existing reader.

## Reference Test (Reftest)

Reftest files live outside the repo — ask the user for directory locations. Use the `run-reftest` and `archive-reftest` skills for execution and promotion.

`identical` = safe. `DIFFERS` = investigate. Uses `same_random = .true.` for determinism.

## Architecture

**Simulation modes** (`simulation_mode` namelist):
- `'chamber'` (default) — ODT turbulence, Dirichlet BCs, nondimensional scalars, adaptive dt.
- `'parcel'` — LEM turbulence, periodic BCs, dimensional scalars, fixed dt. Adiabatic ascent driven by `parcel_file` (piecewise-constant velocity). Particles pre-loaded from aerosol distribution at init (not injected over time), equilibrated to Köhler. Positions wrap periodically (no fallout). Optional entrainment via blob method (`do_entrainment`). `pressure_limit` stops simulation at a target pressure.

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
| **Special effects** | Sidewalls, stochastic fallout (chamber only in output attrs) | Not used |

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
- `ode_integrators.f90` — Cash-Karp RK4(5) explicit and ROS3 Rosenbrock implicit adaptive ODE integrators
- `DGM.f90` — droplet growth model RHS, Jacobian, and Rosenbrock ODE driver
- `special_effects.f90` — sidewall nudging, stochastic fallout
- `writeout.f90` — buffered NetCDF output
- `write_particle.f90` — particle trajectory NetCDF output
- `parcel.f90` — reads piecewise-constant velocity from NetCDF, applies adiabatic forcing (dT, dp)
- `initialize.f90` — namelist I/O, output directory setup, stdout redirect, domain setup, pointer assignment

**Dependency chain:** `globals` → `microphysics` → `particle_types` → `droplets` → `DGM`. `ode_integrators` → `DGM`. `collection_efficiency` → `collision_coalescence` → `droplets`. ODT/LEM use `globals`, `microphysics`, `droplets`, `writeout`.

## Fortran Conventions

- Precision kinds in `globals.f90`: `dp` (double), `sp` (single), `i4` (32-bit int), `i1`/`i2` for smaller.
- Physical constants are `parameter` values in `globals.f90` — do not duplicate.
- Wrap NetCDF calls with `nc_verify()`.
- `implicit none` everywhere.

## Interface Contract with codt_tools

Executable invocation: `codt <NAMELIST_PATH>`. Relative path `aerosol_file` resolves from namelist's parent directory. `output_directory` can be absolute or relative (resolved from cwd); parent directory must exist. No argument → error + `stop 1`. Corresponding spec in `~/dev/CODT_tools/CLAUDE.md`.

## Known Fragilities

**CC ↔ fallout interface:** `collision_coalescence_step` signals fallout by setting `position = -1.0`, then `verify_particle_fallout` detects this and does bookkeeping (`%fellout`, array compaction, `total_n_fellout`). Fragile because: (1) nothing between CC writeback and `verify_particle_fallout` may read `%position` or `%fellout`, (2) `do_random_fallout` can silently recycle CC-removed particles, (3) fellout counting is split across two modules. Future fix: explicit event interface instead of sentinel values.

**Sidewall Ra in parcel mode:** `initialize_special_effects` receives a Rayleigh number from the caller. In chamber mode this is the true Ra = gΔTH³/(T_ref·ν·κ). In parcel mode the LEM Reynolds number is substituted as a placeholder — this needs a proper formulation.

## Planned Modifications

- **Predetermined eddies mode:** Read eddies from `_eddies.bin` instead of Monte Carlo. New namelist flags `use_predetermined_eddies` + `eddy_file`. Mutually exclusive with `write_eddies`.
- **Eddy data reader (codt_tools):** `eddy_io.py` module to read `_eddies.bin`. Integrate as `CODTSimulation.load_eddies()`.

## Output Files

All output to `{output_directory}/`, prefixed with `{simulation_name}`:

| File | Format | Description |
|------|--------|-------------|
| `{name}.nc` | netCDF4 | Profiles + time series (`CODT_output_v1`). Parcel mode adds `parcel_height`, `parcel_pressure`, `parcel_velocity` time series. |
| `{name}.log` | ASCII | Redirected stdout (all initialization and runtime messages) |
| `{name}_particles.nc` | netCDF4 | Particle data (`CODT_particle_output_v1`, if `write_trajectories=.true.`) |
| `{name}_collisions.bin` | Binary stream | Collision/coalescence events (if `write_collisions=.true.`) |
| `{name}_eddies.bin` | Binary stream | Eddy events (if `write_eddies=.true.`) |
| `{name}_DONE` | ASCII | Completion marker with timestamp |

## Cross-Project Sync

When committing CODT changes, update `~/dev/CODT_tools/CLAUDE.md` if any of these change: output formats/structure, namelist parameters, invocation/path resolution, NetCDF variables/dimensions/attributes.
