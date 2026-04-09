# CLAUDE.md

## Habits

BE CONCIOUS OF TOKEN USAGE!

## Build & Run

```bash
source fpm_env          # sets compiler/NetCDF paths (only needed once per shell)
fpm build               # compiles to build/
fpm test                # runs standalone test programs in test/
fpm run -- /path/to/params.nml   # run simulation
```

Tests exercise library routines (paths, file copy, aerosol reader) but do not run full simulations.

## Reference Test (Reftest)

Reftest files live outside the repo — ask the user for directory locations.

```bash
fpm run -- <reftest_input_dir>/params.nml
python3 scripts/compare_reftest_bitwise.py <reftest_output_dir>
```

`identical` = safe. `DIFFERS` = investigate. Uses `same_random = .true.` for determinism.

To promote new reference output (**only when user explicitly requests**):
```bash
cp <reftest_output_dir>/reftest.nc <reftest_output_dir>/reftest_REF.nc
cp <reftest_output_dir>/reftest_particles.nc <reftest_output_dir>/reftest_particles_REF.nc
```

## Architecture

**Simulation modes** (`simulation_mode` namelist):
- `'chamber'` (default) — ODT turbulence, Dirichlet BCs, nondimensional scalars, adaptive dt.
- `'parcel'` — LEM turbulence, periodic BCs, dimensional scalars, fixed dt.

**Turbulence dispatch:** Abstract interfaces in `globals.f90` (`diffuse_iface`, `turbulence_iface`, `sync_iface`). Procedure pointers set at init, called from `main.f90`.

**Fields (chamber):** Dual representation — nondimensional (`T_nd`, `WV_nd`, `Tv_nd`, `W_nd`) for ODT numerics; dimensional (`T`, `WV`, `Tv`, `SS`) for physics. Parcel mode uses dimensional only.

**Time loop** (`app/main.f90`): Each step advances `time` by `dt`:
1. Output (before physics)
2. Diffusion → if fields updated: update_droplets → special_effects → sync
3. Turbulence → if eddy accepted: update_droplets → special_effects → sync

Dim/nondim sync occurs after droplet movement, triggered by diffusion or eddy acceptance. **When modifying the time loop, verify that the diffusion→droplets→sync and turbulence→droplets→sync ordering is preserved.**

**Key modules:**
- `globals.f90` — constants, arrays, triplet map, abstract interfaces, utilities (`nc_verify`, `resolve_path`)
- `ODT.f90` — eddy accept/reject, nondim Crank-Nicolson (Dirichlet), triplet map + addK kernel
- `LEM.f90` — periodic Crank-Nicolson (Sherman-Morrison), -5/3 eddy sampling, periodic triplet map
- `microphysics.f90` — thermodynamic functions, dim/nondim conversion
- `droplets.f90` — `aerosol`→`particle` type hierarchy, Lagrangian tracking
- `DGM.f90` — Cash-Karp RK45 ODE integrator for droplet growth
- `writeout.f90` — buffered NetCDF output
- `write_particle.f90` — particle trajectory NetCDF output
- `initialize.f90` — namelist I/O, domain setup, pointer assignment

**Dependency chain:** `globals` → `microphysics` → `droplets` → `DGM`. ODT/LEM use `globals`, `microphysics`, `droplets`, `writeout`.

## Fortran Conventions

- Precision kinds in `globals.f90`: `dp` (double), `sp` (single), `i4` (32-bit int), `i1`/`i2` for smaller.
- Physical constants are `parameter` values in `globals.f90` — do not duplicate.
- Wrap NetCDF calls with `nc_verify()`.
- `implicit none` everywhere.

## Interface Contract with codt_tools

Executable invocation: `codt <NAMELIST_PATH>`. Relative paths (`aerosol_file`, `bin_data_file`) resolve from namelist's parent directory. `output_directory` must be absolute. No argument → error + `stop 1`. Corresponding spec in `~/dev/CODT_tools/CLAUDE.md`.

## Known Fragilities / Future Cleanup

### Collision-coalescence ↔ droplet fallout interface

`collision_coalescence_step` communicates "this particle fell out" to
`verify_particle_fallout` by setting `lparticles(i)%position = -1.0_dp`
in its writeback loop. `verify_particle_fallout` then detects `position
< 0` and finishes the bookkeeping (sets `%fellout`, compacts the array,
bumps `total_n_fellout`). This works, but is fragile:

- **Order-of-calls coupling.** Nothing between CC's writeback and the
  `verify_particle_fallout` call is allowed to read `%fellout` (still
  `.false.` at that point) or `%position` (sentinel `-1.0` would look
  like a bogus location) without knowing about the contract. Any new
  step inserted in `update_droplets` between CC and
  `verify_particle_fallout` would need to respect or update this.
- **`do_random_fallout` override.** When CC's `handle_fall_event` marks
  a particle dead and the writeback sets `position = -1`,
  `verify_particle_fallout` may instead recycle that particle to the
  top of the domain via `random_fallout` (if `do_random_fallout =
  .true.`). That's the special effect's design, but it means CC's
  decision can be "undone" by a downstream module. CC has no awareness
  of `do_random_fallout`.
- **`total_n_fellout` counter lives in droplets.** CC increments its
  own `collisions_this_step` / `coalescences_this_step`, but the
  fellout counter is bumped by `verify_particle_fallout`. Two
  mechanisms in two modules for closely-related events.

Future cleanup idea: give CC a proper "event → particle state" sink
interface (e.g. a subroutine it calls for each removal/merge) instead
of piggy-backing on `position < 0` sentinels. That would make the
contract explicit and the code easier to read without having to trace
the whole call chain.

## Planned Modifications

### 7. Add predetermined eddies mode

Read eddies from `_eddies.bin` instead of Monte Carlo. New namelist flags `use_predetermined_eddies` + `eddy_file`. Validate header `N` matches. Replace eddy acceptance with sequential reads; implement eddy when simulation time reaches timestamp. Skip `lower_dt`/`raise_dt`. Mutually exclusive with `write_eddies`.

### 8. Add eddy data reader to codt_tools

`eddy_io.py` module: read `_eddies.bin`, return header as dict + eddy records as structured numpy array. Integrate as `CODTSimulation.load_eddies()`.

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `{name}.nc` | netCDF4 | Profiles + time series |
| `{name}_particles.nc` | netCDF4 | Particle data (if `write_trajectories=.true.`) |
| `{name}_eddies.bin` | Binary stream | Eddy events (if `write_eddies=.true.`) |
| `{name}.nml` | ASCII | Namelist copy |

All output to `{output_directory}/{simulation_name}/`.

## Cross-Project Sync

When committing CODT changes, update `~/dev/CODT_tools/CLAUDE.md` if any of these change: output formats/structure, namelist parameters, invocation/path resolution, NetCDF variables/dimensions/attributes, or Fortran modification status.
