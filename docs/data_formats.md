# CODT Data Formats

This document describes the **data files** CODT reads and writes — the NetCDF inputs, the NetCDF outputs, and the raw binary event streams. For the namelist configuration file, see [Input Parameters](input_parameters.md).

Each NetCDF file carries a `conventions` global attribute (e.g. `CODT_output_v1`) that identifies its schema version. Readers should check it; it changes only when the file layout changes in a way that would break an existing reader.

## File overview

| File | Direction | Format | Condition |
|------|-----------|--------|-----------|
| `aerosol_input.nc` (name set by `aerosol_file`) | input | NetCDF4 | when `do_microphysics` |
| `parcel_input.nc` (name set by `parcel_file`) | input | NetCDF4 | `simulation_mode = 'parcel'` |
| Mie absorption table (name set by `mie_data_file`) | input | text | when `do_radiation` (experimental) |
| `{sim_name}.nc` | output | NetCDF4 | always |
| `{sim_name}.log` | output | text | always |
| `{sim_name}_particles.nc` | output | NetCDF4 | `write_trajectories=.true.` |
| `{sim_name}_collisions.bin` | output | binary stream | `write_collisions=.true.` |
| `{sim_name}_eddies.bin` | output | binary stream | `write_eddies=.true.` |
| `{sim_name}_DONE` | output | text | on successful completion |

`{sim_name}` is `simulation_name`; all outputs go to `output_directory/`. Input paths given relative to the namelist are resolved from the namelist's parent directory.

---

# Input files

## Aerosol input — `CODT_aerosol_input_v1`

The aerosol size distribution and solute properties. Required when `do_microphysics = .true.`

**Dimensions**

| Dimension | Description |
|-----------|-------------|
| `aerosol_type` | Number of aerosol species (rows of the composition table) |
| `edge` | Number of background bins + 1 |
| `bin` | Number of background bins |
| `dsd_edge` | Number of DSD output-histogram bins + 1 (independent of `edge`) |
| `time` | Injection time steps |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `n_ions` | (aerosol_type) | — | Van't Hoff factor |
| `molar_mass` | (aerosol_type) | kg/mol | Solute molar mass |
| `solute_density` | (aerosol_type) | kg/m³ | Solute density |
| `category` | (bin) | — | Aerosol category index (output label) |
| `bin_type` | (bin) | — | *Optional.* Composition row each bin is made of. Absent ⇒ all bins are type 1 |
| `edge_radii` | (edge) | nm | Bin edge radii |
| `dsd_bin_edges` | (dsd_edge) | µm | DSD output-histogram bin edges. Uses its own `dsd_edge` dim (output resolution), independent of the aerosol `edge` dim |
| `cumulative_frequency` | (time, bin) | — | Cumulative size distribution sampled for injection |
| `injection_time` | (time) | s | Injection times |
| `injection_rate` | (time) | m⁻³ s⁻¹ | Injection rate (per-volume per-time) at each injection time |

**Global attributes:** `conventions = "CODT_aerosol_input_v1"`, `aerosol_name`

> In **parcel** mode, particles are pre-loaded at initialization from this distribution (using `aerosol_concentration`), not injected over time; `injection_time`/`injection_rate` apply to chamber mode.

### The three per-bin labels

A **bin** is the atom of the size distribution: one sampleable kind of dry aerosol. Sampling draws a bin, and two labels come along with it. They are easy to conflate, so:

| Label | What it is | What it affects |
|-------|-----------|-----------------|
| `bin` | A dry radius that can be drawn | The size that gets sampled |
| `category(bin)` | An **output label only** | Which per-category DSD (`DSD_1`, `DSD_2`, …) the particle is counted in, and its `aerosol_category` tag in trajectory output |
| `bin_type(bin)` | A row of the composition table | The particle's solute physics (`n_ions`, `molar_mass`, `solute_density`), and whether `seed_hydration` applies |

The CDF is over **bins only** — neither `category` nor `bin_type` subdivides it. So you can categorize a background aerosol by size (categories 1–6, say) while every bin shares one composition row, or carry several background materials (types 1–3) under a single category. The two labels are independent.

## Aerosol seeding — the seed group

Seeding introduces a **second aerosol population** with its own bins, its own size distribution, and its own release schedule, kept separate from the background so the file states each explicitly. The group is **optional and all-or-nothing**: include every variable below, or none.

`&MICROPHYSICS do_seeding` is the sole controller. With `do_seeding = .true.` the file **must** contain a seed group (its absence is a fatal error). With `do_seeding = .false.` the seed group is **ignored** — never read, and with no effect on the composition table, DSD categories, or any output — so a file may carry a dormant seed group and still run unseeded; CODT prints a warning ("seeding input detected but do_seeding is false; the seed group is ignored") when it detects one so the unused data is not a surprise. This is what lets one aerosol file serve both a seeded and an unseeded run. (The bundled `input/aerosol_input.nc` ships **with** a seed group for exactly this reason.)

**Dimensions**

| Dimension | Description |
|-----------|-------------|
| `seed_bin` | Number of seed bins |
| `seed_edge` | `seed_bin` + 1 |
| `seed_event` | Number of seeding events |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `seed_edge_radii` | (seed_edge) | nm | Seed bin edge radii |
| `seed_category` | (seed_bin) | — | Output category for each seed bin |
| `seed_bin_type` | (seed_bin) | — | Composition row each seed bin is made of. **Required** |
| `seed_frequency` | (seed_event, seed_bin) | — | Cumulative size distribution, per event. Each event's row runs to 1.0 |
| `seed_coord` | (seed_event) | s or m or Pa | Point on the schedule axis that triggers the event |
| `seed_concentration` | (seed_event) | cm⁻³ | Concentration released by the event |

### Seeding events

An event says *"release this concentration, with this size distribution, when the run reaches this point"*. It fires **once**, the first time the run reaches `seed_coord`, and never again.

`seed_coord` is interpreted per mode:

| Mode | Axis | Meaning |
|------|------|---------|
| chamber | time [s] | "Seed at t = 300 s" |
| parcel, `vertical_axis = 'height'` | height [m] | "Seed at z = 600 m" |
| parcel, `vertical_axis = 'pressure'` | pressure [Pa] | "Seed at p = 85000 Pa" |

Events need **no ordering**, and a parcel that re-crosses a seeded level does **not** seed again — the release is a discrete burst, not an ambient concentration the parcel keeps sweeping up. Two events at the same `seed_coord` are rejected; write one larger event instead. Each event has its own `seed_frequency` row, so successive releases may differ in size distribution.

> **A chemically identical seed still needs its own composition row.** What marks material as seed is that `seed_bin_type` points at it — so seeding NaCl into an NaCl background means duplicating the NaCl row as type 2 and pointing `seed_bin_type` there. That duplicate is what earns the seed its own `seed_hydration` treatment and its own DSD category. A type referenced by both `bin_type` and `seed_bin_type` is rejected, since `seed_hydration` would have no answer for it.

### Hydration at release

`&MICROPHYSICS seed_hydration` sets the wet radius seed particles are born with. It applies **only** to seed material; background particles keep the `initial_wet_radius` multiple of their dry radius.

| Value | Wet radius at release |
|-------|----------------------|
| `equilibrium` (default) | Köhler equilibrium at the local RH (capped at 0.99 to stay on the stable branch), bounded by what the droplet could actually grow to in `seed_growth_time` seconds |
| `double_growth` | Twice the dry radius, regardless of humidity |
| `dry` | The bare dry radius; the growth model wets it from there |

The growth bound on `equilibrium` exists for GCCN: their equilibrium radius is tens of microns, and starting them there would condense water they would really need minutes to collect.

### What a seed becomes when droplets coalesce

With `do_coalescence`, two droplets merge into one, and the survivor takes a single category and material. The rules:

- **Both from the same population:** the **larger** droplet's category and material survive. (Collisions only occur when the faster-falling droplet catches a slower one, and fall speed grows with radius, so the survivor is the larger one.)
- **A seed and a background droplet:** the **seed** survives, regardless of which was larger. A seed collected by a big background droplet still counts as seed, so seeded runs don't lose track of their seed material through collection.

> **Caveat — mass is conserved, dissolved ions are not.** A particle carries one material, so merging droplets of different composition adds the solute masses but treats the total as the survivor's material: its density sets the merged dry radius, and its `n_ions`/`molar_mass` set the Köhler curve. Seeding a material chemically unlike your background makes this approximation bite on every collection event. Seeding a chemical twin of the background (the duplicate-composition-row recipe above) avoids it entirely, since both rows describe the same substance.

## Parcel input — `CODT_parcel_input_v3`

Drives the parcel trajectory when `simulation_mode = 'parcel'` (set by `parcel_file`). The trajectory is an **ordered sequence of waypoint legs**: each leg says "proceed to this level at this signed velocity," and the active leg advances when its target is reached. Because the lookup key is the leg counter (not the current position), trajectories may revisit levels — e.g. ascend to 4 km at 2 m/s, descend to 2 km at 1 m/s, ascend again. **Completing the final leg ends the simulation** (like `pressure_limit`); the `_DONE` marker is written normally.

> **Removed formats:** the time-based `CODT_parcel_input_v1`/`v2` schemas are no longer read (last supported by CODT 1.x). Regenerate old inputs in the v3 waypoint format.

**Dimensions**

| Dimension | Description |
|-----------|-------------|
| `segment` | Number of trajectory legs (≥ 1) |
| `level` | Sounding levels (only when the sounding is present) |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `segment_coord` | (segment) | m or Pa | Leg **target** levels, per `&PARCEL vertical_axis` (`'height'` or `'pressure'`). Consecutive targets must differ; no monotonicity requirement. |
| `velocity` | (segment) | m/s | Signed velocity of each leg (nonzero). Must point from the previous level (the launch level for leg 1: `initial_height`, or the initial `pres` on the pressure axis) toward the leg's target — validated at read. |
| `env_height` | (level) | m | *Sounding:* height (monotonically increasing) |
| `env_pressure` | (level) | Pa | *Sounding:* pressure (monotonically decreasing) |
| `env_temperature` | (level) | K | *Sounding:* temperature |
| `env_RH` | (level) | — | *Sounding:* relative humidity (0–1) |
| `ent_rate` | (segment) | **1/km** | *Optional:* fractional entrainment rate per leg (> 0) |
| `n_blob` | (segment) | — | *Optional:* blobs per event per leg (integer ≥ 1; stored as `int`) |
| `psigma` | (segment) | — | *Optional:* **total** domain fraction replaced per event, per leg (0–1) |

**The sounding is optional**: it is required when `do_entrainment = .true.` or `pressure_mode = 'environment'`; a bare adiabatic hydrostatic parcel runs without it (the `parcel_height_env` diagnostic is then skipped). The four `env_*` variables come as a set.

The three entrainment variables are present together or not at all. When present, they **override** the constant `ent_rate`/`n_blob`/`psigma` from the `&ENTRAINMENT` namelist per leg (`random_entrainment` is still taken from the namelist) — so entrainment during a descent leg can differ from the ascent leg through the same heights. When absent, the namelist constants apply.

**Global attributes:** `conventions = "CODT_parcel_input_v3"`

> Entrainment event *timing* stays in the time domain (interval from `ent_rate`, blob geometry, and |velocity|); only the parameter values are keyed to the active leg. In a run with a per-leg schedule, the active `ent_rate`/`n_blob`/`psigma` are also written to the main output file as time series (see below).

## Mie absorption table (radiation)

> ⚠️ **Experimental.** The radiation code (and this input) is under active development; the format may change.

Plain-text table of Mie absorption data, read by the radiation module when `do_radiation = .true.` (path set by `mie_data_file`). Treat as preliminary.

---

# Output files

## Main output — `CODT_output_v1` (`{sim_name}.nc`)

Profiles and time series. Always written.

**Global attributes**
- `conventions` — schema string (`CODT_output_v1`)
- `code_version` — version from `git describe` at build time (e.g. `v1.0.0`, or `v1.0.0-5-g1a2b3c4-dirty` for dev builds)
- `git_commit` — short commit hash (`-dirty` if the tree had uncommitted changes)
- The namelist parameters, namespaced as `PARAMETERS.N`, `MICROPHYSICS.write_trajectories`, etc. (logicals stored as 0/1; mode-specific groups only present for the relevant mode)
- **Parcel only, derived (not namelist inputs):** `LEM.actual_kolmogorov_scale`,
  `LEM.grid_eddy_scale`, `LEM.diffusivity_length_scale`, `LEM.smallest_eddy_scale`
  (all m), `LEM.smallest_eddy_gridpoints` (int, cells), and
  `LEM.diffusivity_enhancement` (unitless, ≥ 1). Added in v3.0.0; additive, so
  `CODT_output_v1` is unchanged — detect by attribute presence.
  `TURBULENCE_LEM.kolmogorov_length_scale` was **removed** in the same release
  along with the namelist parameter it echoed. `smallest_eddy_scale` is the
  governing scale; `actual_kolmogorov_scale` is reported for reference and never
  governs. Consistency check: `smallest_eddy_scale == smallest_eddy_gridpoints · H/N`.

**Dimensions:** `time` (unlimited), `z`, plus `radius`/`radius_edges` when microphysics is on.

**Always present**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `z` | (z) | m | Vertical coordinate |
| `time` | (time) | s | Output times |
| `T` | (time, z) | °C | Temperature |
| `QV` | (time, z) | kg/kg | Water vapor mixing ratio |
| `Tv` | (time, z) | °C | Virtual temperature |
| `S` | (time, z) | % | Supersaturation |

> **Output temperatures are °C, not K.** The `units` attributes on the file
> (`"celsius"`) are authoritative. Note the asymmetry with the namelist echo: the
> `PARAMETERS.Tref` global attribute is stored in **K** even though `Tref` is
> supplied in °C. Lengths are likewise mixed — the DSD radii are µm while the
> per-particle `solute_radius` is m. Read the `units` attribute rather than
> assuming SI.

**Microphysics variables** (when `do_microphysics`)

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `radius` | (radius) | µm | DSD bin centers |
| `radius_edges` | (radius_edges) | µm | DSD bin edges |
| `DSD` | (time, radius) | count | Droplet size distribution |
| `DSD_1`, `DSD_2` | (time, radius) | count | DSD per aerosol category |
| `Np` | (time) | count | Total particles |
| `Nact` | (time) | count | Activated droplets |
| `Nun` | (time) | count | Unactivated droplets |
| `Ravg` | (time) | µm | Mean radius (wet) |
| `LWC` | (time) | g/m³ | Liquid water content |
| `N_collisions` | (time) | count | Collisions in interval |
| `N_coalescences` | (time) | count | Coalescences in interval |

**Budget variables** (time dim, double, accumulated per write interval)

`budget_inject_solute_mass`, `budget_inject_liquid_mass` (kg); `budget_fallout_liquid_mass`, `budget_fallout_solute_mass` (kg); `budget_condensation` (kg); `budget_dgm_delta_T` (K); `budget_diffusion_delta_T` (K), `budget_diffusion_delta_WV` (kg/kg); `budget_sidewall_delta_T` (K), `budget_sidewall_delta_WV` (kg/kg); `budget_n_injected`, `budget_n_fellout`, `budget_n_coalesced` (counts stored as double).

**Entrainment budget variables** (only when `do_entrainment`)

`budget_detrain_liquid_mass`, `budget_detrain_solute_mass` (kg); `budget_entrain_liquid_mass`, `budget_entrain_solute_mass` (kg); `budget_n_detrained`, `budget_n_entrained` (counts as double).

**Radiation variables** (only when `do_radiation`; both modes)

> ⚠️ **Experimental**, like the radiation code itself.

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `rad_F_net` | (time, z) | W/m² | Net radiative flux |
| `rad_heating_rate` | (time, z) | K/s | Radiative heating rate |
| `budget_radiation_delta_T` | (time) | K | Domain-sum T change from radiation (double) |

**Varying entrainment series** (only with a parcel input carrying a per-leg entrainment schedule)

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `ent_rate` | (time) | 1/km | Active entrainment rate at each output step |
| `n_blob` | (time) | — | Active blob count (`int`) |
| `psigma` | (time) | — | Active total entrained fraction |

With a constant schedule these are not written; the constant values remain available as the `PARCEL.ent_rate` (1/km)/`PARCEL.n_blob`/`PARCEL.psigma` global attributes.

**Additions when `simulation_mode = 'parcel'`**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `parcel_height` | (time) | m | Parcel height (integrated from velocity) |
| `parcel_pressure` | (time) | mb | Parcel pressure |
| `parcel_velocity` | (time) | m/s | Ascent velocity |
| `parcel_height_env` | (time) | m | Environment height at the parcel's pressure (`pressure_mode = 'hydrostatic'` with a sounding only). Its drift from `parcel_height` measures how far the self-integrated pressure has left the sounding's p(z). |

The `PARCEL.pressure_mode`, `PARCEL.vertical_axis`, and `PARCEL.initial_height` global attributes record the trajectory configuration.

## Particle output — `CODT_particle_output_v1` (`{sim_name}_particles.nc`)

Per-particle trajectories, written when `write_trajectories = .true.` (over the `trajectory_start`–`trajectory_end` window at `trajectory_timer` cadence). Stored as a **CF contiguous ragged array**, because the number of live particles changes from step to step (injection, fallout, coalescence).

**Global attributes:** `conventions`, `code_version`, `git_commit`

**Dimensions:** `record` (unlimited), `time_step` (unlimited)

**Per-time-step variables** (indexed by `time_step`)

| Variable | Units | Description |
|----------|-------|-------------|
| `time` | s | Time of each output step |
| `row_sizes` | — | Number of particle records belonging to each step (`cf_role = "ragged_row_sizes"`) |

**Per-record variables** (indexed by `record` — one entry per particle per step)

| Variable | Units | Description |
|----------|-------|-------------|
| `particle_id` | — | Unique particle identifier (stable across steps) |
| `aerosol_id` | — | Solute/aerosol species id |
| `gridcell` | — | Grid cell index the particle occupies |
| `position` | m | Vertical position |
| `temperature` | °C | Local temperature |
| `water_vapor` | g/kg | Local water-vapor mixing ratio |
| `supersaturation` | % | Local supersaturation |
| `radius` | µm | Droplet radius |
| `solute_radius` | m | Dry solute radius |
| `activated` | — | 1 if activated, else 0 |
| `aerosol_category` | — | Aerosol category index |
| `n_collisions` | — | Collisions so far (only when `do_collisions`) |
| `n_coalescences` | — | Coalescences so far (only when `do_collisions`) |
| `radius_before_coalescence` | µm | Radius prior to last coalescence (only when `do_collisions`) |

### How the ragged array is laid out

Every live particle is written at every output step, all concatenated along the single `record` dimension **grouped by time step and in time order**. The `time_step` dimension holds one entry per step:

- `time[k]` — the time of step `k`
- `row_sizes[k]` — how many particle records belong to step `k`

Because the blocks are contiguous and ordered, the records for step `k` occupy a single slice. With 0-based indexing, define the offsets as the cumulative sum of `row_sizes`:

```
offsets = [0, row_sizes[0], row_sizes[0]+row_sizes[1], ...]   # length time_step+1
step k  -> records[offsets[k] : offsets[k+1]]                 # length row_sizes[k]
```

Within a block, records follow CODT's internal particle ordering for that step. **There is no fixed slot per particle across steps** — a given `particle_id` can sit at a different offset in each block, or be absent (a coalesced particle's id disappears; the survivor keeps its id; injected particles appear partway through). So you always locate a particle by matching `particle_id`, never by a constant index.

### Recipe 1 — all particles at one time step (fast, contiguous slice)

```python
import numpy as np
from netCDF4 import Dataset

ds = Dataset("sim_particles.nc")
row_sizes = ds["row_sizes"][:]
offsets = np.concatenate([[0], np.cumsum(row_sizes)])

k = 10                                   # the step you want (time = ds["time"][k])
s, e = offsets[k], offsets[k + 1]
radius_k = ds["radius"][s:e]             # every particle's radius at step k
id_k     = ds["particle_id"][s:e]
```

### Recipe 2 — one particle's trajectory (gather by id across all steps)

The clean vectorized way: expand `time` to one value per record (each step's time repeated `row_sizes[k]` times), then mask by `particle_id`.

```python
ids_all   = ds["particle_id"][:]
radius_all = ds["radius"][:]
time_per_record = np.repeat(ds["time"][:], row_sizes)   # aligns 1:1 with records

target_id = 42
mask  = ids_all == target_id
times = time_per_record[mask]            # the steps where this particle was alive
radii = radius_all[mask]                 # its radius at each of those times
```

`times`/`radii` come out in time order and skip any steps where the particle didn't exist.

> The **CODT_tools** Python package provides readers that do this indexing for you (per-step slices and per-particle trajectories), so you usually don't have to handle the offsets by hand.

## Eddy binary (`{sim_name}_eddies.bin`)

Unformatted Fortran **stream** (no record markers). Mode-aware header, then one record per accepted eddy. Byte sizes: `i1` = 1, `i4` = 4, `f8` = 8.

**Header**
1. `mode_flag` (i1) — `0` = chamber, `1` = parcel
2. `N` (i4), `H` (f8)
3. Mode-specific `f8` fields:
   - **Chamber:** `C2`, `ZC2`, `Tdiff`, `Tref` (4 values)
   - **Parcel:** `integral_length_scale`, `smallest_eddy_scale`, `dissipation_rate` (3 values)

   > **Changed in v3.0.0:** the second parcel value was
   > `kolmogorov_length_scale` (a namelist input). It is now
   > `smallest_eddy_scale`, a *derived* quantity —
   > `max((ν³/ε)^(1/4), 6·dz)`. Field count and types are unchanged, so
   > readers will not fail; they will silently read a different quantity.
   > Distinguish by the file's `code_version`.

**Per-eddy record:** `location` (i4), `length` (i4), `time` (f8)

## Collision binary (`{sim_name}_collisions.bin`)

Unformatted Fortran **stream**, written when `write_collisions = .true.`

**Header:** `N` (i4), `H` (f8), `domain_width` (f8), `volume_scaling` (f8)

**Per-event record:** `id_keep` (i4), `id_kill` (i4), `r_keep` (f8), `r_kill` (f8), `r_after` (f8), `position` (f8), `time` (f8), `flag` (i1)

`flag` is `1` if the event was a coalescence (droplets merged; `r_after` is the merged radius) and `0` if it was a collision without coalescence (`r_after = 0`).

`time` is **absolute simulation time in seconds**, on the same axis as the `time` coordinate of the main netCDF.

> **Files written by CODT 3.0.0 and earlier are different.** In those versions this field held the event's time *within* the current collision-coalescence window (0 → `delta_time`), not absolute time, so values were small (order 1e-8 – 1e-2 s) and reset every window, making the stream non-monotonic. The record layout is **unchanged** — the field was already `f8`, so only its meaning differs and old readers still parse new files correctly. To tell them apart, check the `git_commit` / `code_version` global attribute on the run's netCDF, or simply test whether `time` in the binary ever exceeds one `delta_time`. For recovering absolute time from a pre-3.0.1 file, see `docs/known_issues.md`.

## Log (`{sim_name}.log`)

Plain text — CODT's stdout (all initialization and runtime messages) is redirected here.

## Completion marker (`{sim_name}_DONE`)

Written only on successful completion, with a timestamp. Useful for batch workflows to detect finished runs.

---

## See also

- [Input Parameters](input_parameters.md) — the namelist configuration
- `codt-io` skill — the canonical input/output specification kept in sync with the code
