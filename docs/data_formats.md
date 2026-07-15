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
| `aerosol_type` | Number of aerosol species |
| `edge` | Number of bins + 1 |
| `bin` | Number of bins |
| `time` | Injection time steps |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `n_ions` | (aerosol_type) | — | Van't Hoff factor |
| `molar_mass` | (aerosol_type) | kg/mol | Solute molar mass |
| `solute_density` | (aerosol_type) | kg/m³ | Solute density |
| `category` | (aerosol_type) | — | Aerosol category index |
| `edge_radii` | (edge) | nm | Bin edge radii |
| `dsd_bin_edges` | (edge) | m | DSD bin edges used for output histograms |
| `cumulative_frequency` | (bin, time) | — | Cumulative size distribution sampled for injection |
| `injection_time` | (time) | s | Injection times |
| `injection_rate` | (time) | 1/s | Injection rate at each injection time |

**Global attributes:** `conventions = "CODT_aerosol_input_v1"`, `aerosol_name`

> In **parcel** mode, particles are pre-loaded at initialization from this distribution (using `aerosol_concentration`), not injected over time; `injection_time`/`injection_rate` apply to chamber mode.

## Parcel input — `CODT_parcel_input_v1` / `v2` / `v3`

Drives adiabatic ascent when `simulation_mode = 'parcel'` (set by `parcel_file`). The velocity profile is piecewise-constant in time.

### v1 (basic)

**Dimensions**

| Dimension | Description |
|-----------|-------------|
| `segment` | Number of velocity segments |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `time` | (segment) | s | Start time of each segment (first must be `0`, monotonically increasing) |
| `velocity` | (segment) | m/s | Vertical velocity for each segment (piecewise-constant) |

**Global attributes:** `conventions = "CODT_parcel_input_v1"`, `initial_pressure` (Pa)

### v2 (with environmental profile)

Required when `do_entrainment = .true.` Adds an environmental sounding that the blob method mixes in.

Adds to v1:

**Dimensions**

| Dimension | Description |
|-----------|-------------|
| `level` | Number of environmental sounding levels |

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `env_pressure` | (level) | Pa | Environmental pressure (monotonically decreasing) |
| `env_temperature` | (level) | K | Environmental temperature |
| `env_RH` | (level) | — | Environmental relative humidity (0–1) |

**Global attributes:** `conventions = "CODT_parcel_input_v2"`

### v3 (with time-varying entrainment schedule)

A superset of v2: keeps the environmental sounding and adds a per-segment entrainment schedule on the **same `segment` time axis as `velocity`**, so the entrainment parameters step in time alongside the ascent velocity. When these variables are present, they **override** the constant `ent_rate`/`n_blob`/`psigma` from the `&ENTRAINMENT` namelist (`random_entrainment` is still taken from the namelist). Lookup is piecewise-constant, like `velocity`.

Adds to v2:

**Variables**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `ent_rate` | (segment) | 1/m | Fractional entrainment rate per segment (> 0) |
| `n_blob` | (segment) | — | Blobs per entrainment event per segment (integer ≥ 1; stored as `int`) |
| `psigma` | (segment) | — | Blob fraction of the domain per segment (0–1, with `psigma * n_blob < 1`) |

**Global attributes:** `conventions = "CODT_parcel_input_v3"`

> Readers accept v1/v2/v3. With `do_entrainment = .true.`, the file must be v2 or v3. In a v3 run, the active `ent_rate`/`n_blob`/`psigma` are also written to the main output file as time series (see below).

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

**Dimensions:** `time` (unlimited), `z`, plus `radius`/`radius_edges` when microphysics is on.

**Always present**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `z` | (z) | m | Vertical coordinate |
| `time` | (time) | s | Output times |
| `T` | (time, z) | K | Temperature |
| `QV` | (time, z) | kg/kg | Water vapor mixing ratio |
| `Tv` | (time, z) | K | Virtual temperature |
| `S` | (time, z) | — | Supersaturation |

**Microphysics variables** (when `do_microphysics`)

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `radius` | (radius) | m | DSD bin centers |
| `radius_edges` | (radius_edges) | m | DSD bin edges |
| `DSD` | (time, radius) | count | Droplet size distribution |
| `DSD_1`, `DSD_2` | (time, radius) | count | DSD per aerosol category |
| `Np` | (time) | count | Total particles |
| `Nact` | (time) | count | Activated droplets |
| `Nun` | (time) | count | Unactivated droplets |
| `Ravg` | (time) | m | Mean radius |
| `LWC` | (time) | kg/m³ | Liquid water content |
| `N_collisions` | (time) | count | Collisions in interval |
| `N_coalescences` | (time) | count | Coalescences in interval |

**Budget variables** (time dim, double, accumulated per write interval)

`budget_inject_solute_mass`, `budget_inject_liquid_mass` (kg); `budget_fallout_liquid_mass`, `budget_fallout_solute_mass` (kg); `budget_condensation` (kg); `budget_dgm_delta_T` (K); `budget_diffusion_delta_T` (K), `budget_diffusion_delta_WV` (kg/kg); `budget_sidewall_delta_T` (K), `budget_sidewall_delta_WV` (kg/kg); `budget_n_injected`, `budget_n_fellout`, `budget_n_coalesced` (counts stored as double).

**Entrainment budget variables** (only when `do_entrainment`)

`budget_detrain_liquid_mass`, `budget_detrain_solute_mass` (kg); `budget_entrain_liquid_mass`, `budget_entrain_solute_mass` (kg); `budget_n_detrained`, `budget_n_entrained` (counts as double).

**Time-varying entrainment series** (only with a v3 parcel input, i.e. a time-varying schedule)

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `ent_rate` | (time) | 1/m | Active entrainment rate at each output step |
| `n_blob` | (time) | — | Active blob count (`int`) |
| `psigma` | (time) | — | Active blob fraction |

With a constant (v2) schedule these are not written; the constant values remain available as the `PARCEL.ent_rate`/`PARCEL.n_blob`/`PARCEL.psigma` global attributes.

**Additions when `simulation_mode = 'parcel'`**

| Variable | Dims | Units | Description |
|----------|------|-------|-------------|
| `parcel_height` | (time) | m | Parcel height |
| `parcel_pressure` | (time) | mb | Parcel pressure |
| `parcel_velocity` | (time) | m/s | Ascent velocity |

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
| `supersaturation` | — | Local supersaturation |
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
   - **Parcel:** `integral_length_scale`, `kolmogorov_length_scale`, `dissipation_rate` (3 values)

**Per-eddy record:** `location` (i4), `length` (i4), `time` (f8)

## Collision binary (`{sim_name}_collisions.bin`)

Unformatted Fortran **stream**, written when `write_collisions = .true.`

**Header:** `N` (i4), `H` (f8), `domain_width` (f8), `volume_scaling` (f8)

**Per-event record:** `id_keep` (i4), `id_kill` (i4), `r_keep` (f8), `r_kill` (f8), `r_after` (f8), `position` (f8), `time` (f8), `flag` (i1)

`flag` is `1` if the event was a coalescence (droplets merged; `r_after` is the merged radius) and `0` if it was a collision without coalescence (`r_after = 0`).

## Log (`{sim_name}.log`)

Plain text — CODT's stdout (all initialization and runtime messages) is redirected here.

## Completion marker (`{sim_name}_DONE`)

Written only on successful completion, with a timestamp. Useful for batch workflows to detect finished runs.

---

## See also

- [Input Parameters](input_parameters.md) — the namelist configuration
- `codt-io` skill — the canonical input/output specification kept in sync with the code
