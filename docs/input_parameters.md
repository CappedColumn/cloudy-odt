# CODT Input Parameters

CODT is configured entirely through a single Fortran **namelist file**, whose path is passed as the only command-line argument:

```bash
codt path/to/params.nml
```

The namelist is divided into groups (`&PARAMETERS`, `&MICROPHYSICS`, etc.). Each group is opened with `&NAME` and closed with `/`. Only the groups relevant to your run need to be present — see [Which groups do I need?](#which-groups-do-i-need) below. A complete example lives at [`input/params.nml`](../input/params.nml).

> **Defaults vs. the template.** The tables below list the *code* default — the value used if the parameter is omitted from the namelist. The shipped `input/params.nml` sets example values that may differ. Parameters marked **— (required)** have no default and must be supplied.

> **Units note.** A `—` in the Units column means the parameter is dimensionless or non-numeric (a flag, name, count, fraction, or probability). Temperatures `Tref` and `T_sw` are given in **°C** (converted to K internally); the radiation temperatures `sky_temp` and `T_side` are given in **K**. Pressures are **Pa**.

---

## Which groups do I need?

| Group | Required when | Mode |
|-------|---------------|------|
| `&PARAMETERS` | always | both |
| `&MICROPHYSICS` | `do_microphysics = .true.` (default) | both |
| `&TURBULENCE_ODT` | `simulation_mode = 'chamber'` and `do_turbulence = .true.` | chamber |
| `&TURBULENCE_LEM` | `simulation_mode = 'parcel'` and `do_turbulence = .true.` | parcel |
| `&PARCEL` | `simulation_mode = 'parcel'` | parcel |
| `&ENTRAINMENT` | `simulation_mode = 'parcel'` and `do_entrainment = .true.` | parcel |
| `&RADIATION` | `do_radiation = .true.` | both |
| `&SPECIALEFFECTS` | `simulation_mode = 'chamber'` and `do_special_effects = .true.` | chamber |

---

## `&PARAMETERS` — core run control

Always required. Sets the domain, mode, timing, output, and the master on/off switches for optional physics.

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `N` | integer | — | `2000` | Number of grid cells (domain resolution along the 1-D column). |
| `tmax` | real | s | `100` | Maximum simulation time. |
| `H` | real | m | `1.` | Domain height. |
| `volume_scaling` | real | — | `10` | Cross-sectional area scaling that sets the effective 3-D domain volume (used for number concentrations and collision rates). |
| `Tref` | real | °C | `20.` | Reference temperature. **Chamber:** bottom-boundary temperature. **Parcel:** uniform initial temperature. |
| `pres` | real | Pa | `1.00e5` | Pressure. **Chamber:** constant reference pressure. **Parcel:** initial pressure, evolving per `pressure_mode` (in `'environment'` mode it is replaced by the sounding pressure at `initial_height`). |
| `simulation_mode` | string | — | `'chamber'` | `'chamber'` (ODT, Dirichlet BCs) or `'parcel'` (LEM, periodic BCs, adiabatic ascent). |
| `same_random` | logical | — | `.false.` | Seed the RNG from a fixed state for reproducible/deterministic runs (used by the reftest). |
| `simulation_name` | string | — | **— (required)** | Output file prefix (e.g. `{sim_name}.nc`, `{sim_name}.log`). |
| `output_directory` | string | — | **— (required)** | Directory for all output (absolute, or relative to cwd; parent must exist). |
| `overwrite` | logical | — | `.false.` | Allow overwriting existing output files. |
| `write_timer` | real | s | **— (required)** | Write profile/time-series output every X seconds. |
| `write_buffer` | integer | — | **— (required)** | Number of write-steps buffered in memory before flushing to NetCDF. |
| `write_eddies` | logical | — | `.false.` | Write the eddy event stream to `{sim_name}_eddies.bin`. |
| `do_turbulence` | logical | — | `.true.` | Enable turbulence (ODT in chamber, LEM in parcel). |
| `do_microphysics` | logical | — | `.true.` | Enable aerosol/droplet processes (requires `&MICROPHYSICS`). |
| `do_special_effects` | logical | — | `.false.` | Enable sidewalls / stochastic fallout (requires `&SPECIALEFFECTS`; chamber only). |
| `do_radiation` | logical | — | `.false.` | Enable radiative transfer (requires `&RADIATION`). |
| `do_entrainment` | logical | — | `.false.` | Enable blob entrainment (requires `&ENTRAINMENT`; parcel only). |

---

## `&MICROPHYSICS` — aerosols, droplets, collisions, trajectory output

Required when `do_microphysics = .true.` (the default).

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `aerosol_file` | string | — | **— (required)** | Path to the NetCDF aerosol distribution input (relative paths resolve from the namelist's directory). |
| `init_drop_each_gridpoint` | logical | — | `.false.` | Initialize a droplet at every grid point to reduce spin-up time. |
| `expected_Ndrops_per_gridpoint` | real | — | `1` | Expected droplets per grid point; sizes the initial particle array to avoid reallocation. |
| `initial_wet_radius` | real | × dry radius | **— (required)** | Injected wet radius as a multiple of the dry radius (must be > 1). |
| `aerosol_concentration` | real | cm⁻³ | `0.0` | Initial aerosol number concentration. **Parcel mode only** (chamber injects from `aerosol_file`). |
| `do_seeding` | logical | — | `.false.` | Release seed aerosol at the events defined in `aerosol_file`. Absolute controller: `.true.` **requires** a seed group in the file (absent ⇒ fatal); `.false.` **ignores** any seed group present (not read, no effect on outputs), warning if one is found. Works in both modes. |
| `seed_hydration` | string | — | `'equilibrium'` | Wet radius seed particles are born with: `'equilibrium'` (Köhler solve at local RH, capped at 0.99 and bounded by `seed_growth_time`), `'double_growth'` (2× dry radius), or `'dry'`. Applies to seed material only — background keeps `initial_wet_radius`. |
| `seed_growth_time` | real | s | `5.0` | Growth time a seed is allowed at release, bounding the `'equilibrium'` radius. Keeps GCCN from being born at a radius they would need minutes to grow to. Ignored by the other hydration modes. |
| `do_collisions` | logical | — | `.false.` | Enable collision detection (event-driven 1-D collision-coalescence). |
| `do_coalescence` | logical | — | `.false.` | Merge droplets on collision (requires `do_collisions`). |
| `coalescence_kernel` | string | — | `'hall'` | Collection-efficiency kernel: `'hall'` (Hall 1980), `'long'` (Long 1974), or `'unity'`. |
| `wmax_collision` | real | m/s | `10.0` | Cap on terminal velocity used in collision calculations. |
| `write_collisions` | logical | — | `.false.` | Write the collision/coalescence event stream to `{sim_name}_collisions.bin`. |
| `write_trajectories` | logical | — | `.false.` | Write per-particle trajectories to `{sim_name}_particles.nc`. |
| `trajectory_start` | real | s | `0.` | Start time for trajectory output. |
| `trajectory_end` | real | s | `0.` | End time for trajectory output. |
| `trajectory_timer` | real | s | `1.` | Trajectory write interval. |

---

## `&TURBULENCE_ODT` — chamber turbulence

Required when `simulation_mode = 'chamber'` and `do_turbulence = .true.` Controls the ODT eddy sampling and the convective forcing.

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `Tdiff` | real | °C | `10.` | Top-to-bottom temperature difference driving convection (ΔT). |
| `Lmin` | integer | grid cells | `6` | Minimum eddy size. **Validated at startup (v3.0.0):** must be ≥ 6 and a multiple of 3, or the run aborts. Below 6 each of the triplet map's three segments is a single cell and the eddy carries no sub-eddy structure; the floor is `globals::min_eddy_gridpoints`, the same constant that sets parcel mode's `grid_eddy_scale`. |
| `Lprob` | integer | grid cells | `18` | Eddy-length PDF shape parameter; the length distribution decays as `exp(−2·Lprob/L)`, so larger values favor larger eddies. |
| `max_accept_prob` | real | — | `0.1` | Maximum eddy acceptance probability (caps the eddy event rate). |
| `C2` | real | — | `1.5e3` | Eddy-rate coefficient scaling the buoyant energy available to eddies. |
| `ZC2` | real | — | `1.0e5` | Viscous/length energy penalty; an eddy is accepted only if its available energy exceeds this threshold. |

---

## `&TURBULENCE_LEM` — parcel turbulence

Required when `simulation_mode = 'parcel'` and `do_turbulence = .true.` Sets the Linear Eddy Model eddy sampling from a −5/3 inertial range.

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `integral_length_scale` | real | m | `0.01` | Integral (largest eddy) length scale. |
| `dissipation_rate` | real | m²/s³ | `0.01` | Turbulent kinetic energy dissipation rate. |

### Removed: `kolmogorov_length_scale` (breaking, v3.0.0)

`kolmogorov_length_scale` was **removed** from `&TURBULENCE_LEM`. A namelist that
still declares it is a **fatal read error** — delete the line from existing
`.nml` files.

The smallest turbulence scale is now derived rather than supplied, because a
user-supplied value could contradict the grid. Four scales are computed and
reported in the run log:

| derived quantity | formula | meaning |
|---|---|---|
| `actual_kolmogorov_scale` | `(ν³/ε)^(1/4)` | the physical dissipation scale — **diagnostic only** |
| `grid_eddy_scale` | `6·dz` | the smallest eddy the triplet map can represent at all |
| `diffusivity_length_scale` | `(max(kT,Dv) / (0.1·ε^(1/3)))^(3/4)` | where turbulent diffusivity falls to molecular |
| `smallest_eddy_scale` | `max` of the latter two, rounded **up** to a multiple of 3 cells | **the governing scale** |

`smallest_eddy_scale` sets the eddy-sampler lower bound, the sampler's gridpoint
floor, the Reynolds number `(L/smallest_eddy_scale)^(4/3)`, and the LEM
diffusivity enhancement — so those cannot drift apart. It is also reported in
cells as `smallest_eddy_gridpoints` (always a multiple of 3, always ≥ 6).

Why `actual_kolmogorov_scale` is diagnostic: `0.1·ε^(1/3)·η^(4/3)` reduces to
`0.1·ν` identically (ε cancels), so `diffusivity_length_scale / η = (10·D/ν)^(3/4)
≈ 7.7` for both `kT` and `Dv`. The physical η is always ~8× below the scale at
which this closure's turbulence actually hands off to molecular transport, so it
can never win the `max`. It is reported for reference.

### The LEM diffusivity enhancement

```
diffusivity_enhancement = (smallest_eddy_scale / diffusivity_length_scale)^(4/3)   ≥ 1
```

Two regimes, one formula:

- **`6·dz > diffusivity_length_scale`** (grid-limited) — the grid is coarser than
  the handoff, so eddies between the two are real but unrepresentable. Their
  stirring is absorbed by scaling both LEM diffusivities up by this factor. Here
  `6·dz` is exactly two 3-cell quanta, so the rounding never overshoots:
  `smallest_eddy_gridpoints = 6` and `f = (6·dz / l_D)^(4/3)` exactly. `f` grows
  without bound as the grid coarsens, which is the intent.
- **`6·dz < diffusivity_length_scale`** (diffusivity-limited) — the grid is finer
  than the handoff, so eddies below it are meaningless. The smallest eddy is
  raised to the handoff instead.

Rounding **up** to a multiple of 3 cells is what guarantees the factor is never
below 1, so **LEM diffusion is never slower than molecular**:

| | applies to | value |
|---|---|---|
| `thermal_diffusivity` | LEM diffusion of `T` and the diffusion stability step | `kT · f` |
| `vapor_diffusivity` | LEM diffusion of `WV` and the diffusion stability step | `Dv · f` |

Both are scaled by the *same* `f`, so `Pr` and `Sc` are preserved exactly.

**`f` is a step function of the grid, not a smooth one.** In the
diffusivity-limited regime the quantum `3·dz` can be a large fraction of `l_D` —
up to `l_D/2` right at the crossover — so rounding up can overshoot by as much as
50%, giving `f` up to `(3/2)^(4/3) = 1.717`. It approaches 1 only as the grid
refines. Two nearby values of `N` can therefore give noticeably different
diffusivities: at `ε = 0.01, H = 1`, `N = 1025` gives `f = 1.001` while
`N = 1045` gives `f = 1.675`. This is self-consistent — the model's smallest eddy
really is 9 cells rather than 6 in the second case, and the diffusivity matches
that scale — but it means **`f` should be read from the run log or the
`LEM.diffusivity_enhancement` attribute rather than assumed to be ≈1** whenever
`6·dz` is close to `l_D`.

**Droplet growth is never enhanced.** The DGM computes its own temperature- and
pressure-dependent thermal conductivity and vapor diffusivity internally
(`src/DGM.f90:204-205`, Rogers & Yau Table 7.1) and does not read these values.
Chamber mode is likewise unaffected — ODT diffuses via `Pr`/`Sc` built from the
unmodified `kT`/`Dv` globals.

### Validation

Fatal at startup:

- `smallest_eddy_gridpoints > N` — the domain cannot contain the smallest eddy.
- `smallest_eddy_scale ≥ integral_length_scale` — no inertial range. The −5/3
  sampler would otherwise raise a negative quantity to a fractional power and
  produce silent `NaN`.

Warning (non-fatal): `integral_length_scale / smallest_eddy_scale < 3` (`Re < 4.3`)
— the inertial range is nearly absent and `maps_per_event` will be small. **This
fires on the bundled `input/params.nml`**, where `L = 0.01` and the derived
smallest eddy is 6 mm.

**Consequence for existing runs:** results change, and `maps_per_event` can change
by orders of magnitude where the old input η was far below `6·dz`. To influence
the model's smallest eddy, change `N`/`H` (which moves `6·dz`) or
`dissipation_rate` (which moves `diffusivity_length_scale` as `ε^(-1/4)`).

---

## `&PARCEL` — adiabatic ascent

Required in **parcel** mode (`simulation_mode = 'parcel'`).

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `parcel_file` | string | — | **— (required)** | v3 NetCDF file of waypoint trajectory legs (target level + signed velocity per leg), optionally with an environmental sounding and per-leg entrainment schedule. Completing the last leg ends the run. |
| `initial_RH` | real | fraction (0–1) | `1.0` | Initial relative humidity; sets the initial water-vapor field. |
| `pressure_limit` | real | Pa | `0.0` | Stop the simulation when pressure reaches this target. `0.0` disables the limit. |
| `initial_height` | real | m | `0.0` | Parcel launch height. Leg-1 direction is validated against it (pressure axis: against the initial `pres` instead). In `'environment'` mode the initial pressure is taken from the sounding at this height. |
| `vertical_axis` | string | — | `'height'` | Whether the file's `segment_coord` leg targets are heights (`'height'`, m) or pressures (`'pressure'`, Pa). |
| `pressure_mode` | string | — | `'hydrostatic'` | Parcel pressure evolution. `'hydrostatic'`: self-integrate dp = −ρ·g·w·dt. `'environment'`: follow the sounding's p(z) at the parcel's height, with adiabatic dT from the actual dp (requires the sounding). |

---

## `&ENTRAINMENT` — blob entrainment

Required when `simulation_mode = 'parcel'` and `do_entrainment = .true.` Mixes environmental air into the column via the blob method (environmental profile read from `parcel_file`).

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `ent_rate` | real | **km⁻¹** | `2.0` | Fractional entrainment rate (converted to 1/m internally). |
| `n_blob` | integer | — | `1` | Number of blobs per entrainment event. |
| `psigma` | real | fraction | `0.1` | Blob size as a fraction of the domain (per blob). |
| `random_entrainment` | logical | — | `.true.` | Poisson-randomize entrainment event timing (vs. regular intervals). |

> **Unit change:** `ent_rate` was in 1/m before v3; it is now specified in **1/km** everywhere at the interface (namelist, v3 parcel file, output).

`ent_rate`, `n_blob`, and `psigma` can instead be made **per-leg** by providing per-segment arrays in the `parcel_file` (one value per trajectory leg). When present, those arrays override the constant values here; `random_entrainment` always comes from this namelist. See [Data Formats](data_formats.md#parcel-input--codt_parcel_input_v3).

---

## `&RADIATION` — radiative transfer

> ⚠️ **Experimental.** The radiation code is still under active development and has not been fully validated. Treat results as preliminary and expect parameters and behavior to change.

Required when `do_radiation = .true.`

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `radiation_method` | string | — | `'1d'` | `'1d'` (two-stream) or `'3d'` (Monte Carlo). |
| `mie_data_file` | string | — | **— (required)** | Path to the Mie absorption data table. |
| `eps_top` | real | — | `1.0` | Top-boundary emissivity. |
| `eps_bot` | real | — | `1.0` | Bottom-boundary emissivity. |
| `sky_cooling_flag` | logical | — | `.false.` | Override the top boundary with a radiating sky at `sky_temp`. |
| `sky_temp` | real | K | `263.15` | Sky temperature for cooling. |
| `rad_call_interval` | real | s | `0.0` | Interval between radiation updates; `0.0` calls every physics step. |
| `max_droplets_per_cell` | integer | — | `20` | Cap on droplets per cell included in the radiation calculation. |
| `nPhotons` | integer | — | `700000` | Number of Monte Carlo photons. **`'3d'` only.** |
| `nBins` | integer | — | `30` | Number of spectral bins. **`'3d'` only.** |
| `Lx_rad` | real | m | `2.0` | Horizontal domain extent in x for the 3-D solver. **`'3d'` only.** |
| `Ly_rad` | real | m | `2.0` | Horizontal domain extent in y for the 3-D solver. **`'3d'` only.** |
| `T_side` | real | K | `293.15` | Side-wall temperature for the 3-D solver. **`'3d'` only.** |

---

## `&SPECIALEFFECTS` — sidewalls & stochastic fallout

Required when `simulation_mode = 'chamber'` and `do_special_effects = .true.`

| Parameter | Type | Units | Default | Description |
|-----------|------|-------|---------|-------------|
| `do_sidewalls` | logical | — | `.false.` | Enable sidewall nudging of the scalar fields. |
| `area_sw` | real | m² | **— (required\*)** | Sidewall area. |
| `area_bot` | real | m² | **— (required\*)** | Bottom area. |
| `C_sw` | real | — | **— (required\*)** | Sidewall eddy-velocity coefficient (`velocity_sw = C_sw · velocity_bot`). |
| `sw_nudging_time` | real | s | **— (required\*)** | Sidewall nudging interval. |
| `T_sw` | real | °C | **— (required\*)** | Sidewall temperature. |
| `RH_sw` | real | fraction | **— (required\*)** | Sidewall relative humidity. |
| `P_sw` | real | — | **— (required\*)** | Sidewall tuning parameter. |
| `do_random_fallout` | logical | — | `.false.` | Enable stochastic fallout (droplet "perceived" height at the bottom boundary). |
| `random_fallout_rate` | real | — | `1.` | Rate parameter for stochastic fallout. |

\* The sidewall parameters have no code default; supply them whenever `do_sidewalls = .true.`

---

## See also

- [`input/params.nml`](../input/params.nml) — complete example namelist
- Output file formats — see the `codt-io` reference
- `codt --help` — prints the namelist groups and usage at runtime
