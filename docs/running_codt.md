# Running CODT Manually

This guide covers running the **CODT** Fortran executable directly — getting the
environment in place, preparing inputs, launching a run, and confirming it finished.
It focuses on the *concepts and contract* of the simulation, not any particular
wrapper. (The Python `codt_tools` package automates much of this, but everything
below works with just a shell, a namelist, and the binary.)

---

## 1. Set up the environment

You need a `CODT` executable on your `$PATH` (or an explicit path to one). Two
common routes:

1. **Build from source** with the Fortran Package Manager:

   ```bash
   source fpm_env                            # compiler + NetCDF paths (see fpm_env.template)
   fpm build
   fpm install --prefix <install-dir>        # puts CODT in <install-dir>/bin
   export PATH="<install-dir>/bin:$PATH"
   ```

2. **Use a pre-built executable** if one has been provided to you (e.g. inside a
   shared conda environment). Activate/locate it, then confirm it resolves:

   ```bash
   which CODT          # case-sensitive!
   CODT --help         # prints usage + namelist group summary
   ```

`--help`/`-h` and `--version`/`-v` short-circuit *before* any simulation modules
load, so they're a cheap way to confirm the binary is healthy.

> **NetCDF / shared libraries:** if you see `error while loading shared libraries:
> libnetcdf...` at launch, the executable can't find its NetCDF libraries — make the
> same NetCDF available at runtime that it was built against (e.g. `module load` the
> matching netcdf module, or use an executable whose NetCDF was built with rpath) in
> the same shell/job that runs CODT.

---

## 2. The executable contract

With the executable in your `$PATH` you can run `CODT` directly. If not, you can still run it by providing the full path to the executable: `/path/to/CODT [arguments]`

```
CODT <NAMELIST_PATH>      # run a simulation
CODT --help               # usage and namelist group summary
CODT --version            # "CODT <version> (<git_commit>)"
```

A few hard rules baked into `app/main.f90`:

- **The namelist path must contain a `/`.** CODT rejects a bare filename. Use
  `./params.nml` for a file in the current directory. (This is so it can reliably
  resolve input files relative to the namelist's directory.)
- **Input files are resolved relative to the namelist's parent directory.** Inside
  the namelist, `aerosol_file` and `parcel_file` are interpreted relative to wherever
  the `.nml` lives — *not* relative to your current working directory.
- **`output_directory` is resolved from the current working directory** (it may be
  absolute or relative), and **its parent must already exist**. CODT writes files
  *into* that directory; it does not create a nested per-run subdirectory.
- Exit code **0 = success, 1 = error**. Errors are written to **stderr**, so you see
  them on the terminal even though normal stdout is redirected to the log (below).

---

## 3. Check your inputs

Before launching, walk the namelist top-to-bottom and confirm the pieces line up.

**`&PARAMETERS` (always required):**

- `N`, `tmax`, `H`, `Tref`, `pres`, `volume_scaling` are range-checked; out-of-range
  values are rejected up front (and *all* problems are reported at once, not one per
  run).
- `simulation_mode` is `'chamber'` (default) or `'parcel'`. This choice changes which
  other groups matter — see the mode table below.
- `simulation_name` becomes the **prefix on every output file**. `output_directory` is
  where they land.
- The physics switches — `do_turbulence`, `do_microphysics`, `do_radiation`,
  `do_entrainment`, `do_special_effects` — gate whole feature groups.

**Mode-dependent groups:**

| Group | When it applies |
|-------|-----------------|
| `&TURBULENCE_ODT` | chamber mode turbulence (ODT) |
| `&TURBULENCE_LEM` | parcel mode turbulence (LEM) |
| `&PARCEL` | parcel mode only (`parcel_file`, `initial_RH`, `pressure_limit`) |
| `&ENTRAINMENT` | only when `do_entrainment=.true.` |
| `&RADIATION` | only when `do_radiation=.true.` (`radiation_method` is `'1d'` or `'3d'`) |
| `&SPECIALEFFECTS` | chamber mode only (sidewalls, stochastic fallout) |
| `&MICROPHYSICS` | when `do_microphysics=.true.` (aerosol, collisions, trajectories) |

**Referenced files must exist** (relative to the namelist dir):

- `aerosol_file` (e.g. `aerosol_input.nc`) when microphysics is on.
- `parcel_file` (e.g. `parcel_input.nc`) in parcel mode — **required** for the
  waypoint trajectory (`CODT_parcel_input_v3` schema; v1/v2 are no longer read).
  With `do_entrainment` or `pressure_mode = 'environment'`, the file must also
  carry the environmental sounding — see [Data Formats](data_formats.md).
- `mie_data_file` when radiation is on.

A quick pre-flight check:

```bash
NML=/path/to/params.nml
ls -l "$NML"                                   # namelist exists
grep -E 'aerosol_file|parcel_file|mie_data' "$NML"   # see what it references
ls -l "$(dirname "$NML")"/                      # confirm those files sit next to it
```

---

## 4. Set up the output directory

CODT writes **into** `output_directory`, and its parent must exist. The simplest,
most predictable layout is a per-run folder you create yourself:

```bash
mkdir -p ~/runs/exp01/output     # 'output' (or its parent) must exist before launch
```

Then in the namelist:

```fortran
simulation_name  = "exp01"
output_directory = "output/"     # relative to where you launch CODT, or use an absolute path
overwrite        = .false.       # default: refuse to clobber an existing {name}.nc
```

- With `overwrite=.false.` (the default), CODT **aborts if `{simulation_name}.nc` already
  exists** — a guard against silently destroying a previous run. Set `.true.` only
  when you intend to replace it.
- Input files (the namelist, aerosol NC, parcel NC) are **not** copied into the
  output directory. If you want a self-contained record of a run, be sure to copy them to a safe location.

---

## 5. Run it locally

Launch from a directory where your relative `output_directory` makes sense (or
to the absolute directory you specified in `output_directory`):

```bash
cd ~/runs/exp01
CODT ./params.nml &     # to run in the background
```

Making use of the `nohup` or `screen` programs can keep a local simulation running even if your terminal session is interrupted or exited.

What happens:

- **Normal stdout** (all init banners + runtime progress) is redirected by CODT into
  a single log file: **`{output_directory}/{simulation_name}.log`**. Tail it to watch
  progress:

  ```bash
  tail -f output/exp01.log
  ```

- **Errors go to stderr**, so a fatal problem prints on your terminal even though the
  normal log is redirected.

For anything longer than a quick test, run on a compute node using the batch submission process of your HPC cluster (load the same
modules in the batch script, `cd` to the run dir, then call `CODT ./params.nml`).

---

## 6. Confirm completion

The success signal is a **marker file**, not just the program exiting:

```
{output_directory}/{simulation_name}_DONE
```

This `_DONE` file (with a completion timestamp) is written **only on a successful
run**. The reliable completion check is:

```bash
ls output/exp01_DONE && echo "FINISHED" || echo "not done / failed"
```

On success you'll find, in `output_directory/`:

| File | When |
|------|------|
| `{name}.nc` | always — main profiles + time series |
| `{name}.log` | always — redirected stdout |
| `{name}_DONE` | always on success — completion marker |
| `{name}_particles.nc` | if `write_trajectories=.true.` |
| `{name}_collisions.bin` | if `write_collisions=.true.` |
| `{name}_eddies.bin` | if `write_eddies=.true.` |

If `{name}_DONE` is **absent** but the process has exited, the run failed — read the
tail of `{name}.log` and the stderr you captured. Exit code `1` confirms failure.

```bash
ncdump -h output/exp01.nc | head -40    # quick peek at the main output header
```

---

## 7. Chamber vs. parcel — what changes

The mode flips several inputs' meanings. The common traps:

| Parameter | Chamber | Parcel |
|-----------|---------|--------|
| Turbulence group | `&TURBULENCE_ODT` | `&TURBULENCE_LEM` |
| `Tref` | bottom boundary T (°C→K) | uniform initial T (°C→K) |
| `pres` | constant reference pressure | initial pressure, evolves hydrostatically |
| `Tdiff` | top–bottom ΔT (drives convection) | unused |
| `initial_RH` | ignored | sets initial water-vapor field |
| `aerosol_concentration` | ignored | droplet number concentration (cm⁻³) |
| `parcel_file` | unused | **required** (waypoint trajectory) |
| Particle init | injected over time (aerosol schedule) | pre-loaded + Köhler-equilibrated |
| `&SPECIALEFFECTS` | available (sidewalls, fallout) | not applicable |
| Entrainment | not yet implemented | blob method (`do_entrainment` + `&ENTRAINMENT`) |

---

## 8. Gotchas

A checklist of the things that most often bite when running by hand:

1. **Namelist path needs a slash.** `CODT params.nml` errors out; use
   `CODT ./params.nml`.
2. **Relative input paths are relative to the namelist, not your cwd.**
   `aerosol_file` / `parcel_file` / `mie_data_file` resolve from the `.nml`'s parent
   directory. Keep inputs next to the namelist. One simulation should have one input folder.
3. **`output_directory`'s parent must already exist.** CODT writes into the directory
   but won't create a missing parent — `mkdir -p` first.
4. **`overwrite=.false.` aborts on an existing `{simulation_name}.nc`.** Either bump
   `simulation_name`, point at a fresh `output_directory`, or set `overwrite=.true.`
   deliberately.
5. **No `_DONE` file = the run did not succeed.** Don't infer success from the `.nc`
   existing — a partial/aborted run can leave one behind. Check for `{name}_DONE`. If not present, check your terminal or the location of the redirected `stderr` for an error message.
6. **Errors are on stderr, progress is in the log.** If you redirect stdout to
   `/dev/null` you'll still see fatal errors on the terminal, but you lose the
   progress log — redirect to a file if you want both.
7. **Parameter validation happens before output init, and reports *all* problems at
   once.** If the run dies immediately, read every error line, not just the first.
8. **Mode mismatches are silent in spirit, loud in effect.** Setting `Tdiff` in parcel
   mode (ignored) or forgetting `parcel_file` in parcel mode (fatal) are classic
   mistakes — confirm `simulation_mode` matches the groups you filled in.
9. **A parameter in the wrong namelist group is rejected** ("Invalid parameter in
   `&GROUP`"). Notably: `do_entrainment` lives in `&PARAMETERS`, `pressure_limit`/
   `pressure_mode`/`vertical_axis` in `&PARCEL`, and `ent_rate`/`n_blob`/`psigma`/
   `random_entrainment` in the standalone `&ENTRAINMENT` (`ent_rate` in **1/km**).
10. **`radiation_method` is `'1d'` or `'3d'`** (two-stream vs. Monte Carlo) — not
    `'two_stream'`.
11. **Cross-setting warnings won't stop the run but signal a likely mistake:**
    `do_radiation` without `do_microphysics`, `write_eddies` without `do_turbulence`,
    `do_entrainment` in chamber mode. Treat them as a prompt to re-check intent.
12. **Compiler/NetCDF mismatch at runtime.** If you see `error while loading shared
    libraries: libnetcdf...`, the right module isn't loaded — re-evaluate how you set
    up the environment in the same shell/job that runs CODT.
13. **Login-node etiquette.** Long runs belong on compute/owner nodes via batch
    processing. Respect your HPC center's rules.

---
