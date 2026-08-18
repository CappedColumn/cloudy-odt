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

1. **Build from source** with the Fortran Package Manager. Build flags live in
   `fpm.toml` profiles (`release`/`debug`/`profiled`); `fpm_env` loads the
   compiler module and the netCDF paths live in the `netcdf-local` feature of
   your local `fpm.toml` (see `fpm.toml.template` / `fpm_env.template`):

   ```bash
   source fpm_env                                              # gfortran (default); `source fpm_env nvfortran` to switch
   fpm build   --profile release --compiler gfortran           # optimized build
   fpm install --profile release --compiler gfortran --prefix <install-dir>   # puts CODT in <install-dir>/bin
   export PATH="<install-dir>/bin:$PATH"
   ```

   A bare `fpm build` (no `--profile`) is unoptimized — always pass
   `--profile release` for production, or use `./build.sh` (defaults to release
   and stamps the version into output metadata).

2. **Use a pre-built executable** if one has been provided to you (e.g. inside a
   shared conda environment). Activate/locate it, then confirm it resolves:

   ```bash
   which CODT          # case-sensitive!
   CODT --help         # prints usage + namelist group summary
   ```

`--help`/`-h` and `--version`/`-v` short-circuit *before* any simulation modules
load, so they're a cheap way to confirm the binary is healthy.

> **NetCDF / shared libraries:** CODT's manifest bakes an rpath to its NetCDF
> libraries into the binary (the `netcdf-local` feature's `-Wl,-rpath=...`), so a
> normally-built executable is self-contained and runs without loading any module
> — important for batch/`nohup` jobs. If you *do* see `error while loading shared
> libraries: libnetcdf...` at launch (e.g. an executable built without the rpath),
> make the same NetCDF it was built against available at runtime — `module load`
> the matching netcdf module, or set `LD_LIBRARY_PATH` — in the same shell/job that
> runs CODT.

### 1a. Building for a specific target architecture (CHPC)

CODT's default `release`/`debug`/`profiled`/`vec` profiles target a portable
x86-64 baseline (`sandybridge`), safe to run on any CHPC node. Building for a
*specific* microarchitecture (e.g. `-march=znver2` for AMD Rome nodes on
notchpeak, used by the `notchpeak-rome` profile) can unlock real speedups —
e.g. FMA3 fused multiply-add, which the sandybridge baseline predates — but
the netCDF the compiler links against must match the compiler *and* the
target CPU family, and CHPC doesn't build that combination for every
compiler version.

**CHPC keeps several parallel generations of Spack-built software**, not one
rolling tree. As of 2026-08, under `/uufs/chpc.utah.edu/sys/spack/`:

| Root | Last modified | Notes |
|------|---------------|-------|
| top-level `spack/` | 2026-05-28 | oldest netCDF-Fortran (4.5.3) but broadest compiler coverage (gcc/intel/nvhpc) |
| `v020/` | 2023-07-02 | effectively retired |
| `v019/` | 2025-06-23 | has the only working **gcc-11.2.0 + zen2** netCDF-Fortran build (what `notchpeak-rome` uses) |
| `v11/` | 2026-05-28 | newer generation; only nehalem/gcc-8.5.0 has netCDF-Fortran so far |
| `v10/` | 2026-08-11 (newest) | has gcc-13.4.0/15.1.0 and more arches, but **no netCDF-Fortran build at all yet** — netCDF-C only |

**The newest generation is not necessarily the most complete one.** A newer
root can have a compiler or arch tree with no matching netCDF-Fortran build
yet — check before assuming "newest = best available." As found in 2026-08,
netCDF-Fortran only exists for these arch/compiler pairs, system-wide:

| Arch target | Compiler | Root |
|---|---|---|
| nehalem (portable) | gcc/8.5.0 | `spack/`, `v019/`, `v11/` |
| nehalem (portable) | gcc/11.2.0 | `v019/` |
| nehalem (portable) | intel/18.0.5, 2021.4.0, 2021.7.1 | `spack/`, `v019/` |
| nehalem (portable) | nvhpc/21.5, 21.7 | `spack/` |
| sandybridge | gcc/8.5.0, nvhpc/21.5 | `spack/`, `v019/` |
| **zen2 (AMD Rome match)** | **gcc/11.2.0** | **`v019/`** |

Everything else — gcc 13.x/15.x, nvhpc 20.x/23.x/24.x/25.x, Intel oneAPI
2022+/2025, skylake/cascade-lake-specific builds — exists as a compiler or
netCDF-C-only tree, but has no matching netCDF-Fortran anywhere, and would
need to be built from source (or via user-space Spack against CHPC's
upstream — see `chpc.utah.edu/documentation/software/spack.php`) before it's
usable by CODT.

**Adding a new arch-specific fpm profile:** find the matching
arch/compiler/netCDF-Fortran triple (searching all spack roots, not just the
newest), add a `(compiler, arch)` case to `fpm_env`, and a matching
`[features.optimized-<arch>.*]` / `[features.netcdf-<arch>.*]` pair plus
profile in `fpm.toml`. If `fpm.toml`'s global `link = ["netcdf", "netcdff"]`
is in play, remember netCDF-C and netCDF-Fortran are sometimes **separate**
spack packages (as in the `v019`/zen2 build) — both lib directories need
`-L`/`-rpath`, or `-lnetcdf` silently falls back to a stale system library.

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
| `seed_coord` (in `aerosol_file`) | a **time** [s] | a **height** [m] or **pressure** [Pa], per `vertical_axis` |
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
12. **`do_seeding` is the absolute controller of seeding.** With `do_seeding=.true.`
    the `aerosol_file` **must** contain a seed group — its absence is fatal. With
    `do_seeding=.false.` any seed group present is simply **ignored** (not read, no
    effect on output); CODT prints a warning so you know the seed data is dormant,
    but does not abort. The bundled `input/aerosol_input.nc` ships **with** a seed
    group, so it runs seeded when you set `do_seeding=.true.` and unseeded (with the
    warning) otherwise.
13. **Seeding events fire once, on first arrival.** A parcel that re-crosses a seeded
    level does *not* seed again, and events are keyed to position, not to trajectory
    legs. If you want a level seeded on each pass, that is not currently expressible.
14. **A seed chemically identical to your background still needs its own composition
    row.** What makes material "seed" is that `seed_bin_type` points at it, so seeding
    NaCl into NaCl means duplicating the row. Reusing the background's row is rejected.
15. **Compiler/NetCDF mismatch at runtime.** If you see `error while loading shared
    libraries: libnetcdf...`, the right module isn't loaded — re-evaluate how you set
    up the environment in the same shell/job that runs CODT.
16. **Login-node etiquette.** Long runs belong on compute/owner nodes via batch
    processing. Respect your HPC center's rules.

---
