# Known Issues

Running list of CODT behaviors that need investigation and possibly a source
change. This is a working list, not a bug tracker: items here are things we have
*observed*, with whatever evidence we had at the time. Some have a confirmed
cause, some are still open questions.

**If you are about to change model physics, numerics, or output, read the
relevant section first.** Several of these are latent — they do not crash, they
quietly corrupt results.

## Conventions

Each item is `### Short title` with:

- **Symptom** — what is actually observed, with numbers where we have them.
- **Evidence** — the run, file, or measurement it came from.
- **Suspected location** — best guess at the source file/routine. Marked
  *unverified* when nobody has confirmed it by reading the code.
- **Found in** — the CODT commit the observation came from: the `git_commit` /
  `code_version` global attribute of the simulation output for run-discovered
  items, or the branch/worktree for items found while editing code. Without this
  an item cannot be reproduced, and "still broken?" is unanswerable.
- **Status** — `open`, `cause known / unfixed`, `investigating`, or `fixed in <commit>`.

When an item is fixed, strike it in the same commit as the fix and note the
commit hash. Do not delete it — the history of what was wrong is useful.

Sections: **Physics**, **Input/Output**, **Coding**, **Unknown**. Put an item in
*Unknown* when you cannot yet tell whether it is a physics assumption, a config
trap, or a code defect; move it once you know.

---

## Physics

### Runaway collision-coalescence collapse in parcel mode

- **Symptom** — With collisions and coalescence on and low or no entrainment,
  the droplet population collapses to a single large drop. In a 6-member
  entrainment sweep, `Np` fell 3000 -> 1 for `ent_rate = 0` and `0.1` /km, and
  3000 -> 500 at `0.2`; runs at `0.5`-`1.0` /km retained ~2000-3000. The
  collapse in the undiluted run happened between t ~ 1400 s and ~1800 s, ending
  in one ~150 um drop. Once a single drop is the only condensation sink, the
  domain-mean supersaturation runs away (S -> ~31%), so every thermodynamic
  field after the collapse is meaningless.
- **Cause** — **There is no fallout/sedimentation in the model.** Drops that
  grow large enough to precipitate are retained in the domain and keep
  colliding, so the coalescence cascade has no sink and runs to completion.
  Entrainment only masks this by suppressing collision rates.
- **Consequence** — Any parcel run with `do_coalescence = .true.` and weak
  entrainment is unusable past the onset of collapse. Check `Np(t)` before
  trusting late-time output.
- **Evidence** — `projects/madison/`, SF1 entrainment sweep, runs completed
  2026-07-22 (`madison_SF1_ent{0,0.1,0.2,0.5,0.8,1}`).
- **Found in** — `527cd7f` (`v2.0.0-1-g527cd7f`), `main`; from the output
  `git_commit` attribute.
- **Status** — cause known / unfixed. Needs a sedimentation or large-drop
  removal mechanism in parcel mode.

### Domain liquid budget is hostage to a few high-mass droplets

- **Symptom** — In a small domain, once collision-coalescence concentrates the
  condensate into a handful of large drops, those few particles carry nearly all
  the liquid mass, and any process that removes one of them removes most of the
  domain's water in a single step. Observed in a 3 cm^3 parcel domain: one 121 um
  drop held **96% of the liquid mass** (mass proxy 1.77e6 of 1.84e6 total; 2.01 of
  2.70 g/m3 LWC) while ~245 other droplets held the rest. A single entrainment
  blob event then deleted it, and in one output interval LWC fell 2.698 -> 0.080
  g/m3, the largest drop went 121 -> 13.6 um, and the mass fraction in drops
  r > 40 um went ~96% -> 0%.
- **The tell** — `Np` *rises* across the event (246 -> 510) even though
  coalescence can only lower it: the blob replaced ~10% of the domain
  (`n_blob = 1`, `psigma = 0.1` -> 1 m of a 10 m domain = 0.3 cm^3) and carried
  in ~300 fresh aerosol at 1000 cm^-3. Removing 10% of the domain *volume*
  removed 96% of its liquid *mass*.
- **Why it matters** — This is a sampling problem, not physics. With few
  droplets, per-simulation results become strongly stochastic: the liquid budget
  depends on whether one particle happens to fall inside the next entrainment
  blob, which is roughly a `psigma` coin flip per event. Runs with weak or no
  entrainment keep their large drop and report ~100% large-drop mass; runs with
  entrainment lose it at an essentially random height. Neither number is
  meaningful, and **single-realization results should not be trusted** in this
  regime — quantities dominated by the tail of the mass distribution (LWC,
  radar reflectivity, precipitation-sized mass fraction, effective radius) need
  either a larger domain / higher droplet count or an ensemble.
- **How to detect it** — Check whether a small number of bins dominate
  `sum(DSD * r^3)`. If the top few particles carry most of the mass proxy, the
  domain is undersampled for any mass-weighted diagnostic. A discontinuous jump
  in `Np` marks an entrainment event.
- **Related** — Strongly coupled to the runaway collapse above (no fallout is
  what concentrates the mass into one drop in the first place), but distinct:
  this failure would occur in any small domain whose mass distribution is
  tail-dominated, even with a working sedimentation sink.
- **Evidence** — `projects/madison/`, `madison_SF1_ent0.2`, output indices
  828 -> 829 (z 6822 -> 6826 m). Same signature in `ent0.5` and `ent0.8`;
  `ent0` and `ent0.1` never lose their drop.
- **Found in** — `527cd7f` (`v2.0.0-1-g527cd7f`), `main`; from the output
  `git_commit` attribute.
- **Status** — cause known / unfixed. Needs guidance on minimum droplet count
  (or `volume_scaling`) for mass-weighted diagnostics, and probably an ensemble
  convention for entraining parcel runs.

---

## Input/Output

### `pressure_mode = "environment"` silently launches the parcel at the surface

- **Symptom** — Switching `pressure_mode` from `hydrostatic` to `environment`
  takes the initial pressure from the sounding at `initial_height`, which
  defaults to `0.0`. If `initial_height` is not also set, the parcel launches at
  the surface rather than at the intended level, with no warning.
- **Why it bites** — Studies that launch at the LCL set `pres` directly and
  assume it is honored. Under `environment` mode it is not; `pres` is
  effectively ignored in favor of the sounding lookup.
- **Suspected location** — `parcel.f90:189` (as noted in
  `projects/madison/config.py`) — *unverified against current HEAD*.
- **Found in** — noted while configuring against `527cd7f`, `main`. Not observed
  in a run: the madison sweep deliberately used `pressure_mode = "hydrostatic"`
  to avoid it.
- **Status** — open. At minimum this deserves a warning or a validation error
  when `pressure_mode = "environment"` and `initial_height` is left at default.
  Being examined as part of a planned hydrostatic-vs-environment comparison.

---

## Coding

### `parcel_height_env` is 0 in the first output record

- **Symptom** — `parcel_height_env` (environment height at the parcel pressure)
  writes as exactly `0.0` at t = 0, then jumps to its true value at the first
  step. In the madison SF1 runs it goes `0.0 -> 3510.57 -> 3514.58 ...` for a
  parcel launched near 3.5 km MSL. `parcel_pressure` is correct at t = 0
  (671.0 mb), so the pressure is known — only the derived height is not.
- **Consequence** — Latent and silent. Any budget using the first record picks
  up a spurious ~3.5 km of `gz`. This made an adiabatic parcel appear to *gain*
  34 kJ/kg of moist static energy over its ascent; discarding the first record
  gave -0.05 kJ/kg (conserved, as expected). Anything integrating or
  differencing from t = 0 is affected.
- **Open question** — **Why is it initialized to 0?** Not yet investigated.
  Possibilities: the diagnostic is computed in the step routine and never
  evaluated before the first write; or it is default-initialized and the
  initial-write path skips it. Needs someone to read the write path.
- **Workaround** — Slice `[1:]` when reading it. See
  `projects/madison/analyze.py`.
- **Suspected location** — parcel step/diagnostic path and the initial record in
  the writeout routine — *unverified*.
- **Found in** — `527cd7f` (`v2.0.0-1-g527cd7f`), `main`; from the output
  `git_commit` attribute of the madison SF1 runs.
- **Status** — open, cause not yet investigated.

### `rdx` out-of-bounds access in the tridiagonal solver

- **Symptom** — The tridiagonal solve indexes `rdx` as `rdx(N)` where it should
  be `rdx(N+1)`. Since the commit that introduced the ghost point, the array
  needs the extra element, so the solver reads one element short of what the
  grid now requires.
- **Consequence** — Latent in all versions since the ghost-point change. Will
  not necessarily crash; expect subtly wrong diffusion at the domain edge.
- **Suspected location** — tridiagonal solver routine — *exact file/line
  unverified; needs to be pinned down before fixing.*
- **Found in** — reported from an earlier debugging session on `main`; the
  originating commit/worktree was not recorded. The defect was introduced by the
  ghost-point work, `9a454ff` ("Added ghost points 0:N+1 to nondimensional
  arrays") and `0e070e2` ("Added (ghost) point to top of nondimensional scalar
  arrays"), so it is present on every commit from `0e070e2` forward — including
  current `527cd7f`. Not tied to a specific simulation run.
- **Status** — open. Reported from an earlier debugging session; the fix looks
  mechanical but should be confirmed against the ghost-point commit and covered
  by a reftest before being applied.

### ~~`move_particles_in_eddy` does full-domain work per eddy event (86% of runtime)~~ — fixed in `ab89628`

- **Symptom** — `droplets::move_particles_in_eddy` is **~86% of total runtime**
  in parcel mode, dwarfing everything else. It is a bookkeeping routine, not
  physics. Flat profiles of the madison SF1 runs:

  | routine | ent0 | ent1 |
  |---|---|---|
  | `droplets::move_particles_in_eddy` | **86.65%** (2328 s) | **86.31%** (2983 s) |
  | all DGM / growth / ODE routines combined | ~6.5% | ~7.6% |
  | `lem::diffuse_scalar_periodic` | 2.36% | 1.82% |
  | `collision_coalescence::push_fall_event` | 1.16% | 1.36% |

  Total sampled time 2687 s (ent0) and 3456 s (ent1). Call counts 463,222,331
  and 458,104,324 — roughly 5,700 eddy events per turbulence step across 81,368
  steps. Measured cost is **5.0 us per call**.
- **Cause** — Two full-domain operations are performed on *every* eddy event,
  regardless of eddy size (`droplets.f90:575`):

  ```fortran
  real(dp) :: mapped_z(N)         ! N = 6000
  mapped_z = z                    ! 48 KB copy, every call
  call triplet_map(L, M, mapped_z)

  do i = 1, current_n_particles   ! ~6000 particles, every call
      gc = lparticles(i)%gridcell
      if (mapped_z(gc) /= z(gc)) then ...
  ```

  Cost per eddy is `O(N + n_particles)` when the eddy itself spans only `L`
  cells. `triplet_map` is *not* the problem — it is `O(eddy_length)` and costs
  only 1.38% / 1.05% despite 1.39 billion calls.
- **Why the eddies are small** — with `H = 10`, `N = 6000` the grid spacing is
  `dz = 1.67 mm`, and `kolmogorov_length_scale = 0.001` puts most sampled eddies
  at the 3-gridpoint floor (`LEM.f90:140`, quantized to multiples of 3). So
  typical `L` is ~3-30 against `N = 6000`: **two to three orders of magnitude**
  more work than the eddy requires.
- **Note on an earlier assumption** — this is why removing the eddy acceptance
  method did not speed runs up appreciably. Neither acceptance sampling nor DGM
  was ever the bottleneck; the cost was always this routine.
- **Fix direction (not yet attempted)** — only cells in `[M, M+L-1] (mod N)` can
  move, and only particles whose `gridcell` lies in that window. The `mapped_z`
  copy can be reduced to `O(L)` outright and looks straightforward. The particle
  scan is the harder half: the particle list is not indexed by gridcell, so
  getting below `O(n_particles)` needs a per-cell index or a sorted list, which
  is a real design change. Any fix must be bit-reproducible against a reftest —
  the routine sets particle positions, so an error here corrupts all
  microphysics downstream.
- **Evidence** — gprof flat profiles and call graphs from the original madison
  SF1 runs, which were themselves built with `-pg`:
  `/scratch/general/vast/u1342804/madison/madison_SF1_ent{0,1}/inputs/gmon.out`
  (2026-07-21). Symbols resolved against a rebuild at the same commit and flags
  (gfortran `-O2 -pg`); the original executables were not retained.
- **Suspected location** — `src/droplets.f90:575` (`move_particles_in_eddy`),
  called from `src/LEM.f90:151` and `src/ODT.f90:585`. Verified by reading the
  source, not just inferred from the profile.
- **Found in** — `527cd7f` (`v2.0.0-1-g527cd7f`), `main`; from the output
  `git_commit` attribute of the madison SF1 runs.
- **Status** — **fixed in `ab89628`** (released as `v2.0.1`). Performance only —
  no evidence of a correctness defect. `move_particles_in_eddy` now looks up the
  mapped source cell one at a time via `globals::triplet_map_cell` (an O(1)
  inverse of `triplet_map`) instead of materializing an N-length copy of the
  grid. Exact / bit-identical: both reftests unchanged and
  `test/test_remap_compare` shows 0 mismatches between the old and new remapping
  for a mid-domain and a wrapping eddy. Only the full-domain `mapped_z` copy was
  removed; the per-particle `O(n_particles)` scan remains and is a possible
  future enhancement (would require indexing particles by gridcell — see item 2
  of the profiling discussion, deferred as non-critical).

---

## Unknown

*(No items yet. Use this section for behavior that has not been classified —
move it to Physics, Input/Output, or Coding once the cause is understood.)*
