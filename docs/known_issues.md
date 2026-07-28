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

### Eddy-particle movement (`move_particles_in_eddy`) is ~90% of parcel runtime

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
  > **Superseded in v3.0.0** (`feature/lem-empm-port`, slice 2):
  > `kolmogorov_length_scale` no longer exists. The sampler's lower bound is now
  > `smallest_eddy_scale = max((nu^3/eps)^(1/4), 6*dz)`, so it can never fall
  > below the grid and the pile-up at the floor is gone. The eddy-size
  > *distribution* changes; the `O(N + n_particles)`-per-eddy cost described here
  > was already addressed separately by the unified mover in `a760c51`.
- **Note on an earlier assumption** — this is why removing the eddy acceptance
  method did not speed runs up appreciably. Neither acceptance sampling nor DGM
  was ever the bottleneck; the cost was always this routine.
- **Reframing (2026-07-23) — the 86% attribution above is misleading.** That
  number came from `-pg`/gprof, which lumps the whole routine together. The
  `mapped_z = z` copy is **not** wasted work: it is a per-eddy *cache*. The
  triplet map is computed once per eddy (`O(L)`) into `mapped_z`, and every
  particle then does a cheap array lookup `mapped_z(gc)`. The intrinsic cost is
  the sheer number of eddy-events × particles (a genuinely hot path), not the
  copy.
- **Attempted fix `ab89628` (shipped as v2.0.1) — REVERTED, it was a ~39%
  REGRESSION.** It removed `mapped_z` and had each particle re-derive its mapped
  cell via a new `globals::triplet_map_cell` (integer `modulo` arithmetic),
  called ~2.69e9 times. That function does **not** inline even at `-O2`, so it
  replaced one cheap bulk copy with billions of expensive non-inlined modulo
  calls. Clean 40-core benchmark (ent1, tmax=10, `perf` on release builds):

  | build | wall | hotspot |
  |---|---|---|
  | v2.0.0 (copy+cache) | **49 s** | `move_particles_in_eddy` 90% |
  | v2.0.1 (per-cell)   | **68 s** | `triplet_map_cell` 61% + `move_particles` 32% |

  Lesson: the per-cell recompute loses whenever `n_particles ≈ N` (madison:
  both ≈ 6000). Bit-identical output, but slower — a pessimization.
- **Correct fix direction (planned, not yet done)** — keep the cache idea but
  shrink it: only cells in `[M, M+L-1] (mod N)` move, so build a small
  **eddy-local** mapped array of length `L` once per eddy (`O(L)`, no full-N
  copy, no per-particle modulo), then per-particle array lookup indexed by
  offset within the eddy. This beats *both* v2.0.0 (drops the `O(N)` copy) and
  v2.0.1 (drops the per-particle modulo). Needs planning + a reftest gate
  (the routine sets particle positions; an error corrupts all microphysics).
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
- **Status** — **resolved on `feature/lem-empm-port`** (slice 1). The routine was
  replaced by the composed-map path described in the next entry, which touches
  each droplet once per turbulence step instead of once per map: madison ent1
  went **48.1 s -> 4.0 s (~12x)**, and the mover dropped to 0.2% of runtime with
  the hotspot moving onto DGM/ODE physics.
- **Correction to an earlier claim in this entry** — a previous revision stated
  *"Performance only — no correctness defect in either version."* **That was
  wrong.** The routine had two independent correctness defects, documented in
  the next entry. The performance framing is what kept them unexamined for so
  long: the routine was read as bookkeeping, so nobody checked its arithmetic.

### Droplet transport through triplet maps used the wrong map direction

Two independent defects in `droplets::move_particles_in_eddy`, both silent —
droplet count was conserved, positions stayed in `[0,H)` and indices stayed
valid, so bulk statistics looked plausible. What was corrupted was the pairing
between each droplet and the thermodynamic history of the air around it.

- **Defect A — wrong map direction (affected both parcel and chamber).**
  `triplet_map` (`src/globals.f90`) is written as a *gather*: after the call,
  `field(k)` holds the value pulled from `source_index`, i.e. "what arrived in
  cell k". That receiver's-side view is exactly right for scalars — every cell
  is filled with whatever landed in it, and `T`/`WV` were always correct.

  A droplet asks the opposite question — "where did *my* fluid go?" — which is
  the sender's view, the inverse permutation. `move_particles_in_eddy` indexed
  the gather array at the particle's own cell, so droplets were advected
  *backwards* through the rearrangement, landing in a parcel they had no history
  with.

- **Why it went unnoticed for so long** — the two tables coincide exactly when
  the permutation is its own inverse (all cycles of length 2). That is true for
  eddy lengths **3 and 6**, and false for **9 and above**:

  | eddy length | cells where gather and forward tables differ (N=36) |
  |---|---|
  | 3, 6 | 0 — identical |
  | 9 | 4 |
  | 12 | 10 |
  | 18 | 16 |
  | 36 | 30 |

  Every case in `test/test_move_particles.f90` used `L=6`, except one `L=9` case
  whose particle sat on a fixed point of both maps. The suite was therefore
  direction-blind: no test ever asked a droplet to travel a cycle longer than 2.

- **Defect B — stale `gridcell` across composed maps (parcel only).**
  `move_particles_in_eddy` updated `position` but never `gridcell`. With one
  eddy per step that is harmless, because `move_particles_by_gravity` /
  the collision branch refresh the index every step
  (`src/droplets.f90:143-144`, `352-354`). LEM applies `maps_per_event` maps
  (madison ~5693) between refreshes, so every map after the first looked up its
  displacement against the cell the droplet had already left. Note that moving
  per eddy is *not* wrong in itself — it is wrong only without carrying the
  index forward.

- **Resolution** — a single eddy-sequence API in `src/globals.f90`, used
  identically by both backends:

  | call | does |
  |---|---|
  | `begin_eddy_sequence()` | cell-label tracer to the identity |
  | `accumulate_eddy(L, M)` | folds one eddy in; sits beside the scalar `triplet_map` calls so the two stay in lockstep |
  | `finalize_eddy_sequence()` | inverts the composed tracer into `destination_cell` — the gather->forward flip, and the *only* place it happens |

  `droplets::move_particles_by_cellmap` then displaces each droplet once, setting
  `position` and `gridcell` together. ODT is the one-eddy case
  (`src/ODT.f90`: `begin` in `odt_turbulence_step`, `accumulate` in
  `implement_eddy`, `finalize` + move at the end); LEM is the many-eddy case
  (`src/LEM.f90:lem_turbulence_step`). `move_particles_in_eddy` was deleted.

  Composition order needs no special handling: applying the maps in sequence to
  the tracer yields the correct composed gather automatically. Written out
  explicitly the nesting runs in reverse (eddies 1,2,3 compose as
  `s1(s2(s3(z)))`), which is easy to get backwards — so the tracer approach
  should not be "optimized" into hand-rolled index arithmetic.

- **Consequence for existing results** — parcel results change (this was already
  expected; slice 1 is a correction, not a refactor). **Chamber results also
  change**, for any eddy with `L >= 9`. Chamber was previously believed
  unaffected; it is not. Any completed chamber run with droplets has backwards
  droplet transport for most of its eddies.
- **Verification** — `test/test_eddy_cellmap_standalone.f90` validates the scheme
  using no production code but `triplet_map`: it carries a labelled scalar field
  through the same eddies and asserts each droplet ends in the parcel it started
  in, over 300 random multi-eddy sequences. `test/test_move_particles.f90` adds a
  positive equivalence check (compose-then-move-once == per-eddy moves with the
  index carried forward) and a stale-index regression guard.
- **Found in** — `afc798f`, branch `feature/lem-empm-port`; found by reading the
  code and by an equivalence test, not from a run.
- **Status** — **fixed on `feature/lem-empm-port`** (uncommitted at time of
  writing; update this line with the commit hash when slice 1 lands).

---

## Unknown

### Physics chain may run twice per iteration when an eddy is accepted

- **Naming hazard (read this first)** — four similar names are involved and are
  easy to conflate. In particular `diffusion_step` and `diffusion_timestep`
  differ by one word and are **different quantities**:

  | name | scope | meaning |
  |---|---|---|
  | `dt` | global | the model timestep; `time = time + dt` each iteration (`CODT.f90:69`) |
  | `delta_time` | global | time since the last physics update, `time - last_time_updated` (`CODT.f90:73`) |
  | `diffusion_step` | global | the backstop *threshold* compared against `delta_time` (`CODT.f90:82`); set to `dt` by both backends (`LEM.f90:101`, `ODT.f90:85-86`) |
  | `diffusion_timestep` | LEM local | the diffusive *stability limit* `0.2*dz^2/D` (`LEM.f90:82`); used only to derive `dt`, `maps_per_event`, `steps_between_events` (`LEM.f90:91-99`) |

  `diffusion_step == diffusion_timestep` only in the
  `diffusion_timestep >= convection_timestep` branch, where `dt` is set from
  `diffusion_timestep` (`LEM.f90:92`). In the other branch `dt` comes from
  `convection_timestep` (`:97`) and the two differ. Note also that
  `turbulence_step` receives both as dummies named `ldt` and `ldelta_time`
  (`LEM.f90:141`, `ODT.f90:576-579`).

- **Symptom** — Both backends set `diffusion_step = dt` and the main loop
  advances `time = time + dt`, so `delta_time` is `dt` on essentially every
  iteration and the diffusion backstop at `CODT.f90:82` fires every step.
  `delta_time` is computed once at `CODT.f90:73` and is **not** recomputed after
  that block, so when an eddy is also accepted the chain at `CODT.f90:96-100`
  re-runs `diffuse_step`, `advance_droplets`, special effects and radiation with
  the full `dt` again.
- **Scope** — In ODT eddy acceptance is stochastic and rare, so this is an
  occasional extra `dt`. In LEM with `steps_between_events = 1` (the
  `diffusion_timestep >= convection_timestep` branch, `LEM.f90:91-94` — the
  branch madison takes) `leddy_accepted` is set unconditionally, so every
  iteration would double-apply.
- **Why this may not be a defect** — it may be the intended design: an accepted
  eddy is a distinct physical event, and re-running the chain after the
  rearrangement lets droplets and scalars respond to the new field
  configuration. This is filed here rather than under Coding because nobody has
  established which reading is correct.
- **How to settle it** — instrument one short parcel run: count `diffuse_step`
  calls per iteration and compare integrated diffusion time against elapsed
  `time`. If integrated physics time is ~2x elapsed time in LEM, it is a defect.
- **Note** — floating point softens it: `(time + dt) - time` is not exactly `dt`
  at large `time`, so the backstop occasionally skips and the next iteration
  carries a larger `delta_time`. Any over-application is therefore irregular
  rather than exactly 2x.
- **Found in** — `afc798f`, branch `feature/lem-empm-port`; found by reading the
  code, not from a run.
- **Status** — open question, not yet classified.

### Parcel-mode `reynolds_number` drives `special_effects` sidewall forcing

- **Symptom** — Not yet observed in a run. `LEM::reynolds_number` is public and
  `initialize.f90:58` passes it to `initialize_special_effects` as the *Rayleigh
  number* argument in parcel mode. It lands in `special_effects.f90:84-90` as
  `Ra`, then propagates to `Nuss`, `velocity_bot` and `tau_sw`. An LEM
  turbulence quantity is therefore driving a chamber-derived sidewall nudging
  parameterization through an argument named for a different dimensionless
  group.
- **Why it matters now** — the slice-2 scale rework changes `Re` substantially
  (4.98 -> 1.98 on the bundled `input/params.nml`, because the smallest eddy is
  now derived rather than taken from the removed `kolmogorov_length_scale`). Any
  parcel run with `do_special_effects = .true.` therefore sees changed sidewall
  forcing as a side effect of a turbulence-scale change. Inert at defaults:
  `do_special_effects = .false.`. Chamber mode is unaffected — it passes a real
  Rayleigh number from `initialize.f90:55`.
- **Evidence** — code reading only, during the slice-2 diffusivity audit. No run
  has exercised `do_special_effects = .true.` in parcel mode.
- **Suspected location** — `initialize.f90:58` (the call site) and
  `special_effects.f90:78-90` (the consumer). Both verified by reading.
- **Open question** — whether passing `Re` here was deliberate (as a generic
  "vigour of convection" proxy) or carried over from the chamber path. If
  deliberate the dummy argument should be renamed and the choice documented; if
  not, parcel mode needs its own sidewall closure or should disable it.
- **Found in** — branch `feature/lem-empm-port`, working tree during slice 2.
- **Status** — open.
