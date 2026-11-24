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

### `n_blob` scales the entrained volume per event, not the blob subdivision

- **Symptom** — `psigma` is the volume of **one** blob, so an entrainment event
  replaces `n_blob * psigma` of the domain, and the interval between events
  scales with `n_blob` to compensate. Sweeping `n_blob` at fixed `psigma`
  therefore does **not** subdivide a fixed entrained volume into smaller
  blobs, which is the natural reading of the parameter name and the reading
  under which such a sweep would vary homogeneous vs inhomogeneous mixing.
  Blob *size* is `psigma` alone and is untouched by `n_blob`.
- **Evidence** — `entrainment.f90 place_blobs`: `blob_size = int(psigma * N)`,
  with no `n_blob` term, then `n_blob` such blobs are placed. The validation
  `psigma * n_blob < 1` (`entrainment.f90:144`) only makes sense if the total
  is `n_blob * psigma`. `compute_dt_entm` gives
  `dt = (n_blob / ent_rate) * (psigma / (1 - psigma)) / |vel|`, so the
  entrained fraction per unit height is `ent_rate * (1 - psigma)` and is
  independent of `n_blob`.
- **Not a port error** — EMPM does the same thing, and CODT matches it
  exactly. `EMPM.f90:1594-1595` is explicit:
  `volume_entm = psigma*volume ! volume of ONE entrained blob` followed by
  `volume_entm = volume_entm * n_blob ! fix for n_blob > 1`. `domain.f90:51`
  sizes each blob `psigma*ngrid` with no `n_blob` term, and `fdtime.f90:36`
  has the identical `delta_z = (n_blob/ent_rate)*(psigma/(1.0-psigma))`.
  Logged as a configuration hazard, not a code defect.
- **Consequence** — At `psigma = 0.1`, `n_blob = 5` replaces **50% of the
  domain in a single event**. Combined with a low entrainment rate the event
  spacing can exceed the whole ascent: in EXP001_sf01_seeding at
  `ent_rate = 0.087` /km, `dz` between events is 1274 m (`n_blob=1`), 2549 m
  (`n_blob=2`) and **6372 m (`n_blob=5`) against a ~4070 m ascent** — so those
  members see zero or one entrainment event, and if one fires it is a
  half-domain shock rather than gradual dilution. Realizations intended as
  equivalent stochastic replicates are then qualitatively different
  experiments at the low-entrainment corner of a sweep.
- **Observed impact** — In EXP001_sf01_seeding, all 7 unseeded runs that
  depleted below 100 cm-3 while under 6000 m carried `n_blob >= 2`, and none
  of the `n_blob = 1` members did. Interacts directly with
  [Runaway collision-coalescence collapse](#runaway-collision-coalescence-collapse-in-parcel-mode).
- **Found in** — `e1c03be` (`v3.0.0`), from the EXP001_sf01_seeding
  configuration; behaviour traced in the source and confirmed against the
  run output attributes.
- **Resolution** — **fixed in 3.1.0 by redefining the parameters** (breaking).
  `psigma` is now the *total* fraction replaced per event and `n_blob` only
  subdivides that fixed volume into evenly sized chunks:
  `place_blobs` sizes chunks as `int(psigma*N)` split `n_blob` ways (the integer
  remainder spread one cell at a time, so the total does not drift with
  `n_blob`), and `compute_dt_entm` no longer carries the `n_blob` factor, so
  event timing depends on `psigma` alone. `n_blob` is now a genuine
  homogeneous/inhomogeneous mixing axis. The old `psigma * n_blob < 1` check is
  replaced by `int(psigma*N) >= n_blob` (every chunk needs a cell).

  This deliberately **departs from EMPM**, which keeps the per-blob reading
  (`EMPM.f90:1594-1595`). The departure is intentional: the CODT parameter names
  are now self-consistent, and holding entrained volume fixed while varying
  blob count is the sweep the science needs.

  Validated with two runs identical but for `n_blob` (1 vs 5; 10 m domain,
  3000 droplets, `ent_rate = 10` /km, 150 m ascent): both fired **13** events at
  **identical heights** (11.08 m spacing vs 11.1 m predicted) and entrained
  **exactly 300 particles per event** in both, confirming an identical replaced
  cell count. Detrained particle totals differ slightly (3765 vs 3825) because
  five scattered chunks intersect a different set of droplets than one large
  blob — the intended effect; that difference is far smaller than the
  event-to-event spread within either run (239-354).
- **Impact on existing results** — `n_blob = 1` runs are bit-identical.
  `n_blob > 1` runs are **not comparable** across the change (an old
  `n_blob = 5, psigma = 0.1` event replaced 50% of the domain; a new one
  replaces 10%). The affected ensembles, including EXP001_sf01_seeding, need
  re-running before their realization spread can be read as stochastic. To
  reproduce old behaviour exactly, set `psigma` to the old `n_blob * psigma`
  with `n_blob = 1`. Output carries the same variable/attribute names with the
  new meaning, so distinguish by `code_version`.
- **Status** — **fixed** in 3.1.0 (breaking parameter redefinition).

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

### Beard fall speed: step at 19 um, flat above 1.07 mm, its own viscosity

Three properties of `fall_speed_scheme = 'beard'` (`particle_types.f90
beard_fall_speed`) that are choices or properties of the formula, not coding
errors. None applies to the default `'stokes'`.

- **Step at 19 um diameter** — Beard's two fits do not join exactly. The speed
  steps *down* by 0.17-0.23% going from the slip-corrected Stokes range to the
  Reynolds-number fit (0.225% at 283 K, 1000 hPa), so in a band about 0.1% wide
  in radius (roughly 9.49-9.51 um) a smaller droplet falls faster than a larger
  one. The collision step chooses collector and survivor by radius, so this
  only changes which pairs close on each other inside the band.
  The reference implementation (`test/cluscoll/reference/sdm_motion.f90`) has
  the same step. `test_fall_speed` pins it (exactly one downward step between
  9.4 and 9.6 um).
- **Flat above 1.07 mm** — Beard's raindrop range needs surface tension and was
  left out. Every drop larger than 1.07 mm diameter gets the 1.07 mm speed
  (4.24 m/s at 283 K, 1000 hPa): such drops fall too slowly, by up to about a
  factor of two for the largest raindrops, and cannot collide with each other
  because their fall speeds are equal.
- **Viscosity differs from the rest of the model** — Beard uses the
  temperature-dependent fit of Pruppacher & Klett (1997), 1.766e-5 kg/m/s at
  283 K. Stokes uses the constant `nu * rho_air`, 1.832e-5 at 283 K and
  1000 hPa. The two formulas therefore do not agree exactly even for the
  smallest droplets.
- **Also** — a droplet of exactly 9.5 um radius takes the Reynolds-number fit,
  where the reference takes Stokes, because the 19 um limit is a
  single-precision literal slightly below 19e-6. `beard_fall_speed` guards only
  against a non-positive radius; no path in `src/` was found that passes it an
  unphysical temperature or pressure.
- **Evidence** — independent audit of the uncommitted Beard work, 2026-11-16:
  scan of 200 001 radii at four (T, p) states; formula, units and coefficients
  checked against the reference files.
- **Found in** — branch `feature/beard-fallspeed`, on top of `bee14d6`.
- **Status** — open (documented behaviour).

### ~~Collection efficiency above 1 is truncated to 1~~

- **Symptom** — the coalescence probability per crossing was
  `p_geo * min(E, 1)`, with `p_geo = pi (r_i + r_j)^2 / A`. The Hall table
  exceeds 1 for near-equal sizes (up to 4.0 at size ratio 1): the collector's
  wake captures droplets outside its geometric swept area. Those captures were
  lost, so coalescence among near-equal droplets was undercounted.
- **Evidence** — `test/cluscoll` (local, untracked) Hall pair cases: pooled
  coalescence z over 48 cases was -1.36 with the truncation.
- **Suspected location** — `src/collision_coalescence.f90`
  `handle_pair_event`: it accepted a geometric collision first and only then
  drew against E. Verified by source read.
- **Found in** — `eed495a` (`v3.2.1`), `main`. Also noted while planning the
  Pinsky kernel (2026-08-27).
- **Status** — **fixed in `264fc09`**: one draw `u` per crossing; `u < E p_geo`
  coalesces, `u < p_geo` without coalescing is a geometric collision. A wake
  capture (coalesced with `u >= p_geo`) is counted and recorded as a collision
  with `flag = 1`, so collision counts follow `max(1, E) K`. `E_max` (set by
  `set_kernel_selector`: max of the Hall table, 1 for Long and unity) bounds
  the draw so E is evaluated only when `u < E_max p_geo`. Pair cases after the
  fix, 4 realizations x 12 cases: 122 415 coalescences / 122 570 expected
  (0.999), 230 960 collisions / 231 007 expected (1.000), max |z| 2.9. Runs with
  `do_coalescence = .true.` are not bit-identical to 3.2.1 (the random-number
  sequence changes). Whether the Hall table values themselves are right (the
  CODT "SKK-corrected" table differs from the SDM reference) is a separate,
  open question.

---

### ~~Droplet growth treats every aerosol as NaCl, whatever the composition table says~~

- **Symptom** — The growth ODE uses NaCl solute constants for every droplet.
  `single_droplet_growth` calls `set_aerosol_properties(1, ...)` with the
  species index hard-coded, and `set_aerosol_properties` then takes `n_ions`,
  molar mass, and the `c7` density correction from its own internal
  `select case` table (1 = NaCl), not from the particle's `solute_type`. The
  composition table read from the aerosol file (`n_ions`, `molar_mass`) never
  reaches the growth model. Only `solute_density` does, indirectly, through
  `solute_gross_mass` and `solute_radius`.
- **Why it bites** — The Köhler critical radius and critical supersaturation
  (`particle_critical_kohler`) *do* use the per-type file values, so in a run
  with a non-NaCl type the activation diagnostic and the growth physics
  disagree about the same droplet. For ammonium sulfate the Raoult coefficient
  `n_ions * Mw / M_s` is 2/58.44 instead of 3/132.14 per kg of solute, about
  51% too large by that ratio. Single-species NaCl runs are unaffected, as long
  as the file's type-1 values match the hard-coded ones. The seeding paths
  (`grown_radius`, the equilibration loop near `droplets.f90:1144`) go through
  the same call.
- **Evidence** — source read only; not observed in a run. `git log -S` shows
  the literal `1` has been there since the initial commit.
- **Suspected location** — `src/droplets.f90:848` (the call) and
  `src/DGM.f90:86` (the internal species table). Confirmed by source read; this
  is the only call to `set_aerosol_properties` in `src/`.
- **Found in** — `8d4d6fb` (`v3.1.2`), `main`; found while tracing the growth
  path for the SoA exploration (`SOA_EXPLORATION.md`).
- **Status** — **fixed in `a3f1a71`**.
  `set_aerosol_properties` takes `n_ions`, `molar_mass` and `solute_density`
  from the droplet's `solute_type` (the aerosol file), the same values the
  activation diagnostic uses; the internal species table is gone. The
  solution-density factor is derived: `c7 = 1 - rho_w/solute_density`
  (solute and water volumes additive, so `rho_sol*V - m_s` is the water mass).
  The former constants were exactly that for 2165 (NaCl) and 1774
  ((NH4)2SO4) kg/m3; fumaric acid borrowed the (NH4)2SO4 value (its density
  ~1635 gives 0.388, not 0.436). `build_aerosol_table` now rejects
  `n_ions < 1`, `molar_mass <= 0` or `solute_density <= 0`, and warns for
  `molar_mass > 1 kg/mol` (likely g/mol).
  - **Output changes.** Non-NaCl material changes substantially ((NH4)2SO4
    Raoult ~1.5x smaller). NaCl runs change slightly: the molar mass is now
    the file's double rather than a single-precision literal, and c7 follows
    the file's density. The bundled `input/aerosol_input.nc` was set to NaCl
    2165 kg/m3 in v3.3.0 (the density the former constant encoded; c7 differs
    from it by 3.2e-8). The reftest aerosol files still give 2163, so there c7
    moves 0.5381 -> 0.5377 (8e-4). Not bit-identical to earlier versions.
  - **Coalescence of mixed material** — the survivor's `solute_type` now also
    sets the merged droplet's growth Raoult term (before, it set only the
    activation diagnostic and dry radius), so the documented "survivor's
    material" approximation (`data_formats.md`, `collision_coalescence.f90`)
    has a larger effect than before when unlike materials coalesce.
  - `test_dgm_composition`: Köhler equilibrium written from the volumes
    (independent of DGM's c7 algebra) is a growth equilibrium for three
    compositions to ~1e-13; (NH4)2SO4 and NaCl equilibria differ; densities
    2165/1774 reproduce the former constants.
  - **A/B reftest** `abbce87` vs `a3f1a71` (jobs 17265266 for the `abbce87`
    NaCl legs, 17266084/17266085 for the rest; np-v3, skl; the `a3f1a71`
    legs ran on Gold 6140, the baselines on Gold 6130, same ISA). Single
    realizations; window means.
    - *NaCl chamber*: S fields differ from t = 1 s, then the turbulent run
      is a different realization. 900-1800 s: every quantity within ±2%
      with mixed signs (S +1.7%, Nact -1.7%, Nun -1.8%, LWC -1.6%, mean r
      0.0%). A 0.08% c7 change cannot cause 2%, so this is the noise floor
      of comparing two single chamber runs.
    - *NaCl parcel*: before the first divergence of Np (an entrainment event
      at t = 56 s) Nact is identical, LWC within 0.15%, mean S within
      0.003 percentage points: the systematic effect is negligible. After
      that the runs differ by where blobs land (see "Parcel reftest droplet
      count is dominated by where entrainment blobs land"); cloud-phase
      peaks agree (max Nact 1000 vs 1012, max LWC 0.141 vs 0.144 g/m3, max
      mean S 0.445% both).
    - *(NH4)2SO4 chamber* (reftest chamber case with only the composition
      rows of the aerosol file changed to n_ions 3, 0.1321395 kg/mol,
      1769 kg/m3): clear signal above that noise floor, 900-1800 s: Nun
      +30% (14.9 -> 19.4), fraction of droplets below 1 um +18%, size
      dispersion +1.8%; Nact -1.5%, S +1.1%, LWC -0.4% (within noise).
      As expected: the correct Raoult term is ~2/3 of the NaCl one, so haze
      equilibrates smaller and needs more S to activate.
  - **Code-equivalence proof** (job 17269755, chamber, np-v3, all Gold 6130):
    old code `abbce87` vs new code `07cab53` with an NaCl aerosol file whose
    values reproduce the old single-precision constants exactly in double
    (`solute_density = 2164.999847728025`, so 1 - 1000/rho_s =
    0.5381062030792236 = old c7; `molar_mass = 0.05844280123710632` = old
    58.4428e-3): **bit-identical** (all 32 profile and 16 particle variables,
    collisions.bin). So the file -> `solute_type` -> growth path computes
    exactly what the internal table did. With a plain 2165 / 0.0584428 file
    (c7 +3.2e-8, molar mass -1.2e-9 vs the old constants) the runs diverge
    from t = 8 s and the 900-1800 s window means agree within 1.4%, inside
    the chamber noise floor. The composition table each run used is printed
    in its log (`07cab53`).

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

### Aerosol species are not self-descriptive — only type 1 has a real name

- **Symptom** — The aerosol input carries one global `aerosol_name` attribute
  for the whole composition table. Type 1 gets that name; every other type is
  named `<aerosol_name>_typeN`. A NaCl + ammonium sulfate file with
  `aerosol_name = "NaCl"` names the sulfate `NaCl_type2`, and nothing in the
  input or output says what type 2 actually is. The species identity lives only
  in its `n_ions` / `molar_mass` / `solute_density` row.
- **Why it bites** — Multi-species runs (several background types via
  `bin_type`, or a seed type) cannot be read back without outside notes. The
  obvious workaround of a combined name (`"NaCl+AmSul"`) makes it worse:
  type 1 is then called `NaCl+AmSul` and type 2 `NaCl+AmSul_type2`.
- **Two silent traps on the way in** —
  - `codt_tools` coerces `aerosol_name` with `str(value)`
    (`codt_tools/case/aerosol.py:617`), so passing a list writes the literal
    string `"['NaCl', 'AmSul']"` with no error.
  - `aerosol_name` is `character(20)` (`particle_types.f90:17`), so a long base
    name plus the `_typeN` suffix is truncated without warning.
- **Suspected location** — `droplets.f90 build_aerosol_table` (the naming loop)
  and the attribute read in `read_aerosol_netcdf`. Confirmed by source read.
  The comment there states the cause: "the v1 schema carries no per-type names."
- **Desired behavior** — Each composition row should carry its own name, e.g. a
  per-type `aerosol_name(aerosol_type)` string variable, written through to the
  output so every species is identifiable from the files alone. This is a
  schema change (`CODT_aerosol_input_v1` → v2) and needs a matching
  `codt_tools` `Aerosol` change; keep the global attribute readable for old
  files.
- **Found in** — `8d4d6fb` (`v3.1.2`), `main`; found by source read while
  configuring a two-species injection. Not observed in a run.
- **Status** — open.

## Coding

### Collision scheduler queues pairs that are moving apart, so collisions are undercounted

- **Symptom** — The number of collisions depends on the length of the
  `collision_coalescence_step` window, and is far too low unless the window is
  very short. Gravity-only box (`test/cluscoll`, smoke case: 5001 droplets,
  8-24 um, 20 s, 163.2 collisions expected from the gravitational kernel):

  | window (s) | counted | ratio |
  |---|---|---|
  | 1.0 | 7 | 0.04 |
  | 0.1 | 24 | 0.15 |
  | 0.01 | 89 | 0.55 |
  | 0.001 | 154 | 0.94 |

- **Cause** — The adjacency list runs bottom to top (`head` is the lowest
  droplet, `next(i)` is the droplet above `i`), and `w_fall` is a positive
  downward speed. `push_pair_for_i` schedules pair `(i, j = next(i))` only when
  `vrel = w_fall(i) - w_fall(j) > 0`, that is, when the **lower** droplet falls
  faster. That pair is separating. The pair that actually meets, upper droplet
  faster, is never queued. `handle_pair_event` has the same sign in its guard.
  The event time `dz / vrel` is therefore the gap divided by the rate at which
  the gap is growing.
- **Why short windows hide it** — Positions are advanced correctly and the list
  is re-sorted at the start of every window. With a short window nearly all
  true crossings are between adjacent droplets, and the wrongly queued
  separating pairs have the same statistics of gap and speed difference as the
  approaching ones, so the count comes out about right for the wrong pairs.
  With a long window most true crossings need the cascade of rescheduling after
  each swap, and the wrong sign breaks the cascade (after a swap the new pair
  fails the test).
- **Evidence** — Source read, then an A/B run on a scratch export of `8d4d6fb`
  (`src/` in the repository untouched). Flipping the sign in both places, and
  making `j` the survivor, removes the window dependence:

  | window (s) | counted | ratio | z |
  |---|---|---|---|
  | 1.0 | 170 | 1.04 | 0.5 |
  | 0.1 | 176 | 1.08 | 1.0 |
  | 0.01 | 171 | 1.05 | 0.6 |
  | 0.001 | 178 | 1.09 | 1.2 |

  The four patched runs share one RNG seed, so they are not independent. All
  four are above 1 but each is within 1.2 sigma; a small residual bias is not
  ruled out and needs more seeds or a larger case.
- **Consequence** — In the model the window is the main-loop time step
  (`droplets.f90:139`), so the collision rate of every run with
  `do_collisions = .true.` depends on the time step and is biased low. The
  size of the bias in a production run has not been measured. The wrong pairs
  also collide: the recorded collisions are between droplets that were not
  meeting, so pair statistics in `collisions.bin` (who collected whom) are
  suspect, not only the count. Results that rest on collision rates, including
  [Runaway collision-coalescence collapse](#runaway-collision-coalescence-collapse-in-parcel-mode),
  need to be re-examined after a fix.
- **Related, not yet tested** —
  - The list is not periodic. In parcel mode a droplet that wraps through the
    bottom within a window cannot meet the droplets it re-enters above; the
    `dz < 0 -> dz + H` branch in `push_pair_for_i` cannot fire on a sorted list.
    **Measured 2026-11-16** (after the fix, `test/cluscoll` driver, one
    10 + 40 um pair, Beard fall speeds, `dw` = 16.3 cm/s, 1 m column, 20 000
    collisions expected): the fraction of crossings lost per call is
    `dw * window / (2 H)`. Windows of 0.50 and 0.33 s lost 4.2% and 2.7%
    (4.1% and 2.7% predicted); 0.091 and 0.010 s lost nothing measurable
    (noise 0.7%). The loss resets every call because the list is re-sorted from
    wrapped positions. Chamber mode is not periodic and is not affected. For a
    1 m parcel domain and a 0.01 s step the loss is about 0.05%.
    **Fixed 2026-11-16** on `feature/beard-fallspeed`: in parcel mode the
    neighbour list is now a ring (the highest droplet's upper neighbour is the
    lowest, one domain height up) and the gap is taken modulo `H`. The same
    case then gives 20 063 collisions at every window length, including one
    25 s call in which droplets lap the column up to four times.
    `test_collision_coalescence` gained a boundary-pair case and a two-droplet
    lapping case. Chamber mode keeps the two-ended list and is unchanged.
  - The comment in `handle_pair_event` says the survivor `i` is the larger
    droplet. With the current sign `i` is the faster one of a separating pair;
    with the sign corrected the faster droplet is `j`, so `keep`/`kill` must
    change with it.
- **A sign flip alone is not a fix** — three more places in `handle_pair_event`
  assume `i` is the larger droplet and must change in the same commit:
  - `collection_efficiency_E(r_i, r_j)` takes (collector, collectee) and does
    not order its arguments. Left as is after the flip it is called with
    (small, large). Measured on the scratch copy: coalescences fell to about a
    quarter of sum(E) (4 of 16.1-18.9 expected, z = -3.3 to -3.7; 7 of 17.8 at
    the 0.1 s window). With the arguments swapped: 17/17.1, 15/18.3, 18/18.8,
    20/19.0, |z| <= 0.85. The runs in the table above had this error; it does
    not affect the collision count.
  - `keep = i; kill = j` must become `keep = j; kill = i`, or the smaller
    droplet's id, category and `n_coalescences` survive.
  - The two `write_collision` calls for a collision without coalescence pass
    `i` then `j`. `docs/data_formats.md` names the fields `id_keep`, `id_kill`,
    `r_keep`, `r_kill`, so `j` must come first to keep "larger droplet first"
    in every record.
- **Suspected location** — `src/collision_coalescence.f90:462`
  (`push_pair_for_i`) and `:216` (`handle_pair_event`). Verified by source read
  and by the A/B run.
- **Found in** — `8d4d6fb` (`v3.1.2`), `main`; from the `git_commit` attribute
  of the `test/cluscoll` output. The line dates from the commit that added the
  module, `ba68801` (2026-04-08).
- **Resolution** — fixed in four commits, one per change above: `d62c592`
  (closing speed in `push_pair_for_i` and `handle_pair_event`), `8a1474f`
  (survivor is `j`), `8c10a96` (efficiency argument order), `661e47a` (record
  order for collisions without coalescence). The first three are each wrong
  without the others; do not cherry-pick one alone.

  Verified with the `test/cluscoll` smoke case, three RNG seeds at four window
  lengths. Collisions counted / expected from the gravitational kernel:

  | window (s) | seed 12345 | seed 777 | seed 2024 |
  |---|---|---|---|
  | 1.0 | 168 / 163.3 | 154 / 161.8 | 154 / 162.0 |
  | 0.1 | 171 / 163.4 | 160 / 161.8 | 160 / 162.0 |
  | 0.01 | 173 / 163.4 | 159 / 161.7 | 163 / 161.8 |
  | 0.001 | 172 / 163.3 | 179 / 161.9 | 154 / 161.9 |

  All twelve pass the kernel check (|z| <= 1.4) and the efficiency check
  (|z| <= 1.5). Per-seed totals over the four windows are 1.047, 1.007 and
  0.974 of expected, and all twelve runs together are 1967 / 1948.3 =
  1.010 +/- 0.023, so the 1.04-1.09 seen with a single seed is not an
  overcount. This table shows window independence only to about +/- 5% per
  window length; the exact demonstration came with the periodic ring (below).

  `test/test_collision_coalescence.f90` had the error built in: its two-droplet
  case put the large droplet below the small one and expected them to merge.
  The geometry is corrected, and the old geometry is kept as a new subtest that
  must produce no collision (25 of 25 checks pass).
- **Impact on existing results** — Every run with `do_collisions = .true.` made
  before this fix undercounts collisions by an amount that depends on its time
  step, and its `collisions.bin` records pairs that were not meeting. Such runs
  are not comparable with runs after the fix; distinguish by `code_version`.
  The non-periodic list (first item under "Related") is still open.
- **Status** — **fixed** in `d62c592`..`661e47a`.

### Collision step: open items from the audit of the scheduler fix

Found by an independent audit of `d62c592`..`bee14d6` on 2026-11-16, which
checked the event logic with a line-by-line Python port (not by running the
Fortran). None changes collision counts measurably at realistic settings.

- **The faster droplet is assumed to be the larger one** — `handle_pair_event`
  takes the upper, faster droplet `j` as collector and survivor. Fall speed
  also depends on the droplet's own `virt_temp` (Stokes speed goes as
  `virt_temp * r^2`), so `j` can be faster yet smaller when the radii are
  within about 0.17% per kelvin of `virt_temp` difference; with
  `fall_speed_scheme = 'beard'` the step at 19 um diameter is a second route.
  Then `collection_efficiency_E` is called with (smaller, larger), the smaller
  droplet's id and category survive, and the `collisions.bin` record lists the
  smaller droplet first. The efficiency error is small because the radii are
  nearly equal (up to 1-2% of E for 40-70 um collectors under Hall). Not seen
  in 47 276 records of the `test/cluscoll` pair test, which cannot reach
  either route. **Status** — **fixed** in the commit after `9931d16`
  ("Choose collector and survivor by radius"): `handle_pair_event` now takes
  the larger radius as collector and survivor (equal radii fall back to the
  faster droplet), for the efficiency call, the merge and the record order.
  `test_collision_coalescence` has a case with a smaller, warmer, faster
  droplet.
- **Pair events are not invalidated when a droplet's speed changes** — after a
  coalescence the survivor falls faster, but a pair event already queued for
  it keeps the time computed from its old speed, passes the adjacency check,
  and can fire while the two droplets are still apart. The pair is then
  tested (and merged or swapped) early; the count is unchanged to first order
  but the event time and position are wrong, and `zcur(keep) = zcur(i)` moves
  the survivor to the other droplet's position. In an extreme port case
  (collision probability 0.25, E = 1) about 1-2% of pair events fired early;
  at realistic probabilities (about 1e-4) coalescences are rare relative to
  crossings and the effect is negligible. The mechanism predates the fix.
  An early swap also leaves the list out of position order, so the next gap
  for that pair is read as negative and wrapped to nearly a full lap
  (`modulo(dz, H)` in parcel mode, `dz + H` in chamber mode), losing that
  meeting. **Status** — **fixed** after `16fc44f`: each droplet carries a
  fall-speed version, bumped when a coalescence changes its speed; each pair
  event records both versions when it is queued, and `is_pair_valid` drops
  it if either has changed. This also covers a pair that stops being
  adjacent and later becomes adjacent again. `test_collision_coalescence`
  has a three-droplet case that fails without the check (2 coalescences
  instead of 1). The `test/cluscoll` Hall pair cases are unchanged at a
  0.1 s window (no stale event fired at realistic probabilities). Runs with
  coalescence are not bit-identical to earlier versions: stale events used
  to consume random draws.
- **Equal-size droplets in different cells collide** — the same `virt_temp`
  dependence gives equal radii a non-zero closing speed, and the Hall table's
  size-ratio-1 column reaches 4.0. Predates the fix. **Status** — accepted,
  not fixing: radii are double precision and grow independently, so exactly
  equal radii are not expected in practice.
- **`collisions.bin` `position` is not wrapped** — it is the working position
  of the lower droplet `i`, which starts each window in [0, `H`) and only moves
  down, so in parcel mode it can be slightly negative (at most `w_fall(i)`
  times the window). It never reaches `H`: the earlier note here that it could
  be "offset by a multiple of `H`" was wrong. **Status** — **fixed** after
  `a03d461`: `write_collision` wraps the value into [0, `H`) when the column is
  periodic. Files from 3.2.0 and earlier still need `modulo(position, H)`.
- **`test/scripts/test_cc_invariants.py` reads 48-byte records** — the record
  has been 49 bytes since the coalesced flag was added, so the script misreads
  current files. It also imported classes that CODT_tools 1.0 no longer has,
  and treated every record as a coalescence. **Status** — **fixed** after
  `2e3e3aa`: it reads the log with CODT_tools `Simulation.load_collisions`,
  checks coalescence records by the flag, and skips the per-particle check
  when there are no trajectory records.
- **Unit test gaps** — `test_collision_coalescence` would catch the sign error
  and a reverted survivor, but not a reverted efficiency argument order (it
  uses the `unity` kernel) or record order (it does not write
  `collisions.bin`); those are covered only by the `test/cluscoll` driver. It
  has no chamber-mode or three-droplet case. The test fails at the
  intermediate commits `d62c592`..`661e47a` and passes from `c87b7e0`, which
  matters for `git bisect`. **Status** — open.
- **Found in** — `bee14d6` (`v3.1.2-6-gbee14d6`), branch
  `fix/cc-pair-scheduling`.

### Parcel reftest droplet count is dominated by where entrainment blobs land

- **Symptom** — Two parcel reftest runs that differ in anything touching the
  random stream end with very different droplet counts, for reasons unrelated
  to the change under test. Comparing `b7e6f3d` (baseline) with `bee14d6` (the
  collision scheduler fix): final `Np` 852 vs 1217. Injection, seeding and
  entrainment (15 events of 100, same times) are identical; coalescence removed
  1 vs 3 droplets; the whole gap is detrainment, 1797 vs 1430. One event at
  t = 278 s detrained 651 in the baseline and 212 in the new run.
- **Cause** — `random_entrainment = .false.` fixes only the *timing* of
  entrainment events. The blob *position* is drawn from the global random
  stream in `place_blobs` (`src/entrainment.f90`), and the collision step draws
  from the same stream for every pair event it processes
  (`collision_coalescence.f90`, the `random_number` call in
  `handle_pair_event`). Droplets are strongly clustered during the ascent (607
  of 1038 in one tenth of the domain at t = 275 s in the baseline), so the
  number detrained per event ranges from 0 to 651 depending on where the blob
  falls. The cause of the clustering was not investigated.
- **Consequence** — The parcel reftest cannot isolate the effect of a change on
  droplet number unless the change leaves the random stream untouched. Fields
  already differ at t = 1 s. To isolate a collision change, run with
  `do_entrainment = .false.`, or compare ensembles of seeds; separate random
  streams for collisions and blob placement would remove the coupling but is a
  code change.
- **Evidence** — `devCODT/reftest/parcel/`, `reftest.nc` vs `reftest_REF.nc`
  and the two `collisions.bin` files, SLURM job 17126431; budget from the
  per-interval `budget_n_*` variables summed over time (both budgets close
  exactly). Analysis done by a subagent on 2026-11-16; the numbers were not
  independently rechecked.
- **Found in** — `bee14d6` (`v3.1.2-6-gbee14d6`) against `b7e6f3d`.
- **Status** — open.

### `test_seeding` fails and then segfaults

- **Symptom** — `fpm test test_seeding` reports `FAIL: n seed bins` (expected 2,
  got 0) and then dies with SIGSEGV (exit 139). The reader prints "seeding
  input detected but do_seeding is false; the seed group is ignored" first.
- **Suspected location** — `test/test_seeding.f90`: it reads its own seeded
  aerosol file without turning `do_seeding` on, so the seed group is skipped
  and the later checks index seed arrays that were never filled. *Unverified*
  beyond the log message; not debugged.
- **Evidence** — `release` build, gcc 8.5.0, 2026-11-16. Seen while running the
  unit tests for the collision scheduler fix. The test touches no collision
  code, so it is taken to predate that fix, but it was not rerun on `main`.
- **Found in** — branch `fix/cc-pair-scheduling` at `661e47a` (`v3.1.2` plus the
  collision fix).
- **Status** — open.

### `collisions.bin` `time` is not simulation time

- **Symptom** — The per-event `time` field spans only `7.6e-08 -> 1.20e-02` s in
  a run whose simulation time covers `0.0 -> 1954.3` s, and it is **not
  monotonic** across the event stream (first five events: 2.95e-4, 3.64e-3,
  1.17e-3, 1.04e-2, 2.22e-3). The magnitudes look like a time *within* the
  current timestep; no step index or absolute time is written alongside it.
- **Consequence** — Per-event absolute time is unrecoverable from the binary,
  and with it per-event parcel altitude. Any analysis that wants "when did
  collision-coalescence start" or "at what height did this collision happen"
  cannot use this field. Silent: the values are plausible small floats, so a
  naive onset-time or time-height plot produces a figure rather than an error.
- **Workaround** — Use `N_collisions(time)` / `N_coalescences(time)` from the
  main netCDF; they are per write interval, correctly timed, and their totals
  match the binary exactly (25025 / 3253 on the run checked). Map to altitude
  through `parcel_height_env` as usual. The binary remains authoritative for
  per-event *sizes* (`r_keep`, `r_kill`, `r_after`, `flag`), which check out:
  `r_after = 0` for every non-coalescence, and `r_after^3 = r_keep^3 + r_kill^3`
  to 2e-15 relative error for every coalescence.
- **Recovering per-event time anyway** — the binary is a Fortran stream written
  in event order, so the events partition into write intervals by the
  cumulative `N_collisions`: the k-th block of `N_collisions[k]` events belongs
  to interval k. That restores per-event time, altitude and parcel state
  without the `time` field. Validated on the 250-run EXP001_sf01_seeding
  ensemble: under this mapping **zero** of 32043 events carry a collector
  radius exceeding the largest drop present in the DSD at the assigned time,
  versus ~37-41% for a shuffled control, and mean collector radius rises
  monotonically with interval (r = 0.92). Used by
  `projects/SF01_Seeding/analyze.py`. Note this assumes the writer never
  buffers events out of order across intervals.
- **Cause** — confirmed by source read. `current_time` is a *local* in
  `collision_coalescence_step`, reset to `0.0` on every call and thereafter only
  assigned from the event heap key, which is an event time *within* the window
  `[0, ldt]`. It was passed straight to `write_collision` as `time_s`. The module
  never imported the global `time`, so the accumulated clock was simply
  unavailable to the writer — the second of the two options above. The magnitudes
  are bounded by `ldt`; the non-monotonicity is the per-step clock resetting, so
  the stream is a sawtooth of within-window times.
- **Fix** — `collision_coalescence_step` now takes `ltime` (absolute time at the
  window END; the main loop advances `time` before physics, `CODT.f90:69`) and
  stamps events with `ltime - ldt + current_time`. The field was already
  `real(dp)`, so **the binary layout is unchanged** and old readers parse new
  files — only the meaning of the field changed. Validated on a 100 s parcel
  ascent with large aerosol: 499 events spanning `0.0100 -> 76.54` s, strictly
  monotonic (0 backward steps), totals still matching `N_collisions`, volume
  conservation unchanged at 2.2e-15.
- **Residual error** (before the double-physics fix on `fix/double-physics`) —
  accuracy was capped by the separate "Physics chain ran twice per iteration
  when an eddy was accepted" issue: `delta_time` is not
  recomputed between the backstop and eddy branches (`CODT.f90:73` vs `:88`), so
  both CC calls receive the same `ltime`/`ldt` and the second call's events are
  backdated by up to one `delta_time`. On the validation run this showed up as 4
  of 101 write intervals having events swapped with an adjacent interval (totals
  preserved). Judged rare enough to accept; the interval-block reconstruction
  below is unaffected and remains exact.
- **Found in** — `e1c03be` (`v3.0.0`), from the `git_commit` attribute of the
  EXP001_sf01_seeding runs.
- **Status** — **fixed** (offset added). Pre-fix files still need the
  interval-block reconstruction above; detect them from the run's
  `git_commit` / `code_version` attribute, or by testing whether any `time`
  in the binary exceeds one `delta_time`. The fix is `c67809c` (3.0.1).
  Confirmed still correct on 2026-11-16: in the parcel reftest files from
  `b7e6f3d` and `bee14d6`, record times span the run (18-537 s) and agree with
  `N_collisions(time)` in the netCDF.

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
- **Cause** — `parcel_height_env` starts at `0.0` (`parcel.f90:61`) and was
  computed only in `apply_adiabatic_forcing`, the per-step routine. The first
  output record is written in `initialize_simulation` (`initialize.f90:71`)
  before any step, so it held the default.
- **Status** — **fixed in `a9f98c7`**:
  `initialize_parcel` now evaluates it for the initial pressure as soon as
  `write_height_env` is set. Only the t = 0 value of `parcel_height_env`
  changes (hydrostatic runs with a sounding); `[1:]` slicing is no longer
  needed for output from this version on. `test_parcel_v3_ascent` has a check
  that fails without the fix.

### `rdx` out-of-bounds access in the tridiagonal solver

- **Symptom** — The tridiagonal solve indexes `rdx` as `rdx(N)` where it should
  be `rdx(N+1)`. Since the commit that introduced the ghost point, the array
  needs the extra element, so the solver reads one element short of what the
  grid now requires.
- **Consequence** — Latent in all versions since the ghost-point change. Will
  not necessarily crash; expect subtly wrong diffusion at the domain edge.
  Confirmed to be an out-of-bounds **write**, not just a short read. `rdx(N)`
  is an automatic array sized by the run-time `N`, which gfortran allocates on
  the **heap** (no `-fstack-arrays` in the profiles), so in a release build the
  write lands just past that allocation (see Status for what that means).
- **Location** — *pinned* (`e1c03be`): `src/ODT.f90:556` declares
  `real(dp) :: b, rdx(N)` while the loop at `:561-562` runs `do k = 2, N+1` and
  assigns `rdx(k)`. Reproduced under the `debug` profile from a chamber run on
  the bundled `input/params.nml`:
  `Index '2001' of dimension 1 of array 'rdx' above upper bound of 2000`,
  via `tridiagonal` (`:562`) <- `diffuse_scalar` (`:544`) <- `diffusion` (`:142`)
  <- `odt_diffuse_step` (`:588`) <- `run_simulation` (`CODT.f90:83`). Because it
  aborts on the first diffusion call, **any chamber run under `debug` currently
  dies immediately** — which is why this is easy to miss under `release`.
- **Found in** — reported from an earlier debugging session on `main`; the
  originating commit/worktree was not recorded. The defect was introduced by the
  ghost-point work, `9a454ff` ("Added ghost points 0:N+1 to nondimensional
  arrays") and `0e070e2` ("Added (ghost) point to top of nondimensional scalar
  arrays"), so it is present on every commit from `0e070e2` forward — including
  current `527cd7f`. Not tied to a specific simulation run.
- **Status** — **fixed in `4b80979`**: `rdx` is
  declared `rdx(N+1)`. Checked against the caller: `diffuse_scalar` builds an
  (N+1)-point system (`l`, `d`, `r`, `xsc` are all `N+1`), the forward loop
  writes `rdx(2:N+1)` and the back substitution reads `rdx(N+1)` for `ui(N)`,
  the cell below the top boundary. Measured: A/B reftest `b1c3a0d` vs
  `4b80979` (job 17264748, np-v3, Gold 6130, N = 2000) is bit-identical in
  chamber and parcel mode, so past release output at that N is unaffected.
  *Inferred, not observed* (from glibc malloc chunk rounding, per the
  2026-11-19 audit; the commit message's "onto the stack" was wrong): for even
  N the 8N-byte request gets 8 bytes of rounding slack, so the stray element
  lands in slack and reads back intact; for odd N there is no slack and the
  write overwrites the next heap chunk's header, which could corrupt the heap
  or crash at `free`. Runs with an odd `N` before this fix deserve suspicion.
  The fix removes the out-of-bounds write and lets chamber runs proceed under
  `debug`.

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
- **Status** — **resolved in `a760c51`** (LEM/EMPM slice 1, released in
  v3.0.0). The routine was
  replaced by the composed-map path described in the next entry, which touches
  each droplet once per turbulence step instead of once per map: madison ent1
  went **48.1 s -> 4.0 s (~12x)**, and the mover dropped to 0.2% of runtime with
  the hotspot moving onto DGM/ODE physics.
- **Correction to an earlier claim in this entry** — a previous revision stated
  *"Performance only — no correctness defect in either version."* **That was
  wrong.** The routine had two independent correctness defects, documented in
  the next entry. The performance framing is what kept them unexamined for so
  long: the routine was read as bookkeeping, so nobody checked its arithmetic.

### `ros3_integrate` heap-allocates its work arrays on every droplet call

- **Symptom** — `malloc`/`free` are **29% of parcel runtime and 13% of chamber
  runtime** (`perf`, `np-v3` build; 27% and 12% on the gcc 8.5 `release`
  build). gprof never showed this because libc is not instrumented.
- **Cause** — `ros3_integrate` declares ten local arrays sized by its `nvar`
  dummy argument (`J_mat`, `W_matrix`, `f_current`, `f_stage2`, `k1`-`k3`,
  `y_trial`, `y_err`, `y_stage`). gfortran puts runtime-sized local arrays on
  the heap unless `-fstack-arrays` is set, so each of the ~4.3e8 droplet calls
  in the parcel profiling case allocates and frees them. `nvar` is always 3.
  Confirmed by `perf --call-graph dwarf`: the `malloc` samples sit directly
  under `ros3_integrate`.
- **Fix direction** — fixed-size (`parameter`) work arrays.
- **Evidence** — `devCODT/profiling/archive/8d4d6fb/soa_step0/`
  (`npv3_*`, `rel8_*`, `stacks_*`); notchpeak Xeon Gold 6130, 2026-11-13.
  Summary in `SOA_RESULTS.md`.
- **Suspected location** — `src/ode_integrators.f90:233`. Verified by source
  read and call stacks.
- **Found in** — `8d4d6fb` (`v3.1.2`), built from a `git archive` export
  (binaries report `3.1.2-src`).
- **Status** — **fixed in `507d660`** (cherry-picked from `84411dd` on
  `explore/soa-particles`). The ten work arrays
  are sized by the constant `ROS3_NVAR = 3`; `ros3_integrate` returns
  `ierr = -2` for any other `nvar`. Measured on the SoA branch (np-v3, Gold
  6130): bit-identical, user CPU 590.9 -> 389.9 s parcel (1.52x) and
  94.7 -> 80.2 s chamber (1.18x). On `main`, together with `deaf0a6`,
  `fbde630` and `b1c3a0d` (the other bit-identical DGM savings): A/B reftest
  v3.2.2 (`cb0c5a1`, job 17260660) vs `b1c3a0d` (job 17264748), np-v3, Gold
  6130, is bit-identical in chamber and parcel mode; wall time 328 -> 259 s
  chamber (1.27x) and 1299 -> 760 s parcel (1.71x), single runs per leg.

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
- **Status** — **fixed in `a760c51`** (LEM/EMPM slice 1, "Unify droplet eddy
  transport for LEM and ODT; fix map direction"; released in v3.0.0).

### ~~DGM growth equations are written out three times~~

- **Symptom** — `growth_rhs`, `growth_jacobian` and `growth_rhs_jac` in
  `src/DGM.f90` each compute the same thermodynamic intermediates (Lv, K, Dv,
  `esat`, jump lengths, kinetic factors, solution density, Köhler terms,
  diffusion denominator), and the two Jacobian routines repeat the same
  derivative code. Any physics change has to be made three times. Adding
  `do_kinetic_correction` showed the risk: the first draft had two errors in
  `growth_rhs_jac` only (kinetic factors set to 0 instead of 1, and
  `d_vvp_dr` assigned where `d_vvp_dt` was meant). `test_kinetic_correction`
  now checks all three routines against finite differences.
- **Cause** — `growth_rhs_jac` was added in `16a6cd4` (2026-05-22) as a fused
  copy of the other two, so ROS3 stage 1 computes the intermediates once.
  `growth_rhs` is still needed for the RHS-only stages. `growth_jacobian` was
  left behind: no model code calls it; only `test/test_jac_magnitudes.f90`
  (unregistered) and `test/test_kinetic_correction.f90` use it.
- **Proposed fix** — (1) delete `growth_jacobian` and point its tests at
  `growth_rhs_jac`. (2) Move the shared code into one new `pure` module
  procedure in `DGM` that both remaining routines call. It must stay in the
  same module so gfortran can still inline it on this hot path. Keep the
  operation order so output stays bit-identical (check with an A/B reftest),
  and compare runtimes.
- **Found in** — `cb0c5a1` (`v3.2.2`), branch
  `feature/kinetic-correction-switch`; found by source read.
- **Status** — **fixed** on `main` (the commit that adds this line).
  `growth_jacobian` is deleted. `growth_terms` (private, pure) computes the
  shared terms once into a `growth_terms_t`, and `growth_tendencies` builds
  dydt; `growth_rhs` and `growth_rhs_jac` both call them, and the derivative
  code in `growth_rhs_jac` is unchanged (it reads the shared terms through
  `associate`). Checked bit-identical before/after on 542 states (grid of
  radius, S, T, dry size, kinetic on/off, plus the evaporation limiter and the
  guarded state) under the debug and release profiles, and by A/B reftest
  `29c6b95` vs `abbce87` (job 17265266, np-v3, Gold 6130) in chamber and
  parcel mode. gfortran does **not** inline `growth_terms` (release: called
  from 3 sites, `growth_tendencies` inlined; np-v3: neither inlined), but the
  extra call costs nothing measurable: 256 -> 252 s chamber, 755 -> 765 s
  parcel (single runs). Fallback if that ever changes: a textual `include` of
  the shared block in both routines.
  `test_kinetic_correction` now checks that `growth_rhs_jac`'s dydt equals
  `growth_rhs` bit for bit.

### A failed growth RHS silently truncates the droplet's growth step

- **Symptom** — `growth_terms` returns `ierr = 1` for an unphysical state
  (`qv < 0` or T outside 193-320 K). `ros3_integrate` returns that code
  immediately (`ierr = rhs_err; return`), including when it comes from a
  *trial* stage state `y + k1`, not an accepted one. `integrate_ODE` aborts
  only for `ierr < 0`, and the model's only call (`single_droplet_growth`,
  `droplets.f90:~863`) passes no `stat`. So a positive code ends that
  droplet's integration at its last accepted state with t < t_end, with no
  message: the droplet gets a truncated growth step and the budgets book the
  partial change.
- **Plausible trigger** — strong evaporation in a very dry cell (e.g. a dry
  entrained blob), where a stage-2 trial state overshoots to `qv < 0`.
  Not observed in a run; frequency unknown.
- **Suspected location** — `src/DGM.f90` `integrate_ODE` (the `ierr < 0`
  test) and `src/ode_integrators.f90` `ros3_integrate` (stage RHS failures
  returned rather than treated as a rejected step).
- **Fix options** — (a) in `ros3_integrate`, treat a stage RHS failure as a
  rejected step (shrink `h`, retry); (b) in `integrate_ODE`, abort on any
  `ierr /= 0`; (c) at least count and report occurrences. (a) changes results
  only where the truncation happens today; decide before changing.
- **Found in** — 2026-11-19 audit of `main` at `abbce87`
  (`v3.2.2-11-gabbce87`); dates from `5dffb4c`. Source read only.
- **Status** — open, now **reported** (the commit that adds this line):
  `integrate_ODE` calls `warn_rhs_failure` for `ierr > 0`, which writes a
  `WARNING: droplet growth step cut short` line to stdout (the run's log) with
  the code, the requested interval and the state the droplet keeps; the first
  10 occurrences in full, then one line every 1000 with the running count.
  Behaviour and output are unchanged. Next step: check production logs for
  the warning before choosing a fix.
- **Related minor notes from the same audit** (no action taken):
  `growth_rhs_jac` column 1 is non-zero when `y(1) < r_floor` although the
  RHS does not depend on r there (`max(y(1), r_floor)`), costing only W-method
  efficiency; `n_ions` is an integer in the aerosol schema, so a fractional
  effective van't Hoff factor cannot be entered.

### ~~Physics chain ran twice per iteration when an eddy was accepted~~

- **Naming hazard (read this first)** — four similar names are involved and are
  easy to conflate. In particular `diffusion_step` and `diffusion_timestep`
  differ by one word and are **different quantities**:

  | name | scope | meaning |
  |---|---|---|
  | `dt` | global | the model timestep; `time = time + dt` each iteration (`CODT.f90`) |
  | `delta_time` | global | time since the last physics update, `time - last_time_updated` |
  | `diffusion_step` | global | the backstop *threshold* compared against `delta_time`; LEM sets it to `dt` (`LEM.f90`), ODT to `H^2/(nu N^2)` and starts `dt` there (`ODT.f90`) |
  | `diffusion_timestep` | LEM local | the diffusive *stability limit* `0.2*dz^2/D`; used only to derive `dt`, `maps_per_event`, `steps_between_events` |

- **Symptom** — the main loop had two copies of the physics chain (diffuse,
  droplets, special effects, radiation, sync): a diffusion backstop
  (`delta_time >= diffusion_step`) and an eddy branch (`eddy_accepted`).
  `delta_time` was computed once per iteration. The backstop reset
  `last_time_updated` but not `delta_time`, so when both fired in one
  iteration the eddy branch re-ran the whole chain over the **same**
  `delta_time`. Droplet growth, collision-coalescence, diffusion and radiation
  then advanced `2*delta_time` while adiabatic forcing and entrainment
  advanced `dt`.
- **Why it was erratic** — in LEM with `steps_between_events = 1` an eddy is
  accepted every iteration and `diffusion_step = dt`, so whether the backstop
  *also* fired was decided by floating-point rounding of `(time + dt) - time`
  against `dt`. The rounding is fixed within each binade of `time`
  ([2^k, 2^(k+1)), where the spacing of doubles is constant), so the doubling
  is **all-or-nothing per binade**: either every iteration in that range ran
  the chain twice or none did. For the reftest `dt` = 2.1286e-3 s the doubled
  ranges are 8-16 s, 64-128 s and 256-512 s; 0-8, 16-64, 128-256 and >= 512 s
  are clean. Which binades double depends only on `dt`, so any past run's
  doubled time ranges can be recomputed from its `dt` without rerunning (a
  Python replica of the loop reproduces the Fortran counts exactly). Overall
  physics/elapsed ratio: 1.06-1.83 for `dt` between 1e-3 and 1e-2 s over a
  1500 s run (2.13e-3 s: 1.22 over 1500 s, 1.44 over 100 s). Ensemble members
  with different `dt` (domain size, dissipation rate) carried different biases,
  in different phases of the run.
- **Measured** — counters in the main loop (diagnostic worktree), reftest
  configurations, `release` build on kingspeak:

  | case | iterations | backstop | eddy chains | both fired | physics / elapsed |
  |---|---|---|---|---|---|
  | parcel (LEM), 100 s, `dt` = 2.13e-3 s | 46 979 | 20 795 | 46 979 | 20 795 (44%) | **1.442645** (144.27 s / 100.00 s) |
  | chamber (ODT), 1800 s, `diffusion_step` = 1.68e-2 s | 197 166 661 | 65 290 | 102 394 | 39 | 1.000359 (+0.65 s) |

  In ODT the eddy acceptance rate is low (~5e-4 per iteration) and adaptive
  `dt` is ~1e-5 s, so coincidences are rare and chamber results are
  effectively unaffected. The parcel case is badly affected: supersaturation
  is depleted by too much condensation per unit cooling, and
  collision-coalescence runs ~1.2-1.8x too long.
- **Related consequences** — the doubled call is also the "adjacent-interval
  swap" noted under `collisions.bin` `time` (the second CC call in an iteration
  reused the same `ltime`/`ldt`, backdating its events) and the -1/0
  oscillation in the `collisions.bin` short-records entry. Both should
  disappear with the fix (not re-checked).
- **Found in** — `afc798f`, branch `feature/lem-empm-port` (code read);
  quantified 2026-11-24 on `74ef16c` (`v3.3.0`).
- **Status** — **fixed** on branch `fix/double-physics` (the commit that
  strikes this entry). One physics block runs after the turbulence step when
  `eddy_accepted .or. delta_time >= diffusion_step`, so the chain runs at most
  once per iteration over `delta_time`. `eddy_accepted` is now reset to
  `.false.` every iteration (it was uninitialized when `do_turbulence` is off).
  Re-measured with the same counters: parcel physics/elapsed **1.000000**
  (100.0013 s / 100.0013 s, wall time 82 -> 58 s); chamber 0.999995 (the
  0.0085 s remainder is the final partial `diffusion_step` at loop exit).
  In LEM the chain now always runs *after* the maps; before, it ran before
  them on backstop iterations. In ODT an accepted eddy now acts on fields up to
  one `diffusion_step` old, as non-backstop iterations already did.
- **A/B reftest** — `74ef16c` (v3.3.0) vs `2912f7d` (fix), np-v3, skl
  (job 17456438; outputs in `devCODT/reftest_ab/{74ef16c,2912f7d}`, staged in
  `reftest_ab/compare/`). Self-comparison control passed in both modes. The A
  legs are bit-identical to the v3.3.0 `_REF` files (the A chamber leg ran on
  a Gold 6140, `_REF` on a Gold 6130, so CPU model within `skl` does not
  matter).
  - *Parcel* — identical through t = 8 s; parcel heights and the 15
    entrainment event times identical. Activation (t ~ 8-16 s) falls in a
    doubled binade, so in v3.3.0 droplets grew twice as fast there: with the
    fix the domain-mean supersaturation peak is higher and later (0.555% at
    t = 20 s vs 0.445% at 10 s; domain max 0.65% vs 0.50%), LWC lags (B/A =
    0.60 at 10 s, 0.77 at 20 s) and mean radius is smaller (1.53 vs 1.82 um
    at 10 s). The runs re-converge by t ~ 42 s (16-64 s is a clean binade;
    LWC B/A = 1.00, mean S within 0.01%). After the first entrainment event
    (56 s) droplet number is dominated by blob placement (see "Parcel reftest
    droplet count is dominated by where entrainment blobs land"), so later
    Np/Nact/LWC differences (final Np 885 vs 869) are not attributable.
    Collisions are too few to compare (67 vs 47).
  - *Chamber* — domain-mean shifts are at the level of random-stream
    divergence: S max shift 0.21%, T 0.19 C, Np/Nact mean rel. diff 2.5%,
    LWC 4.4%, final mean radius 5.89 vs 5.90 um.
- **Impact on existing results** — every parcel run before the fix is biased by
  its own `dt`-dependent pattern of doubled time ranges; not comparable with
  later runs. The effect is largest where condensation is fast relative to
  forcing (activation, and wherever S is far from quasi-steady) and falls in a
  doubled binade.
  Chamber runs change only on the ~40 coincident iterations per 1800 s, but are
  not bit-identical (fewer CC random draws).

---

## Unknown

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

### `collisions.bin` can be short a few records in the final write interval

- **Symptom** — On an EXP002 smoke run (`c529ae9`, v3.1.0, seeded + entraining,
  `n_blob=5`), `sum(N_collisions)` = 4997 but the binary holds **4994** records
  with **0** trailing bytes. `N_coalescences` matched the `flag==1` count exactly
  (159 / 159), so only non-coalescing collisions were short.
- **Localized to the tail** — the running `nc`-minus-binary difference oscillates
  between -1 and 0 for the whole ascent (that oscillation is the accepted
  adjacent-interval swap from the "physics chain twice" issue, since fixed — totals
  preserved), then steps to +3 **only in the last write interval**
  (t = 1960.9 s of a 1960.9 s run). Zero events follow it. Reads as events
  counted into `N_collisions` but never flushed to the stream before
  termination.
- **Consequence** — small and benign *here*: the interval-block reconstruction
  in `projects/SF01_Seeding/analyze.py:372` partitions by `cumsum(N_collisions)`
  positionally, so a deficit at the very end leaves every preceding event's
  interval assignment untouched, and the existing `np.clip` absorbs the tail. It
  would matter if the shortfall ever landed mid-stream: every subsequent event
  would shift by that many positions. Worth a guard that compares
  `sum(N_collisions)` to the record count and reports where they diverge, rather
  than assuming the exact match documented above.
- **Not universal** — an adhoc `c529ae9` run matched exactly (9526 / 9526), as
  did the pre-fix run in the entry above (25025 / 3253). Intermittent, so
  probably conditioned on how the run terminates (this one ended on
  `pressure_limit`) or on activity in the final partial interval.
- **Evidence** — one run, arithmetic on the file and the netCDF counters only.
  No source read yet; the writer/flush path has not been inspected.
- **Found in** — `c529ae9` (`v3.1.0`), EXP002 smoke run.
- **Status** — open, not investigated.
