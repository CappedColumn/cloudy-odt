#!/usr/bin/env python3
"""Collision-coalescence invariant tests on finished CODT runs.

Reads an existing simulation directory (run with ``write_collisions = .true.``)
and verifies:

  T1. Per-event water-volume conservation over the coalescence records
  T2. Collision bookkeeping: Σ N_collisions == #records,
      Σ N_coalescences == #coalescence records (flag 1), and every surviving
      particle has n_coalescences <= n_collisions
  T3. Particle-count invariant: Np_final = Np_0 + injected - fellout - coalesced
  T4. Determinism: two runs produce identical collision logs (optional)
  T5. Global water-mass closure via budgets (sanity)

The collision log is read with ``Simulation.load_collisions`` from
codt_tools, so the record layout lives in one place. Every collision is
logged, coalescing or not; the ``coalesced`` flag tells them apart.

Usage:
    conda run -n research python test/scripts/test_cc_invariants.py \\
        <sim_dir> [--sim-b <second_sim_dir>]
"""
import argparse
import pathlib
import sys

import numpy as np

from codt_tools import Simulation



PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"


# ---------------------------------------------------------------------------
# Collision log
# ---------------------------------------------------------------------------

def coalescence_records(log):
    """Records of events where the droplets merged (flag 1)."""
    ev = log["events"]
    return ev[ev["coalesced"] == 1]


# ---------------------------------------------------------------------------
# Invariant checks
# ---------------------------------------------------------------------------

def t1_water_volume_per_event(log, rtol=1e-10):
    rec = coalescence_records(log)
    if len(rec) == 0:
        return False, "no coalescence events in log — test inconclusive"
    v_pre = rec["r_keep"] ** 3 + rec["r_kill"] ** 3
    v_post = rec["r_after"] ** 3
    err = np.abs(v_post - v_pre) / np.maximum(v_pre, 1e-30)
    worst = float(err.max())
    return worst < rtol, f"max rel err = {worst:.2e} over {len(rec)} events (tol {rtol:.0e})"


def t2_bookkeeping(sim: Simulation, log):
    n_log = int(len(log["events"]))
    n_log_coal = int(len(coalescence_records(log)))
    ts_coll = int(sim.N_collisions.sum())
    ts_coal = int(sim.N_coalescences.sum())
    ok_coll_log = (ts_coll == n_log)
    ok_coal_log = (ts_coal == n_log_coal)
    ok_coal_le_coll = (ts_coal <= ts_coll)

    # Surviving-particle check on the final snapshot:
    # each surviving particle's n_coalescences <= n_collisions.
    # Needs trajectory records (write_trajectories = .true. over the run's end);
    # skipped without them.
    pds = sim.load_trajectories() if sim.has_trajectories else None
    if pds is not None and pds["row_sizes"].size > 0:
        row_sizes = pds["row_sizes"].values.astype(int)
        end = int(row_sizes.sum())
        start = end - int(row_sizes[-1])
        p_coll = pds["n_collisions"].values[start:end].astype(int)
        p_coal = pds["n_coalescences"].values[start:end].astype(int)
        ok_particles = bool(np.all(p_coal <= p_coll))
    else:
        ok_particles = "skipped (no trajectory records)"

    ok = ok_coll_log and ok_coal_log and ok_coal_le_coll and ok_particles is not False
    msg = (
        f"Σ N_collisions={ts_coll}, log records={n_log}; "
        f"Σ N_coalescences={ts_coal}, coalescence records={n_log_coal}; "
        f"coll==log: {ok_coll_log}, coal==log: {ok_coal_log}, "
        f"coal<=coll: {ok_coal_le_coll}, "
        f"per-particle coal<=coll: {ok_particles}"
    )
    return ok, msg


def t3_particle_count(sim: Simulation, log):
    Np = sim.Np.values.astype(int)
    n_inj = int(sim.budget_n_injected.values.sum())
    n_fel = int(sim.budget_n_fellout.values.sum())
    n_coal_log = int(len(coalescence_records(log)))
    n_coal_budget = int(sim.budget_n_coalesced.values.sum())

    Np_0 = int(Np[0])
    Np_f = int(Np[-1])
    expected = Np_0 + n_inj - n_fel - n_coal_budget
    ok = (Np_f == expected) and (n_coal_budget == n_coal_log)
    msg = (
        f"Np_final={Np_f}, expected={expected} "
        f"(Np_0={Np_0} + inj={n_inj} - fell={n_fel} - coal={n_coal_budget})"
        f"; coal_budget==coal_log: {n_coal_budget == n_coal_log}"
    )
    return ok, msg


def t4_determinism(log_a: pathlib.Path, log_b: pathlib.Path):
    ba = log_a.read_bytes()
    bb = log_b.read_bytes()
    return (ba == bb), f"sizes {len(ba)} vs {len(bb)}"


def t5_water_mass_closure(sim: Simulation, tol=1e-6):
    """Liquid water mass closure: QC(final) = QC(init) + injected - fallout + condensation.

    All terms in kg. domain_volume = volume_scaling * domain_width^2 * H,
    where domain_width = 0.001 m (hardcoded in globals.f90).
    """
    DOMAIN_WIDTH = 0.001  # m, from globals.f90
    ds = sim._ds
    H = float(ds.attrs.get("PARAMETERS.H", ds.attrs.get("H")))
    vs = float(ds.attrs.get("PARAMETERS.volume_scaling", ds.attrs.get("volume_scaling")))
    domain_volume = vs * DOMAIN_WIDTH**2 * H

    LWC = sim.LWC.values
    qc_init = float(LWC[0]) * 1e-3 * domain_volume   # g/m3 -> kg
    qc_final = float(LWC[-1]) * 1e-3 * domain_volume

    inj = float(sim.budget_inject_liquid_mass.values.sum())
    fel = float(sim.budget_fallout_liquid_mass.values.sum())
    cond = float(sim.budget_condensation.values.sum())

    delta_qc = qc_final - qc_init
    delta_budget = inj - fel + cond
    residual = delta_qc - delta_budget
    scale = abs(inj) + abs(fel) + abs(cond) + 1e-30
    rel = abs(residual) / scale
    ok = rel < tol
    msg = (
        f"ΔQC={delta_qc:.3e} kg, Σ(inj-fel+cond)={delta_budget:.3e} kg, "
        f"rel_err={rel:.2e} (tol {tol:.0e})"
    )
    return ok, msg


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "sim_dir",
        help="Path to an existing simulation output directory",
    )
    ap.add_argument(
        "--sim-b",
        default=None,
        help="Path to a second simulation for determinism check (T4). "
             "If omitted, T4 is skipped.",
    )
    args = ap.parse_args()

    sim_dir = pathlib.Path(args.sim_dir)
    sim_a = Simulation(sim_dir)
    log_path_a = sim_dir / f"{sim_a.name}_collisions.bin"

    if not sim_a.has_collisions:
        print(f"ERROR: collision log not found: {log_path_a}")
        sys.exit(1)

    log_a = sim_a.load_collisions()

    print(f"Simulation: {sim_dir}")
    print(f"Collision log: {len(log_a['events'])} collisions, "
          f"{len(coalescence_records(log_a))} coalescences")
    print("=" * 60)

    results = []
    results.append(("T1 per-event water volume",       *t1_water_volume_per_event(log_a)))
    results.append(("T2 bookkeeping",                   *t2_bookkeeping(sim_a, log_a)))
    results.append(("T3 particle-count invariant",      *t3_particle_count(sim_a, log_a)))

    if args.sim_b is not None:
        sim_b_dir = pathlib.Path(args.sim_b)
        log_path_b = sim_b_dir / f"{Simulation(sim_b_dir).name}_collisions.bin"
        results.append(("T4 determinism (collision log)",   *t4_determinism(log_path_a, log_path_b)))
    else:
        results.append(("T4 determinism (collision log)",   True, "skipped (no --sim-b provided)"))

    results.append(("T5 water-mass closure (proxy)",    *t5_water_mass_closure(sim_a)))

    print()
    n_fail = 0
    for name, ok, msg in results:
        tag = PASS if ok else FAIL
        print(f"  {tag}  {name}")
        print(f"         {msg}")
        if not ok:
            n_fail += 1

    print()
    if n_fail:
        print(f"{FAIL}: {n_fail} of {len(results)} tests failed")
        sys.exit(1)
    print(f"{PASS}: all {len(results)} tests passed")


if __name__ == "__main__":
    main()
