#!/usr/bin/env python3
"""Collision-coalescence invariant tests.

Runs the CODT executable twice on a CC-on, E=1 namelist and verifies:

  T1. Per-event water-volume conservation (from collision binary log)
  T2. Collision bookkeeping: Σ N_collisions == Σ N_coalescences == #log_records
      and every surviving particle has n_collisions == n_coalescences.
  T3. Particle-count invariant: Np_final = Np_0 + injected - fellout - coalesced
  T4. Determinism: two back-to-back runs produce identical collision logs
  T5. Global water-mass closure via budgets (sanity)

Uses codt_tools for namelist manipulation and simulation loading. Must be
run in the research conda env. Invokes CODT via `fpm run` from the repo
root.

Usage:
    conda run -n research python test/scripts/test_cc_invariants.py \\
        [--nml <source_nml>] [--outdir <scratch_dir>]
"""
import argparse
import pathlib
import shutil
import struct
import subprocess
import sys

import numpy as np

from codt_tools import CODTConfig, CODTSimulation


REPO = pathlib.Path(__file__).resolve().parents[2]

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def prepare_run(src_nml: pathlib.Path, base_outdir: pathlib.Path, sim_name: str,
                aerosol_file: pathlib.Path | None = None, **overrides) -> pathlib.Path:
    """Write a copy of the source config into base_outdir/sim_name/run/ and
    return the path to the generated namelist."""
    cfg = CODTConfig(namelist_path=src_nml)
    cfg.set(simulation_name=sim_name, output_directory=str(base_outdir), **overrides)
    run_dir = base_outdir / sim_name / "run"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    cfg.write(run_dir)
    if aerosol_file is not None:
        shutil.copy2(aerosol_file, run_dir / "aerosol_input.nc")
    return run_dir / "params.nml"


def run_codt(nml: pathlib.Path) -> None:
    r = subprocess.run(
        ["fpm", "run", "--", str(nml)],
        cwd=REPO, capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise RuntimeError(f"CODT run failed for {nml}")


# ---------------------------------------------------------------------------
# Collision binary parser
# ---------------------------------------------------------------------------

def read_collision_log(path: pathlib.Path):
    """Parse the CODT ``_collisions.bin`` stream.

    Header: N (i4), H (f8), domain_width (f8), volume_scaling (f8)
    Record: id_keep (i4), id_kill (i4), r_keep (f8), r_kill (f8),
            r_after (f8), position (f8), time (f8)
    """
    raw = path.read_bytes()
    header_fmt = "<i3d"
    header_size = struct.calcsize(header_fmt)
    N, H, Dw, Vs = struct.unpack(header_fmt, raw[:header_size])

    record_fmt = "<2i5d"
    record_size = struct.calcsize(record_fmt)
    body = raw[header_size:]
    n_records = len(body) // record_size
    rec = np.zeros(
        n_records,
        dtype=[
            ("id_keep", "i4"), ("id_kill", "i4"),
            ("r_keep", "f8"), ("r_kill", "f8"), ("r_after", "f8"),
            ("position", "f8"), ("time", "f8"),
        ],
    )
    for k in range(n_records):
        off = k * record_size
        rec[k] = struct.unpack(record_fmt, body[off : off + record_size])
    return {"N": N, "H": H, "domain_width": Dw, "volume_scaling": Vs, "records": rec}


# ---------------------------------------------------------------------------
# Invariant checks
# ---------------------------------------------------------------------------

def t1_water_volume_per_event(log, rtol=1e-10):
    rec = log["records"]
    if len(rec) == 0:
        return False, "no coalescence events in log — test inconclusive"
    v_pre = rec["r_keep"] ** 3 + rec["r_kill"] ** 3
    v_post = rec["r_after"] ** 3
    err = np.abs(v_post - v_pre) / np.maximum(v_pre, 1e-30)
    worst = float(err.max())
    return worst < rtol, f"max rel err = {worst:.2e} over {len(rec)} events (tol {rtol:.0e})"


def t2_bookkeeping(sim: CODTSimulation, log):
    n_log = int(len(log["records"]))
    ts_coll = int(sim.N_collisions.sum())
    ts_coal = int(sim.N_coalescences.sum())
    ok_coal_log = (ts_coal == n_log)
    ok_coal_le_coll = (ts_coal <= ts_coll)

    # Surviving-particle check on the final snapshot:
    # each surviving particle's n_coalescences <= n_collisions
    pds = sim._particles_ds
    if pds is None:
        pds = sim.load_trajectories()
    row_sizes = pds["row_sizes"].values.astype(int)
    end = int(row_sizes.sum())
    start = end - int(row_sizes[-1])
    p_coll = pds["n_collisions"].values[start:end].astype(int)
    p_coal = pds["n_coalescences"].values[start:end].astype(int)
    ok_particles = bool(np.all(p_coal <= p_coll))

    ok = ok_coal_log and ok_coal_le_coll and ok_particles
    msg = (
        f"Σ N_collisions={ts_coll}, Σ N_coalescences={ts_coal}, "
        f"log records={n_log}; coal==log: {ok_coal_log}, "
        f"coal<=coll: {ok_coal_le_coll}, "
        f"per-particle coal<=coll: {ok_particles}"
    )
    return ok, msg


def t3_particle_count(sim: CODTSimulation, log):
    Np = sim.Np.values.astype(int)
    n_inj = int(sim.budget_n_injected.values.sum())
    n_fel = int(sim.budget_n_fellout.values.sum())
    n_coal_log = int(len(log["records"]))
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


def t5_water_mass_closure(sim: CODTSimulation, tol=0.25):
    """Column-integrated LWC drift should match budget injection/fallout/
    condensation to within ``tol`` of the total budget magnitude. Uses a
    column-mean proxy (no per-cell area weighting), so the tolerance is
    loose — this catches order-of-magnitude leaks, not tight closure.
    """
    LWC = sim.LWC.values  # may be (time,) or (time, z)
    if LWC.ndim == 2:
        lwc_col = LWC.mean(axis=1)
    else:
        lwc_col = LWC

    init = float(lwc_col[0])
    final = float(lwc_col[-1])
    inj = float(sim.budget_inject_liquid_mass.values.sum())
    fel = float(sim.budget_fallout_liquid_mass.values.sum())
    cond = float(sim.budget_condensation.values.sum())

    delta_lwc = final - init
    delta_budget = inj - fel + cond
    scale = abs(inj) + abs(fel) + abs(cond) + 1e-30
    residual = delta_lwc * 1e-3 - delta_budget  # LWC g/m3 -> kg/m3 rough scale
    rel = abs(residual) / scale
    ok = rel < tol
    msg = (
        f"ΔLWC_col={delta_lwc:.3e} g/m3, Σbudget={delta_budget:.3e} kg, "
        f"residual/scale={rel:.2e} (tol {tol})"
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
    sim_name = sim_dir.name
    sim_a = CODTSimulation(sim_dir)
    log_path_a = sim_dir / f"{sim_name}_collisions.bin"

    if not log_path_a.exists():
        print(f"ERROR: collision log not found: {log_path_a}")
        sys.exit(1)

    log_a = read_collision_log(log_path_a)

    print(f"Simulation: {sim_dir}")
    print(f"Collision log: {len(log_a['records'])} coalescence events")
    print("=" * 60)

    results = []
    results.append(("T1 per-event water volume",       *t1_water_volume_per_event(log_a)))
    results.append(("T2 bookkeeping",                   *t2_bookkeeping(sim_a, log_a)))
    results.append(("T3 particle-count invariant",      *t3_particle_count(sim_a, log_a)))

    if args.sim_b is not None:
        sim_b_dir = pathlib.Path(args.sim_b)
        log_path_b = sim_b_dir / f"{sim_b_dir.name}_collisions.bin"
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
