#!/usr/bin/env python3
"""Physics-level comparison of reftest output against reference.

Unlike compare_reftest.py (bitwise), this script evaluates whether
differences are physically meaningful by reporting:
  - Time series alignment and drift
  - Profile-level statistics (RMSE, max abs diff, relative diff)
  - Domain-mean evolution comparison
  - Particle count and size distribution differences

Usage:
    python3 compare_reftest_physics.py <reftest_output_dir>
"""
import sys
import os
import numpy as np
import netCDF4 as nc


def fmt(val, units=""):
    """Format a number with appropriate precision."""
    if isinstance(val, (int, np.integer)):
        return f"{val}{' ' + units if units else ''}"
    if abs(val) < 1e-6:
        return f"{val:.2e}{' ' + units if units else ''}"
    return f"{val:.6g}{' ' + units if units else ''}"


def compare_time_series(ref_time, test_time):
    """Compare time axes for alignment and drift."""
    print("  Time axis:")
    print(f"    REF steps: {len(ref_time)},  Test steps: {len(test_time)}")

    n = min(len(ref_time), len(test_time))
    rt, tt = ref_time[:n], test_time[:n]
    drift = tt - rt

    print(f"    Comparable steps: {n}")
    print(f"    REF time range:  [{rt[0]:.4f}, {rt[-1]:.4f}] s")
    print(f"    Test time range: [{tt[0]:.4f}, {tt[-1]:.4f}] s")
    print(f"    Max time drift:  {fmt(np.max(np.abs(drift)), 's')}")
    if n > 1:
        ref_dt = np.diff(rt)
        test_dt = np.diff(tt)
        print(f"    REF  mean dt:  {fmt(np.mean(ref_dt), 's')}")
        print(f"    Test mean dt:  {fmt(np.mean(test_dt), 's')}")

    return n


def compare_profiles(ref_ds, test_ds, varname, n_times, units=""):
    """Compare a 2D (time, z) profile variable."""
    r = ref_ds[varname][:n_times, :]
    t = test_ds[varname][:n_times, :]

    if np.array_equal(r, t):
        print(f"  {varname}: identical")
        return

    diff = t - r
    abs_diff = np.abs(diff)
    rmse = np.sqrt(np.mean(diff**2))

    # Relative difference (avoid division by zero)
    r_range = np.max(r) - np.min(r)
    rel_rmse_pct = (rmse / r_range * 100) if r_range > 0 else float('inf')

    # Domain-mean evolution
    ref_mean = np.mean(r, axis=1)
    test_mean = np.mean(t, axis=1)
    mean_diff = test_mean - ref_mean
    max_mean_diff = np.max(np.abs(mean_diff))

    print(f"  {varname} [{units}]:")
    print(f"    RMSE:           {fmt(rmse, units)}")
    print(f"    RMSE / range:   {rel_rmse_pct:.4f}%")
    print(f"    Max abs diff:   {fmt(np.max(abs_diff), units)}")
    print(f"    Mean abs diff:  {fmt(np.mean(abs_diff), units)}")
    print(f"    Max domain-mean shift: {fmt(max_mean_diff, units)}")

    # Check if differences grow in time
    rmse_per_step = np.sqrt(np.mean(diff**2, axis=1))
    if len(rmse_per_step) > 3:
        early = np.mean(rmse_per_step[:3])
        late = np.mean(rmse_per_step[-3:])
        if early > 0:
            growth = late / early
            trend = "growing" if growth > 2 else "stable" if growth < 1.5 else "mildly growing"
            print(f"    Temporal trend:  {trend} (late/early RMSE ratio: {growth:.1f}x)")


def compare_scalars(ref_ds, test_ds, varname, n_times, units=""):
    """Compare a 1D (time,) scalar variable."""
    r = ref_ds[varname][:n_times]
    t = test_ds[varname][:n_times]

    if np.array_equal(r, t):
        print(f"  {varname}: identical")
        return

    diff = t.astype(float) - r.astype(float)
    abs_diff = np.abs(diff)

    ref_mean = np.mean(np.abs(r.astype(float)))
    rel_pct = (np.mean(abs_diff) / ref_mean * 100) if ref_mean > 0 else float('inf')

    print(f"  {varname} [{units}]:")
    print(f"    Max abs diff:     {fmt(np.max(abs_diff))}")
    print(f"    Mean abs diff:    {fmt(np.mean(abs_diff))}")
    print(f"    Mean rel diff:    {rel_pct:.4f}%")
    print(f"    REF final value:  {fmt(r[-1])}")
    print(f"    Test final value: {fmt(t[-1])}")


def compare_particles(output_dir):
    """Compare particle output files at a physics level."""
    ref_path = os.path.join(output_dir, "reftest_particles_REF.nc")
    test_path = os.path.join(output_dir, "reftest_particles.nc")

    if not os.path.exists(ref_path):
        print("  SKIP: reference particle file not found")
        return
    if not os.path.exists(test_path):
        print("  FAIL: test particle file not found")
        return

    ref = nc.Dataset(ref_path)
    test = nc.Dataset(test_path)

    ref_rows = ref["row_sizes"][:]
    test_rows = test["row_sizes"][:]
    n_steps = min(len(ref_rows), len(test_rows))

    print(f"  Time steps: REF={len(ref_rows)}, Test={len(test_rows)}")
    print(f"  Total records: REF={ref.dimensions['record'].size}, "
          f"Test={test.dimensions['record'].size}")

    # Compare particle counts per time step
    count_diff = test_rows[:n_steps].astype(int) - ref_rows[:n_steps].astype(int)
    print(f"  Particle count diff (per step): "
          f"max={np.max(np.abs(count_diff))}, mean={np.mean(np.abs(count_diff)):.1f}")

    # Compare aggregate radius distributions at last common time step
    ref_offset = int(np.sum(ref_rows[:n_steps-1]))
    ref_count = int(ref_rows[n_steps-1])
    test_offset = int(np.sum(test_rows[:n_steps-1]))
    test_count = int(test_rows[n_steps-1])

    ref_radii = ref["radius"][ref_offset:ref_offset+ref_count]
    test_radii = test["radius"][test_offset:test_offset+test_count]

    print(f"  Final step particle count: REF={ref_count}, Test={test_count}")
    print(f"  Final step mean radius: REF={fmt(np.mean(ref_radii), 'um')}, "
          f"Test={fmt(np.mean(test_radii), 'um')}")
    print(f"  Final step std radius:  REF={fmt(np.std(ref_radii), 'um')}, "
          f"Test={fmt(np.std(test_radii), 'um')}")

    ref.close()
    test.close()


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)

    output_dir = sys.argv[1]
    ref_path = os.path.join(output_dir, "reftest_REF.nc")
    test_path = os.path.join(output_dir, "reftest.nc")

    if not os.path.exists(ref_path):
        print(f"Reference file not found: {ref_path}")
        sys.exit(1)
    if not os.path.exists(test_path):
        print(f"Test file not found: {test_path}")
        sys.exit(1)

    ref = nc.Dataset(ref_path)
    test = nc.Dataset(test_path)

    # Time alignment
    print("=" * 60)
    print("TIME ALIGNMENT")
    print("=" * 60)
    n = compare_time_series(ref["time"][:], test["time"][:])

    # Profile variables (time, z)
    print()
    print("=" * 60)
    print("PROFILE VARIABLES")
    print("=" * 60)
    profile_vars = {
        "T": "C", "QV": "g/kg", "Tv": "C", "S": "%"
    }
    for var, units in profile_vars.items():
        if var in ref.variables and var in test.variables:
            compare_profiles(ref, test, var, n, units)

    # Scalar time series
    print()
    print("=" * 60)
    print("SCALAR TIME SERIES")
    print("=" * 60)
    scalar_vars = {
        "Np": "#", "Nact": "#", "Nun": "#",
        "Ravg": "um", "LWC": "g/m3"
    }
    for var, units in scalar_vars.items():
        if var in ref.variables and var in test.variables:
            compare_scalars(ref, test, var, n, units)

    # DSD comparison
    print()
    print("=" * 60)
    print("DROPLET SIZE DISTRIBUTION")
    print("=" * 60)
    for dsd_var in ["DSD", "DSD_1", "DSD_2"]:
        if dsd_var in ref.variables and dsd_var in test.variables:
            r = ref[dsd_var][:n, :]
            t = test[dsd_var][:n, :]
            if np.array_equal(r, t):
                print(f"  {dsd_var}: identical")
            else:
                diff = t.astype(float) - r.astype(float)
                total_ref = np.sum(r, axis=1).astype(float)
                total_test = np.sum(t, axis=1).astype(float)
                count_diff = total_test - total_ref
                print(f"  {dsd_var}:")
                print(f"    Max bin count diff:   {int(np.max(np.abs(diff)))}")
                print(f"    Mean bin count diff:  {np.mean(np.abs(diff)):.2f}")
                print(f"    Total count diff (per step): "
                      f"max={int(np.max(np.abs(count_diff)))}, "
                      f"mean={np.mean(np.abs(count_diff)):.1f}")

    ref.close()
    test.close()

    # Particle-level comparison
    print()
    print("=" * 60)
    print("PARTICLE TRAJECTORIES")
    print("=" * 60)
    compare_particles(output_dir)


if __name__ == "__main__":
    main()
