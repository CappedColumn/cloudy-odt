#!/usr/bin/env python3
"""Compare reftest simulation output against reference files bitwise.

Usage:
    python3 compare_reftest.py <reftest_output_dir>

The output directory should contain:
    reftest.nc              — current run output
    reftest_REF.nc          — reference profiles
    reftest_particles.nc    — current run particle output (optional)
    reftest_particles_REF.nc — reference particles (optional)

Exit code 0 if all variables are identical, 1 if any differ.
"""
import sys
import os
import numpy as np
import netCDF4 as nc


def compare_files(test_path, ref_path, label=""):
    if not os.path.exists(ref_path):
        print(f"  SKIP: reference file not found: {ref_path}")
        return True
    if not os.path.exists(test_path):
        print(f"  FAIL: test file not found: {test_path}")
        return False

    ref = nc.Dataset(ref_path)
    test = nc.Dataset(test_path)
    all_match = True

    # Compare variables common to both files
    common_vars = set(ref.variables) & set(test.variables)
    ref_only = set(ref.variables) - set(test.variables)
    test_only = set(test.variables) - set(ref.variables)

    for var in sorted(common_vars):
        r = ref[var][:]
        t = test[var][:]
        if r.shape != t.shape:
            print(f"  {label}{var}: SHAPE MISMATCH  ref={r.shape} test={t.shape}")
            all_match = False
            continue
        if np.array_equal(r, t):
            print(f"  {label}{var}: identical")
        else:
            diff = np.abs(r - t)
            print(f"  {label}{var}: DIFFERS  max_diff={np.nanmax(diff)}")
            all_match = False

    for var in sorted(ref_only):
        print(f"  {label}{var}: SKIP (ref only)")
    for var in sorted(test_only):
        print(f"  {label}{var}: NEW (test only)")

    ref.close()
    test.close()
    return all_match


def compare_binary(test_path, ref_path):
    if not os.path.exists(ref_path):
        print("  SKIP: reference file not found")
        return True
    if not os.path.exists(test_path):
        print("  SKIP: test file not found")
        return True

    with open(ref_path, 'rb') as f:
        ref_bytes = f.read()
    with open(test_path, 'rb') as f:
        test_bytes = f.read()

    if ref_bytes == test_bytes:
        print("  identical")
        return True
    else:
        print(f"  DIFFERS  ref_size={len(ref_bytes)} test_size={len(test_bytes)}")
        return False


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)

    output_dir = sys.argv[1]

    print("Profiles:")
    profiles_ok = compare_files(
        os.path.join(output_dir, "reftest.nc"),
        os.path.join(output_dir, "reftest_REF.nc"),
    )

    print("\nParticles:")
    particles_ok = compare_files(
        os.path.join(output_dir, "reftest_particles.nc"),
        os.path.join(output_dir, "reftest_particles_REF.nc"),
    )

    print("\nCollisions:")
    collisions_ok = compare_binary(
        os.path.join(output_dir, "reftest_collisions.bin"),
        os.path.join(output_dir, "reftest_collisions_REF.bin"),
    )

    if profiles_ok and particles_ok and collisions_ok:
        print("\nAll variables identical.")
        sys.exit(0)
    else:
        print("\nDifferences detected!")
        sys.exit(1)


if __name__ == "__main__":
    main()
