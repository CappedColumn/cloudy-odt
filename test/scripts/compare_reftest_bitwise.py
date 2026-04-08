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

    for var in ref.variables:
        r = ref[var][:]
        t = test[var][:]
        if np.array_equal(r, t):
            print(f"  {label}{var}: identical")
        else:
            diff = np.abs(r - t)
            print(f"  {label}{var}: DIFFERS  max_diff={np.nanmax(diff)}")
            all_match = False

    ref.close()
    test.close()
    return all_match


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

    if profiles_ok and particles_ok:
        print("\nAll variables identical.")
        sys.exit(0)
    else:
        print("\nDifferences detected!")
        sys.exit(1)


if __name__ == "__main__":
    main()
