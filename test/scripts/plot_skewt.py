"""Plot CODT parcel ascent/descent on a Skew-T diagram using MetPy.

Usage:
    python plot_skewt.py <path_to_output.nc> [output_image.png]

Overlays the domain-mean T and Td profiles (ascent in solid, descent in dashed)
on a Skew-T log-P diagram. Prints T, QV, and S at 10 mb intervals.
"""
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units

def load_parcel_data(ncpath):
    ds = nc.Dataset(ncpath, 'r')
    pres = ds['parcel_pressure'][:]  # mb
    T_mean = np.mean(ds['T'][:], axis=1)  # domain-mean T (C)
    QV_mean = np.mean(ds['QV'][:], axis=1)  # domain-mean QV (g/kg)
    S_mean = np.mean(ds['S'][:], axis=1)  # domain-mean SS (%)
    ds.close()
    return pres, T_mean, QV_mean, S_mean


def split_ascent_descent(pres):
    peak_idx = np.argmin(pres)
    return slice(0, peak_idx + 1), slice(peak_idx, len(pres))


def print_table(pres, T, QV, S, label, targets):
    print(f'\n=== {label} ===')
    print(f'{"P (mb)":>8} {"T (C)":>8} {"QV (g/kg)":>10} {"S (%)":>8}')
    print('-' * 40)
    for tgt in targets:
        idx = np.argmin(np.abs(pres - tgt))
        if np.abs(pres[idx] - tgt) < 5:
            print(f'{pres[idx]:8.1f} {T[idx]:8.2f} {QV[idx]:10.3f} {S[idx]:8.3f}')


def dewpoint_from_qv(qv_kgkg, pres_mb):
    e = (qv_kgkg * pres_mb * 100.0) / (0.622 + qv_kgkg)  # vapor pressure in Pa
    e_unit = (e / 100.0) * units.mbar
    return mpcalc.dewpoint(e_unit).magnitude  # C


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    ncpath = sys.argv[1]
    outpath = sys.argv[2] if len(sys.argv) > 2 else ncpath.replace('.nc', '_skewt.png')

    pres, T, QV, S = load_parcel_data(ncpath)
    asc, desc = split_ascent_descent(pres)

    # Print tables at 10 mb intervals
    p_min = int(np.ceil(np.min(pres) / 10.0) * 10)
    targets_up = np.arange(950, p_min - 1, -10)
    targets_down = np.arange(p_min, 951, 10)

    print_table(pres[asc], T[asc], QV[asc], S[asc], 'ASCENT', targets_up)
    print_table(pres[desc], T[desc], QV[desc], S[desc], 'DESCENT', targets_down)

    # Compute dewpoint
    Td_asc = np.array([dewpoint_from_qv(q, p) for q, p in zip(QV[asc], pres[asc])])
    Td_desc = np.array([dewpoint_from_qv(q, p) for q, p in zip(QV[desc], pres[desc])])

    # Plot
    fig = plt.figure(figsize=(10, 12))
    skew = SkewT(fig, rotation=45)

    # Ascent (solid)
    skew.plot(pres[asc] * units.mbar, T[asc] * units.degC, 'r-', linewidth=2, label='T ascent')
    skew.plot(pres[asc] * units.mbar, Td_asc * units.degC, 'g-', linewidth=2, label='Td ascent')

    # Descent (dashed)
    skew.plot(pres[desc] * units.mbar, T[desc] * units.degC, 'r--', linewidth=1.5, label='T descent')
    skew.plot(pres[desc] * units.mbar, Td_desc * units.degC, 'g--', linewidth=1.5, label='Td descent')

    # Moist adiabat from LCL (where T ≈ Td on ascent)
    lcl_idx = np.argmin(np.abs(T[asc] - Td_asc))
    lcl_p = pres[asc][lcl_idx] * units.mbar
    lcl_T = T[asc][lcl_idx] * units.degC
    p_profile = np.arange(float(lcl_p.magnitude), float(max(p_min - 20, 800)) - 1, -1) * units.mbar
    moist_T = mpcalc.moist_lapse(p_profile, lcl_T, lcl_p)
    skew.plot(p_profile, moist_T, 'k--', linewidth=0.8, label='Moist adiabat (LCL)')

    # Reference lines
    skew.plot_dry_adiabats(alpha=0.3)
    skew.plot_moist_adiabats(alpha=0.3)
    skew.plot_mixing_lines(alpha=0.3)

    skew.ax.set_ylim(960, max(p_min - 20, 800))
    skew.ax.set_xlim(14, 18)
    skew.ax.legend(loc='upper left')
    skew.ax.set_title('CODT Parcel Ascent/Descent')

    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f'\nSaved: {outpath}')


if __name__ == '__main__':
    main()
