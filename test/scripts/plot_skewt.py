"""Plot CODT parcel ascent trajectories on a Skew-T diagram using MetPy.

Usage:
    python plot_skewt.py <nc_file1> [nc_file2 ...] [-o output_image.png] [-e parcel_input.nc]

Overlays domain-mean T and Td ascent profiles from multiple simulations.
Each simulation gets a unique color. Optionally plots the environmental
sounding and a moist adiabat from the parcel starting conditions.
Prints T, QV, and S at 10 mb intervals.
"""
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units
from pathlib import Path

COLORS = ['r', 'b', 'darkorange', 'purple', 'green', 'brown', 'teal', 'magenta']


def load_parcel_data(ncpath):
    ds = nc.Dataset(ncpath, 'r')
    pres = ds['parcel_pressure'][:]
    T_mean = np.mean(ds['T'][:], axis=1)
    QV_mean = np.mean(ds['QV'][:], axis=1)
    S_mean = np.mean(ds['S'][:], axis=1)
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
    e = (qv_kgkg * pres_mb * 100.0) / (0.622 + qv_kgkg)
    e_unit = (e / 100.0) * units.mbar
    return mpcalc.dewpoint(e_unit).magnitude


def dewpoint_from_rh(T_K, rh):
    T_C = T_K - 273.15
    e_sat = mpcalc.saturation_vapor_pressure(T_C * units.degC)
    e = rh * e_sat
    return mpcalc.dewpoint(e).magnitude


def label_from_path(ncpath):
    return Path(ncpath).stem


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    args = sys.argv[1:]
    outpath = None
    env_path = None

    if '-o' in args:
        oi = args.index('-o')
        outpath = args[oi + 1]
        args = args[:oi] + args[oi + 2:]

    if '-e' in args:
        ei = args.index('-e')
        env_path = args[ei + 1]
        args = args[:ei] + args[ei + 2:]

    nc_files = args
    if not outpath:
        outpath = str(Path(__file__).resolve().parent.parent / 'data' / 'figs' / 'skewt_comparison.png')

    fig = plt.figure(figsize=(10, 12))
    skew = SkewT(fig, rotation=45)

    all_pmin = []
    all_T = []
    all_Td = []
    all_P = []

    for i, ncpath in enumerate(nc_files):
        color = COLORS[i % len(COLORS)]
        label = label_from_path(ncpath)

        pres, T, QV, S = load_parcel_data(ncpath)
        asc, desc = split_ascent_descent(pres)
        all_pmin.append(np.min(pres))

        p_min = int(np.ceil(np.min(pres) / 10.0) * 10)
        targets_up = np.arange(950, p_min - 1, -10)
        print_table(pres[asc], T[asc], QV[asc], S[asc], f'{label} ASCENT', targets_up)

        Td_asc = np.array([dewpoint_from_qv(q, p) for q, p in zip(QV[asc], pres[asc])])

        all_T.extend(T[asc])
        all_Td.extend(Td_asc)
        all_P.extend(pres[asc])

        skew.plot(pres[asc] * units.mbar, T[asc] * units.degC, '-',
                  color=color, linewidth=2, label=f'{label} T')
        skew.plot(pres[asc] * units.mbar, Td_asc * units.degC, '--',
                  color=color, linewidth=1.5, label=f'{label} Td')

    # Environmental sounding and moist adiabat
    if env_path:
        ds_env = nc.Dataset(env_path, 'r')
        env_p = ds_env['env_pressure'][:] / 100.0  # Pa -> mb
        env_T_K = ds_env['env_temperature'][:]
        env_T_C = env_T_K - 273.15
        env_rh = ds_env['env_RH'][:]
        ds_env.close()

        env_Td_C = np.array([dewpoint_from_rh(tk, rh) for tk, rh in zip(env_T_K, env_rh)])

        skew.plot(env_p * units.mbar, env_T_C * units.degC, 'k-',
                  linewidth=1.5, label='Env T')
        skew.plot(env_p * units.mbar, env_Td_C * units.degC, 'k--',
                  linewidth=1.5, label='Env Td')

        p_hi = max(all_P)
        p_lo = min(all_pmin)
        mask = (env_p <= p_hi + 10) & (env_p >= p_lo - 10)
        all_T.extend(env_T_C[mask])
        all_Td.extend(env_Td_C[mask])
        all_P.extend(env_p[mask])

    # Moist adiabat from parcel starting conditions
    p0 = all_P[0]
    T0 = all_T[0]
    p_min_all = int(np.ceil(min(all_pmin) / 10.0) * 10)
    p_profile = np.arange(p0, max(p_min_all - 20, 700) - 1, -1) * units.mbar
    moist_T = mpcalc.moist_lapse(p_profile, T0 * units.degC, p0 * units.mbar)
    skew.plot(p_profile, moist_T, 'k:', linewidth=1.0, label='Moist adiabat')

    skew.plot_dry_adiabats(alpha=0.3)
    skew.plot_moist_adiabats(alpha=0.3)
    skew.plot_mixing_lines(alpha=0.3)

    skew.ax.set_ylim(960, max(p_min_all - 20, 700))

    # Skew-aware x-limits
    p_bot = skew.ax.get_ylim()[0]
    rot = 45.0
    all_T = np.array(all_T)
    all_Td = np.array(all_Td)
    all_P = np.array(all_P)
    shift = rot * np.log10(p_bot / all_P)
    x_t = all_T + shift
    x_td = all_Td + shift
    x_hi = max(np.max(x_t), np.max(x_td))
    x_lo = min(np.min(all_T), np.min(all_Td))
    pad = 5.0
    skew.ax.set_xlim(x_lo, x_hi + pad)

        skew.ax.legend(loc='upper left', fontsize=8)
    skew.ax.set_title('CODT Parcel Entrainment Comparison')

    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f'\nSaved: {outpath}')


if __name__ == '__main__':
    main()
