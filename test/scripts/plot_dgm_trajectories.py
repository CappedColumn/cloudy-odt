"""Plot DGM solver trajectories (RK45 vs ROS3) in S-r space with Köhler curves."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from io import StringIO

# --- Köhler curve parameters for (NH4)2SO4 ---
rho_solute = 2165.0
molar_mass_solute = 132.1395e-3
n_ions = 3.0
Mw = 18.015e-3
rho_l = 997.0
sigma_w = 7.392730e-2
Rv = 461.5
c7 = 0.4363021
pi_43 = 4.0 / 3.0 * np.pi
T = 280.0


def kohler_curve(r, r_solute):
    """Equilibrium supersaturation (fraction) vs radius."""
    sol_mass = pi_43 * r_solute**3 * rho_solute
    solute_c7 = sol_mass * c7
    raoult_coeff = n_ions * (Mw / molar_mass_solute) * sol_mass

    rho_sol = (r**3 * pi_43 * rho_l + solute_c7) / (r**3 * pi_43)
    denom = pi_43 * r**3 * rho_sol - sol_mass
    denom = np.where(denom > 0, denom, np.nan)

    kelvin = 2.0 * sigma_w / (Rv * T * rho_sol * r)
    raoult = raoult_coeff / denom

    return kelvin - raoult


repo = Path(__file__).resolve().parents[2]
csv_path = repo / "test" / "data" / "dgm_solver_trajectories.csv"

# --- Run the test and capture CSV ---
result = subprocess.run(
    ["bash", "-c", "source fpm_env && fpm test test_dgm_solvers 2>/dev/null"],
    capture_output=True, text=True, cwd=str(repo),
)
lines = [l for l in result.stdout.strip().split("\n") if "," in l]
csv_text = "\n".join(lines)

csv_path.parent.mkdir(parents=True, exist_ok=True)
csv_path.write_text(csv_text)
print(f"Saved CSV to {csv_path}")

df = pd.read_csv(StringIO(csv_text), skipinitialspace=True)

# --- Mappings ---
aero_colors = {
    "rs5nm": "#1f77b4",
    "rs50nm": "#ff7f0e",
    "rs200nm": "#2ca02c",
    "rs1000nm": "#d62728",
}
aero_nice = {
    "rs5nm": "5 nm",
    "rs50nm": "50 nm",
    "rs200nm": "200 nm",
    "rs1000nm": "1000 nm",
}
r_solute_map = {
    "rs5nm": 5.0e-9,
    "rs50nm": 50.0e-9,
    "rs200nm": 200.0e-9,
    "rs1000nm": 1000.0e-9,
}
solver_markers = {"RK45": "o", "ROS3": "x"}
solver_sizes = {"RK45": 5, "ROS3": 6}
solver_ls = {"RK45": "-", "ROS3": "--"}

aerosols = df["aerosol"].unique()
ics = df["ic"].unique()
solvers = df["solver"].unique()

# --- Figure 1: Trajectories in S-r space ---
fig, ax = plt.subplots(figsize=(12, 8))

for aero in aerosols:
    color = aero_colors[aero]
    rs = r_solute_map[aero]
    r_floor = rs * 1.01
    r_max = df.loc[df["aerosol"] == aero, "radius"].max() * 1.5
    r_kohler = np.geomspace(r_floor * 1.1, max(30e-6, r_max), 500)
    ss_kohler = kohler_curve(r_kohler, rs)
    ax.plot(r_kohler * 1e6, ss_kohler * 100, color=color, lw=2,
            alpha=0.35, label=f"Köhler {aero_nice[aero]}", zorder=1)

for aero in aerosols:
    color = aero_colors[aero]
    for ic_name in ics:
        for solver in solvers:
            sub = df[(df["aerosol"] == aero) & (df["ic"] == ic_name) & (df["solver"] == solver)]
            sub = sub[sub["step"] <= 20]  # fine steps only, skip coarse
            if sub.empty:
                continue
            r_um = sub["radius"].values * 1e6
            ss_pct = sub["SS"].values * 100
            mk = solver_markers.get(solver, "x")
            fillstyle = "none" if solver == "RK45" else "full"
            ms = solver_sizes.get(solver, 5)
            ax.plot(r_um, ss_pct, color=color,
                    marker=mk, fillstyle=fillstyle,
                    linestyle=solver_ls.get(solver, "-"),
                    markersize=ms, lw=1.0, alpha=0.8, zorder=2)
            ax.scatter(r_um[0], ss_pct[0], s=50, marker="o",
                       facecolors=color, edgecolors="black", linewidths=1.2, zorder=3)
            end_color = "saddlebrown" if solver == "RK45" else "deeppink"
            ax.scatter(r_um[-1], ss_pct[-1], s=50, marker="o",
                       facecolors=color, edgecolors=end_color, linewidths=1.2, zorder=3)

from matplotlib.lines import Line2D
handles = []
for aero in aerosols:
    handles.append(Line2D([0], [0], color=aero_colors[aero], lw=2,
                          label=f"{aero_nice[aero]} aerosol"))
handles.append(Line2D([0], [0], color="gray", lw=2, alpha=0.35, label="Köhler curve"))
for solver in solvers:
    fs = "none" if solver == "RK45" else "full"
    handles.append(Line2D([0], [0], color="k",
                          marker=solver_markers[solver], fillstyle=fs,
                          linestyle=solver_ls[solver],
                          markersize=solver_sizes[solver], lw=1.0, label=solver))
handles.append(Line2D([0], [0], color="gray", marker="o", linestyle="None",
                      markeredgecolor="black", markeredgewidth=1.2,
                      markersize=7, label="Start"))
handles.append(Line2D([0], [0], color="gray", marker="o", linestyle="None",
                      markeredgecolor="deeppink", markeredgewidth=1.2,
                      markersize=7, label="End (ROS3)"))
handles.append(Line2D([0], [0], color="gray", marker="o", linestyle="None",
                      markeredgecolor="saddlebrown", markeredgewidth=1.2,
                      markersize=7, label="End (RK45)"))
ax.legend(handles=handles, fontsize=8, ncol=2, loc="best")

ax.set_xlabel("Radius (µm)")
ax.set_ylabel("Supersaturation (%)")
ax.set_title("DGM Solver Trajectories: RK45 vs ROS3\n4 aerosol sizes × 4 initial conditions")
ax.axhline(0, color="k", lw=0.5, ls="--")
ax.set_xscale("log")
ax.set_ylim(-71, 1)

fig.tight_layout()
png_path = repo / "test" / "data" / "dgm_solver_trajectories.png"
fig.savefig(png_path, dpi=150)
print(f"Saved trajectory plot to {png_path}")

# --- Figure 2: ROS3 - RK45 differences ---
rk = df[df["solver"] == "RK45"].set_index(["aerosol", "ic", "step"])
ros = df[df["solver"] == "ROS3"].set_index(["aerosol", "ic", "step"])

diff = rk[["time", "radius"]].copy()
diff.rename(columns={"radius": "r_rk45"}, inplace=True)
diff["r_ros3"] = ros["radius"]
diff["dr_ros3"] = ros["radius"].values - rk["radius"].values
diff["rel_dr"] = diff["dr_ros3"] / rk["radius"].values
diff = diff.reset_index()

diff_path = repo / "test" / "data" / "dgm_solver_diffs.csv"
diff.to_csv(diff_path, index=False, float_format="%.6e")
print(f"Saved difference table to {diff_path}")

# Compute critical radius for sub/super-critical classification
def critical_radius(r_solute):
    rs_arr = np.geomspace(r_solute * 1.02, r_solute * 1e4, 5000)
    ss_eq = kohler_curve(rs_arr, r_solute)
    return rs_arr[np.nanargmax(ss_eq)]

r_crit_map = {aero: critical_radius(r_solute_map[aero]) for aero in aerosols}

ic_markers = {"haze": "o", "critical": "s", "growing": "^", "large": "D"}
ic_ls_map = {"haze": "-", "critical": "--", "growing": "-.", "large": ":"}

fig2, (ax_abs, ax_rel) = plt.subplots(1, 2, figsize=(14, 6))

for aero in aerosols:
    color = aero_colors[aero]
    for ic_name in ics:
        sub = diff[(diff["aerosol"] == aero) & (diff["ic"] == ic_name)]
        fine = sub[sub["step"] <= 20]
        mk = ic_markers[ic_name]
        ls_style = ic_ls_map[ic_name]

        label = f"{aero_nice[aero]} {ic_name}"
        ax_abs.plot(fine["step"].values, fine["dr_ros3"].values * 1e6,
                    color=color, marker=mk, linestyle=ls_style,
                    markersize=4, lw=1.0, label=label)
        ax_rel.plot(fine["step"].values, fine["rel_dr"].values,
                    color=color, marker=mk, linestyle=ls_style,
                    markersize=4, lw=1.0, label=label)

ax_abs.axhline(0, color="k", lw=0.5, ls="--")
ax_abs.set_xlabel("Step")
ax_abs.set_ylabel("Δr = r_ROS3 − r_RK45  [µm]")
ax_abs.set_title("Absolute radius difference")
ax_abs.ticklabel_format(axis="y", style="scientific", scilimits=(-2, 2))
ax_abs.legend(fontsize=6, loc="best", ncol=2)

ax_rel.axhline(0, color="k", lw=0.5, ls="--")
ax_rel.set_xlabel("Step")
ax_rel.set_ylabel("Relative Δr/r")
ax_rel.set_title("Relative radius difference")
ax_rel.ticklabel_format(axis="y", style="scientific", scilimits=(-2, 2))

fig2.suptitle("ROS3 − RK45 Radius Differences (dt=0.01s steps)", fontsize=13)
fig2.tight_layout(rect=(0, 0, 1, 0.95))

diff_png = repo / "test" / "data" / "dgm_solver_diffs.png"
fig2.savefig(diff_png, dpi=150)
print(f"Saved difference plot to {diff_png}")

# Print summary
print("\n" + "=" * 70)
print("ROS3 vs RK45 radius differences")
print("=" * 70)
for aero in aerosols:
    for ic_name in ics:
        sub = diff[(diff["aerosol"] == aero) & (diff["ic"] == ic_name)]
        print(f"\n{aero_nice[aero]}  |  {ic_name}  |  r0 = {sub['r_rk45'].iloc[0]:.3e} m")
        print(f"  {'step':>4s}  {'time':>8s}  {'r_rk45':>12s}  {'dr_ros3':>12s}  {'rel_dr':>12s}")
        print(f"  {'----':>4s}  {'--------':>8s}  {'------------':>12s}  {'------------':>12s}  {'------------':>12s}")
        for _, row in sub.iterrows():
            print(f"  {int(row['step']):4d}  {row['time']:8.4f}  {row['r_rk45']:12.6e}  {row['dr_ros3']:+12.4e}  {row['rel_dr']:+12.4e}")
