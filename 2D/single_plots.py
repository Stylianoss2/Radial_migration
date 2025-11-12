#!/usr/bin/env python3
import os
import argparse
import math
import matplotlib.pyplot as plt

def load_cols(path, t_col=0, y_col=1, max_points=None):
    t, y = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) <= max(t_col, y_col):
                continue
            try:
                t.append(float(parts[t_col]))
                y.append(float(parts[y_col]))
            except ValueError:
                continue
    if max_points and len(t) > max_points:
        step = max(1, len(t) // max_points)
        t = t[::step]
        y = y[::step]
    return t, y

def main():
    p = argparse.ArgumentParser(description="Plot Lz, Energy, and Jacobi vs time with perturbation ramp.")
    p.add_argument("--labels", nargs="+", default=["L0","L1","L2","L3"],
                   help="Star labels to plot (reads orbit_<label>_parameters.dat)")
    p.add_argument("--pert-file", default="pert_strength_time_series.dat",
                   help="File with time and s(t) in cols 1–2")
    p.add_argument("--out-lz", default="Lz_vs_perturbation.png", help="Output image for Lz plot")
    p.add_argument("--out-energy", default="Energy_vs_perturbation.png", help="Output image for Energy plot")
    p.add_argument("--out-jacobi", default="Jacobi_vs_perturbation.png", help="Output image for Jacobi plot")
    p.add_argument("--title-lz", default="Lz of selected stars vs perturbation")
    p.add_argument("--title-energy", default="Energy vs perturbation")
    p.add_argument("--title-jacobi", default="Jacobi integral vs perturbation")
    p.add_argument("--max-points", type=int, default=5000, help="Downsample per series")
    args = p.parse_args()

    # Load Lz series
    series_Lz = []
    for lab in args.labels:
        fname = f"orbit_{lab}_parameters.dat"
        if not os.path.exists(fname):
            print(f"[warn] missing {fname}, skipping")
            continue
        t, Lz = load_cols(fname, 0, 1, max_points=args.max_points)   # col1 = L_z
        if t:
            series_Lz.append((lab, t, Lz))
    if not series_Lz:
        raise SystemExit("No Lz data loaded. Check files like orbit_L0_parameters.dat.")

    # Load perturbation ramp s(t)
    s_t, s_val = [], []
    if os.path.exists(args.pert_file):
        s_t, s_val = load_cols(args.pert_file, 0, 1, max_points=args.max_points)  # take s(t)
    else:
        print(f"[warn] {args.pert_file} not found; right axis will be empty")

    # ---------- Figure 1: Lz ----------
    fig1, axL1 = plt.subplots(figsize=(9, 5.5))
    for lab, t, Lz in series_Lz:
        axL1.plot(t, Lz, lw=1.2, label=f"{lab}")
    axL1.set_xlabel("Time [Myr]")
    axL1.set_ylabel(r"$L_z$  [kpc$^2$ Myr$^{-1}$]")
    axL1.grid(True, alpha=0.25)

    axR1 = axL1.twinx()
    if s_t and s_val:
        axR1.plot(s_t, s_val, lw=1.8, linestyle="--", color="k", label="Perturbation")
        axR1.set_ylabel("Perturbation strength")
        axR1.set_ylim(0.0, 1.05)
        linesL, labelsL = axL1.get_legend_handles_labels()
        linesR, labelsR = axR1.get_legend_handles_labels()
        axL1.legend(linesL + linesR, labelsL + labelsR, loc="upper left", ncols=2, frameon=False)
    else:
        axL1.legend(loc="upper left", ncols=2, frameon=False)
    plt.title(args.title_lz)
    plt.tight_layout()
    plt.savefig(args.out_lz, dpi=200)
    print(f"Wrote {args.out_lz}")

    # ---------- Figure 2: Energy E ----------
    fig2, axL2 = plt.subplots(figsize=(9, 5.5))
    for lab in args.labels:
        fname = f"orbit_{lab}_parameters.dat"
        if not os.path.exists(fname):
            continue
        tE, E = load_cols(fname, 0, 4, max_points=args.max_points)   # col4 = E
        if tE:
            axL2.plot(tE, E, lw=1.2, label=f"{lab}")
    axL2.set_xlabel("Time [Myr]")
    axL2.set_ylabel(r"Energy $E$  [kpc$^2$ Myr$^{-2}$]")
    axL2.grid(True, alpha=0.25)

    axR2 = axL2.twinx()
    if s_t and s_val:
        axR2.plot(s_t, s_val, lw=1.8, linestyle="--", color="k", label="Perturbation")
        axR2.set_ylabel("Perturbation strength")
        axR2.set_ylim(0.0, 1.05)
        linesL, labelsL = axL2.get_legend_handles_labels()
        linesR, labelsR = axR2.get_legend_handles_labels()
        axL2.legend(linesL + linesR, labelsL + labelsR, loc="upper left", ncols=2, frameon=False)
    else:
        axL2.legend(loc="upper left", ncols=2, frameon=False)
    plt.title(args.title_energy)
    plt.tight_layout()
    plt.savefig(args.out_energy, dpi=200)
    print(f"Wrote {args.out_energy}")

    # ---------- Figure 3: Jacobi E_J ----------
    fig3, axL3 = plt.subplots(figsize=(9, 5.5))
    for lab in args.labels:
        fname = f"orbit_{lab}_parameters.dat"
        if not os.path.exists(fname):
            continue
        tJ, EJ = load_cols(fname, 0, 5, max_points=args.max_points)  # col5 = E_J
        if tJ:
            axL3.plot(tJ, EJ, lw=1.2, label=f"{lab}")
    axL3.set_xlabel("Time [Myr]")
    axL3.set_ylabel(r"Jacobi integral $E_J$  [kpc$^2$ Myr$^{-2}$]")
    axL3.grid(True, alpha=0.25)

    axR3 = axL3.twinx()
    if s_t and s_val:
        axR3.plot(s_t, s_val, lw=1.8, linestyle="--", color="k", label="Perturbation")
        axR3.set_ylabel("Perturbation strength")
        axR3.set_ylim(0.0, 1.05)
        linesL, labelsL = axL3.get_legend_handles_labels()
        linesR, labelsR = axR3.get_legend_handles_labels()
        axL3.legend(linesL + linesR, labelsL + labelsR, loc="upper left", ncols=2, frameon=False)
    else:
        axL3.legend(loc="upper left", ncols=2, frameon=False)
    plt.title(args.title_jacobi)
    plt.tight_layout()
    plt.savefig(args.out_jacobi, dpi=200)
    print(f"Wrote {args.out_jacobi}")

if __name__ == "__main__":
    main()
