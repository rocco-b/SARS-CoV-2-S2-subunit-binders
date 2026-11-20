#!/usr/bin/env python3
"""Core-distance analysis of the SMD replicates (revision, Table S2 columns and Table S4 rows on core distances).

Reads data/04_steered_md/core_distance_analysis/<system>_core_analysis_<snapshot>ns.tab (written by
macros/md_analyze_local_S2_CORE.mcr, one frame per ns) and writes, in results/tables/:
  core_distance_per_replicate.csv   per replicate: mean 0-25 ns of the Calpha distances 1-339, 5-335, 10-330, 20-320,
                                    the values at t = 0, and the end-to-end distances of residues 1-10 and 330-339
  core_distance_windows.csv         mean +- SD (n - 1) over the three replicates of each system, for every distance and
                                    window, with the replicate values and whether all replicates of the system lie
                                    below all apo-like replicates
Usage (from the repository root):  python scripts/core_distance_analysis.py
"""
import os

import numpy as np
import pandas as pd

SYSTEMS = [("apo-like", "apo"), ("ripretinib (1830)", "1830lig"), ("fumagillin (1864)", "1864lig")]
SNAPSHOTS = (250, 275, 300)
DIST = {"D(1-339)": "NCdist_A", "D(5-335)": "D_Ca5_335", "D(10-330)": "D_Ca10_330", "D(20-320)": "D_Ca20_320"}
WINDOWS = [(0, 25), (0, 50), (0, 100), (25, 50), (50, 100)]
IN = os.path.join("data", "04_steered_md", "core_distance_analysis")
OUT = os.path.join("results", "tables")


def load(prefix, snap):
    d = pd.read_csv(os.path.join(IN, f"{prefix}_core_analysis_{snap}ns.tab"), sep=r"\s+", engine="python")
    for c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["Time[ns]"])


def window(d, lo, hi, col):
    return d[(d["Time[ns]"] >= lo) & (d["Time[ns]"] <= hi)][col].mean()


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {(s, k): load(p, k) for s, p in SYSTEMS for k in SNAPSHOTS}
    rows = []
    for s, _ in SYSTEMS:
        for k in SNAPSHOTS:
            d = data[(s, k)]
            r = {"system": s, "snapshot_ns": k}
            for name, col in DIST.items():
                r[f"mean_{name}_0_25ns"] = round(window(d, 0, 25, col), 2)
                r[f"{name}_t0"] = round(d[col].iloc[0], 2)
            r["mean_d(1-10)_0_25ns"] = round(window(d, 0, 25, "NTerm_1_10"), 2)
            r["d(1-10)_t0"] = round(d["NTerm_1_10"].iloc[0], 2)
            r["mean_d(330-339)_0_25ns"] = round(window(d, 0, 25, "CTerm_330339"), 2)
            rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "core_distance_per_replicate.csv"), index=False)

    out = []
    for name, col in DIST.items():
        for lo, hi in WINDOWS:
            per = {s: np.array([window(data[(s, k)], lo, hi, col) for k in SNAPSHOTS]) for s, _ in SYSTEMS}
            apo_min = per["apo-like"].min()
            for s, _ in SYSTEMS:
                v = per[s]
                out.append({
                    "distance": name, "window_ns": f"{lo}-{hi}", "system": s,
                    "mean": round(v.mean(), 1), "sd": round(v.std(ddof=1), 1),
                    "replicates": " / ".join(f"{x:.1f}" for x in v),
                    "all_below_all_apo": "" if s == "apo-like" else bool(v.max() < apo_min),
                    "difference_from_apo": "" if s == "apo-like" else round(per["apo-like"].mean() - v.mean(), 1),
                })
    pd.DataFrame(out).to_csv(os.path.join(OUT, "core_distance_windows.csv"), index=False)
    print("written:", os.path.join(OUT, "core_distance_per_replicate.csv"), os.path.join(OUT, "core_distance_windows.csv"))


if __name__ == "__main__":
    main()
