"""
SMD mechanical readout along the pulling coordinate itself: D(t), the
N-/C-terminal C-alpha distance (column NCdist of a dense per-frame
md_analyze table), plus F(t) and cumulative work vs. extension derived from
it. Complements analizza_smd.py, whose protein RMSD is not the variable
conjugate to the applied SMD force.

Usage:
    ./my_venv/Scripts/python.exe utilities/analizza_smd_meccanica.py --paper <slug> [--out DIR]

F(t) is NOT read from a logged column but derived analytically from D(t) with
the deterministic force law of the SMD macro used for the runs. Two laws are
supported (config smd_mechanics.force_law):
  "extension" (default): F = k * (D - D0)  (force zero at t = 0)
  "absolute":            F = k * D         (Reilly-like custom bond force)
  F [pN] = k [kJ/mol/nm^2] * x [A] * 1e23/N_A, x = extension or D, N_A =
  Avogadro's constant. With a force that is a function of D only, F(t) and the
  cumulative work are fully determined by D(t): they add no independent
  information, they only express the load actually applied.
  Valid only below the macro's acceleration cap (accmax; 0 or absent = no cap,
  check skipped). If the macro's force law differs from both, do not use it.

Requires a dense per-frame table containing NCdist. The sparse
"new-maximum-only" logs written by some SMD macros are not suitable.

Dataset-specific values come from papers/<slug>/config/analysis_params.json,
section "smd_mechanics":
    systems            ordered ids (defined in the top-level "systems")
    data_dir           folder with {system id}_analysis_{replicate id}.tab
    replicates         list of {id, label, linestyle}
    bondforce_const_kJ_mol_nm2, accmax_pm_ps2 (0/absent = no cap)
    force_law            "extension" (default) or "absolute"
    ncdist_column        name of the D column (default "NCdist"; "NCdist_A" for monomer tables)
    sanity_ligand_mass_gmol   approximate mass for the cap check (default 36.033)
    pull_duration_ns, rolling_frames

Output (in --out, default papers/<slug>/bozze/figures):
    smd_fig_forceext.png, smd_table_forceext_summary.csv
"""
import os
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from _common import load_paper, apply_style, read_analysis_tab, roll

P = load_paper(__doc__.split('\n\n')[0])
plt = apply_style()
cfg = P.section('smd_mechanics')
DATA = P.path(cfg['data_dir'])
OUT = P.out

SYSD = P.systems(cfg['systems'])
SYS = [s['id'] for s in SYSD]
LABEL = {s['id']: s['label'] for s in SYSD}
COL = {s['id']: s['color'] for s in SYSD}
REPS = [r['id'] for r in cfg['replicates']]
REP_LABEL = {r['id']: r['label'] for r in cfg['replicates']}
REP_STYLE = {r['id']: r.get('linestyle', '-') for r in cfg['replicates']}
DUR = cfg['pull_duration_ns']
ROLL_N = cfg.get('rolling_frames', 9)

BONDFORCE_CONST = cfg['bondforce_const_kJ_mol_nm2']
ACCMAX_PM_PS2 = cfg.get('accmax_pm_ps2', 0)
FORCE_LAW = cfg.get('force_law', 'extension')
DCOL = cfg.get('ncdist_column', 'NCdist')
if FORCE_LAW not in ('extension', 'absolute'):
    raise SystemExit("force_law must be 'extension' or 'absolute'")
LIGMASS_APPROX = cfg.get('sanity_ligand_mass_gmol', 3 * 12.011)
AVOGADRO = 6.02214076e23
PN_PER_KJMOLNM2_A = 1e23 / AVOGADRO   # F[pN] = BONDFORCE_CONST * extension[A] * this

# ---------------------------------------------------------------- load
runs = {}
for s in SYS:
    for rep in REPS:
        d = read_analysis_tab(os.path.join(DATA, f'{s}_analysis_{rep}.tab'),
                              drop_summary_rows=True)
        d = d.sort_values('Time[ns]').reset_index(drop=True)
        d['NCdist'] = d[DCOL]
        d['Extension'] = d['NCdist'] - d['NCdist'].iloc[0]
        load = d['Extension'] if FORCE_LAW == 'extension' else d['NCdist']
        d['Force_pN'] = BONDFORCE_CONST * load * PN_PER_KJMOLNM2_A
        # trapezoidal cumulative work, Force[pN] * dExtension[A]
        dd = np.diff(d['Extension'].values, prepend=d['Extension'].values[0])
        favg = np.concatenate([[0], (d['Force_pN'].values[1:] + d['Force_pN'].values[:-1]) / 2])
        d['CumWork_pN_A'] = np.cumsum(favg * dd)
        runs[(s, rep)] = d

# ---------------------------------------------------------------- table: summary
recs = []
for s in SYS:
    for rep in REPS:
        d = runs[(s, rep)]
        recs.append(dict(
            system=LABEL[s], replicate=REP_LABEL[rep], n_frames=len(d),
            ncdist_start_A=d['NCdist'].iloc[0], ncdist_end_A=d['NCdist'].iloc[-1],
            extension_end_A=d['Extension'].iloc[-1],
            extension_max_A=d['Extension'].max(),
            force_start_pN=d['Force_pN'].iloc[0],
            force_end_pN=d['Force_pN'].iloc[-1],
            force_max_pN=d['Force_pN'].max(),
            cumwork_end_pN_A=d['CumWork_pN_A'].iloc[-1],
            cumwork_end_kcal_mol=d['CumWork_pN_A'].iloc[-1] * 0.014393,
        ))
tab = pd.DataFrame(recs)
tab.to_csv(os.path.join(OUT, 'smd_table_forceext_summary.csv'), index=False)
print(tab.round(3).to_string(index=False))

# Sanity check: distance of the observed extensions from the acceleration cap.
if ACCMAX_PM_PS2 > 0:
    key = 'Extension' if FORCE_LAW == 'extension' else 'NCdist'
    max_ext = max(runs[(s, rep)][key].max() for s in SYS for rep in REPS)
    accel_at_max = BONDFORCE_CONST * max_ext * 100 / LIGMASS_APPROX
    print(f'\nMax observed load coordinate across all runs: {max_ext:.2f} A -> '
          f'approx. acceleration {accel_at_max:.1f} pm/ps^2 '
          f'(cap is {ACCMAX_PM_PS2:.0f} pm/ps^2, '
          f'{"so the cap is not reached" if accel_at_max < ACCMAX_PM_PS2 else "CAP REACHED: F(t) formula not valid"}).')
else:
    print('\nNo acceleration cap configured (accmax_pm_ps2 = 0/absent): cap check skipped.')

# ---------------------------------------------------------------- fig: D(t), F(t), work vs extension
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
ax_d, ax_f, ax_w = axes

for s in SYS:
    for rep in REPS:
        d = runs[(s, rep)]
        ls = REP_STYLE[rep]
        ax_d.plot(d['Time[ns]'], roll(d['NCdist'], ROLL_N), color=COL[s], lw=1.4, ls=ls)
        ax_f.plot(d['Time[ns]'], roll(d['Force_pN'], ROLL_N), color=COL[s], lw=1.4, ls=ls)
        ax_w.plot(roll(d['Extension'], ROLL_N), roll(d['CumWork_pN_A'], ROLL_N),
                  color=COL[s], lw=1.4, ls=ls)

ax_d.set_xlabel('Pulling time (ns)')
ax_d.set_ylabel(r'N-/C-terminal C$\alpha$ distance, $D$ ($\mathrm{\AA}$)')
ax_d.set_title('Force-conjugate extension, D(t)')
ax_d.set_xlim(0, DUR)

ax_f.set_xlabel('Pulling time (ns)')
ax_f.set_ylabel('Applied force, F (pN)')
ax_f.set_title('Applied force, F(t) (derived from D(t))')
ax_f.set_xlim(0, DUR)

ax_w.set_xlabel(r'Extension from run start, $D - D_0$ ($\mathrm{\AA}$)')
ax_w.set_ylabel('Cumulative work (pN·Å)')
ax_w.set_title('Work vs. extension')

handles = [Line2D([], [], color=COL[s], lw=1.6, label=LABEL[s]) for s in SYS] + \
          [Line2D([], [], color='black', lw=1.4, ls=REP_STYLE[r], label=REP_LABEL[r]) for r in REPS]
fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.08), fontsize=8)
fig.tight_layout(rect=[0, 0.08, 1, 1])
fig.savefig(os.path.join(OUT, 'smd_fig_forceext.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

print('\nfigure and table written to', OUT)
