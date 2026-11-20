"""
Classical-MD stability analysis (YASARA md_analyze output), any number of
systems (apo/control + N ligand complexes).

Usage:
    ./my_venv/Scripts/python.exe utilities/analizza_stabilita_md.py --paper <slug> [--out DIR]

Covers: global protein stability (TotalEnergy, RMSDCa/Bb/All), ligand pose
stability (RMSDLigMove/RMSDLigConf, systems with has_ligand=true), and
binding-pocket residue flexibility (RMSF, bound chain vs free chains within
the same trimer, paired Wilcoxon per residue). No comparison against the
apo/control system is drawn: each system is analyzed for its own internal
stability.

All dataset-specific values come from papers/<slug>/config/analysis_params.json:
    systems            per-system label/code/color/linestyle/has_ligand
    md.systems         ordered ids to analyze
    md.data_dir        folder with {id}_analysis.tab and {id}_rmsfres.tab
    md.windows_ns      id -> [start, end] equilibration window (researcher-defined,
                       never auto-detected)
    md.duration_ns     total simulated time (x-axis limit, "full" window)
    md.plateau_threshold_A   a metric has plateaued in its window when the
                       linear-trend change over that window is below this (A)
    md.rolling_frames  smoothing window for the plotted curves
    pocket.residues / numbering_offset / bound_chain / free_chains

Output (in --out, default papers/<slug>/bozze/figures):
    fig_energy.png, fig_rmsd.png, fig_ligand_rmsd.png, fig_pocket_rmsf.png
    table_stability_summary.csv, table_ligand_pose.csv,
    table_pocket_rmsf.csv, table_pocket_rmsf_test.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats as spstats

from _common import (load_paper, apply_style, read_analysis_tab, read_rmsf,
                     roll, bar_layout)

P = load_paper(__doc__.split('\n\n')[0])
plt = apply_style()
cfg = P.section('md')
DATA = P.path(cfg['data_dir'])
OUT = P.out

SYSD = P.systems(cfg['systems'])
SYS = [s['id'] for s in SYSD]
LABEL = {s['id']: s['label'] for s in SYSD}
CODE = {s['id']: s['code'] for s in SYSD}
COL = {s['id']: s['color'] for s in SYSD}
STYLE = {s['id']: s['linestyle'] for s in SYSD}
LIG = [s['id'] for s in SYSD if s['has_ligand']]
WIN = {k: tuple(v) for k, v in cfg['windows_ns'].items()}
DUR = cfg['duration_ns']
ROLL_N = cfg.get('rolling_frames', 21)
PLATEAU_THRESHOLD_A = cfg.get('plateau_threshold_A', 1.0)
POCKET, OFFSET, pk_cfg = P.pocket()
BOUND = pk_cfg.get('bound_chain', 'A')
FREE = pk_cfg.get('free_chains', ['B', 'C'])
POCKET_NUMS = [n for n, _ in POCKET]

# ---------------------------------------------------------------- load
A, Rf = {}, {}
for s in SYS:
    A[s] = read_analysis_tab(os.path.join(DATA, f'{s}_analysis.tab'), drop_summary_rows=True)
    Rf[s] = read_rmsf(os.path.join(DATA, f'{s}_rmsfres.tab'))


def series(s, k, mask_energy=True):
    y = A[s][k].copy()
    if k == 'TotalEnergy' and mask_energy:
        y = y.mask(y.abs() > 1e7)  # drop non-physical startup spikes
    return A[s]['Time[ns]'].values, y.values


def eqmask(s):
    t = A[s]['Time[ns]'].values
    a, b = WIN[s]
    return (t >= a) & (t <= b)


def trend(t, y):
    """Linear trend over a window: slope (unit/ns), total change over the
    window, and the p-value of the slope (OLS, scipy.stats.linregress)."""
    ok = ~np.isnan(y)
    t, y = t[ok], y[ok]
    if len(t) < 3:
        return np.nan, np.nan, np.nan
    r = spstats.linregress(t, y)
    return r.slope, r.slope * (t.max() - t.min()), r.pvalue


def windows(s):
    return (('0-%g ns (full)' % DUR, np.ones(len(A[s]), bool)),
            (f'{WIN[s][0]}-{WIN[s][1]} ns (eq.)', eqmask(s)))


# ---------------------------------------------------------------- table 1: global protein stability
metrics = ['TotalEnergy', 'RMSDCa', 'RMSDBb', 'RMSDAll']
recs = []
for s in SYS:
    for tag, m in windows(s):
        r = {'system': LABEL[s], 'code': CODE[s], 'window': tag, 'n_frames': int(m.sum())}
        for k in metrics:
            t, y = series(s, k)
            v = pd.Series(y[m]).dropna()
            r[f'{k}_mean'] = v.mean()
            r[f'{k}_sd'] = v.std()
            slope, dtot, p = trend(t[m], y[m])
            r[f'{k}_trend_change'] = dtot
            r[f'{k}_trend_p'] = p
            if k != 'TotalEnergy':
                r[f'{k}_plateau(<{PLATEAU_THRESHOLD_A}A)'] = abs(dtot) < PLATEAU_THRESHOLD_A
        recs.append(r)
tab = pd.DataFrame(recs)
tab.to_csv(os.path.join(OUT, 'table_stability_summary.csv'), index=False)

print(tab[['system', 'window', 'n_frames', 'TotalEnergy_mean', 'TotalEnergy_sd',
           'RMSDCa_mean', 'RMSDCa_trend_change', 'RMSDBb_trend_change',
           'RMSDAll_trend_change']].to_string(index=False))

# ---------------------------------------------------------------- table 2: ligand pose stability
recs2 = []
for s in LIG:
    for tag, m in windows(s):
        r = {'system': LABEL[s], 'code': CODE[s], 'window': tag, 'n_frames': int(m.sum())}
        for k in ['RMSDLigMove', 'RMSDLigConf']:
            t, y = series(s, k)
            v = pd.Series(y[m]).dropna()
            r[f'{k}_mean'] = v.mean()
            r[f'{k}_sd'] = v.std()
            slope, dtot, p = trend(t[m], y[m])
            r[f'{k}_trend_change'] = dtot
            r[f'{k}_trend_p'] = p
            r[f'{k}_plateau(<{PLATEAU_THRESHOLD_A}A)'] = abs(dtot) < PLATEAU_THRESHOLD_A
        recs2.append(r)
tab_pose = pd.DataFrame(recs2)
tab_pose.to_csv(os.path.join(OUT, 'table_ligand_pose.csv'), index=False)
print()
print(tab_pose.to_string(index=False))

# ---------------------------------------------------------------- table 3+4: pocket RMSF
pkd = {}
for s in SYS:
    a = Rf[s][(Rf[s].chain == BOUND) & (Rf[s].num_i.isin(POCKET_NUMS))]
    pkd[LABEL[s]] = a.set_index('num_i')['rmsf']
tab_pocket = pd.DataFrame(pkd)
tab_pocket.insert(0, 'residue_model', [f'{aa}{n}' for n, aa in POCKET])
tab_pocket.insert(1, 'residue_spike', [f'{aa}{n + OFFSET}' for n, aa in POCKET])
tab_pocket.to_csv(os.path.join(OUT, 'table_pocket_rmsf.csv'), index=False)

rec3 = []
for s in SYS:
    d = Rf[s]
    pockA = d[(d.chain == BOUND) & (d.num_i.isin(POCKET_NUMS))].set_index('num_i')['rmsf']
    free = [d[(d.chain == c) & (d.num_i.isin(POCKET_NUMS))].set_index('num_i')['rmsf']
            for c in FREE]
    pockBC = pd.concat(free, axis=1).mean(axis=1)
    common = pockA.index.intersection(pockBC.index)
    wstat, wp = spstats.wilcoxon(pockA.loc[common], pockBC.loc[common])
    rec3.append(dict(system=LABEL[s], code=CODE[s], n_residues=len(common),
                     chainA_pocket_mean=pockA.loc[common].mean(),
                     chainBC_pocket_mean=pockBC.loc[common].mean(),
                     median_diff_A_minus_BC=(pockA.loc[common] - pockBC.loc[common]).median(),
                     wilcoxon_stat=wstat, wilcoxon_p=wp))
tab_rmsf_test = pd.DataFrame(rec3)
tab_rmsf_test.to_csv(os.path.join(OUT, 'table_pocket_rmsf_test.csv'), index=False)
print()
print(f'Wilcoxon signed-rank test, pocket-residue RMSF: chain {BOUND} (bound in ligand '
      f'systems) vs mean(chains {", ".join(FREE)}), same trimer, paired by residue:')
print(tab_rmsf_test.round(4).to_string(index=False))

# ---------------------------------------------------------------- fig 1: TotalEnergy
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for s in SYS:
    t, y = series(s, 'TotalEnergy')
    ax.plot(t, y, color=COL[s], lw=0.35, alpha=0.25)
    ax.plot(t, roll(y, ROLL_N), color=COL[s], lw=1.4, ls=STYLE[s], label=LABEL[s])
    ax.axvline(WIN[s][0], color=COL[s], lw=0.6, ls=':', alpha=0.5)
ax.set_xlabel('Time (ns)')
ax.set_ylabel(r'Total potential energy (kcal mol$^{-1}$)')
ax.set_title('Total potential energy vs. simulation time')
ax.set_xlim(0, DUR)
ax.legend(loc='best', fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_energy.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# ---------------------------------------------------------------- fig 2: RMSD (Ca, Bb, All)
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=False)
for ax, k, ttl in zip(axes, ['RMSDCa', 'RMSDBb', 'RMSDAll'],
                       [r'C$\alpha$ RMSD', 'Backbone RMSD', 'All-atom RMSD']):
    for s in SYS:
        t, y = series(s, k)
        ax.plot(t, y, color=COL[s], lw=0.35, alpha=0.22)
        ax.plot(t, roll(y, ROLL_N), color=COL[s], lw=1.4, ls=STYLE[s], label=LABEL[s])
        ax.axvline(WIN[s][0], color=COL[s], lw=0.6, ls=':', alpha=0.5)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'RMSD ($\mathrm{\AA}$)')
    ax.set_title(ttl)
    ax.set_xlim(0, DUR)
axes[0].legend(loc='upper left', fontsize=7.5, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_rmsd.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

# ---------------------------------------------------------------- fig 3: ligand pose stability
if LIG:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    for ax, k, ttl in zip(
            axes, ['RMSDLigMove', 'RMSDLigConf'],
            ['Ligand positional RMSD', 'Ligand internal (conformational) RMSD']):
        for s in LIG:
            t, y = series(s, k)
            ax.plot(t, y, color=COL[s], lw=0.35, alpha=0.25)
            ax.plot(t, roll(y, ROLL_N), color=COL[s], lw=1.4, label=LABEL[s])
            ax.axvline(WIN[s][0], color=COL[s], lw=0.6, ls=':', alpha=0.5)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel(r'RMSD ($\mathrm{\AA}$)')
        ax.set_title(ttl)
        ax.set_xlim(0, DUR)
    axes[0].legend(loc='upper left', fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_ligand_rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

# ---------------------------------------------------------------- fig 4: pocket RMSF
fig, ax = plt.subplots(figsize=(9.5, 4.2))
x = np.arange(len(POCKET))
w, offs = bar_layout(len(SYS))
for s, off in zip(SYS, offs):
    vals = [tab_pocket.loc[tab_pocket.residue_model == f'{aa}{n}', LABEL[s]].iloc[0]
            for n, aa in POCKET]
    ax.bar(x + off, vals, w, color=COL[s], label=LABEL[s],
           edgecolor='white', linewidth=0.3)
ax.set_xticks(x)
ax.set_xticklabels([f'{aa}{n}' for n, aa in POCKET], rotation=60, ha='right')
ax.set_ylabel(rf'RMSF ($\mathrm{{\AA}}$), chain {BOUND}')
ax.set_title(f'Binding-pocket residue fluctuation, chain {BOUND} (bound in ligand systems)')
ax.legend(loc='upper right', fontsize=8, frameon=False, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_pocket_rmsf.png'), dpi=300, bbox_inches='tight')
plt.close(fig)

print('\nfigures and tables written to', OUT)
