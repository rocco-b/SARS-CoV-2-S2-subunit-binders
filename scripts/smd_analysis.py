#!/usr/bin/env python3
"""Tables and figures for the steered-MD (SMD) analysis of the apo-like model, ripretinib, and fumagillin.

Run from the repository root:
    python scripts/smd_analysis.py

Inputs (this repository)
    data/04_steered_md/analysis_tables/<system>_analysis_<snapshot>ns.tab
        md_analyze output of each SMD replicate (every 0.25 ns): N-C distance (NCdist_A), Rg, RMSD, energies, cell lengths.
    data/04_steered_md/pocket_image_analysis/<system>_pocket_analysis_<snapshot>ns.tab
        md_analyze_local_S2_POCKET output (every 1 ns): pocket contacts, ligand-pocket distances, terminal segments, TrueExt/TrueGap.
    data/04_steered_md/pocket_image_analysis/<system>_image_dist_<snapshot>ns.csv
        shortest protein-to-periodic-image distance per frame (scripts/image_distance.py).
Systems: apo, 1830lig (ripretinib), 1864lig (fumagillin); snapshots 250, 275, 300 ns of the classical MD.

Outputs
    results/tables/*.tsv, *.csv   (Table 4 rows, Tables S3-S6, per-replicate metrics, image-distance summary)
    results/figures/*.png         (Figure 4, Figure S6, Figure S7)
Requires numpy, pandas, matplotlib.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
SMD = REPO / 'data' / '04_steered_md' / 'analysis_tables'
POCKET = REPO / 'data' / '04_steered_md' / 'pocket_image_analysis'
TAB = REPO / 'results' / 'tables'
FIG = REPO / 'results' / 'figures'
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SYSTEMS = [('apo', 'apo-like', '#4d4d4d'), ('1830lig', 'ripretinib (1830)', '#1f78b4'), ('1864lig', 'fumagillin (1864)', '#d95f02')]
SNAPS = [(250, ':', 'o'), (275, '--', 's'), (300, '-', '^')]
K_KJ_MOL_NM2 = 5.0                    # force constant of the SMD macro, |F| = k * D
PN_PER_UNIT = 0.1660539               # pN per (kJ mol-1 nm-2 * Angstrom)
KCAL_PER_PN_A = 0.014393              # kcal/mol per (pN * Angstrom)
CUTOFF_A = 8.0                        # non-bonded cutoff
WINDOW_PRIMARY = 25.0                 # ns, primary comparison window
POCKET_RES = {'C83': 'Phe83', 'C84': 'Ala84', 'C87': 'Met87', 'C101': 'Leu101', 'C105': 'Gln105', 'C107': 'Leu107', 'C108': 'Ile108',
              'C109': 'Ala109', 'C111': 'Gln111', 'C112': 'Phe112', 'C235': 'Met235', 'C250': 'Val250', 'C251': 'Thr251',
              'C252': 'Tyr252', 'C255': 'Ala255', 'C294': 'Phe294'}
# contacts listed for the docked poses (Table 3 of the manuscript), used only to mark Figure S7
DOCKED_CONTACTS = {
    '1830lig': {'Gln111', 'Thr251', 'Val250', 'Met235', 'Met87', 'Phe112', 'Ala109', 'Leu101', 'Phe83', 'Tyr252'},
    '1864lig': {'Gln105', 'Ala84', 'Phe112', 'Thr251', 'Phe83', 'Leu101', 'Met87', 'Tyr252', 'Ile108', 'Leu107', 'Phe294'}}


def read_tab(path):
    d = pd.read_csv(path, sep=r'\s+')
    d = d[pd.to_numeric(d['Time[ps]'], errors='coerce').notna()].reset_index(drop=True)
    return d.apply(pd.to_numeric, errors='coerce')


def roll(x, n=9):
    return pd.Series(x).rolling(n, center=True, min_periods=1).mean().values


def ms(x, nd=1):
    x = np.asarray(x, float)
    return f'{x.mean():.{nd}f} ± {x.std(ddof=1):.{nd}f}' if len(x) > 1 else f'{x[0]:.{nd}f}'


def rng(x, nd=1):
    x = np.asarray(x, float)
    return f'{x.min():.{nd}f}–{x.max():.{nd}f}' if len(x) > 1 else f'{x[0]:.{nd}f}'


def load_all():
    smd, pocket, image = {}, {}, {}
    for s, _, _ in SYSTEMS:
        for n, _, _ in SNAPS:
            d = read_tab(SMD / f'{s}_analysis_{n}ns.tab')
            smd[(s, n)] = d[d['Time[ns]'] <= 100.0001].reset_index(drop=True)
            pocket[(s, n)] = read_tab(POCKET / f'{s}_pocket_analysis_{n}ns.tab')
            image[(s, n)] = pd.read_csv(POCKET / f'{s}_image_dist_{n}ns.csv')
    return smd, pocket, image


def per_replicate(smd, pocket, image):
    rows = []
    for s, lab, _ in SYSTEMS:
        for n, _, _ in SNAPS:
            d, im = smd[(s, n)], image[(s, n)]
            t, D = d['Time[ns]'].values, d.NCdist_A.values
            F = K_KJ_MOL_NM2 * PN_PER_UNIT * D
            sm = roll(D)
            rg = roll(d.RadGyration.values)
            m25, m50, m80 = t <= 25, t <= 50, t >= 80
            below = np.where(im.min_image_dist.values < CUTOFF_A)[0]
            first = im.time_ns.values[below[0]] if len(below) else np.inf
            V = d.CellLengthX * d.CellLengthY * d.CellLengthZ
            ok = d.Coulomb.abs() <= 1e6
            E = d.TotalEnergy.where(ok)
            rows.append(dict(
                system=s, label=lab, snapshot=n, D0=D[0], F0_pN=F[0], Fmean_0_50_pN=F[m50].mean(), Fmax_pN=F.max(),
                work_kcal_mol=np.trapezoid(F, D) * KCAL_PER_PN_A, rate_0_50_A_ns=(sm[m50][-1] - sm[0]) / 50, rate_0_100_A_ns=(sm[-1] - sm[0]) / 100,
                D_0_25=D[m25].mean(), D_0_50=D[m50].mean(), D_0_100=D.mean(), D_80_100=D[m80].mean(),
                ext_0_25=(D[m25] - D[0]).mean(), ext_0_50=(D[m50] - D[0]).mean(), ext_0_100=(D - D[0]).mean(),
                Rg0=d.RadGyration.iloc[0], Rg_mean_0_25=d.RadGyration[m25].mean(), Rg_max_0_25=rg[m25].max(), Rg_max_0_100=rg.max(),
                cell_x=d.CellLengthX.iloc[0], cell_y=d.CellLengthY.iloc[0], cell_z=d.CellLengthZ.iloc[0], frames=len(d),
                coulomb_artifact_frames=int((~ok).sum()), cell_volume_change_first_ns_pct=100 * (V[(t > 0) & (t <= 1)].mean() / V.iloc[0] - 1),
                energy_change_1_10_vs_90_100_pct=100 * (E[t >= 90].mean() / E[(t > 1) & (t <= 10)].mean() - 1),
                min_image_dist_A=im.min_image_dist.min(), frames_image_lt8=int((im.min_image_dist < CUTOFF_A).sum()), frames_image=len(im),
                first_image_contact_ns=first, clean_0_25=first > 25, clean_0_50=first > 50))
    P = pd.DataFrame(rows)
    P.to_csv(TAB / 'per_replicate_metrics.csv', index=False)
    P[['system', 'snapshot', 'min_image_dist_A', 'frames_image_lt8', 'frames_image', 'first_image_contact_ns']].to_csv(TAB / 'image_distance_summary.csv', index=False)
    return P


def table_s2(P):
    L = ['System\tSnapshot (ns)\tD(0) (Å)\tF(0) (pN)\tMean F, 0–50 ns (pN)\tMaximum F (pN)\tWork along D, 0–100 ns (kcal/mol)\t'
         'Mean extension rate, 0–50 ns (Å/ns)\tMean extension rate, 0–100 ns (Å/ns)']
    for _, r in P.iterrows():
        L.append('\t'.join([r.label, str(r.snapshot), f'{r.D0:.1f}', f'{r.F0_pN:.1f}', f'{r.Fmean_0_50_pN:.1f}', f'{r.Fmax_pN:.1f}',
                            f'{r.work_kcal_mol:.1f}', f'{r.rate_0_50_A_ns:.2f}', f'{r.rate_0_100_A_ns:.2f}']))
    for _, lab, _ in SYSTEMS:
        q = P[P.label == lab]
        L.append('\t'.join([lab + ', mean ± SD', '', ms(q.D0), ms(q.F0_pN), ms(q.Fmean_0_50_pN), ms(q.Fmax_pN), ms(q.work_kcal_mol), ms(q.rate_0_50_A_ns, 2), ms(q.rate_0_100_A_ns, 2)]))
    (TAB / 'TableS3_loading.tsv').write_text('\n'.join(L) + '\n', encoding='utf-8')


def table_s3(P):
    L = ['Window (ns)\tReplicates included\tSystem\tn\tMean D (Å), mean ± SD\tMean D (Å), replicate range\tMean extension D − D(0) (Å), mean ± SD\tMean extension (Å), replicate range']

    def add(window, incl, cd, ce, mask):
        for _, lab, _ in SYSTEMS:
            q = P[(P.label == lab) & mask(P)]
            L.append('\t'.join([window, incl, lab, str(len(q)), ms(q[cd]), rng(q[cd]), ms(q[ce]), rng(q[ce])]))
    allm = lambda p: np.ones(len(p), bool)  # noqa: E731
    add('0–25', 'all replicates', 'D_0_25', 'ext_0_25', allm)
    add('0–25', 'replicates without image contacts (< 8 Å) before 25 ns', 'D_0_25', 'ext_0_25', lambda p: p.clean_0_25)
    add('0–50', 'all replicates', 'D_0_50', 'ext_0_50', allm)
    add('0–50', 'replicates without image contacts before 50 ns', 'D_0_50', 'ext_0_50', lambda p: p.clean_0_50)
    add('0–100', 'all replicates (five of nine with image contacts before 100 ns)', 'D_0_100', 'ext_0_100', allm)
    (TAB / 'TableS4_windows.tsv').write_text('\n'.join(L) + '\n', encoding='utf-8')
    # D(0) dependence in the apo-like replicates (0-25 ns)
    a = P[P.system == 'apo']
    b, c0 = np.polyfit(a.D0, a.ext_0_25, 1)
    bD, _ = np.polyfit(a.D0, a.D_0_25, 1)
    lines = [f'apo-like, 0–25 ns: extension = {b:.2f}·D(0) + {c0:.1f} Å; slope of the mean D on D(0) = {bD:.2f} Å per Å']
    for s, lab, _ in SYSTEMS[1:]:
        q = P[P.system == s]
        res = q.ext_0_25 - (b * q.D0 + c0)
        lines.append(f'{lab}: difference from the apo-like relation (snapshots {", ".join(map(str, q.snapshot))}): ' + ', '.join(f'{x:+.1f}' for x in res) + ' Å')
    (TAB / 'TableS4_D0_dependence_note.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def table_s4(pocket):
    """Table S5: one row per replicate, ligand retention in the pocket and terminal segments."""
    names = list(POCKET_RES.values())
    cols = list(POCKET_RES)
    L = ['System\tSnapshot (ns)\tLigand heavy atoms within 4.5 Å of the pocket, median (minimum) of 33\tLigand–pocket centroid distance, median (maximum) (Å)\t'
         'Pocket residues in contact in at least 70% of the frames\tEnd-to-end distance, residues 1–10, median (range) (Å)\tEnd-to-end distance, residues 330–339, median (range) (Å)']
    occ = {}
    for s, lab, _ in SYSTEMS:
        for n, _, _ in SNAPS:
            t = pocket[(s, n)]
            n1, n2 = t.NTerm_1_10.values, t.CTerm_330339.values
            term = [f'{np.median(n1):.1f} ({n1.min():.1f}–{n1.max():.1f})', f'{np.median(n2):.1f} ({n2.min():.1f}–{n2.max():.1f})']
            if s == 'apo':
                L.append('\t'.join([lab, str(n), '–', '–', '–'] + term))
                continue
            o = (t[cols] > 0).mean()
            occ[(s, n)] = o
            heavy = int(t.LigHeavy.iloc[0])
            near = t.LigInPocket
            res70 = ', '.join(nm for c, nm in POCKET_RES.items() if o[c] >= 0.70)
            L.append('\t'.join([lab, str(n), f'{near.median():.0f} ({near.min():.0f}) of {heavy}', f'{t.LigPockDist.median():.1f} ({t.LigPockDist.max():.1f})', res70] + term))
    (TAB / 'TableS5_ligand_retention_terminals.tsv').write_text('\n'.join(L) + '\n', encoding='utf-8')
    O = pd.DataFrame({f'{s}_{n}': v.rename(POCKET_RES) for (s, n), v in occ.items()})
    O.to_csv(TAB / 'contact_occupancy_per_replicate.csv')
    return O


def table_s5(P):
    L = ['System\tSnapshot (ns)\tCell X (Å)\tCell Y (Å)\tCell Z (Å)\tFrames\tFrames with non-physical Coulomb energy\tCell volume change, first ns (%)\t'
         'Energy change, 1–10 ns vs 90–100 ns (%)\tMinimum distance between protein and periodic images (Å)\tFrames (1 ns) with this distance < 8 Å\tFirst such frame (ns)']
    for _, r in P.iterrows():
        first = f'{r.first_image_contact_ns:.0f}' if np.isfinite(r.first_image_contact_ns) else '–'
        L.append('\t'.join([r.label, str(r.snapshot), f'{r.cell_x:.1f}', f'{r.cell_y:.1f}', f'{r.cell_z:.1f}', str(r.frames), str(r.coulomb_artifact_frames),
                            f'{r.cell_volume_change_first_ns_pct:.1f}', f'{r.energy_change_1_10_vs_90_100_pct:.2f}', f'{r.min_image_dist_A:.1f}',
                            f'{r.frames_image_lt8}/{r.frames_image}', first]))
    (TAB / 'TableS6_cell_images.tsv').write_text('\n'.join(L) + '\n', encoding='utf-8')


def table4_rows(P):
    L = ['Metric\tapo-like\tripretinib (1830)\tfumagillin (1864)']
    f = lambda col, sel=None, nd=1: [ms((P[P.label == lab] if sel is None else P[(P.label == lab) & sel(P)])[col], nd) for _, lab, _ in SYSTEMS]  # noqa: E731
    L.append('F(0) (pN)\t' + '\t'.join(f('F0_pN')))
    L.append('Mean D, 0–25 ns (Å)\t' + '\t'.join(f('D_0_25')))
    clean = lambda p: p.clean_0_25  # noqa: E731
    ns = [str(int(P[(P.label == lab) & P.clean_0_25].shape[0])) for _, lab, _ in SYSTEMS]
    L.append('Mean D, 0–25 ns, replicates without image contacts (Å)\t' + '\t'.join(f'{v} (n = {n})' for v, n in zip(f('D_0_25', clean), ns)))
    L.append('Mean extension, 0–25 ns (Å)\t' + '\t'.join(f('ext_0_25')))
    L.append('Mean D, 0–50 ns (Å)\t' + '\t'.join(f('D_0_50')))
    (TAB / 'Table4_new_rows.tsv').write_text('\n'.join(L) + '\n', encoding='utf-8')


def figure4(smd, image, P):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.6))
    (a, b), (c, d) = ax
    first = {(r.system, r.snapshot): r.first_image_contact_ns for r in P.itertuples()}
    for s, lab, col in SYSTEMS:
        L = min(len(smd[(s, n)]) for n, _, _ in SNAPS)
        for key, axx in [('NCdist_A', a), ('RadGyration', b)]:
            ys = []
            for n, ls, _ in SNAPS:
                t = smd[(s, n)]
                y = roll(t[key].values)
                axx.plot(t['Time[ns]'], y, color=col, ls=ls, lw=0.9, alpha=0.7)
                ys.append(y[:L])
                if np.isfinite(first[(s, n)]):
                    i = int(np.argmin(np.abs(t['Time[ns]'].values - first[(s, n)])))
                    axx.plot(first[(s, n)], y[i], 'x', color='k', ms=6, mew=1.3, zorder=5)
            axx.plot(smd[(s, 250)]['Time[ns]'].values[:L], np.mean(ys, axis=0), color=col, lw=2.4, label=lab)
    for axx in (a, b):
        axx.axvspan(0, WINDOW_PRIMARY, color='0.9', zorder=0)
        axx.set_xlim(0, 100)
        axx.set_xlabel('Pulling time (ns)')
    a.set_ylabel('N–C Cα distance, D (Å)')
    b.set_ylabel('Radius of gyration, Rg (Å)')
    a.plot([], [], 'x', color='k', label='first image contact (< 8 Å)')
    a.legend(frameon=False, fontsize=8, loc='lower right')
    for axx, key, yl, ttl in [(c, 'NCdist_A', 'D (Å)', 'Filled: mean D, 0–25 ns; open: D(0)'), (d, 'RadGyration', 'Rg (Å)', 'Filled: mean Rg, 0–25 ns; open: Rg(0)')]:
        for i, (s, lab, col) in enumerate(SYSTEMS):
            m25 = []
            for j, (n, _, mk) in enumerate(SNAPS):
                t = smd[(s, n)]
                v = t[key][t['Time[ns]'] <= WINDOW_PRIMARY].mean()
                m25.append(v)
                x = i + (j - 1) * 0.12
                affected = first[(s, n)] <= WINDOW_PRIMARY
                axx.plot(x, v, mk, color=col, ms=6.5, mec='k' if affected else col, mew=1.6 if affected else 1)
                axx.plot(x, t[key].iloc[0], mk, mfc='white', color=col, ms=6.5)
                axx.plot([x, x], [t[key].iloc[0], v], color=col, lw=0.6, alpha=0.6)
                if affected:
                    axx.annotate('†', (x, v), xytext=(5, 3), textcoords='offset points', fontsize=11)
            axx.hlines(np.mean(m25), i - 0.3, i + 0.3, color=col, lw=2)
        axx.set_xticks(range(3))
        axx.set_xticklabels([x[1] for x in SYSTEMS], fontsize=8)
        axx.set_ylabel(yl)
        axx.set_title(ttl, fontsize=8)
    for axx, lt in zip([a, b, c, d], 'abcd'):
        axx.text(-0.12, 1.05, f'({lt})', transform=axx.transAxes, fontweight='bold')
    fig.tight_layout()
    fig.savefig(FIG / 'Figure4_SMD_distance_Rg.png', dpi=300)
    plt.close(fig)


def figure_s5(smd):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    ax_d, ax_f, ax_w = axes
    for s, lab, col in SYSTEMS:
        for n, ls, _ in SNAPS:
            d = smd[(s, n)]
            D = d.NCdist_A.values
            F = K_KJ_MOL_NM2 * PN_PER_UNIT * D
            W = np.concatenate([[0], np.cumsum((F[1:] + F[:-1]) / 2 * np.diff(D))])
            ax_d.plot(d['Time[ns]'], roll(D), color=col, ls=ls, lw=1.4)
            ax_f.plot(d['Time[ns]'], roll(F), color=col, ls=ls, lw=1.4)
            ax_w.plot(roll(D - D[0]), roll(W), color=col, ls=ls, lw=1.4)
    ax_d.set_ylabel('N–C Cα distance, D (Å)')
    ax_f.set_ylabel('Applied force, F = k·D (pN)')
    ax_w.set_ylabel('Cumulative work along D (pN·Å)')
    ax_w.set_xlabel('Extension from run start, D − D(0) (Å)')
    for axx in (ax_d, ax_f):
        axx.set_xlabel('Pulling time (ns)')
        axx.set_xlim(0, 100)
    handles = [plt.Line2D([], [], color=c, lw=1.6, label=lab) for _, lab, c in SYSTEMS] + \
              [plt.Line2D([], [], color='black', lw=1.4, ls=ls, label=f'{n} ns snapshot') for n, ls, _ in SNAPS]
    fig.legend(handles=handles, loc='lower center', ncol=6, bbox_to_anchor=(0.5, -0.06), fontsize=8, frameon=False)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(FIG / 'FigureS6_force_work.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def figure_s6(O):
    """Heat map of the contact occupancy; residues listed for the docked pose are marked with an asterisk."""
    resnames = list(POCKET_RES.values())
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 6.4), sharey=False, gridspec_kw={'wspace': 0.42})
    for axx, (s, title) in zip(axes, [('1830lig', 'ripretinib (1830)'), ('1864lig', 'fumagillin (1864)')]):
        cols = [f'{s}_{n}' for n, _, _ in SNAPS]
        M = np.array([[O[c][r] for c in cols] for r in resnames])
        im = axx.imshow(M, cmap='Blues', vmin=0, vmax=1, aspect='auto')
        axx.set_xticks(range(3))
        axx.set_xticklabels([f'{n} ns' for n, _, _ in SNAPS], fontsize=8)
        axx.set_yticks(range(len(resnames)))
        axx.set_yticklabels([r + (' *' if r in DOCKED_CONTACTS[s] else '') for r in resnames], fontsize=8)
        for i in range(len(resnames)):
            for j in range(3):
                axx.text(j, i, f'{100 * M[i, j]:.0f}', ha='center', va='center', fontsize=7.5, color='white' if M[i, j] > 0.6 else 'black')
        axx.set_title(title, fontsize=9)
        axx.set_xlabel('Snapshot of the classical MD', fontsize=8)
    cb = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.03)
    cb.set_label('Fraction of frames in contact (heavy atoms within 4.5 Å)')
    fig.suptitle('Contact occupancy of the pocket residues (percent of frames; * = residue listed for the docked pose, Table 3)', fontsize=8.5)
    fig.savefig(FIG / 'FigureS7_contact_occupancy.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    smd, pocket, image = load_all()
    P = per_replicate(smd, pocket, image)
    table_s2(P)
    table_s3(P)
    O = table_s4(pocket)
    table_s5(P)
    table4_rows(P)
    figure4(smd, image, P)
    figure_s5(smd)
    figure_s6(O)
    print('tables ->', TAB.relative_to(REPO), '| figures ->', FIG.relative_to(REPO))


if __name__ == '__main__':
    main()
