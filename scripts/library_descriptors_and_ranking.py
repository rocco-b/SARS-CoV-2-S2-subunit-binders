#!/usr/bin/env python3
"""Size and flexibility of the docked library, and robustness of the hit list to the scoring criterion.

Run from the repository root:
    python scripts/library_descriptors_and_ranking.py

Question addressed (Reviewer 3): the three reference compounds of Reilly et al. fail CNN_VS >= 6.3 because of a low
CNN pose-confidence score rather than a low CNN affinity, which could mean that the filter is biased against large,
flexible scaffolds; would the prioritized set survive a ranking by the affinity component alone?

Inputs (all in data/02_virtual_screening/):
  docking_scores_all_2342.csv                    top-pose scores of the 2,342 docked records (molecule_id = position in the SDF)
  ligands_2342.sdf                               the docking input (the 3D structures, with explicit hydrogens)
  primary_hits_30.csv, candidates_after_in_situ_evaluation.csv
  reference_compounds_redocking/*_log.txt        GNINA logs of bemcentinib, WYS-633, WYS-694 (mode 1) and *_ligand.sdf
  excluded_compound_odanacatib/..._excerpt.tsv   scores of the neutral form of odanacatib (compound 2006), mode 1

Odanacatib (2006): the CSV holds the scores of the wrongly protonated form (CNN_VS 6.51); here the record is replaced
by the neutral re-docking (CNN score 0.8375, CNN affinity 6.546, CNN_VS 5.48), so that it is not a hit.

Descriptors (RDKit; hydrogens not counted):
  heavy atoms                  number of non-hydrogen atoms of the record
  rotatable bonds              rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True), the default
                               definition (terminal and triple bonds, amide C-N bonds and other restricted bonds not
                               counted); the strict=False count (every non-ring single bond between two non-terminal
                               heavy atoms) is reported as a sensitivity analysis.
  68 records fail default RDKit sanitization (non-standard valences); for them the descriptors are computed after a
  partial sanitization (aromaticity, ring perception, hybridization). Every statistic is also given without them.

Correlations are Spearman rank correlations (rho, two-sided p, 95% bootstrap interval from 2,000 resamples, seed 1).
Ranks use method 'min' (tied values share the best rank); "top x%" is rank / number of records.

Writes:
  results/tables/library_descriptors_2342.csv             one row per record
  results/tables/TableS1_descriptor_columns.tsv           heavy atoms and rotatable bonds of the 30 primary hits (extra columns for Table S1)
  results/tables/library_descriptors_and_ranking.txt      every number quoted in Section 2.1 and in the response
  results/figures/FigureS1_CNN_score_distribution.png     six panels (a-f); a and b are the original panels
Requires numpy, pandas, scipy, matplotlib and rdkit.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from scipy.stats import spearmanr

RDLogger.DisableLog('rdApp.*')
REPO = Path(__file__).resolve().parents[1]
VS = REPO / 'data' / '02_virtual_screening'
REF = VS / 'reference_compounds_redocking'
OUT_T = REPO / 'results' / 'tables'
OUT_F = REPO / 'results' / 'figures'
THRESHOLD = 6.3
ODANACATIB = 2006
ADVANCED = [1830, 1864, 2168, 2195]          # the four compounds of Table 2
NAMES = {1830: 'ripretinib', 1864: 'fumagillin', 2168: 'afoxolaner', 2195: 'dexamethasone palmitate'}
REFERENCES = {'bemcentinib': 'Bemcentinib', 'WYS-633': 'WYS-633', 'WYS-694': 'WYS-694'}
rng = np.random.default_rng(1)
OUT = []


def say(*a):
    line = ' '.join(str(x) for x in a)
    print(line)
    OUT.append(line)


# ---------------------------------------------------------------- inputs
def mol_descriptors(m):
    heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
    ok = True
    try:
        mm = Chem.Mol(m)
        Chem.SanitizeMol(mm)
    except Exception:
        ok = False
        mm = Chem.Mol(m)
        mm.UpdatePropertyCache(strict=False)
        ops = (Chem.SANITIZE_SETAROMATICITY | Chem.SANITIZE_SETCONJUGATION | Chem.SANITIZE_SETHYBRIDIZATION |
               Chem.SANITIZE_SYMMRINGS | Chem.SANITIZE_FINDRADICALS)
        Chem.SanitizeMol(mm, sanitizeOps=ops, catchErrors=True)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mm, strict=True)
    rb_loose = rdMolDescriptors.CalcNumRotatableBonds(mm, strict=False)
    return heavy, rb, rb_loose, ok


def read_reference(name):
    """Heavy atoms, rotatable bonds and mode-1 scores of a reference compound."""
    m = next(iter(Chem.SDMolSupplier(str(REF / f'{REFERENCES[name]}_ligand.sdf'), sanitize=False, removeHs=False)))
    heavy, rb, rb_loose, ok = mol_descriptors(m)
    with open(REF / f'{REFERENCES[name]}_log.txt') as f:
        for line in f:
            mt = re.match(r'\s*1\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d\.\d+)\s+(\d+\.\d+)', line)
            if mt:
                score, aff = float(mt.group(3)), float(mt.group(4))
                break
    return dict(name=name, heavy=heavy, rb=rb, rb_loose=rb_loose, score=score, aff=aff)


def load():
    df = pd.read_csv(VS / 'docking_scores_all_2342.csv')
    sup = Chem.SDMolSupplier(str(VS / 'ligands_2342.sdf'), sanitize=False, removeHs=False)
    rows = [mol_descriptors(m) for m in sup]
    assert len(rows) == len(df), 'SDF and score table differ in length'
    df['heavy'], df['rb'], df['rb_loose'], df['sanitized'] = zip(*rows)
    # odanacatib: neutral form (see the folder excluded_compound_odanacatib)
    ex = pd.read_csv(VS / 'excluded_compound_odanacatib' / 'odanacatib_neutral_redocking_gnina_excerpt.tsv', sep='\t')
    first = ex.iloc[0]
    i = df.index[df.molecule_id == ODANACATIB][0]
    say(f'odanacatib (2006) as in the CSV: CNN score {df.loc[i, "cnn_pose_score"]}, CNN affinity {df.loc[i, "cnn_affinity"]}, '
        f'CNN_VS {df.loc[i, "cnn_vs"]} (protonated form)')
    df.loc[i, ['cnn_pose_score', 'cnn_affinity']] = [first.cnn_pose_score, first.cnn_affinity]
    df['cnn_vs'] = df.cnn_pose_score * df.cnn_affinity
    say(f'odanacatib (2006) replaced by the neutral form: CNN score {first.cnn_pose_score}, CNN affinity {first.cnn_affinity}, '
        f'CNN_VS {df.loc[i, "cnn_vs"]:.2f}')
    df['hit'] = df.cnn_vs >= THRESHOLD - 1e-9
    df['advanced'] = df.molecule_id.isin(ADVANCED)
    return df


# ---------------------------------------------------------------- statistics
def rho_ci(x, y, n=2000):
    r, p = spearmanr(x, y)
    idx = np.arange(len(x))
    bs = []
    for _ in range(n):
        s = rng.choice(idx, len(idx))
        bs.append(spearmanr(x[s], y[s])[0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return r, p, lo, hi


def fmt_p(p):
    return '<1e-300' if p == 0 else f'{p:.1e}'


def top_pct(values, x, higher_is_better=True):
    v = np.asarray(values)
    return 100.0 * ((v >= x).sum() if higher_is_better else (v <= x).sum()) / len(v)


def describe(s):
    return f'median {s.median():g} (IQR {s.quantile(.25):g}-{s.quantile(.75):g}; range {s.min():g}-{s.max():g})'


def main():
    OUT_T.mkdir(parents=True, exist_ok=True)
    OUT_F.mkdir(parents=True, exist_ok=True)
    df = load()
    refs = [read_reference(n) for n in REFERENCES]
    n = len(df)
    hits = df[df.hit]
    adv = df[df.advanced]

    say('\n== 1. Library and hit list ==')
    say(f'records: {n}; hits (CNN_VS >= {THRESHOLD}): {len(hits)}; records failing default sanitization: {(~df.sanitized).sum()}')
    say(f'CNN_VS median {df.cnn_vs.median():.2f}, 95th percentile {np.percentile(df.cnn_vs, 95):.2f}; '
        f'threshold at the {100 - 100 * len(hits) / n:.1f}th percentile; CNNaffinity >= {THRESHOLD}: {(df.cnn_affinity >= THRESHOLD).sum()} '
        f'records ({100 * (df.cnn_affinity >= THRESHOLD).mean():.1f}%)')
    prim = pd.read_csv(VS / 'primary_hits_30.csv')
    say('hit list identical to primary_hits_30.csv:', set(prim.molecule_id) == set(hits.molecule_id))

    say('\n== 2. Reference compounds of Reilly et al. (mode 1 of the redocking) ==')
    for r in refs:
        say(f"{r['name']:12s} heavy atoms {r['heavy']}, rotatable bonds {r['rb']} (strict) / {r['rb_loose']} (non-strict); "
            f"CNN score {r['score']:.2f}, CNN affinity {r['aff']:.2f}, CNN_VS {r['score'] * r['aff']:.2f}; "
            f"affinity top {top_pct(df.cnn_affinity, r['aff']):.1f}%, pose confidence bottom {top_pct(df.cnn_pose_score, r['score'], False):.1f}%")
        r['ha_pct'] = 100 * (df.heavy < r['heavy']).mean()
        say(f"{'':12s} larger than {r['ha_pct']:.0f}% of the library and than {100 * (hits.heavy < r['heavy']).mean():.0f}% of the hits (heavy atoms); "
            f"rotatable bonds: more than {100 * (df.rb < r['rb']).mean():.0f}% of the library, more than {100 * (hits.rb < r['rb']).mean():.0f}% of the hits")
    say(f"lowest CNN score among the hits: {hits.cnn_pose_score.min():.2f} (median {hits.cnn_pose_score.median():.2f})")

    say('\n== 3. Size and flexibility of the library and of the hits ==')
    for col, lab in [('heavy', 'heavy atoms'), ('rb', 'rotatable bonds (strict)'), ('rb_loose', 'rotatable bonds (non-strict)')]:
        say(f'{lab}: library {describe(df[col])}')
        say(f'{" " * len(lab)}  30 hits {describe(hits[col])}; hit median at the {100 * (df[col] < hits[col].median()).mean():.0f}th '
            f'percentile of the library; {100 * (hits[col] > df[col].median()).mean():.0f}% of the hits above the library median')
        say(f'{" " * len(lab)}  library without the {(~df.sanitized).sum()} unsanitized records: {describe(df[df.sanitized][col])}')
        say(f'{" " * len(lab)}  four prioritized compounds: ' + ', '.join(f"{NAMES[r.molecule_id]} {r[col]}" for _, r in adv.iterrows()))

    say('\n== 4. Spearman correlations over the 2,342 records ==')
    for col, lab in [('heavy', 'heavy atoms'), ('rb', 'rotatable bonds (strict)'), ('rb_loose', 'rotatable bonds (non-strict)')]:
        for sc, sl in [('cnn_pose_score', 'CNN score'), ('cnn_affinity', 'CNN affinity'), ('cnn_vs', 'CNN_VS')]:
            r, p, lo, hi = rho_ci(df[sc].to_numpy(), df[col].to_numpy())
            say(f'{sl:13s} vs {lab:30s} rho = {r:+.2f} (95% CI {lo:+.2f} to {hi:+.2f}), p = {fmt_p(p)}')
    d2 = df[df.sanitized]
    for col in ['heavy', 'rb']:
        for sc in ['cnn_pose_score', 'cnn_affinity']:
            r, p = spearmanr(d2[sc], d2[col])
            say(f'  without unsanitized records: {sc} vs {col}: rho = {r:+.2f}')
    r, p = spearmanr(df.heavy, df.rb)
    say(f'heavy atoms vs rotatable bonds (strict): rho = {r:+.2f}')

    say('\n== 5. Fraction of records passing CNN_VS >= 6.3, by size class (quartiles of the library) ==')
    for col, lab in [('heavy', 'heavy atoms'), ('rb', 'rotatable bonds (strict)')]:
        q = pd.qcut(df[col], 4, duplicates='drop')
        g = df.assign(q=q).groupby('q', observed=True)
        t = g.agg(n=('hit', 'size'), hits=('hit', 'sum'), pass_pct=('hit', lambda s: 100 * s.mean()),
                  median_score=('cnn_pose_score', 'median'), median_affinity=('cnn_affinity', 'median'))
        say(f'{lab}:')
        say(t.round(2).to_string())

    say('\n== 6. Ranking of the library by alternative criteria ==')
    df['rank_vs'] = df.cnn_vs.rank(ascending=False, method='min')
    df['rank_aff'] = df.cnn_affinity.rank(ascending=False, method='min')
    df['rank_score'] = df.cnn_pose_score.rank(ascending=False, method='min')
    df['aff_per_heavy'] = df.cnn_affinity / df.heavy
    df['rank_aff_per_heavy'] = df.aff_per_heavy.rank(ascending=False, method='min')
    slope, icpt = np.polyfit(np.log(df.heavy), df.cnn_affinity, 1)
    df['aff_resid'] = df.cnn_affinity - (slope * np.log(df.heavy) + icpt)
    df['rank_aff_resid'] = df.aff_resid.rank(ascending=False, method='min')
    hits, adv = df[df.hit], df[df.advanced]          # refresh the subsets with the rank columns
    cols = [('rank_vs', 'CNN_VS (criterion used)'), ('rank_aff', 'CNN affinity alone'), ('rank_score', 'CNN score alone'),
            ('rank_aff_per_heavy', 'CNN affinity per heavy atom'), ('rank_aff_resid', 'CNN affinity residual from ln(heavy atoms)')]
    say('position (top %) of the four prioritized compounds:')
    for rc, lab in cols:
        say(f'  {lab:45s} ' + '; '.join(f"{NAMES[r.molecule_id]} {int(r[rc])} ({100 * r[rc] / n:.1f}%)" for _, r in adv.iterrows()))
    say('overlap of the 30 primary hits with the best-ranked records under each alternative criterion:')
    hit_ids = set(hits.molecule_id)
    for rc, lab in cols[1:]:
        s = [f'top {k}: {len(hit_ids & set(df.loc[df[rc] <= k, "molecule_id"]))}' for k in (30, 100, 617)]
        say(f'  {lab:45s} ' + '; '.join(s) + f'; ties at rank 30: {(df[rc] <= 30).sum()} records')
    top30 = df[df.rank_aff <= 30]
    say(f'30 best by CNN affinity: heavy atoms {describe(top30.heavy)}, rotatable bonds {describe(top30.rb)}, CNN score {describe(top30.cnn_pose_score)}')
    say('the four prioritized compounds are in the top 617 by CNN affinity by construction (CNN_VS >= 6.3 implies CNN affinity >= 6.3).')
    say('where the reference compounds would rank by CNN affinity among the 2,342 records: ' +
        '; '.join(f"{r['name']} {1 + int((df.cnn_affinity > r['aff']).sum())}" for r in refs))

    # ---- tables
    df.to_csv(OUT_T / 'library_descriptors_2342.csv', index=False,
              columns=['molecule_id', 'cnn_pose_score', 'cnn_affinity', 'cnn_vs', 'hit', 'heavy', 'rb', 'rb_loose', 'sanitized',
                       'rank_vs', 'rank_aff', 'rank_score', 'rank_aff_per_heavy', 'rank_aff_resid'], float_format='%.4g')
    t1 = hits[['molecule_id', 'heavy', 'rb', 'rb_loose']].rename(columns={
        'molecule_id': 'Molecule ID', 'heavy': 'Ligand heavy atoms', 'rb': 'Rotatable bonds', 'rb_loose': 'Rotatable bonds (non-strict)'})
    t1.to_csv(OUT_T / 'TableS1_descriptor_columns.tsv', sep='\t', index=False)

    # ---- figure
    figure(df, hits, adv, refs)
    (OUT_T / 'library_descriptors_and_ranking.txt').write_text('\n'.join(OUT) + '\n', encoding='utf-8')


def figure(df, hits, adv, refs):
    BLUE, ORANGE, GRAY, INK = '#5b7a9d', '#c98a5a', '#aeb6bf', '#222222'
    DARK_BLUE, DARK_ORANGE = '#2f5d8a', '#c4661f'
    plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(2, 3, figsize=(10, 6.7))
    n = len(df)
    npass = int(df.hit.sum())
    # (a) CNN_VS
    a = ax[0, 0]
    a.hist(df.cnn_vs, bins=60, color=BLUE, edgecolor='none')
    a.axvline(THRESHOLD, color='k', ls='--', lw=1.1)
    a.annotate(f'CNN_VS = {THRESHOLD:g}\n{npass}/{n} ({100 * npass / n:.1f}%)', xy=(THRESHOLD, a.get_ylim()[1] * 0.75),
               xytext=(5, 0), textcoords='offset points', fontsize=7)
    a.set_xlabel('CNN_VS'); a.set_ylabel('compounds')
    # (b) CNN score
    b = ax[0, 1]
    b.hist(df.cnn_pose_score, bins=60, color=ORANGE, edgecolor='none')
    ymax = b.get_ylim()[1]
    for r in refs:
        b.axvline(r['score'], color='k', lw=0.9)
        b.annotate(r['name'], xy=(r['score'], ymax * 0.85), rotation=90, fontsize=6.5, ha='right', va='top')
    b.set_xlabel('CNN pose-confidence score'); b.set_ylabel('compounds')

    def scatter(axis, ycol, ylabel, label_refs):
        axis.scatter(df.heavy, df[ycol], s=5, c=GRAY, alpha=0.45, linewidths=0, rasterized=True, label='library (2,342)')
        axis.scatter(hits.heavy, hits[ycol], s=16, c=DARK_BLUE, linewidths=0, label='primary hits (30)', zorder=3)
        axis.scatter(adv.heavy, adv[ycol], s=60, facecolors='none', edgecolors=DARK_ORANGE, linewidths=1.3,
                     label='four prioritized', zorder=4)
        key = 'score' if ycol == 'cnn_pose_score' else 'aff'
        axis.scatter([r['heavy'] for r in refs], [r[key] for r in refs], s=34, marker='D', c='k', zorder=5,
                     label='reference compounds')
        if label_refs:
            for r in refs:
                axis.annotate(r['name'], (r['heavy'], r[key]), xytext=(5, -2), textcoords='offset points', fontsize=6.5)
        axis.set_xlabel('ligand heavy atoms'); axis.set_ylabel(ylabel)
        rho = spearmanr(df.heavy, df[ycol])[0]
        axis.text(0.97, 0.97 if ycol == 'cnn_affinity' else 0.03, f'Spearman ρ = {rho:+.2f}', transform=axis.transAxes,
                  ha='right', va='top' if ycol == 'cnn_affinity' else 'bottom', fontsize=7.5)

    scatter(ax[0, 2], 'cnn_pose_score', 'CNN pose-confidence score', True)
    scatter(ax[1, 0], 'cnn_affinity', 'CNN affinity (pK)', False)

    def ecdf(axis, col, xlabel):
        for d, c, lab in [(df, GRAY, 'library (2,342)'), (hits, DARK_BLUE, 'primary hits (30)')]:
            x = np.sort(d[col].to_numpy())
            axis.step(x, np.arange(1, len(x) + 1) / len(x), where='post', color=c, lw=1.6, label=lab)
        axis.plot([r[{'heavy': 'heavy', 'rb': 'rb'}[col]] for r in refs], [0.02] * len(refs), 'kD', ms=4.5,
                  clip_on=False, label='reference compounds', zorder=5)
        axis.set_xlabel(xlabel); axis.set_ylabel('cumulative fraction')
        axis.set_ylim(0, 1.02)

    ecdf(ax[1, 1], 'heavy', 'ligand heavy atoms')
    ecdf(ax[1, 2], 'rb', 'rotatable bonds (RDKit, strict)')
    for axis, lab in zip(ax.ravel(), 'abcdef'):
        axis.set_title(f'({lab})', loc='left', fontsize=10)
        axis.tick_params(width=0.8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='o', ls='', ms=4, mfc=GRAY, mec=GRAY, label='library (2,342)'),
               Line2D([], [], marker='o', ls='', ms=5, mfc=DARK_BLUE, mec=DARK_BLUE, label='primary hits (30)'),
               Line2D([], [], marker='o', ls='', ms=8, mfc='none', mec=DARK_ORANGE, mew=1.3, label='four prioritized (Table 2)'),
               Line2D([], [], marker='D', ls='', ms=5, mfc='k', mec='k', label='reference compounds (Reilly et al.)')]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=8, handletextpad=0.3, columnspacing=2.0)
    fig.tight_layout(w_pad=1.6, h_pad=1.8, rect=[0, 0.045, 1, 1])
    fig.savefig(OUT_F / 'FigureS1_CNN_score_distribution.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
