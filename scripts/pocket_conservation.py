#!/usr/bin/env python3
"""Conservation of the 16 hinge-pocket residues across SARS-CoV-2 lineages and across betacoronaviruses.

Run from the repository root:
    python scripts/pocket_conservation.py --fetch      # query LAPIS (internet) and rewrite the cached table, then analyze
    python scripts/pocket_conservation.py              # analyze with the cached table in data/05_conservation/

The pocket residues are those of Section 3.4 (model numbering; Spike number = model number + 815):
Phe83, Ala84, Met87, Leu101, Gln105, Leu107, Ile108, Ala109, Gln111, Phe112, Met235, Val250, Thr251, Tyr252, Ala255, Phe294
= Phe898, Ala899, Met902, Leu916, Gln920, Leu922, Ile923, Ala924, Gln926, Phe927, Met1050, Val1065, Thr1066, Tyr1067,
Ala1070, Phe1109 of the Wuhan-Hu-1 spike (UniProt P0DTC2).

(1) SARS-CoV-2 lineages. For each lineage (with its sublineages) the open LAPIS instance of CoV-Spectrum
    (https://lapis.cov-spectrum.org/open/v2, data: Nextstrain open dataset of public GenBank sequences) returns the amino-acid
    substitutions of the spike observed with a proportion >= 0.005% together with the number of sequences that cover the position.
    Positions are in the Wuhan-Hu-1 numbering. For every pocket position the table gives the proportion of sequences with any
    substitution, and the most frequent substitution. Lineages are nested (JN.1* includes KP.3*, XEC* and LP.8.1*, BA.2.86* includes
    JN.1*), and 'all sequences' is the whole open dataset. A substitution seen in few sequences may be a sequencing error.
    The query is cached in data/05_conservation/lapis_pocket_substitutions.tsv together with the date and the data version.
(2) Betacoronaviruses. The spike sequences in data/05_conservation/spike_sequences_uniprot.fasta (UniProt accessions in the headers:
    SARS-CoV-2 P0DTC2; bat sarbecoviruses RaTG13 A0ABF7PLN6, WIV1 U5WI05, RsSHC014 U5WLK5; SARS-CoV P59594; MERS-CoV K9N5Q8;
    HCoV-HKU1 isolate N5 Q0ZME7; HCoV-OC43 P36334) are aligned with FAMSA (pyfamsa, UPGMA guide tree). For every pocket residue the
    residue in each virus is classed as identical, similar (BLOSUM62 score > 0) or different/gap. For context the identity over the
    S2 region of the model (Spike 816-1154) and over the 323 residues of that region outside the pocket are also computed.

Writes (results/tables, results/figures):
  pocket_conservation_lineages.tsv, pocket_conservation_betacoronaviruses.tsv, pocket_conservation_identity.tsv (wide tables),
  pocket_conservation_A_by_position.tsv, pocket_conservation_B_alignment.tsv, pocket_conservation_C_identity.tsv (narrow tables of the SI),
  pocket_conservation_summary.txt, FigureS_pocket_conservation.png
Requires numpy, pandas, requests, matplotlib, biopython and pyfamsa.
"""
import argparse
import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Align import substitution_matrices
from pyfamsa import Aligner, Sequence

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'data' / '05_conservation'
OUT_T = REPO / 'results' / 'tables'
OUT_F = REPO / 'results' / 'figures'
MODEL = [83, 84, 87, 101, 105, 107, 108, 109, 111, 112, 235, 250, 251, 252, 255, 294]
SPIKE = [m + 815 for m in MODEL]
EXPECTED = 'FAMLQLIAQFMVTYAF'            # one-letter code of the 16 residues, in the order above
LAPIS = 'https://lapis.cov-spectrum.org/open/v2/sample/'
MIN_PROPORTION = 0.00005
LINEAGES = [('Alpha', 'B.1.1.7*'), ('Beta', 'B.1.351*'), ('Gamma', 'P.1*'), ('Delta', 'B.1.617.2*'), ('Omicron BA.1', 'BA.1*'),
            ('Omicron BA.2', 'BA.2*'), ('Omicron BA.5', 'BA.5*'), ('XBB.1.5', 'XBB.1.5*'), ('BA.2.86', 'BA.2.86*'), ('JN.1', 'JN.1*'),
            ('KP.3', 'KP.3*'), ('XEC', 'XEC*'), ('LP.8.1', 'LP.8.1*'), ('NB.1.8.1', 'NB.1.8.1*'), ('XFG', 'XFG*'), ('all sequences', None)]
VIRUSES = [('P0DTC2', 'SARS-CoV-2 (Wuhan-Hu-1)'), ('A0ABF7PLN6', 'Bat RaTG13'), ('U5WI05', 'Bat WIV1'), ('U5WLK5', 'Bat RsSHC014'),
           ('P59594', 'SARS-CoV'), ('K9N5Q8', 'MERS-CoV'), ('Q0ZME7', 'HCoV-HKU1 (N5)'), ('P36334', 'HCoV-OC43')]
THREE = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
         'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'}
OUT = []


def say(*a):
    line = ' '.join(str(x) for x in a)
    print(line)
    OUT.append(line)


def lapis(endpoint, **params):
    r = requests.get(LAPIS + endpoint, params=params, timeout=300)
    r.raise_for_status()
    return r.json()


def fetch():
    rows = []
    for label, lin in LINEAGES:
        flt = {'nextcladePangoLineage': lin} if lin else {}
        n = lapis('aggregated', **flt)['data'][0]['count']
        res = lapis('aminoAcidMutations', minProportion=MIN_PROPORTION, limit=1000000, **flt)
        info = res['info']
        for m in res['data']:
            if m['sequenceName'] == 'S' and m['position'] in SPIKE:
                rows.append(dict(lineage=label, query=lin or '(none)', n_sequences=n, position=m['position'], mutation=m['mutation'],
                                 count=m['count'], coverage=m['coverage'], proportion=m['proportion']))
        rows.append(dict(lineage=label, query=lin or '(none)', n_sequences=n, position=0, mutation='-', count=0, coverage=0, proportion=0.0))
        print(label, n)
    df = pd.DataFrame(rows)
    df.attrs['info'] = info
    df.to_csv(DATA / 'lapis_pocket_substitutions.tsv', sep='\t', index=False)
    (DATA / 'lapis_query_info.txt').write_text(
        f"queried on {datetime.date.today()}; LAPIS {info.get('lapisVersion')}; dataVersion {info.get('dataVersion')}; "
        f"{info.get('requestInfo')}; minProportion {MIN_PROPORTION}\n", encoding='utf-8')


def lineage_table():
    df = pd.read_csv(DATA / 'lapis_pocket_substitutions.tsv', sep='\t')
    rows = []
    for label, lin in LINEAGES:
        d = df[df.lineage == label]
        n = int(d.n_sequences.iloc[0])
        row = {'Lineage': label, 'Sequences': n}
        for sp, mp in zip(SPIKE, MODEL):
            x = d[d.position == sp]
            ref = EXPECTED[SPIKE.index(sp)]
            if x.empty:
                row[f'{THREE[ref]}{sp}'] = f'0 (<{100 * MIN_PROPORTION:g}%)'
            else:
                tot = 100 * x.proportion.sum()
                top = x.sort_values('proportion').iloc[-1]
                row[f'{THREE[ref]}{sp}'] = f'{tot:.3f}% ({top.mutation.split(":")[1]})'
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT_T / 'pocket_conservation_lineages.tsv', sep='\t', index=False)
    sub = df[df.position > 0].copy()
    sub['pct'] = 100 * sub.proportion
    return df, sub


def align():
    recs = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(DATA / 'spike_sequences_uniprot.fasta', 'fasta')}
    seqs = [Sequence(acc.encode(), recs[acc].encode()) for acc, _ in VIRUSES]
    msa = Aligner(guide_tree='upgma').align(seqs)
    aligned = {s.id.decode(): s.sequence.decode() for s in msa}
    ref = aligned['P0DTC2']
    # map Wuhan residue number -> alignment column
    col_of, k = {}, 0
    for c, ch in enumerate(ref):
        if ch != '-':
            k += 1
            col_of[k] = c
    assert ''.join(ref[col_of[p]] for p in SPIKE) == EXPECTED, 'pocket residues do not match the Wuhan-Hu-1 sequence'
    return aligned, col_of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fetch', action='store_true')
    a = ap.parse_args()
    OUT_T.mkdir(parents=True, exist_ok=True)
    OUT_F.mkdir(parents=True, exist_ok=True)
    if a.fetch:
        fetch()
    say((DATA / 'lapis_query_info.txt').read_text(encoding='utf-8').strip())

    say('\n== (1) SARS-CoV-2 lineages: substitutions at the 16 pocket positions ==')
    raw, sub = lineage_table()
    for label, _ in LINEAGES:
        d = sub[sub.lineage == label]
        n = int(raw[raw.lineage == label].n_sequences.iloc[0])
        say(f'{label:15s} {n:>9,d} sequences; ' + ('none above 0.005% at any of the 16 positions' if d.empty else
            '; '.join(f'{m.mutation.split(":")[1]} {m.pct:.3f}% ({int(m["count"])}/{int(m.coverage)})' for _, m in d.sort_values("position").iterrows())))
    allseq = sub[sub.lineage == 'all sequences']
    say(f'positions with a substitution in more than 0.1% of the sequences of any lineage listed: ' +
        (', '.join(sorted(set(sub[sub.pct > 0.1].mutation))) if (sub.pct > 0.1).any() else 'none'))
    say(f'largest proportion over all lineages (excluding "all sequences"): {sub[sub.lineage != "all sequences"].pct.max():.2f}% '
        f'({sub.loc[sub[sub.lineage != "all sequences"].pct.idxmax(), "mutation"]} in {sub.loc[sub[sub.lineage != "all sequences"].pct.idxmax(), "lineage"]})')
    say(f'positions with any substitution in all sequences (>= 0.005%): {len(allseq)} substitution types at positions '
        f'{sorted(set(allseq.position))}; invariant at the 0.005% level: ' +
        ', '.join(f'{THREE[EXPECTED[i]]}{p}' for i, p in enumerate(SPIKE) if p not in set(allseq.position)))

    say('\n== (2) Betacoronaviruses: FAMSA alignment of the full-length spikes ==')
    aligned, col_of = align()
    B62 = substitution_matrices.load('BLOSUM62')
    rows, ident = [], []
    s2 = range(816, 1155)
    pocket_set = set(SPIKE)
    outside = [p for p in s2 if p not in pocket_set]
    for acc, name in VIRUSES:
        seq = aligned[acc]
        row = {'Virus': name, 'UniProt': acc}
        n_id = n_sim = n_dif = 0
        for i, p in enumerate(SPIKE):
            ref = EXPECTED[i]
            r = seq[col_of[p]]
            cls = 'identical' if r == ref else ('different (gap)' if r == '-' else ('similar' if B62[ref][r] > 0 else 'different'))
            if cls == 'identical': n_id += 1
            elif cls == 'similar': n_sim += 1
            else: n_dif += 1
            row[f'{THREE[ref]}{p}'] = r
        row['Identical'], row['Similar'], row['Different'] = n_id, n_sim, n_dif
        rows.append(row)
        ident.append(dict(Virus=name, UniProt=acc,
                          pocket_identity_pct=round(100 * n_id / 16, 1),
                          S2_816_1154_identity_pct=round(100 * sum(seq[col_of[p]] == aligned['P0DTC2'][col_of[p]] for p in s2) / len(s2), 1),
                          S2_outside_pocket_identity_pct=round(100 * sum(seq[col_of[p]] == aligned['P0DTC2'][col_of[p]] for p in outside) / len(outside), 1),
                          full_spike_identity_pct=round(100 * sum(1 for c, ch in enumerate(aligned['P0DTC2']) if ch != '-' and seq[c] == ch) /
                                                        sum(1 for ch in aligned['P0DTC2'] if ch != '-'), 1)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT_T / 'pocket_conservation_betacoronaviruses.tsv', sep='\t', index=False)
    ti = pd.DataFrame(ident)
    ti.to_csv(OUT_T / 'pocket_conservation_identity.tsv', sep='\t', index=False)
    for r in rows:
        say(f"{r['Virus']:26s} identical {r['Identical']:2d}/16, similar {r['Similar']:2d}, different {r['Different']:2d}")
    say(ti.to_string(index=False))
    # positions identical in every virus
    allid = [f'{THREE[EXPECTED[i]]}{p}' for i, p in enumerate(SPIKE) if all(aligned[acc][col_of[p]] == EXPECTED[i] for acc, _ in VIRUSES)]
    say(f'identical in all 8 sequences: {len(allid)} positions: ' + ', '.join(allid))
    sarbeco = [a for a, _ in VIRUSES[:5]]
    ids_sarb = [f'{THREE[EXPECTED[i]]}{p}' for i, p in enumerate(SPIKE) if all(aligned[acc][col_of[p]] == EXPECTED[i] for acc in sarbeco)]
    say(f'identical in the 5 sarbecoviruses (SARS-CoV-2, three bat, SARS-CoV): {len(ids_sarb)} positions')
    nonsarb = [a for a, _ in VIRUSES[5:]]
    for acc, name in VIRUSES[5:]:
        diff = [f'{THREE[EXPECTED[i]]}{p}->{aligned[acc][col_of[p]]}' for i, p in enumerate(SPIKE) if aligned[acc][col_of[p]] != EXPECTED[i]]
        say(f'{name}: differences at ' + ', '.join(diff))
    narrow_tables(sub, aligned, col_of, B62, ti, t)
    figure(raw, aligned, col_of, B62)
    (OUT_T / 'pocket_conservation_summary.txt').write_text('\n'.join(OUT) + '\n', encoding='utf-8')


def narrow_tables(sub, aligned, col_of, B62, ti, t):
    """Three narrow tables for the Supporting Information (positions as rows)."""
    rowsA = []
    lin_names = [l for l, _ in LINEAGES if l != 'all sequences']
    for i, p in enumerate(SPIKE):
        ref = EXPECTED[i]
        allx = sub[(sub.lineage == 'all sequences') & (sub.position == p)]
        a_all = '<0.005%' if allx.empty else f'{100 * allx.proportion.sum():.3f}% ({allx.sort_values("proportion").iloc[-1].mutation.split(":")[1]})'
        per = {}
        for l in lin_names:
            x = sub[(sub.lineage == l) & (sub.position == p)]
            if not x.empty:
                per[l] = (100 * x.proportion.sum(), x.sort_values('proportion').iloc[-1].mutation.split(':')[1])
        if per:
            l, (v, m) = max(per.items(), key=lambda kv: kv[1][0])
            hi = f'{v:.3f}% ({l}, {m})'
        else:
            hi = '<0.005%'
        n01 = sum(1 for v, _ in per.values() if v >= 0.1)
        rowsA.append({'Pocket residue (Spike)': f'{THREE[ref]}{p}', 'Model residue': f'{THREE[ref]}{MODEL[i]}',
                      'All sequences: % with a substitution (most frequent)': a_all,
                      'Highest frequency in a lineage (lineage, substitution)': hi,
                      f'Lineages with a substitution in >= 0.1% of the sequences (of {len(lin_names)})': n01})
    pd.DataFrame(rowsA).to_csv(OUT_T / 'pocket_conservation_A_by_position.tsv', sep='	', index=False)
    rowsB = []
    for i, p in enumerate(SPIKE):
        r = {'Pocket residue (Spike)': f'{THREE[EXPECTED[i]]}{p}'}
        for acc, name in VIRUSES:
            r[name] = aligned[acc][col_of[p]]
        rowsB.append(r)
    pd.DataFrame(rowsB).to_csv(OUT_T / 'pocket_conservation_B_alignment.tsv', sep='	', index=False)
    c = t[['Virus', 'UniProt', 'Identical', 'Similar', 'Different']].merge(ti, on=['Virus', 'UniProt'])
    c = c.rename(columns={'pocket_identity_pct': 'Identity, 16 pocket residues (%)', 'S2_816_1154_identity_pct': 'Identity, S2 region 816-1154 (%)',
                          'S2_outside_pocket_identity_pct': 'Identity, S2 region outside the pocket (%)'}).drop(columns=['full_spike_identity_pct'])
    c.to_csv(OUT_T / 'pocket_conservation_C_identity.tsv', sep='	', index=False)


def figure(raw, aligned, col_of, B62):
    plt.rcParams.update({'font.size': 8})
    labels = [f'{THREE[EXPECTED[i]]}{p}' for i, p in enumerate(SPIKE)]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7.6), gridspec_kw={'height_ratios': [len(LINEAGES), len(VIRUSES)]})
    # (a) lineages: proportion of sequences with a substitution (log-ish scale as text)
    P = np.zeros((len(LINEAGES), 16))
    for i, (lab, _) in enumerate(LINEAGES):
        d = raw[(raw.lineage == lab) & (raw.position > 0)]
        for j, p in enumerate(SPIKE):
            P[i, j] = 100 * d[d.position == p].proportion.sum()
    im = a1.imshow(np.minimum(P, 5), cmap='Oranges', vmin=0, vmax=5, aspect='auto')
    for i in range(P.shape[0]):
        for j in range(16):
            a1.text(j, i, EXPECTED[j] if P[i, j] < 0.005 else f'{EXPECTED[j]}\n{P[i, j]:.2f}', ha='center', va='center', fontsize=6.5 if P[i, j] >= 0.005 else 8)
    a1.set_yticks(range(len(LINEAGES)))
    a1.set_yticklabels([f"{l} (n = {int(raw[raw.lineage == l].n_sequences.iloc[0]):,})" for l, _ in LINEAGES])
    a1.set_xticks(range(16)); a1.set_xticklabels(labels, rotation=45, ha='right'); a1.xaxis.tick_top()
    a1.set_title('(a) SARS-CoV-2: % of sequences with a substitution (shading saturates at 5%; value shown when ≥ 0.005%)', loc='left', fontsize=9, pad=34)
    # (b) betacoronaviruses
    colors = {'identical': '#2f5d8a', 'similar': '#a9c4de', 'different': '#c4661f'}
    for i, (acc, name) in enumerate(VIRUSES):
        for j, p in enumerate(SPIKE):
            ref, r = EXPECTED[j], aligned[acc][col_of[p]]
            cls = 'identical' if r == ref else ('similar' if (r != '-' and B62[ref][r] > 0) else 'different')
            a2.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color=colors[cls], ec='white', lw=1))
            a2.text(j, i, r, ha='center', va='center', fontsize=8.5, color='white' if cls != 'similar' else '#222222')
    a2.set_xlim(-.5, 15.5); a2.set_ylim(len(VIRUSES) - .5, -.5)
    a2.set_yticks(range(len(VIRUSES))); a2.set_yticklabels([n for _, n in VIRUSES])
    a2.set_xticks(range(16)); a2.set_xticklabels(labels, rotation=45, ha='right')
    a2.set_title('(b) Betacoronaviruses: residue aligned to each pocket position (FAMSA alignment of the full-length spikes)', loc='left', fontsize=9)
    from matplotlib.patches import Patch
    a2.legend(handles=[Patch(color=colors['identical'], label='identical'), Patch(color=colors['similar'], label='similar (BLOSUM62 > 0)'),
                       Patch(color=colors['different'], label='different or gap')], loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=3, frameon=False)
    for ax in (a1, a2):
        for s in ax.spines.values(): s.set_visible(False)
        ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(OUT_F / 'FigureS_pocket_conservation.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
