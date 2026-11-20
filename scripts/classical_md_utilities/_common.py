"""
Shared helpers for the scripts in utilities/.

Every analysis script is paper-agnostic: it takes --paper <slug>, resolves the
folder papers/<slug>/, and reads all dataset-specific values (system names,
equilibration windows, pocket residues, replicates, paths) from
papers/<slug>/config/analysis_params.json. Nothing about a specific paper
belongs in the scripts themselves.

Schema of analysis_params.json: see .claude/skills/utilities-scripts/SKILL.md
(template: papers/_template/config/analysis_params.json).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

UTIL_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(UTIL_DIR)
PAPERS_DIR = os.path.join(WORKSPACE, 'papers')

_FALLBACK_COLORS = ['#1f78b4', '#d95f02', '#6a3d9a', '#33a02c', '#e31a1c',
                    '#b15928', '#a6761d', '#666666']


class Paper:
    """A resolved papers/<slug>/ folder plus its analysis configuration."""

    def __init__(self, slug, out=None):
        self.slug = slug
        self.dir = os.path.join(PAPERS_DIR, slug)
        if not os.path.isdir(self.dir):
            avail = [d for d in sorted(os.listdir(PAPERS_DIR))
                     if os.path.isdir(os.path.join(PAPERS_DIR, d)) and not d.startswith('_')]
            sys.exit(f"Paper '{slug}' not found in {PAPERS_DIR}. Available: {avail}")
        cfg_path = os.path.join(self.dir, 'config', 'analysis_params.json')
        if not os.path.isfile(cfg_path):
            sys.exit(f'Missing {cfg_path} (copy it from papers/_template/config/).')
        with open(cfg_path, encoding='utf-8') as fh:
            self.cfg = json.load(fh)
        self.out = os.path.abspath(out) if out else os.path.join(self.dir, 'bozze', 'figures')
        os.makedirs(self.out, exist_ok=True)

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def section(self, name):
        if name not in self.cfg:
            sys.exit(f"config/analysis_params.json of '{self.slug}' has no '{name}' section.")
        return self.cfg[name]

    # ---- systems -------------------------------------------------------
    def systems(self, ids):
        """Ordered list of system ids -> merged definition dicts."""
        defs = self.cfg.get('systems', {})
        out = []
        for i, sid in enumerate(ids):
            if sid not in defs:
                sys.exit(f"System '{sid}' is not defined in 'systems' of the config.")
            d = dict(defs[sid])
            d['id'] = sid
            d.setdefault('label', sid)
            d.setdefault('code', sid)
            d.setdefault('linestyle', '-')
            d.setdefault('has_ligand', True)
            d.setdefault('color', _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
            out.append(d)
        return out

    def pocket(self):
        pk = self.section('pocket')
        res = [(int(n), aa) for n, aa in pk['residues']]
        return res, int(pk.get('numbering_offset', 0)), pk


def parse_args(description):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--paper', required=True,
                    help='slug of the folder in papers/ (e.g. locking-the-spike)')
    ap.add_argument('--out', default=None,
                    help='output directory (default: papers/<slug>/bozze/figures)')
    return ap.parse_args()


def load_paper(description):
    args = parse_args(description)
    return Paper(args.paper, args.out)


def apply_style(big=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if big:
        plt.rcParams.update({'font.size': 12, 'axes.labelsize': 14,
                             'axes.titlesize': 16, 'axes.linewidth': 1.5})
    else:
        plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                             'axes.titlesize': 12, 'axes.linewidth': 1.0})
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']})
    return plt


def read_analysis_tab(path, drop_summary_rows=False):
    """YASARA md_analyze *_analysis.tab -> numeric DataFrame.
    drop_summary_rows removes trailing Mean/Min/Max rows (no valid Time)."""
    d = pd.read_csv(path, sep=r'\s+')
    d = d.apply(pd.to_numeric, errors='coerce')
    if drop_summary_rows:
        d = d[d['Time[ns]'].notna()]
    return d


def read_rmsf(path):
    """YASARA md_analyze *_rmsfres.tab -> DataFrame(res, num, chain, rmsf, num_i)."""
    rows = []
    with open(path, encoding='utf-8') as fh:
        for ln in fh.read().splitlines()[1:]:
            p = ln.split()
            if len(p) == 4:
                rows.append((p[0], p[1], p[2], float(p[3])))
    rf = pd.DataFrame(rows, columns=['res', 'num', 'chain', 'rmsf'])
    rf['num_i'] = pd.to_numeric(rf['num'], errors='coerce')
    return rf


def roll(y, n):
    return pd.Series(y).rolling(int(n), center=True, min_periods=1).mean().values


def bar_layout(n_items, total_width=0.8):
    """Bar width and per-item x offsets for grouped bars, any number of items."""
    w = total_width / n_items
    return w, [(i - (n_items - 1) / 2) * w for i in range(n_items)]
