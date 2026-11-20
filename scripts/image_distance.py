#!/usr/bin/env python3
"""Exact distance between the protein and its periodic images in an SMD run (R2-9).

Input (one run = one folder, all files of the same target):
  <target>_pocket_labxyz_<snapshot>.tab   written by md_analyze_local_S2_POCKET.mcr with pdbsaved=1:
                                          x y z of the protein heavy atoms in the SIMULATION-CELL frame
  <target>_pocket_analysis*.tab           same macro: CellLengthX/Y/Z and TrueGapX/Y/Z per analyzed snapshot
Needs numpy and scipy only.

For every snapshot:
  1. rebuilds the chain atom after atom (minimum-image correction) -> true extent,
  2. VALIDATES the coordinates against the macro: the true extent must equal cell - TrueGap of the same
     snapshot (both come from the same PosAtom calls); otherwise the run is flagged and must not be used,
  3. computes the shortest distance between a protein atom and any atom of the 26 neighbouring periodic
     images of the same protein (the quantity that matters: below the 8 A non-bonded cutoff the protein
     interacts with its own image) and the number of atoms closer than 8 A to an image atom.

Usage:
  python rev_image_dist.py FOLDER --tab FOLDER/<target>_pocket_analysis.tab [--out image_dist.csv]
Time column: snapshot number x 0.25 ns (250 ps saving interval).
"""
import argparse
import glob
import itertools
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

CUTOFF = 8.0
DT_NS = 0.25


def read_xyz(path):
    out = []
    for line in open(path):
        t = line.split()
        if len(t) < 3:
            continue
        try:
            out.append([float(t[0]), float(t[1]), float(t[2])])
        except ValueError:
            continue
    return np.array(out)


def read_tab(path):
    hdr, rows = None, []
    for line in open(path):
        t = line.split()
        if not t:
            continue
        if hdr is None and 'CellLengthX' in t:
            hdr = t
            continue
        if hdr is not None:
            try:
                float(t[0])
            except ValueError:
                continue
            if len(t) == len(hdr):
                rows.append([float(x) for x in t])
    d = np.array(rows)
    col = lambda n: d[:, hdr.index(n)]
    cell = np.column_stack([col('CellLengthX'), col('CellLengthY'), col('CellLengthZ')])
    gap = np.column_stack([col('TrueGapX'), col('TrueGapY'), col('TrueGapZ')])
    return cell, gap


def unwrap(xyz, L):
    out = xyz.copy()
    for i in range(1, len(out)):
        d = out[i] - out[i - 1]
        out[i] = out[i - 1] + d - L * np.round(d / L)
    return out


def image_min_distance(xyz, L):
    tree = cKDTree(xyz)
    best, axes, close = np.inf, '', np.zeros(len(xyz), bool)
    for s in itertools.product((-1, 0, 1), repeat=3):
        if s == (0, 0, 0):
            continue
        d, _ = tree.query(xyz + np.array(s) * L)
        close |= d < CUTOFF
        if d.min() < best:
            best, axes = d.min(), ''.join('xyz'[i] for i, v in enumerate(s) if v)
    return best, axes, int(close.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('folder')
    ap.add_argument('--tab', required=True, help='<target>_pocket_analysis.tab of the same run')
    ap.add_argument('--out')
    a = ap.parse_args()
    files = glob.glob(os.path.join(a.folder, '*_pocket_labxyz_*.tab'))
    num = lambda f: int(re.search(r'_labxyz_(\d+)\.tab$', f).group(1))
    files = sorted(files, key=num)
    if not files:
        sys.exit('no *_pocket_labxyz_<n>.tab files in %s (set pdbsaved=1 in the macro)' % a.folder)
    cell, gap = read_tab(a.tab)
    if len(cell) != len(files):
        print('WARNING: %d coordinate files but %d rows in %s: matched by order' % (len(files), len(cell), a.tab))
    rows, bad = [], 0
    for k, f in enumerate(files):
        xyz = read_xyz(f)
        L = cell[k]
        u = unwrap(xyz, L)
        ext = u.max(0) - u.min(0)
        ok = np.abs((L - ext) - gap[k]).max() < 0.5
        bad += (not ok)
        dmin, axes, nclose = image_min_distance(u, L)
        rows.append((num(f) * DT_NS, *ext, *(L - ext), dmin, axes, nclose, int(ok)))
    hdr = 'time_ns,extent_x,extent_y,extent_z,truegap_x,truegap_y,truegap_z,min_image_dist,axes,atoms_within_8A_of_image,matches_macro'
    lines = [hdr] + [','.join(str(round(v, 2)) if isinstance(v, float) else str(v) for v in r) for r in rows]
    if a.out:
        open(a.out, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines[:6] + (['...'] if len(lines) > 6 else [])))
    d = np.array([r[7] for r in rows])
    print('\nVALIDATION (coordinates vs macro TrueGap): %d/%d snapshots consistent -> %s'
          % (len(rows) - bad, len(rows), 'OK' if bad == 0 else 'FAILED, do not use'))
    print('frames: %d | min image distance: min %.1f, median %.1f A | frames < %.0f A: %d | max atoms within cutoff: %d'
          % (len(rows), d.min(), np.median(d), CUTOFF, (d < CUTOFF).sum(), max(r[9] for r in rows)))


if __name__ == '__main__':
    main()
