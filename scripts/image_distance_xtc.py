#!/usr/bin/env python3
"""Shortest distance between the protein and its periodic images, computed directly from an .xtc trajectory.

Needs numpy, scipy, and mdtraj (only for reading the .xtc); no topology file is required.
The atom connectivity is inferred from the interatomic distances of the first analyzed frame (with the minimum-image
convention): X-H bonds < 1.25 A, heavy-atom bonds < 2.1 A. Connected components with more than --min-atoms atoms
(default 1000) are taken as protein chains; water, ions, and ligands are smaller and are ignored. Each chain is
rebuilt across the periodic boundaries along the bond tree (frame by frame), the chains are placed around the first
chain, and the shortest distance between a protein heavy atom (an atom with no neighbor closer than 1.15 A) and any
atom of the 26 neighboring periodic images of the whole protein is computed, together with the number of protein
heavy atoms closer than 8 A to an image atom.

The box of every frame is read from the .xtc. Coordinates must be in the simulation-cell frame, i.e. not superposed
on a reference (a check on the rebuilt chains is printed: the largest distance between bonded atoms after rebuilding
must be below 8 A; transient close contacts that are not bonds may stretch to about 5 A).

Usage (from the repository root):
  python scripts/image_distance_xtc.py data/03_classical_md/trajectories/apo.xtc --stride 10 --out apo_image_dist.csv
  python scripts/image_distance_xtc.py data/04_steered_md/trajectories/apo_250ns.xtc --stride 4 --out x.csv
"""
import argparse
import itertools

import numpy as np
from mdtraj.formats import XTCTrajectoryFile
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components
from scipy.spatial import cKDTree

CUTOFF = 8.0
NM = 10.0  # nm -> A


def mic(d, L):
    return d - L * np.round(d / L)


def build_graph(x0, L):
    """Bonds from the first frame: (i, j) pairs and hydrogen mask.

    A pair closer than 1.15 A is an X-H bond. Which of the two atoms is the hydrogen is decided by the number of
    neighbors within 1.2 A and, for ties (e.g. R-O-H), within 1.7 A: the hydrogen has the fewer neighbors.
    """
    tree = cKDTree(x0 % L, boxsize=L)
    pairs = tree.query_pairs(2.1, output_type="ndarray")
    d = np.linalg.norm(mic(x0[pairs[:, 0]] - x0[pairs[:, 1]], L), axis=1)
    n = len(x0)
    deg12 = np.bincount(pairs[d < 1.2].ravel(), minlength=n)
    deg17 = np.bincount(pairs[d < 1.7].ravel(), minlength=n)
    key = deg12 * 1000 + deg17
    short = d < 1.15
    a, b = pairs[short, 0], pairs[short, 1]
    hmask = np.zeros(n, bool)
    hmask[np.where(key[a] <= key[b], a, b)] = True
    keep = np.where(hmask[pairs[:, 0]] | hmask[pairs[:, 1]], d < 1.25, True)
    return pairs[keep], hmask


def persistent_pairs(f, pairs, hmask, nframes, nsample=8):
    """Keep only the pairs that stay bonded (X-H < 1.35 A, heavy-heavy < 2.3 A) in frames spread over the trajectory;
    close contacts of the first frame (clashes, ions, solvent) do not survive this test."""
    keep = np.ones(len(pairs), bool)
    ish = hmask[pairs[:, 0]] | hmask[pairs[:, 1]]
    for fr in np.unique(np.linspace(0, nframes - 1, nsample).astype(int)):
        f.seek(int(fr))
        xyz, _, _, box = f.read(1)
        x, L = xyz[0] * NM, np.diag(box[0]) * NM
        d = np.linalg.norm(mic(x[pairs[:, 0]] - x[pairs[:, 1]], L), axis=1)
        keep &= np.where(ish, d < 1.35, d < 2.3)
    f.seek(0)
    return pairs[keep]


def tree_levels(pairs, n, members):
    """BFS tree of the chain: order, predecessor, depth."""
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
    g = g + g.T
    order, pred = breadth_first_order(g, members[0], directed=False, return_predecessors=True)
    depth = np.zeros(n, int)
    for i in order[1:]:
        depth[i] = depth[pred[i]] + 1
    return order, pred, depth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xtc")
    ap.add_argument("--stride", type=int, default=1, help="analyze every n-th frame")
    ap.add_argument("--min-atoms", type=int, default=1000)
    ap.add_argument("--out")
    a = ap.parse_args()

    with XTCTrajectoryFile(a.xtc) as f:
        nframes = len(f)
        xyz, time, _, box = f.read(1)
        x0 = xyz[0] * NM
        L0 = np.diag(box[0]) * NM
        pairs, hmask = build_graph(x0, L0)
        pairs = persistent_pairs(f, pairs, hmask, nframes)
        n = len(x0)
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
        ncomp, lab = connected_components(g, directed=False)
        sizes = np.bincount(lab)
        chains = [np.where(lab == c)[0] for c in np.argsort(-sizes) if sizes[c] > a.min_atoms]
        chains.sort(key=lambda m: m[0])
        print("atoms %d, frames %d, protein chains found: %s" % (n, nframes, [len(c) for c in chains]))
        trees = [tree_levels(pairs, n, m) for m in chains]
        heavy = np.concatenate([m[~hmask[m]] for m in chains])
        bonded = pairs  # for the sanity check
        in_chain = np.zeros(n, bool)
        for m in chains:
            in_chain[m] = True
        bp = bonded[in_chain[bonded[:, 0]] & in_chain[bonded[:, 1]]]

        rows, worst_bond = [], 0.0
        f.seek(0)
        k = 0
        while True:
            xyz, time, _, box = f.read(1)
            if len(xyz) == 0:
                break
            if k % a.stride == 0:
                x = xyz[0] * NM
                L = np.diag(box[0]) * NM
                u = np.zeros_like(x)
                for m, (order, pred, depth) in zip(chains, trees):
                    u[m[0]] = x[m[0]]
                    idx_all = order[1:]
                    for dpt in range(1, depth[idx_all].max() + 1):
                        idx = idx_all[depth[idx_all] == dpt]
                        u[idx] = u[pred[idx]] + mic(x[idx] - x[pred[idx]], L)
                ref = u[chains[0]].mean(0)
                for m in chains[1:]:
                    u[m] -= L * np.round((u[m].mean(0) - ref) / L)
                bond = np.linalg.norm(u[bp[:, 0]] - u[bp[:, 1]], axis=1).max()
                worst_bond = max(worst_bond, bond)
                p = u[heavy]
                tree = cKDTree(p)
                best, axes, close = np.inf, "", np.zeros(len(p), bool)
                for s in itertools.product((-1, 0, 1), repeat=3):
                    if s == (0, 0, 0):
                        continue
                    dd, _ = tree.query(p + np.array(s) * L)
                    close |= dd < CUTOFF
                    if dd.min() < best:
                        best, axes = dd.min(), "".join("xyz"[i] for i, v in enumerate(s) if v)
                ext = p.max(0) - p.min(0)
                rows.append((time[0] / 1000.0, *ext, best, axes, int(close.sum()), bond))
            k += 1

    hdr = "time_ns,extent_x,extent_y,extent_z,min_image_dist,axes,atoms_within_8A_of_image,max_bonded_distance"
    lines = [hdr] + [",".join(str(round(v, 2)) if isinstance(v, float) else str(v) for v in r) for r in rows]
    if a.out:
        open(a.out, "w").write("\n".join(lines) + "\n")
    d = np.array([r[4] for r in rows])
    print("\n".join(lines[:6]))
    print("\nCHECK: largest bonded distance after rebuilding the chains = %.2f A -> %s"
          % (worst_bond, "OK" if worst_bond < 8.0 else "FAILED, do not use"))
    print("frames analyzed: %d | min image distance: min %.1f, median %.1f A | frames < %.0f A: %d | max atoms within cutoff: %d"
          % (len(rows), d.min(), np.median(d), CUTOFF, (d < CUTOFF).sum(), max(r[6] for r in rows)))


if __name__ == "__main__":
    main()
