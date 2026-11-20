#!/usr/bin/env python3
"""Check that a reference PDB describes the atoms of an .xtc trajectory: same number of atoms, same atom order.

Needs numpy, scipy, and mdtraj. Usage (from the repository root):
  python scripts/check_reference_pdb.py data/03_classical_md/trajectories/apo.pdb data/03_classical_md/trajectories/apo.xtc

The coordinates of the PDB may differ from those of the first frame (a scene can be saved in another orientation or
before a short minimization), so they are not compared. Instead, the pairs of atoms closer than 1.7 A in the PDB
(covalent bonds and hydrogens) must still be close (< 2.5 A, minimum-image convention) in the first frame of the
.xtc: this holds for more than 99% of the pairs if the two files list the same atoms in the same order, and for far
fewer if the order or the protonation differs.
"""
import sys

import numpy as np
from mdtraj.formats import XTCTrajectoryFile
from scipy.spatial import cKDTree


def read_pdb(path):
    xyz = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(xyz)


def main():
    pdb, xtc = sys.argv[1:3]
    p = read_pdb(pdb)
    with XTCTrajectoryFile(xtc) as f:
        x, _, _, box = f.read(1)
    x = x[0] * 10.0
    L = np.diag(box[0]) * 10.0
    print("atoms: PDB %d, xtc %d" % (len(p), len(x)))
    if len(p) != len(x):
        sys.exit("FAILED: different number of atoms, this PDB does not belong to this trajectory")
    pairs = cKDTree(p).query_pairs(1.7, output_type="ndarray")
    d = x[pairs[:, 0]] - x[pairs[:, 1]]
    d -= L * np.round(d / L)
    ok = np.linalg.norm(d, axis=1) < 2.5
    frac = 100.0 * ok.mean()
    print("bonded pairs of the PDB still close in the first frame of the xtc: %.2f%% of %d" % (frac, len(pairs)))
    print("OK" if frac > 99.0 else "FAILED: the atom order (or the protonation) of the PDB differs from that of the trajectory")


if __name__ == "__main__":
    main()
