#!/usr/bin/env python3
"""Checks reported in Sections 3.1 and 3.2 of the manuscript.

Run from the repository root:
    python scripts/library_and_minimization_checks.py

1. Library re-check (RDKit): records of data/02_virtual_screening/ligands_2342.sdf that fail default sanitization,
   groups of records sharing the same InChIKey, number of unique valid structures, and whether any of the 30 primary hits
   (data/02_virtual_screening/primary_hits_30.csv) is among these records.
2. Effect of the YASARA energy minimization on the AlphaFold3 model: C-alpha RMSD between the original model
   (data/01_structure_preparation/Reilly.alphafold3_694spiketrimershort_model.cif) and the minimized protein
   (data/01_structure_preparation/protein_minimized_ligand_removed.pdb), for the trimer and for the 16 pocket residues.
3. Net charge of the ligands of the selected candidates as present in the docking input (ligands_2342.sdf).
Requires numpy and rdkit. Output is printed to the console and written to results/tables/library_and_minimization_checks.txt.
"""
import csv
import itertools
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog('rdApp.*')
REPO = Path(__file__).resolve().parents[1]
VS = REPO / 'data' / '02_virtual_screening'
PREP = REPO / 'data' / '01_structure_preparation'
POCKET = [83, 84, 87, 101, 105, 107, 108, 109, 111, 112, 235, 250, 251, 252, 255, 294]   # model numbering (Spike = model + 815)
OUT = []


def say(*a):
    line = ' '.join(str(x) for x in a)
    print(line)
    OUT.append(line)


def library_check():
    suppl = Chem.SDMolSupplier(str(VS / 'ligands_2342.sdf'), sanitize=False, removeHs=False)
    keys, fail, formulas = {}, [], {}
    n = 0
    for i, m in enumerate(suppl, 1):
        n += 1
        if m is None:
            fail.append(i)
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            fail.append(i)
            continue
        keys.setdefault(inchi.MolToInchiKey(m), []).append(i)
        formulas[i] = (CalcMolFormula(m), Chem.GetFormalCharge(m))
    dups = {k: v for k, v in keys.items() if len(v) > 1}
    say('records docked:', n)
    say('records failing default RDKit sanitization:', len(fail))
    say('InChIKey groups with more than one record:', len(dups), '| records in these groups:', sum(len(v) for v in dups.values()),
        '| redundant records:', sum(len(v) - 1 for v in dups.values()))
    say('unique valid structures (InChIKey):', len(keys))
    with open(VS / 'primary_hits_30.csv', encoding='utf-8') as f:
        hits = [int(r['molecule_id']) for r in csv.DictReader(f)]
    bad = set(fail) | {i for v in dups.values() for i in v}
    say('primary hits among invalid or duplicated records:', sorted(set(hits) & bad))
    for mid, name in [(1830, 'ripretinib'), (1864, 'fumagillin'), (2168, 'afoxolaner'), (2195, 'dexamethasone palmitate'), (2006, 'odanacatib (excluded)')]:
        if mid in formulas:
            say(f'ligand {mid} ({name}) in the docking input: {formulas[mid][0]}, net charge {formulas[mid][1]:+d}')


def cif_ca(path):
    ca, hdr, inloop = {}, [], False
    for line in open(path, encoding='utf-8', errors='ignore'):
        if line.startswith('_atom_site.'):
            hdr.append(line.strip())
            inloop = True
            continue
        if inloop and line.startswith(('ATOM', 'HETATM')):
            g = dict(zip(hdr, line.split()))
            if g['_atom_site.label_atom_id'] == 'CA' and g['_atom_site.group_PDB'] == 'ATOM':
                ca[(g['_atom_site.auth_asym_id'], int(g['_atom_site.auth_seq_id']))] = [float(g['_atom_site.Cartn_' + c]) for c in 'xyz']
        elif inloop and ca and line.strip() in ('#', ''):
            break
    return ca


def pdb_ca(path):
    ca = {}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            ca[(line[21], int(line[22:26]))] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return ca


def kabsch_rotation(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    return U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt


def minimization_check():
    a = cif_ca(PREP / 'Reilly.alphafold3_694spiketrimershort_model.cif')
    b = pdb_ca(PREP / 'protein_minimized_ligand_removed.pdb')
    keys = sorted(set(a) & set(b))
    P, Q = np.array([a[k] for k in keys]), np.array([b[k] for k in keys])
    R = kabsch_rotation(P, Q)
    rmsd = np.sqrt((((P - P.mean(0)) @ R - (Q - Q.mean(0))) ** 2).sum(1).mean())
    say('C-alpha atoms compared:', len(keys), '| C-alpha RMSD trimer: %.2f A' % rmsd)
    kA = [k for k in keys if k[0] == 'A']
    Pa, Qa = np.array([a[k] for k in kA]), np.array([b[k] for k in kA])
    R = kabsch_rotation(Pa, Qa)
    dev = np.sqrt((((Pa - Pa.mean(0)) @ R - (Qa - Qa.mean(0))) ** 2).sum(1))
    idx = [kA.index(k) for k in kA if k[1] in POCKET]
    say('pocket residues (n=%d) C-alpha RMSD after chain-A superposition: %.2f A (maximum %.2f A)' % (len(idx), np.sqrt((dev[idx] ** 2).mean()), dev[idx].max()))


if __name__ == '__main__':
    library_check()
    minimization_check()
    (REPO / 'results' / 'tables').mkdir(parents=True, exist_ok=True)
    (REPO / 'results' / 'tables' / 'library_and_minimization_checks.txt').write_text('\n'.join(OUT) + '\n', encoding='utf-8')
