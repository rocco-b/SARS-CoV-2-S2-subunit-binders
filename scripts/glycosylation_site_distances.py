#!/usr/bin/env python3
"""Distance of the N-glycosylation sites Asn1098 and Asn1134 (model residues 283 and 319) from the hinge pocket.

Run from the repository root:
    python scripts/glycosylation_site_distances.py

Uses the minimized, ligand-free model (data/01_structure_preparation/protein_minimized_ligand_removed.pdb, chains A-C) and the top
poses of ripretinib (1830) and fumagillin (1864) in data/02_virtual_screening/hit_poses/ (same coordinate frame). The pocket is the set of
16 residues of Section 3.4 of chain A. The model carries no glycans: the distances are those of the asparagine itself. Hydrogens are ignored.
Writes results/tables/glycosylation_site_distances.txt.
"""
import numpy as np
from pathlib import Path
import contextlib, io
buf = io.StringIO()
_print = print
def print(*a, **k):
    _print(*a, **k); _print(*a, file=buf, **k)
P='data/01_structure_preparation/protein_minimized_ligand_removed.pdb'
atoms=[]
for l in open(P):
    if l.startswith(('ATOM','HETATM')):
        el=(l[76:78].strip() or l[12:16].strip()[0])
        if el=='H': continue
        atoms.append((l[21],int(l[22:26]),l[17:20],l[12:16].strip(),np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
chains=sorted({a[0] for a in atoms}); print('chains',chains, 'n heavy', len(atoms))
POCKET=[83,84,87,101,105,107,108,109,111,112,235,250,251,252,255,294]
pocA=[a for a in atoms if a[0]=='A' and a[1] in POCKET]
pocA_ca=np.array([a[4] for a in pocA if a[3]=='CA']); cen=pocA_ca.mean(0)
def ligand(path):
    xyz=[];m=0
    for l in open(path):
        if l.startswith('MODEL'): m+=1
        if m==1 and l.startswith(('ATOM','HETATM')):
            el=(l[76:78].strip() or l[12:16].strip()[0]).upper()
            if el!='H': xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return np.array(xyz)
lig={k:ligand(f'data/02_virtual_screening/hit_poses/molecule_{k}.pdb') for k in (1830,1864)}
for k,v in lig.items(): print('ligand',k,'heavy atoms',len(v),'COM-to-pocket-centroid',round(np.linalg.norm(v.mean(0)-cen),2))
for ch in chains:
    res=[a for a in atoms if a[0]==ch and a[1]==283]
    if not res: print(ch,'no residue 283'); continue
    name=res[0][2]
    nd2=[a for a in res if a[3]=='ND2']; xyz=np.array([a[4] for a in res]); nd2x=nd2[0][4]
    ca=[a[4] for a in res if a[3]=='CA'][0]
    dpoc=min(np.linalg.norm(xyz-a[4],axis=1).min() for a in pocA)
    dpoc_nd2=min(np.linalg.norm(nd2x-a[4]) for a in pocA)
    # nearest pocket residue
    nr=min(pocA,key=lambda a:np.linalg.norm(nd2x-a[4]))
    out=f'{ch} {name}283: Calpha-to-pocket-centroid {np.linalg.norm(ca-cen):.1f} A; ND2-to-centroid {np.linalg.norm(nd2x-cen):.1f}; min heavy-atom dist to any pocket-residue heavy atom {dpoc:.1f} (ND2 only {dpoc_nd2:.1f}, nearest {nr[2]}{nr[1]}); '
    for k,v in lig.items(): out+=f'ND2-ligand {k}: {np.linalg.norm(v-nd2x,axis=1).min():.1f} A; '
    print(out)
# for context: distance of every Asn283/N134 equivalent? also N1134 = model 319
for ch in chains:
    res=[a for a in atoms if a[0]==ch and a[1]==319]
    if res:
        nd2x=[a[4] for a in res if a[3]=='ND2'][0]
        print(ch,res[0][2],'319 (N1134): ND2 to pocket centroid %.1f; min dist to pocket heavy atoms %.1f'%(np.linalg.norm(nd2x-cen),min(np.linalg.norm(nd2x-a[4]) for a in pocA)))

Path('results/tables/glycosylation_site_distances.txt').write_text(buf.getvalue(), encoding='utf-8')
