# YASARA MACRO
# TITLE:       Report protonation states of a loaded scene (read-only)
# DESCRIPTION: Prints the protonation state of histidines, protonated Asp/Glu, neutral Lys
#              and of the ligand N/O atoms for a scene with the protonation already assigned
#              (e.g. the *_water.sce used for the classical MD, or the 0-th snapshot of a run).
#              Nothing is modified or saved: copy the console output (Space extends it).
#
# USE: Options > Macro > Set target -> choose the scene WITHOUT extension
#      (e.g. 1864lig_water), then run this macro. Repeat for apo, 1830lig, 1864lig, 2006lig.

# FindLigand is defined in md_library (as in the SMD macro)
include md_library

Console On
Clear
LoadSce (MacroTarget)

# ---- protein: counts of titratable groups -------------------------------------------------
nasp=CountAtom Element H with bond to Atom OD1 OD2 Protein
nglu=CountAtom Element H with bond to Atom OE1 OE2 Protein
nlys=CountAtom Atom NZ Protein
hlys=CountAtom Element H with bond to Atom NZ Protein
Print ''
Print '=== SCENE: (MacroTarget) ==='
Print 'Protonated Asp, H on OD1/OD2: (nasp)'
Print 'Protonated Glu, H on OE1/OE2: (nglu)'
Print 'Lys: (nlys), H atoms on NZ: (hlys)  [3 H per Lys = all charged; neutral Lys = 3*Lys - H]'

# ---- histidines: H on ND1 / NE2 -------------------------------------------------------------
# both = HIP (+1), ND1 only = HID, NE2 only = HIE (both neutral)
hnum()=ListRes His Protein,Format='RESNUM'
hmol()=ListRes His Protein,Format='MOLNAME'
hnam()=ListRes His Protein,Format='RESNAME'
Print ''
Print 'Histidines: name, number, chain, H on ND1, H on NE2:'
for i=1 to count hnum
  nd1=CountAtom Element H with bond to Atom ND1 Res (hnum(i)) Mol (hmol(i)) Protein
  ne2=CountAtom Element H with bond to Atom NE2 Res (hnum(i)) Mol (hmol(i)) Protein
  state='neutral HID or HIE'
  if nd1 and ne2
    state='doubly protonated, charged HIP'
  Print '  (hnam(i)) (hnum(i)) chain (hmol(i)): ND1-H=(nd1) NE2-H=(ne2) -> (state)'

# ---- ligand -----------------------------------------------------------------------------------
ligsel=FindLigand 0
if ligsel==''
  Print ''
  Print 'No ligand found, apo scene.'
else
  lignam()=ListRes (ligsel),Format='RESNAME'
  ligat=CountAtom (ligsel)
  ligh=CountAtom Element H (ligsel)
  Print ''
  Print 'Ligand: (lignam1), atoms=(ligat), hydrogens=(ligh), heavy atoms=(ligat-ligh)'
  Print 'Ligand N and O atoms: atom name, heavy neighbors, H neighbors:'
  anum()=ListAtom Element N O (ligsel),Format='ATOMNUM'
  anam()=ListAtom Element N O (ligsel),Format='ATOMNAME'
  for j=1 to count anum
    nh=CountAtom Element H with bond to Atom (anum(j))
    nheavy=CountAtom Element !H with bond to Atom (anum(j))
    Print '  (anam(j)) [atom (anum(j))]: heavy neighbors=(nheavy), H=(nh)'
  Print 'Hint: an N with heavy+H = 4 is a protonated, charged amine; an O with 1 heavy neighbor'
  Print 'and 1 H is a neutral acid/alcohol; an O with 1 heavy neighbor and 0 H is carbonyl or carboxylate.'
  Print 'The net charge assigned by AutoSMILES is printed in the MD console log as'
  Print '"Topology for residue ... net charge is X": quote that line for each ligand.'
Print ''
Print '=== END ==='
