#!/usr/bin/env python3
"""Export the virtual-screening workbook (VS.xlsx) to CSV files and package the docking poses.

Run from the repository root:
    python scripts/export_vs_tables.py

Reads  data/02_virtual_screening/VS.xlsx
Writes data/02_virtual_screening/docking_scores_all_2342.csv        (all docked compounds, top pose)
       data/02_virtual_screening/primary_hits_30.csv                 (CNN_VS >= 6.3, odanacatib excluded; see the excluded-compound folder)
       data/02_virtual_screening/candidates_after_in_situ_evaluation.csv   (the four candidates of Table 2)
       data/02_virtual_screening/interaction_profile_candidates.csv  (Table 3)
Requires openpyxl.
"""
import csv
from pathlib import Path

import openpyxl

REPO = Path(__file__).resolve().parents[1]
VS = REPO / 'data' / '02_virtual_screening'
EXCLUDED_ID = 2006   # odanacatib: ligand protonated incorrectly during preparation, see excluded_compound_odanacatib/


def write(path, header, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print('wrote', path.relative_to(REPO), len(rows), 'rows')


def main():
    wb = openpyxl.load_workbook(VS / 'VS.xlsx', read_only=True, data_only=True)
    sheets = {ws.title: [r for r in ws.iter_rows(values_only=True)] for ws in wb}
    # all docked compounds: Molecule (order in ligands_2342.sdf), CNN score, CNN affinity, CNN_VS
    rows = sheets['Dati']
    write(VS / 'docking_scores_all_2342.csv', ['molecule_id', 'cnn_pose_score', 'cnn_affinity', 'cnn_vs'], rows[1:])
    # primary hits (the workbook lists 31 including odanacatib, whose docking used the wrongly protonated ligand)
    hits = [r for r in sheets['Filtrati'][1:] if r[0] is not None and r[0] != EXCLUDED_ID]
    write(VS / 'primary_hits_30.csv', ['molecule_id', 'cnn_pose_score', 'cnn_affinity', 'cnn_vs', 'in_situ'], hits)
    # indications as in Table 2 of the manuscript (the workbook has them in Italian)
    indication = {1830: 'Advanced gastrointestinal stromal tumor', 1864: 'Intestinal microsporidiosis diarrhea',
                  2168: 'Veterinary insecticide and acaricide', 2195: 'Anti-inflammatory and immunosuppressant'}
    cands = [r[:5] + (indication[r[0]],) for r in sheets['Hit'][1:] if r[0] is not None and r[0] != EXCLUDED_ID]
    write(VS / 'candidates_after_in_situ_evaluation.csv', ['molecule_id', 'cnn_pose_score', 'cnn_affinity', 'cnn_vs', 'drug', 'indication'], cands)
    # Table 3: the sheet holds a first draft and, below it, the final table with an English header
    sheet = sheets['Foglio2']
    start = next(i for i, r in enumerate(sheet) if r[0] == 'Molecule ID')
    inter = [r[:5] for r in sheet[start + 1:] if r[0] is not None and r[0] != EXCLUDED_ID]
    write(VS / 'interaction_profile_candidates.csv', ['molecule_id', 'drug', 'hydrogen_bonds', 'hydrophobic_interactions', 'aromatic_interactions'], inter)


if __name__ == '__main__':
    main()
