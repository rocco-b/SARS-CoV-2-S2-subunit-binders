# File formats and columns

All `.tab` files are whitespace-separated text tables written by YASARA macros (one header line, one row per analyzed snapshot, and summary rows at the end of some files).

## Classical MD (`data/03_classical_md/analysis_tables/`)
- `<system>_analysis.tab` — one row per 0.1 ns (3,001 rows). Main columns: `Time[ns]`, cell lengths, `TotalEnergy`, energy components, `RadGyration`, secondary-structure content, `RMSDCa`, `RMSDBb`, `RMSDAll` (solute RMSD from the starting structure), and for the complexes `RMSDLigMove` (ligand heavy-atom RMSD after superposing the receptor C-alpha atoms) and `RMSDLigConf` (ligand RMSD after superposing the ligand on itself).
- `<system>_pca.tab` — projections on the principal components of the backbone motion.
- `<system>_rmsfres.tab` — per-residue root-mean-square fluctuations.
- Systems: `apo` (apo-like trimer), `1830lig` (ripretinib), `1864lig` (fumagillin).

## SMD (`data/04_steered_md/analysis_tables/`)
- `<system>_analysis_<snapshot>ns.tab` — one row per 0.25 ns of pulling (about 400 rows). In addition to the columns above: `NCdist_A` (N–C terminal Cα distance, D), `RgLoc_*` (local radii of gyration; the columns `RgLoc_*` are not used in the manuscript), `RMSDCa`, `RMSDLigMove`, `RMSDLigConf`. Energy values with |Coulomb| > 10⁶ kcal/mol are artifacts of the energy re-evaluation in the analysis (see the manuscript, Section 3.6).
- `<system>_rmsfres_<snapshot>ns.tab` — per-residue fluctuations.
- `<snapshot>` is the snapshot of the classical MD (250, 275, or 300 ns) used as the starting structure.

## Pocket and periodic-image analysis (`data/04_steered_md/pocket_image_analysis/`)
- `<system>_pocket_analysis_<snapshot>ns.tab` — analyzed every 1 ns (every fourth snapshot) with `md_analyze_local_S2_POCKET.mcr`. Additional columns:
  - `C83`, `C84`, … `C294` — number of heavy atoms of each pocket residue (model numbering) within 4.5 Å of a ligand heavy atom.
  - `LigHeavy` — number of ligand heavy atoms; `LigInPocket` — number of ligand heavy atoms within 4.5 Å of a pocket-residue heavy atom.
  - `LigPockDist` — distance between the ligand center of mass and the centroid of the Cα atoms of the 16 pocket residues.
  - `NTerm_1_10`, `CTerm_330339` — Cα end-to-end distance of residues 1–10 and 330–339.
  - `ImageGapX/Y/Z` — cell length minus extent of the raw (possibly wrapped) bounding box of the protein heavy atoms.
  - `TrueExtX/Y/Z`, `TrueGapX/Y/Z` — extent of the protein along the cell axes after reconstructing the molecule across the periodic boundaries, and cell length minus this extent (a lower bound of the distance to the periodic images; it can be negative without any contact).
- `<system>_image_dist_<snapshot>ns.csv` — computed by `scripts/image_distance.py` from the coordinates in the cell frame: `time_ns`, `extent_*`, `truegap_*`, `min_image_dist` (shortest distance between a protein heavy atom and an atom of any of the 26 neighboring periodic images, in Å), `axes` (cell axes along which the shortest distance occurs), `atoms_within_8A_of_image`, `matches_macro` (1 when the reconstructed extent equals the one computed by the macro).

## Results (`results/tables/`)
- `per_replicate_metrics.csv` — one row per SMD replicate with all the per-replicate quantities used in the tables.
- `Table4_new_rows.tsv`, `TableS3_loading.tsv`, `TableS4_windows.tsv`, `TableS5_ligand_retention_terminals.tsv`, `TableS6_cell_images.tsv` — tables of the manuscript and Supporting Information, tab-separated.
- `contact_occupancy_per_replicate.csv` — fraction of frames in which each pocket residue is in contact with the ligand (Figure S7).
- `image_distance_summary.csv` — minimum image distance, number of frames below 8 Å, and first such frame, per replicate.
