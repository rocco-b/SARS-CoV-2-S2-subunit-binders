# SARS-CoV-2 S2 Subunit binders — data, macros, and scripts

Supporting data for the manuscript **"Virtual Screening and Molecular Dynamics Identify Ripretinib and Fumagillin as Candidate Binders of the Transient Hinge Pocket of the SARS-CoV-2 S2 Subunit"**
(R. Buccheri, E. Luca, A. Rescifina; Int. J. Mol. Sci., submitted). *[Citation and DOI of the published article: to be added.]*

The repository contains everything needed to inspect and re-derive the numbers, tables, and figures of the manuscript. The full molecular dynamics trajectories (compressed `.xtc` files, about 3 GB) are too large for GitHub and are deposited on Zenodo: https://zenodo.org/records/24377037 (see "Trajectories").

## Workflow in one paragraph

An AlphaFold3 model of the SARS-CoV-2 S2 trimer (residues 816–1154) with an open hinge pocket, provided by the Reilly and Ingber group, was minimized with YASARA and used as the docking target. A library of 2,342 compounds was docked with GNINA, and compounds were kept when CNN_VS ≥ 6.3 (30 primary hits). Four candidates were selected by visual inspection of the pose (in situ fit); two of them (ripretinib, fumagillin) passed the translational evaluation and were simulated by classical MD (300 ns, pH 4.5, 310 K) together with an apo-like model. Chain A of the apo-like model and of the two complexes, taken from the classical-MD snapshots at 250, 275, and 300 ns, was then pulled for 100 ns at its N- and C-termini (steered MD, SMD, 300 K), giving three replicates per system.

Residue numbering: the models use the numbering 1–339, which corresponds to Spike residues 816–1154 (Spike number = model number + 815).

## Repository layout

- `data/01_structure_preparation/` — the AlphaFold3 model provided by Reilly and Ingber (`.cif`), the protein after energy minimization with the ligand removed (docking target), and the ligand of the model (WYS-694), whose position defined the docking box.
- `data/02_virtual_screening/` — docking input and output.
  - `ligands_2342.sdf` (docked library), `docking_config.txt` (GNINA parameters and grid box), `docking_scores_log_2342.txt` (GNINA scores), `docking_poses_all_2342.zip` (all poses).
  - `VS.xlsx` and its CSV exports: `docking_scores_all_2342.csv`, `primary_hits_30.csv` (Table 1), `candidates_after_in_situ_evaluation.csv` (Table 2), `interaction_profile_candidates.csv` (Table 3).
  - `hit_poses/` — docked poses of the four candidates.
  - `reference_compounds_redocking/` — redocking of bemcentinib, WYS-633, and WYS-694 (Figure 1).
  - `excluded_compound_odanacatib/` — files documenting why odanacatib was excluded (see its README).
- `data/03_classical_md/` — YASARA macros for the classical MD (`md_run_S2.mcr`) and its analysis (`md_analyze.mcr`), an example starting scene, a macro that reports protonation states (`check_protonation.mcr`), and the analysis tables of the apo-like model, ripretinib, and fumagillin (energy, RMSD, ligand RMSD, PCA projections, per-residue RMSF).
- `data/04_steered_md/` — YASARA macros for the SMD (`md_runsteered_S2_REILLY_FIXED.mcr`) and its analyses (`md_analyze_local_S2_POCKET.mcr`), the analysis tables of the nine replicates, the pocket and periodic-image analyses (`pocket_image_analysis/`), and the core-distance analysis (`core_distance_analysis/`, produced by `macros/md_analyze_local_S2_CORE.mcr`: one table per replicate with the Cα distances 1–339, 5–335, 10–330, and 20–320 every 1 ns).
- `data/03_classical_md/trajectories/` and `data/04_steered_md/trajectories/` — placeholders for the compressed trajectories, which are hosted on Zenodo (see "Trajectories").
- `scripts/` — Python scripts (see below).
- `results/tables/`, `results/figures/` — outputs of the scripts that appear in the manuscript and Supporting Information.
- `docs/FILE_FORMATS.md` — description of the columns of the analysis tables.

## Reproducing the results

Python 3.10 or later and the packages in `requirements.txt` are needed (`pip install -r requirements.txt`). Run all commands from the repository root.

- `python scripts/export_vs_tables.py` — converts `VS.xlsx` to the CSV files of the screening.
- `python scripts/library_and_minimization_checks.py` — library re-check (invalid and duplicated records, unique structures), effect of the minimization on the AlphaFold3 model, and net charge of the docked ligands.
- `python scripts/core_distance_analysis.py` — Cα distances 1–339, 5–335, 10–330, and 20–320 per replicate and window (blocks B of Tables S2 and S4).
- `python scripts/smd_analysis.py` — all SMD results: per-replicate metrics, Table 4 rows, Tables S3–S6, Figure 4, Figures S6 and S7.
- `python scripts/image_distance_xtc.py FILE.xtc --stride N` — shortest distance between the protein and its periodic images directly from an `.xtc` (blocks B of Table S6 for the classical MD).
- `python scripts/image_distance.py FOLDER --tab FILE` — shortest distance between the protein and its periodic images from the coordinate dumps written by the pocket analysis macro. The coordinate dumps are not included (see below); the resulting distances are in `data/04_steered_md/pocket_image_analysis/*_image_dist_*.csv`.
- `scripts/classical_md_utilities/` — scripts that produced the stability analysis of the classical MD (Figures S2 and S3) and the force and work of the SMD. They were run in the authors' workspace layout (`papers/<slug>/` with a configuration file); they are provided for transparency, and their outputs are in `results/tables/` (`table_stability_summary.csv`, `table_ligand_pose.csv`).
- `scripts/split_docked_pdb.py`, `scripts/gnina_log_to_table.py` — helpers used to split the docked poses and to convert the GNINA log to a table.

The YASARA macros require YASARA Structure (the SMD macro and the pocket analysis macro are documented in their headers; the analysis macro requires a recent version, see `RequireVersion`).

## Main parameters

- Classical MD: AMBER14 force field, explicit water in a cuboid cell with 5 Å padding, 0.9% NaCl, 310 K, pH 4.5, PME electrostatics, 8 Å cutoff, positional restraint on the C-terminal Cα atom of each subunit, 300 ns.
- SMD: chain A only, force |F| = k·D with D the N–C terminal Cα distance and k = 5.0 kJ mol⁻¹ nm⁻² applied to the two terminal Cα atoms in opposite directions, no acceleration cap, 300 K, 100 ns, snapshots every 250 ps; starting structures are the classical-MD snapshots at 250, 275, and 300 ns (the solvent is regenerated by the macro).
- Docking: GNINA, exhaustiveness 8, `--cnn_scoring rescore --cnn_rotation 15`, grid box 46 × 44 × 40 Å (80,960 Å³) centered at (0.947, −18.434, −11.559).

## Trajectories

The trajectories were converted from the YASARA `.sim` format to GROMACS `.xtc` and are deposited on Zenodo: https://zenodo.org/records/24377037. Download the archives and place the files in the folders listed below.

- Classical MD (files `MD_apo.*`, `MD_1830lig.*` for ripretinib, `MD_1864lig.*` for fumagillin): 300 ns, **without water and ions** (protein and, for the complexes, the ligand only) to limit the file size. For each system the archive holds the trajectory (`.xtc`), a reference structure (`.pdb`, with the same atoms in the same order as the `.xtc`; the coordinates are in the orientation of the YASARA scene and are not superposed on the trajectory), and the YASARA scene (`.sce`).
- SMD (files `SMD_<system>_<snapshot>ns.*` for apo, 1830lig, and 1864lig × 250, 275, and 300 ns): chain A with solvent, 100 ns of pulling, with the trajectory, the reference structure, and the YASARA scene of each of the nine replicates. Every replicate was solvated separately, so each has its own reference structure.
- `python scripts/check_reference_pdb.py FILE.pdb FILE.xtc` verifies that a reference structure matches a trajectory (number of atoms and atom order); the 12 pairs of the archive pass the check.
- Frame interval: 100 ps for the classical MD (3,001 frames, 0–300 ns; box in the `.xtc`: cuboid cell, 5 Å padding at the start) and 250 ps for the SMD (400 frames, 0–100 ns). The cell vectors of every frame are stored in the `.xtc`.
- `scripts/image_distance_xtc.py` computes the shortest distance between the protein and its periodic images directly from an `.xtc` (no topology needed: bonds are inferred from the first frame). It reproduces the distances obtained from the YASARA coordinate dumps (`data/04_steered_md/pocket_image_analysis/*_image_dist_*.csv`) to within 0.01 Å in all 100 frames checked (ripretinib, 275 ns snapshot); the classical-MD results are in `data/03_classical_md/image_distance/`.

Because the classical-MD trajectories contain no water, quantities that depend on the solvent (for example the energy terms) cannot be recomputed from them; the energy and structural time courses of the manuscript were computed with YASARA on the full trajectories and are provided in the analysis tables.

## What is not included

- The coordinate dumps (`*_pocket_labxyz_*.tab`) from which the periodic-image distances were computed (about 7 MB per replicate); available on request.
- The molecular dynamics of odanacatib: this compound was excluded because its amine was protonated incorrectly during ligand preparation (see `data/02_virtual_screening/excluded_compound_odanacatib/README.md`).

## License

MIT License
