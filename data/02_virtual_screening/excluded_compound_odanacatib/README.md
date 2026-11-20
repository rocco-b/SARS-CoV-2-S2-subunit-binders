# Odanacatib (compound 2006): excluded from the study

Odanacatib was among the 31 compounds that met CNN_VS ≥ 6.3 in the first screening and was initially advanced to molecular dynamics. While documenting the protonation states of the ligands, we found that the automatic ligand preparation (rule-based protonation at pH 7.4) had assigned a positive charge to its secondary amine, a group that carries a trifluoromethyl group in β position and is reported to be non-basic. The docking input record (`ligands_2342.sdf`, record 2006) has the formula C25H28F4N3O3S⁺ (net charge +1) instead of the neutral C25H27F4N3O3S.

The compound was re-docked in its neutral form with the same docking parameters (box, exhaustiveness, nine modes, and CNN rescoring). The complete GNINA log is `odanacatib_neutral_redocking_gnina_log.txt`; the scores of the nine poses are also in `odanacatib_neutral_redocking_gnina_excerpt.tsv`. The first pose has CNN pose score 0.8375 and CNN affinity 6.546, so CNN_VS = 0.8375 × 6.546 = 5.48, below the threshold of 6.3. Odanacatib therefore is not a primary hit and was removed from all the analyses of the manuscript (30 primary hits, four candidates, two compounds simulated).

Files in this folder:
- `odanacatib_2006_as_docked_protonated_ligands_record.sdf` — the record of `ligands_2342.sdf` that was docked (protonated form, net charge +1).
- `odanacatib_2006_pose_protonated_form.pdb` — the pose obtained with the protonated form.
- `odanacatib_2006_as_docked_neutral_ligands_record.sdf` — the input of the re-docking (neutral form, C25H27F4N3O3S).
- `odanacatib_2006_pose_neutral_form.pdb` — the poses obtained with the neutral form.
- `odanacatib_neutral_redocking_gnina_log.txt` — the complete GNINA log of the re-docking (file paths shortened).
- `odanacatib_neutral_redocking_gnina_excerpt.tsv` — the GNINA scores of the neutral form.

In `../docking_scores_all_2342.csv` the record of odanacatib is flagged in the column `status`: its scores are those of the protonated form.

The classical MD of odanacatib, which was run with the protonated form, is not part of this repository and is not used in the manuscript.
