# SV10 drug-likeness by structural group (cl7 vs cl4)

Unit: ligand level = one row per ligand of interest whose chains are >= 80 % in one group; entry level = one row per entry (its majority group by the same rule, the ligand held by most of its chains). Continuous: median [IQR], Mann-Whitney; counts: units with >= 1 (share), Fisher exact. q = BH over the parameters within each level. rho = Spearman with heavy atoms (ligand level). Size-adjusted q: group term of OLS rank(parameter) ~ group + rank(heavy atoms), BH within each level (ligand / entry). Label 'adds beyond size' when the size-adjusted q < 0.05 at both levels.

These parameters describe oral drug-likeness, not binding, conformation or activity; the set includes the natural hormone, fragments and tool compounds. LE and LipE use pChEMBL, which mixes assay types (IC50, Ki, Kd from different assays).

| parameter | definition | cl7 n | cl7 median [IQR] or share | cl4 n | cl4 median [IQR] or share | p ligand | q ligand | p entry | q entry | rho heavy atoms | size-adj q ligand / entry | label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mw | Lipinski: molecular weight, Descriptors.MolWt (g/mol) | 192 | 344 [296-411] | 104 | 470 [436-512] | 2.9e-18 | 1.5e-17 | 1.5e-22 | 7.9e-22 | +0.97 | 0.027 / 0.11 | restates size |
| clogp | Lipinski: cLogP, Crippen.MolLogP | 192 | 4.58 [3.61-5.5] | 104 | 5.58 [4.87-5.95] | 1.5e-09 | 4.4e-09 | 3e-15 | 1.1e-14 | +0.59 | 0.69 / 0.14 | restates size |
| hbd | Lipinski: H-bond donors, Lipinski.NumHDonors | 192 | 2 [2-2] | 104 | 2 [1-2] | 2e-06 | 3.4e-06 | 5.1e-10 | 9e-10 | -0.14 | 0.00011 / 8.3e-08 | adds beyond size |
| hba | Lipinski: H-bond acceptors, Lipinski.NumHAcceptors | 192 | 3 [2-5] | 104 | 4 [3-6] | 6.5e-05 | 9.7e-05 | 3.9e-06 | 5.5e-06 | +0.63 | 0.0025 / 0.0027 | adds beyond size |
| lipinski_violations | Lipinski violations: #(MW > 500, cLogP > 5, HBD > 5, HBA > 10) | 192 | 85/192 (44 %) | 104 | 80/104 (77 %) | 5.3e-08 | 1.2e-07 | 5.5e-13 | 1.4e-12 | +0.69 | 0.25 / 0.66 | restates size |
| rotb | Veber: rotatable bonds, rdMolDescriptors.CalcNumRotatableBonds | 192 | 3 [2-5] | 104 | 6 [5-8] | 1.2e-22 | 1.3e-21 | 5.4e-28 | 5.7e-27 | +0.88 | 0.0013 / 1.1e-05 | adds beyond size |
| tpsa | Veber: TPSA, rdMolDescriptors.CalcTPSA (A^2) | 192 | 58.2 [40.5-87] | 104 | 60.8 [53-71.4] | 0.26 | 0.27 | 0.39 | 0.41 | +0.47 | 2.1e-05 / 2.6e-07 | adds beyond size |
| veber_violations | Veber violations: #(rotatable bonds > 10, TPSA > 140) | 192 | 7/192 (4 %) | 104 | 11/104 (11 %) | 0.022 | 0.029 | 0.024 | 0.03 | +0.38 | 0.14 / 0.11 | restates size |
| qed | QED.qed (weighted desirability, 0-1) | 192 | 0.638 [0.489-0.757] | 104 | 0.473 [0.337-0.539] | 3.3e-10 | 1.4e-09 | 4e-14 | 1.2e-13 | -0.77 | 0.066 / 0.33 | restates size |
| qed_mw | QED property MW | 192 | 344 [296-411] | 104 | 470 [436-512] | 2.9e-18 | 1.5e-17 | 1.5e-22 | 7.9e-22 | +0.97 | 0.027 / 0.11 | restates size |
| qed_alogp | QED property ALOGP | 192 | 4.58 [3.61-5.5] | 104 | 5.58 [4.87-5.95] | 1.5e-09 | 4.4e-09 | 3e-15 | 1.1e-14 | +0.59 | 0.69 / 0.14 | restates size |
| qed_hba | QED property HBA (QED acceptor SMARTS) | 192 | 3 [2-5] | 104 | 4 [4-6] | 4.4e-08 | 1.1e-07 | 1.7e-09 | 2.8e-09 | +0.69 | 0.082 / 0.078 | restates size |
| qed_hbd | QED property HBD | 192 | 2 [2-2] | 104 | 2 [1-2] | 2e-06 | 3.4e-06 | 5.1e-10 | 9e-10 | -0.14 | 0.00011 / 8.3e-08 | adds beyond size |
| qed_psa | QED property PSA | 192 | 58.2 [40.5-87] | 104 | 60.8 [53-71.4] | 0.26 | 0.27 | 0.39 | 0.41 | +0.47 | 2.1e-05 / 2.6e-07 | adds beyond size |
| qed_rotb | QED property ROTB | 192 | 3 [2-5] | 104 | 6 [5-8] | 1.2e-22 | 1.3e-21 | 5.4e-28 | 5.7e-27 | +0.88 | 0.0013 / 1.1e-05 | adds beyond size |
| qed_arom | QED property AROM (aromatic rings) | 192 | 3 [2-3] | 104 | 3 [3-3] | 3.9e-07 | 8.1e-07 | 1.1e-11 | 2.6e-11 | +0.54 | 0.86 / 0.4 | restates size |
| qed_alerts | QED property ALERTS (QED structural-alert matches) | 192 | 95/192 (49 %) | 104 | 40/104 (38 %) | 0.087 | 0.11 | 0.91 | 0.91 | +0.12 | 0.00028 / 0.0035 | adds beyond size |
| brenk_n | Brenk alerts: FilterCatalog BRENK match count | 192 | 95/192 (49 %) | 104 | 53/104 (51 %) | 0.9 | 0.9 | 0.31 | 0.36 | +0.37 | 5e-06 / 1.7e-05 | adds beyond size |
| pains_n | PAINS alerts: FilterCatalog PAINS_A + PAINS_B + PAINS_C match count | 192 | 0/192 (0 %) | 104 | 10/104 (10 %) | 2.1e-05 | 3.5e-05 | 9.2e-05 | 0.00012 | +0.14 | 0.00073 / 0.0029 | adds beyond size |
| le | Ligand efficiency LE = 1.37 * pChEMBL / heavy atoms; pChEMBL from results/potency_by_cluster.csv mixes IC50, Ki and Kd from different assays (median per ChEMBL molecule, then per ligand) | 27 | 0.441 [0.394-0.524] | 52 | 0.358 [0.329-0.39] | 0.00043 | 0.0006 | 5.6e-11 | 1.2e-10 | -0.79 | 0.75 / 0.78 | restates size |
| lipe | Lipophilic efficiency LipE = pChEMBL - cLogP; pChEMBL from results/potency_by_cluster.csv mixes IC50, Ki and Kd from different assays (median per ChEMBL molecule, then per ligand) | 27 | 3.04 [2.33-3.92] | 52 | 2.84 [2.19-3.68] | 0.25 | 0.27 | 2.7e-07 | 4e-07 | -0.33 | 0.69 / 0.29 | restates size |

QED properties equal to a Lipinski / Veber row for all 306 parsed ligands (same RDKit function): qed_mw = mw, qed_alogp = clogp, qed_hbd = hbd, qed_psa = tpsa, qed_rotb = rotb; kept as separate rows as the run file lists them. qed_hba equals hba for 250 of 306.
Direction of the size-adjusted group term (cl4 vs cl7 at equal heavy atoms): hbd cl4 lower; hba cl4 lower; rotb cl4 higher; tpsa cl4 lower; qed_hbd cl4 lower; qed_psa cl4 lower; qed_rotb cl4 higher; qed_alerts cl4 lower; brenk_n cl4 lower; pains_n cl4 higher.

Entry level: n and median [IQR] or share per group are in druglikeness_by_group.csv (ent_* columns).
Source: scripts/sv/SV10_druglikeness.py; inputs results/ligands.csv, results/ext/ligand_descriptors.csv, results/potency_by_cluster.csv, results/sv/SV00/chain_table.csv.
