## Supervisor requests ⟦2026-10-28‖FB‖Meeting (\d{4}-\d\d-\d\d)⟧

Answers to the requests of the meeting on ⟦2026-10-28‖FB‖Meeting (\d{4}-\d\d-\d\d)⟧ and to the notes received ⟦2026-11-04‖FB‖received (\d{4}-\d\d-\d\d)⟧,
produced by the unattended run `docs/sv/RUN_SV1028.md` (status `docs/sv/RUN_STATUS.md`). One paragraph and one figure
per step, written from each step's `results/sv/<ID>/SUMMARY.md` and verdict files; every number in this section is
matched to the file it comes from by `scripts/sv/SV09_number_check.py` (`results/sv/SV09/number_check.csv`). Fixed
choices: the H12 window stays ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ for all primary analyses, the superposition core stays v3, and cl7 / cl4
are the S08 clusters at ⟦1.5‖S00sum‖ARI at (1\.5) A⟧ A. Crystal structures show H12 order and position, so results are worded as H12 disorder
(H12 not modelled) or H12 mobility (local B), never as receptor activity or degradation. SV00-SV11 describe ligands by
chemistry only; the action classes of SV12-SV15 are provisional ChEMBL-mechanism labels that the teammate's labels
will replace. New limitations are in `docs/limitations.md`.

### Main figures for the next meeting

| figure | step | what it shows |
|---|---|---|
| `figures/sv/SV01b/SV01b_forest_completeness.png`, `figures/sv/SV01b/SV01b_forest_mobility.png` | SV01b (the SV01 forest plot, redrawn with lattice terms) | odds ratios (complete H12) and coefficients (H12 mobility) per predictor, coloured by block Q / X / B, within cl7 and cl4; chain level coloured, entry level grey; each block's gain entered first and last |
| `figures/sv/SV02/SV02_helix_track.png` | SV02 | pydssp helix fraction per residue for all chains, cl7 and cl4, with the consensus helices and the historical H12 windows |
| `figures/sv/SV02/rama/SV02_rama_537.png` | SV02 | phi / psi of residue ⟦537‖S02_nr‖positions/0⟧ by cluster (named a main figure in the run file) |
| `figures/sv/SV04/h12_map_partial.png` | SV03 with the SV04 overlay | the SV03 H12 map (H12 distance to the H11 axis vs H12-H3 angle) with the partial-H12 chains added, hollow, at their SV04 midpoint estimates |
| `figures/sv/SV06/SV06_outliers_map.png` | SV06 | consensus outliers on the shared MDS map (A) and their RMSD to their own vs the other complete-H12 medoid (B) |
| `figures/sv/SV07/SV07_spike_profile.png` | SV07 | the Y10 cl7-cl4 displacement profile with every spike labelled by its audit verdict, and the share of chains observing each position |
| `figures/sv/SV13/SV13_h12_class_panels.png` | SV13 | H12 and the ligand of one representative per provisional action class, with the H12 traces of all chains of the class, in one saved view |
| `figures/sv/SV13/SV13_forest.png` | SV13 | the planned class comparisons; every row shows its coverage-gate reason in place of an estimate |
| `figures/sv/SV14/SV14_pocket_vs_ligand.png` | SV14 | pocket vs ligand volume per chain as % of the E2 references, cluster as fill, action class as marker |

### Results per step

**SV00 - data restore and shared chain table.** All ⟦526‖S00m‖#count⟧ mmCIF files of the frozen entry list were restored,
⟦526‖S00m‖#count‖sha_matches_frozen=True⟧ of them with the frozen sha256, and the reproduction gate passed (ARI ⟦1.0‖S00g‖ari_by_cut/1.5⟧ at all four cuts, max |dRMSD|
⟦9.2e-15‖S00g‖max_abs_delta_rmsd_A⟧ A, ⟦0‖S00g‖n_differing_chains⟧ differing chains), so every later step works on the same ⟦788‖S00c‖regenerated‖quantity=kept_chains⟧ chains from ⟦351‖S00c‖regenerated‖quantity=kept_entries⟧ entries. The shared
chain table (`results/sv/SV00/chain_table.csv`, ⟦53‖S00t‖#ncols⟧ columns) adds dimer partner, bound peptide and SIFTS
offsets per chain; crystallisation pH is missing for ⟦285‖S00miss‖n_missing‖column=ph⟧ chains, because only the mmCIF pH field and the RCSB Data API were read.
Figure: none.

**SV01a - H12 completeness, first pass.** *Request: "check distribution of completeness, and if completeness correlates
with any experimental conditions/parameters, or is it biology".* H12 ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ is complete in ⟦431‖S01a_h‖complete‖level=chain&group=cl7⟧ of ⟦445‖S01a_h‖n‖level=chain&group=cl7⟧ cl7 chains
(⟦0.969‖S01a_h‖frac_complete‖level=chain&group=cl7⟧) but in ⟦120‖S01a_h‖complete‖level=chain&group=cl4⟧ of ⟦329‖S01a_h‖n‖level=chain&group=cl4⟧ cl4 chains (⟦0.365‖S01a_h‖frac_complete‖level=chain&group=cl4⟧); ⟦207‖S01a_pe‖#count‖side=C-side⟧ of the ⟦212‖S01a_pe‖#nrows⟧ partial H12 lack only
C-terminal residues. In multi-chain entries the H12 label differs between chains in ⟦59‖S01a_cc‖k‖measure=H12 label differs (complete/partial/not modelled)⟧ of ⟦317‖S01a_cc‖n‖measure=H12 label differs (complete/partial/not modelled)⟧ entries, the cluster
in only ⟦5‖S01a_cc‖k‖measure=S08 cluster differs⟧. Logistic models within each cluster, with blocks Q (quality), X (experiment) and B (biology: peptide,
construct mutations, ligand chemistry), give the [DEFAULT] verdict **mixed**: entry level is artifact (Q+X vs B gain
cl7 ⟦0.613‖S01a_vt‖gain_QX‖outcome=complete&level=entry&cluster=7&order=first⟧ vs ⟦0.315‖S01a_vt‖gain_B‖outcome=complete&level=entry&cluster=7&order=first⟧, cl4 ⟦0.346‖S01a_vt‖gain_QX‖outcome=complete&level=entry&cluster=4&order=first⟧ vs ⟦0.117‖S01a_vt‖gain_B‖outcome=complete&level=entry&cluster=4&order=first⟧), but at chain level cl7 is a near tie (Q+X ⟦0.224‖S01a_vt‖gain_QX‖outcome=complete&level=chain&cluster=7&order=first⟧ vs B ⟦0.225‖S01a_vt‖gain_B‖outcome=complete&level=chain&cluster=7&order=first⟧; cl4
⟦0.208‖S01a_vt‖gain_QX‖outcome=complete&level=chain&cluster=4&order=first⟧ vs ⟦0.079‖S01a_vt‖gain_B‖outcome=complete&level=chain&cluster=4&order=first⟧). cl7 has only ⟦14‖S01a_mb‖n_incomplete‖outcome=complete&cluster=7&level=chain&block=Q&order=first⟧ incomplete chains for ⟦16‖S01a_mb‖n_predictors_full‖outcome=complete&cluster=7&level=chain&block=Q&order=first⟧ predictors. Figure (backup):
`figures/sv/SV01a/SV01a_completeness_profile.png`.

**SV01b - with coordinates: lattice contacts and H12 mobility.** Two additions: heavy-atom lattice contacts of residues
⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧ with ASU partners and symmetry mates, added to block X, and a second outcome, local H12 mobility (C-alpha B
z-score of the observed residues of ⟦528‖S01b_sum‖observed residue of (528)-547⟧-⟦547‖S01b_sum‖observed residue of 528-(547)⟧ against the same chain's core v3). The re-issued verdict is
**artifact** for both outcomes at chain and entry level: H12 complete, chain level, Q+X vs B gain cl7 ⟦0.414‖S01b_vt‖gain_QX‖outcome=complete&level=chain&cluster=7&order=first⟧ vs ⟦0.225‖S01b_vt‖gain_B‖outcome=complete&level=chain&cluster=7&order=first⟧ and
cl4 ⟦0.523‖S01b_vt‖gain_QX‖outcome=complete&level=chain&cluster=4&order=first⟧ vs ⟦0.079‖S01b_vt‖gain_B‖outcome=complete&level=chain&cluster=4&order=first⟧; H12 mobility cl7 ⟦0.184‖S01b_vt‖gain_QX‖outcome=mobility&level=chain&cluster=7&order=first⟧ vs ⟦0.076‖S01b_vt‖gain_B‖outcome=mobility&level=chain&cluster=7&order=first⟧ and cl4 ⟦0.210‖S01b_vt‖gain_QX‖outcome=mobility&level=chain&cluster=4&order=first⟧ vs ⟦0.143‖S01b_vt‖gain_B‖outcome=mobility&level=chain&cluster=4&order=first⟧; the
sensitivity runs (fraction observed, entries without TLS refinement) agree. The change from SV01a is the lattice
terms' share: per SD of symmetry-mate contact atoms the odds of a complete H12 rise ⟦20.2‖S01b_or‖estimate‖outcome=complete&role=primary&cluster=4&level=chain&predictor=lattice_sym⟧-fold in cl4 (⟦12.6‖S01b_or‖estimate‖outcome=complete&role=primary&cluster=7&level=chain&predictor=lattice_sym⟧-fold in cl7)
and H12 mobility falls by ⟦0.268‖S01b_or‖estimate‖outcome=mobility&cluster=4&level=chain&predictor=lattice_sym‖abs⟧ z in cl4 (⟦0.17‖S01b_or‖estimate‖outcome=mobility&cluster=7&level=chain&predictor=lattice_sym‖abs⟧ in cl7). Caveat (first line and Reading of the SUMMARY): for completeness
this share is partly mechanical, because contacts are counted on modelled atoms. In cl4, complete chains carry
⟦24.3‖S01b_lg‖mean_sym_atoms‖cluster_group=cl4&h12_status=complete⟧ symmetry-contact atoms on ⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧ on average and partial chains ⟦9.1‖S01b_lg‖mean_sym_atoms‖cluster_group=cl4&h12_status=partial⟧, but ⟦12.3‖S01b_lg‖mean_sym_atoms_546_550‖cluster_group=cl4&h12_status=complete⟧ of the former sit on ⟦546‖S01b_sum‖sit on (546)-550⟧-⟦550‖S01b_sum‖sit on 546-(550)⟧, which
partial chains hardly model; on ⟦538‖S01b_sum‖On (538)-545, which both model⟧-⟦545‖S01b_sum‖On 538-(545), which both model⟧, modelled by both, the means are ⟦8.22‖S01b_lg‖mean_sym_atoms_538_545‖cluster_group=cl4&h12_status=complete⟧ vs ⟦6.66‖S01b_lg‖mean_sym_atoms_538_545‖cluster_group=cl4&h12_status=partial⟧. The data cannot tell whether a
contact holds the C-terminal end in place or exists because the end is modelled; the mobility outcome is not tied to
missing atoms in that way. The medoid clash test was not run (it runs only if the verdict is mixed). Main figure: the
SV01 forest plots `figures/sv/SV01b/SV01b_forest_completeness.png` and `figures/sv/SV01b/SV01b_forest_mobility.png`.

**SV01c - sensitivity to the H12 window.** On the data-derived window ⟦538‖S01c_n‖window_proposed/0⟧-⟦544‖S01c_n‖window_proposed/1⟧ (SV02 proposal), ⟦204‖S01c_st‖n‖cluster_group=all&level=chain&h12_status_538_547=partial&h12_status_538_544=complete⟧ chains that are partial
on ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ become complete, so an incomplete H12 becomes a rare event (cl7 ⟦14‖S01c_n‖n_incomplete_chains/538_547/cl7⟧ -> ⟦10‖S01c_n‖n_incomplete_chains/538_544/cl7⟧ chains, cl4 ⟦209‖S01c_n‖n_incomplete_chains/538_547/cl4⟧ -> ⟦19‖S01c_n‖n_incomplete_chains/538_544/cl4⟧),
mostly H12 not modelled at all. The completeness verdict moves from artifact to **mixed** at chain level only (cl7 Q+X
⟦0.130‖S01c_sb‖gain_QX_538_544‖outcome=complete&level=chain&cluster=7&order=first⟧ vs B ⟦0.386‖S01c_sb‖gain_B_538_544‖outcome=complete&level=chain&cluster=7&order=first⟧, cl4 ⟦0.506‖S01c_sb‖gain_QX_538_544‖outcome=complete&level=chain&cluster=4&order=first⟧ vs ⟦0.216‖S01c_sb‖gain_B_538_544‖outcome=complete&level=chain&cluster=4&order=first⟧); entry level stays artifact but rests on ⟦4‖S01c_n‖n_incomplete_entries/538_544/cl7⟧ (cl7) and ⟦6‖S01c_n‖n_incomplete_entries/538_544/cl4⟧ (cl4)
incomplete entries, and H12 mobility stays artifact. Because the verdict was mixed, the medoid clash test ran: of ⟦31‖S01c_cl‖n_chains‖cluster_group=all&measure=both positions blocked⟧
partial or not-modelled chains on ⟦538‖S01c_n‖window_proposed/0⟧-⟦544‖S01c_n‖window_proposed/1⟧, ⟦3‖S01c_cl‖k‖cluster_group=all&measure=both positions blocked⟧ have both medoid H12 positions blocked by lattice neighbours, ⟦0‖S01c_cl‖k‖cluster_group=all&measure=only cl7 position blocked⟧ only the
cl7 and ⟦0‖S01c_cl‖k‖cluster_group=all&measure=only cl4 position blocked⟧ only the cl4 position, ⟦28‖S01c_cl‖k‖cluster_group=all&measure=neither position blocked⟧ neither, so the lattice does not favour one H12 position. The
completeness verdict therefore depends on the H12 window, which the team has yet to sign off. Figure: none (tables
`results/sv/SV01c/side_by_side.csv` and `results/sv/SV01c/clash_summary.csv`).

**SV02 - consensus secondary structure.** *Request: "consensus alignment of secondary structures".* pydssp on all
⟦788‖S02_nc‖groups/global⟧ chains (fallback for mkdssp, which is not installed; its 'H' covers every helix type) gives ⟦8‖S02_nc‖n_consensus/global⟧ global consensus
helices (at least ⟦4‖S02_hw‖>= (\d+) consecutive positions⟧ consecutive positions helical in at least ⟦80‖S02_hw‖in >= (\d+) % of the observing⟧ % of the observing chains). The data-derived H12
window is ⟦538‖S02_hb‖global_start‖helix_name=H12⟧-⟦544‖S02_hb‖global_end‖helix_name=H12⟧ (cl7 ⟦538‖S02_hb‖cl7_start‖helix_name=H12⟧-⟦545‖S02_hb‖cl7_end‖helix_name=H12⟧, cl4 ⟦537‖S02_hb‖cl4_start‖helix_name=H12⟧-⟦544‖S02_hb‖cl4_end‖helix_name=H12⟧), a proposal only: it matches none of the
historical windows (Jaccard ⟦0.70‖S02_nc‖h12_comparisons/0/jaccard⟧ with each of ⟦538‖S02_nc‖h12_comparisons/0/start⟧-⟦547‖S02_nc‖h12_comparisons/0/end⟧, ⟦537‖S02_nc‖h12_comparisons/1/start⟧-⟦546‖S02_nc‖h12_comparisons/1/end⟧ and ⟦536‖S02_nc‖h12_comparisons/2/start⟧-⟦545‖S02_nc‖h12_comparisons/2/end⟧), and ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ stays
primary. The clusters also differ in where other helices end: H10/11 ends at ⟦530‖S02_hb‖cl7_end‖helix_name=H10/11⟧ in cl7 and ⟦525‖S02_hb‖cl4_end‖helix_name=H10/11⟧ in cl4, H9
starts at ⟦473‖S02_hb‖cl7_start‖helix_name=H9⟧ in cl7 and ⟦466‖S02_hb‖cl4_start‖helix_name=H9⟧ in cl4. Residue ⟦537‖S02_rs‖position‖position=537⟧ alone separates the clusters: circular-mean psi
⟦+156.6‖S02_rs‖psi_mean_cl7‖position=537⟧ deg in cl7 vs ⟦-46.8‖S02_rs‖psi_mean_cl4‖position=537⟧ deg in cl4 (separation ⟦157.2‖S02_rs‖separation_cl7_cl4‖position=537⟧ deg; overlap ⟦0‖S02_rs‖overlap_chain‖position=537⟧ at chain and ⟦0‖S02_rs‖overlap_entry‖position=537⟧ at entry
level), as recorded earlier (⟦+157‖S02_nr‖check_537_psi/cl7/expected_psi⟧ / ⟦-47‖S02_nr‖check_537_psi/cl4/expected_psi⟧). Caveat: pydssp has no chain-break handling, so beta bridges
across breaks are missed (helices are not affected). Main figures: `figures/sv/SV02/SV02_helix_track.png` and
`figures/sv/SV02/rama/SV02_rama_537.png`; the SSE heatmap and the other Ramachandran plots are backup.

**SV03 - helix geometry.** *Request: "we can try geometrical approach since it's more sensitive for things like the
length of the helix".* ⟦16‖S03_nt‖n_tests/(a) chain/size⟧ descriptors per chain (helix lengths, inter-helix angles, H12 distances to the H3 and
H11 axes, H11 kink, H11-H12 hinge) in the core v3 frame. cl7 vs cl4: ⟦12‖S03_nt‖a_separates_both_levels‖‖len⟧ of ⟦16‖S03_nt‖n_tests/(a) entry/size⟧ differ at q < ⟦0.05‖S03_sum‖differ at q < (0\.05)⟧ at chain and
entry level, all ⟦7‖S03_sum‖all (\d+) H12 descriptors⟧ H12 descriptors among them; the largest are the H12-H3 angle (median ⟦79.05‖S03_w‖median0_chain‖descriptor=angle_H12_H3_deg⟧ vs ⟦46.63‖S03_w‖median1_chain‖descriptor=angle_H12_H3_deg⟧ deg)
and the distance from the H12 centroid to the H11 axis (⟦12.48‖S03_w‖median0_chain‖descriptor=dist_H12c_H11axis_A⟧ vs ⟦22.26‖S03_w‖median1_chain‖descriptor=dist_H12c_H11axis_A⟧ A). Within clusters, ⟦37‖S03_nt‖b_meets_default⟧ covariate
comparisons meet the [DEFAULT] rule (|r| >= ⟦0.3‖S03_sum‖\|r\| >= (0\.3)⟧ and q < ⟦0.05‖S03_sum‖\|r\| >= 0\.3 and q < (0\.05)⟧ at both levels; ⟦19‖S03_sum‖\((\d+) in cl7, \d+ in cl4⟧ in cl7, ⟦18‖S03_sum‖\(\d+ in cl7, (\d+) in cl4⟧ in cl4), so by that
rule the descriptors are more sensitive than RMSD for those subgroups. The effects are small (angles ⟦0.5‖S03_sum‖kink and hinge (0\.5)-5\.4 deg⟧-⟦5.4‖S03_sum‖kink and hinge 0\.5-(5\.4) deg⟧ deg,
lengths ⟦1‖S03_sum‖lengths (1)-4 residues⟧-⟦4‖S03_sum‖lengths 1-(4) residues⟧ residues) and ⟦4‖S03_sum‖; (\d+) of them rest on a group of <= 5 entries⟧ of the ⟦37‖S03_nt‖b_meets_default⟧ rest on a group of at most ⟦5‖S03_sum‖rest on a group of <= (5) entries⟧ entries. Caveat: helix ends
come from the pydssp fallback, and cl4 H12 axes rest mostly on ⟦7‖S03_sum‖cl4 H12 axes rest mostly on (\d+) C-alpha⟧ C-alpha, where the axis fit on an ideal helix
is off by ⟦4.5‖S03_ax‖angle_error_deg‖n_ca=7⟧ deg, so H12 angle differences of a few degrees inside a cluster are within the method's range.
Figure (backup): `figures/sv/SV03/SV03_descriptors_by_cluster.png`; the SV03 H12 map is shown with the SV04 overlay.

**SV04 - helix midpoints for incomplete helices.** *Request (Ye): compare incomplete helices by "taking the middle
atoms".* Each helix is reduced to the midpoint of its fitted axis, validated by deleting residues from complete
helices; the H12 curve rests on the ⟦75‖S04_ve‖n_complete_chains‖helix=H12&end=N&k_removed=0⟧ chains whose whole ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ is helical. Clustering chains on midpoints
reproduces S08: ARI ⟦0.954‖S04_ari‖ari‖distance=midpoint_rms&level=chain&K=2⟧ at K = ⟦2‖S04_ari‖K‖distance=midpoint_rms&level=chain&n=774⟧ (⟦774‖S04_ari‖n‖distance=midpoint_rms&level=chain&K=2⟧ cl7 and cl4 chains) and ⟦0.931‖S04_ari‖ari‖distance=midpoint_rms&level=chain&K=7⟧ at K = ⟦7‖S04_ari‖K‖distance=midpoint_rms&level=chain&n=788⟧ (⟦788‖S04_ari‖n‖distance=midpoint_rms&level=chain&K=7⟧ chains);
entry level ⟦0.977‖S04_ari‖ari‖distance=midpoint_rms&level=entry&K=2⟧ / ⟦0.959‖S04_ari‖ari‖distance=midpoint_rms&level=entry&K=7⟧. This mostly restates the H12 split - the cl7 and cl4 H12 midpoints lie a median
⟦12.58‖S04_ds‖cl7_vs_cl4_median‖helix=H12&measure=midpoint_distance_A⟧ A apart, every other helix at most ⟦0.30‖S04_sum‖every other helix <= (0\.30) A⟧ A - and is not independent evidence, since S08 used the
same observed H12 C-alpha. Of ⟦212‖S04_nc‖partial/n⟧ partial-H12 chains, ⟦199‖S04_nc‖partial/h12_eligible⟧ have a usable H12 midpoint, and every one of the ⟦191‖S04_sum‖\((\d+) of 191\)⟧
in cl7 or cl4 lands on its own S08 side. Main figure (with SV03): `figures/sv/SV04/h12_map_partial.png`; panel (a)
places complete chains at SV03 positions and partial chains at SV04 estimates, panel (b) uses the SV04 estimator for all
chains, because midpoint and centroid differ by a cluster-dependent offset in x (complete cl4 ⟦+1.51‖S04_nc‖sv04_minus_sv03_map_x_by_status_group/complete|cl4/median⟧ A, cl7 ⟦-0.68‖S04_nc‖sv04_minus_sv03_map_x_by_status_group/complete|cl7/median⟧ A).

**SV05 - stable core v4 from complete helices.** *Request: "take helices without the problem of completeness, find a
stable core".* From ⟦132‖S05_n‖n_candidates⟧ candidate positions (SV02 consensus helices outside ⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧, observed in at least ⟦98‖S05_sum‖observed in < (98) % of the 788 chains⟧ % of
chains), iterative superposition keeps ⟦128‖S05_n‖n_core_v4⟧ positions with spread at most ⟦0.5‖S05_sum‖spread <= (0\.5) A⟧ A as core v4: it adds H1 and H9,
drops ⟦25‖S05_sum‖v3 positions not in v4: (\d+)⟧ v3 positions at helix ends, and has Jaccard ⟦0.64‖S05_n‖jaccard_v4_v3⟧ with v3. Recomputed with v4, pairwise RMSDs barely
change (Spearman ⟦0.9975‖S05_n‖spearman_pairwise_v4_vs_v3⟧ with the primary matrix; ARI with S08 ⟦0.990‖S05_ari‖ari_v4_vs_s08‖n_clusters_v3_s08=7⟧ at ⟦1.5‖S05_ari‖cut_A‖n_clusters_v3_s08=7⟧ A, ⟦3‖S05_ari‖chains_changed_best_match‖n_clusters_v3_s08=7⟧ chains change cluster), and
⟦98.5‖S05_n‖profile_v4/share_528_550‖‖x100⟧ % of the squared cl7-cl4 difference stays in ⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧ (v3 ⟦98.6‖S05_n‖profile_v3_y10/share_528_550‖‖x100⟧ %). Whether to adopt v4 is a team decision;
SV03, SV04 and SV07 used v3. Figure (backup): `figures/sv/SV05/SV05_core_v4.png`.

**SV06 - outliers on the clustering map.** *Request: "color outliers in clustering (but no pattern is also
possible)".* Of the ⟦22‖S06_n‖n_outliers⟧ consensus outliers (Y15), ⟦2‖S06_n‖outliers_by_group/cl7⟧ are in cl7 (⟦445‖S06_n‖chains_by_group/cl7⟧ chains), ⟦12‖S06_n‖outliers_by_group/cl4⟧ in cl4 (⟦329‖S06_n‖chains_by_group/cl4⟧) and ⟦8‖S06_n‖outliers_by_group/minor⟧ in the minor
clusters (⟦14‖S06_n‖chains_by_group/minor⟧ chains). Between the main clusters this is a pattern (hypergeometric, BH over ⟦8‖S06_h‖#count‖tested=True⟧ tests): ⟦2‖S06_h‖observed‖level=chain&test=cl7 vs cl4 (main only)⟧ in cl7
against ⟦8.0‖S06_h‖expected‖level=chain&test=cl7 vs cl4 (main only)⟧ expected (q ⟦0.0023‖S06_h‖q_bh‖level=chain&test=cl7 vs cl4 (main only)⟧), and at entry level ⟦2‖S06_h‖observed‖level=entry&test=cl7 vs cl4 (main only)⟧ against ⟦7.0‖S06_h‖expected‖level=entry&test=cl7 vs cl4 (main only)⟧ (q ⟦0.0034‖S06_h‖q_bh‖level=entry&test=cl7 vs cl4 (main only)⟧). The minor-cluster
enrichment is expected by construction, because the detectors measure distance from the main clusters. None
(⟦0‖S06_n‖n_main_closer_to_other⟧ of ⟦14‖S06_n‖n_main⟧) of the main-cluster outliers is closer to the other cluster's medoid than to its own (own/other
RMSD ratio ⟦0.21‖S06_n‖ratio_main_min⟧-⟦0.63‖S06_n‖ratio_main_max⟧), so they sit at the edge of their own cluster rather than between the two. Caveat:
panel B uses the complete-H12 medoids `5TM1|A|1` and `7UJ7|C|1` (⟦3.68‖S06_n‖rmsd_medoid_cl7_cl4⟧ A apart), because the S12 cl4 medoid has no
modelled H12. Main figure: `figures/sv/SV06/SV06_outliers_map.png`.

**SV07 - spike audit of the displacement profile.** *Request: "check spikes on figure in slide ⟦5‖FB‖figure in slide (\d)⟧".* ⟦19‖S07_n‖n_spikes⟧ positions
outside ⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧ exceed ⟦1.0‖S07_sum‖C-alpha distance > (1\.0) A⟧ A in the Y10 profile (cl7-cl4 C-alpha distance of the cluster means or of the medoid
pair). ⟦14‖S07_n‖verdict_counts_chain/genuine loop difference⟧ are genuine loop differences (⟦331‖S07_v‖position‖position=331&verdict=genuine loop difference⟧-⟦340‖S07_v‖position‖position=340&verdict=genuine loop difference⟧, ⟦464‖S07_v‖position‖position=464&verdict=genuine loop difference⟧, ⟦466‖S07_v‖position‖position=466&verdict=genuine loop difference⟧-⟦468‖S07_v‖position‖position=468&verdict=genuine loop difference⟧): the distance stays above ⟦1.0‖S07_sum‖stays > (1\.0) A after removing⟧ A
after removing the ⟦3‖S07_sum‖removing the (\d+) most influential⟧ most influential entries, at chain and entry level. ⟦5‖S07_n‖verdict_counts_chain/driven by <= 3 entries⟧ are driven by at most ⟦3‖S07_sum‖driven by <= (\d+) entries⟧ entries:
⟦330‖S07_v‖position‖position=330&exceeds=medoid&verdict=driven by <= 3 entries⟧, ⟦397‖S07_v‖position‖position=397&exceeds=medoid&verdict=driven by <= 3 entries⟧, ⟦420‖S07_v‖position‖position=420&exceeds=medoid&verdict=driven by <= 3 entries⟧ and ⟦459‖S07_v‖position‖position=459&exceeds=medoid&verdict=driven by <= 3 entries⟧ exceed only on the medoid pair, and ⟦462‖S07_v‖position‖position=462&exceeds=mean&verdict=driven by <= 3 entries⟧ drops to ⟦0.96‖S07_s‖d_mean_minus_1_entries_A‖position=462⟧ A without one entry. No
spike is a numbering error. Among the genuine spikes, ligand chemotype is associated (q < ⟦0.05‖S07_sum‖associated \(q < (0\.05) at chain and entry level\)⟧ at both levels) at
⟦331‖S07_v‖position‖position=331&covariates_associated~ligand_chemotype⟧-⟦332‖S07_v‖position‖position=332&covariates_associated~ligand_chemotype⟧ and ⟦334‖S07_v‖position‖position=334&covariates_associated~ligand_chemotype⟧-⟦339‖S07_v‖position‖position=339&covariates_associated~ligand_chemotype⟧, dimer-partner source at ⟦334‖S07_v‖position‖position=334&covariates_associated~dimer_partner_source⟧-⟦340‖S07_v‖position‖position=340&covariates_associated~dimer_partner_source⟧ and ⟦467‖S07_v‖position‖position=467&covariates_associated~dimer_partner_source⟧-⟦468‖S07_v‖position‖position=468&covariates_associated~dimer_partner_source⟧, and the bound peptide at
⟦331‖S07_v‖position‖position=331&covariates_associated~peptide_bound⟧ and ⟦466‖S07_v‖position‖position=466&covariates_associated~peptide_bound⟧-⟦468‖S07_v‖position‖position=468&covariates_associated~peptide_bound⟧; at ⟦335‖S07_ct‖position‖position=335&covariate=ligand_chemotype&level=chain⟧ the cl4-cl7 coefficient falls from ⟦2.74‖S07_ct‖b_cl4_without_cov_A‖position=335&covariate=ligand_chemotype&level=chain⟧ to ⟦0.51‖S07_ct‖b_cl4_with_cov_A‖position=335&covariate=ligand_chemotype&level=chain⟧ A once chemotype is in
the model, but chemotype is unevenly spread over the clusters, so the adjusted value rests on few chains. Partly
BLOCKED: the 'sparse-sample effect' rule is degenerate, because Y10 already averages only the chains that observe a
position (restricted means equal Y10 within ⟦2e-14‖S07_n‖restricted_vs_y10/max_abs_diff_d_mean_A⟧ A). Main figure: `figures/sv/SV07/SV07_spike_profile.png`.

**SV08 - literature: is only H12 different?** *Request: "check literature if only H12 is significant? so far seems
so".* From ⟦14‖S08_sum‖\*\*(\d+) distinct sources opened\*\*⟧ opened sources (⟦25‖S08_sum‖\*\*(\d+) claims\*\*⟧ claims; ⟦8‖S08_sum‖\*\*(\d+) sources\s+could not be opened⟧ sources could not be opened) the answer is no: they report differences
between ERalpha LBD structures in H11 (shifting with ligand chemistry), H3, H7, the H11-H12 loop and the
His524 / Leu525 rotamers, and one deposited structure leaves residues ⟦331‖S08_lit‖correspond to residues (331)-340⟧-⟦340‖S08_lit‖correspond to residues 331-(340)⟧ and ⟦462‖S08_lit‖residues (462)-463 in the loop⟧-⟦463‖S08_lit‖residues 462-(463) in the loop⟧ unmodelled
(Li et al. ⟦2011‖S08_lit‖Li et al\. (2011), J Med Chem⟧) - the two loop regions where the genuine SV07 spikes lie. Y537S, D538G and L536S are documented as
conformation-driving construct mutations, but no opened source gives the residue-level effect of L536S (carried by
⟦243‖S01a_h‖n‖level=chain&group=L536S⟧ of our chains), and none attributes a difference outside H12 to the coactivator peptide. Table:
`docs/sv/literature.md`. Figure: none.

**SV10 - ligand drug-likeness.** Of ⟦21‖S10_n‖n_parameters⟧ drug-likeness parameters of the crystallised ligands, ⟦10‖S10_n‖adds_beyond_size‖‖len⟧ add information
beyond ligand size (⟦7‖S10_n‖n_distinct_beyond_size⟧ distinct quantities: H-bond donors and acceptors, rotatable bonds, TPSA, QED alerts, Brenk and
PAINS alerts) and ⟦11‖S10_n‖restates_size‖‖len⟧ restate size (heavy atoms median cl7 ⟦24‖S10_n‖heavy_atoms_ligand/cl7/median⟧ vs cl4 ⟦34‖S10_n‖heavy_atoms_ligand/cl4/median⟧); table
`results/sv/SV10/druglikeness_by_group.md`.

**SV11 - pocket sites and steroid flag.** Part A: P2Rank, which saw only the receptor chain and its ligand (no peptide,
no partner chain), ranks the ligand pocket first in ⟦788‖S11_n‖label_counts_rank1/ligand pocket⟧ of ⟦788‖S11_n‖n_chains⟧ chains, with no cl7-cl4 difference there. A
pocket in the AF-2 groove is found in ⟦0.53‖S11_st‖value_cl7‖site=AF-2 groove&outcome=share with pocket&level=chain⟧ of cl7 chains (⟦234‖S11_st‖k_cl7‖site=AF-2 groove&outcome=share with pocket&level=chain⟧ of ⟦445‖S11_st‖n_cl7‖site=AF-2 groove&outcome=share with pocket&level=chain⟧) but ⟦0.02‖S11_st‖value_cl4‖site=AF-2 groove&outcome=share with pocket&level=chain⟧ of cl4 chains (⟦5‖S11_st‖k_cl4‖site=AF-2 groove&outcome=share with pocket&level=chain⟧ of ⟦329‖S11_st‖n_cl4‖site=AF-2 groove&outcome=share with pocket&level=chain⟧;
entries ⟦0.72‖S11_st‖value_cl7‖site=AF-2 groove&outcome=share with pocket&level=entry⟧ vs ⟦0.01‖S11_st‖value_cl4‖site=AF-2 groove&outcome=share with pocket&level=entry⟧), and within cl7 only in chains with a bound peptide (⟦234‖S11_sum‖occur in (\d+)/414 peptide-bound⟧ of ⟦414‖S11_sum‖occur in 234/(\d+) peptide-bound⟧, ⟦0‖S11_sum‖and (\d+)/31 without a peptide⟧ of
⟦31‖S11_sum‖and 0/(\d+) without a peptide⟧ without); a dimer-interface pocket in ⟦0.18‖S11_st‖value_cl7‖site=dimer interface&outcome=share with pocket&level=chain⟧ vs ⟦0.31‖S11_st‖value_cl4‖site=dimer interface&outcome=share with pocket&level=chain⟧ of chains (q ⟦0.0034‖S11_st‖q‖site=dimer interface&outcome=share with pocket&level=chain⟧ chain level, ⟦0.022‖S11_st‖q‖site=dimer interface&outcome=share with pocket&level=entry⟧ entry level).
Both are low-probability P2Rank pockets (median probability ⟦0.010‖S11_sum‖median probability (0\.010) and⟧ and ⟦0.013‖S11_sum‖median probability 0\.010 and (0\.013)⟧ vs ⟦0.924‖S11_sum‖and 0\.013 vs (0\.924) for the ligand pocket⟧ for the ligand pocket),
and the groove rank contrasts rest on ⟦5‖S11_st‖n_cl4‖site=AF-2 groove&outcome=best rank&level=chain⟧ cl4 chains. Part B: a steroid ring system is present in ⟦23‖S11_ster‖n_ligands_flagged_all‖level=ligand⟧ of ⟦306‖S11_ster‖n_ligands_all‖level=ligand⟧
ligands, more often in cl7 entries (⟦41‖S11_ster‖steroid_cl7‖level=entry⟧ of ⟦217‖S11_ster‖n_cl7‖level=entry⟧) than in cl4 entries (⟦5‖S11_ster‖steroid_cl4‖level=entry⟧ of ⟦124‖S11_ster‖n_cl4‖level=entry⟧; q ⟦0.00011‖S11_ster‖q‖level=entry⟧), and not
significantly at ligand level (q ⟦0.071‖S11_ster‖q‖level=ligand⟧). Figure (backup): `figures/sv/SV11/SV11_site_share.png`.

#### Addition: notes received ⟦2026-11-04‖FB‖received (\d{4}-\d\d-\d\d)⟧

**SV12 - ligand identifiers and provisional action classes.** *Request: "⟦6‖FB‖(\d) different action class types⟧ different action class types and based
on their DOIs, show in the plot"; "get the SMILES, InChiKey from api and using chembl api".* SMILES and InChIKeys came
from RCSB for all ⟦306‖S12_n‖n_ccd⟧ ligand codes and were matched to ChEMBL. The teammate's label file did not exist at run time, so
classes come from ChEMBL mechanism records only; the DOI route failed its validation gate (⟦1‖S12_n‖gate/n_agree⟧ of ⟦2‖S12_n‖gate/n_compared⟧ rows agree)
and is not used. ⟦51‖S12_sum‖covers (\d+) of 362 rows⟧ of ⟦362‖S12_n‖n_rows⟧ entry x ligand rows are labelled - agonist ⟦27‖S12_n‖final_rows/agonist⟧, SERM ⟦20‖S12_n‖final_rows/SERM⟧, SERD ⟦4‖S12_n‖final_rows/SERD⟧; partial agonist
⟦0‖S12_cov‖n_entries‖set=primary&action_class=partial_agonist&source=any&cluster=all⟧, antagonist ⟦0‖S12_cov‖n_entries‖set=primary&action_class=antagonist&source=any&cluster=all⟧, other ⟦0‖S12_cov‖n_entries‖set=primary&action_class=other&source=any&cluster=all⟧ - and ⟦311‖S12_n‖final_rows/unlabelled⟧ stay unlabelled. Only agonist (⟦27‖S12_cov‖n_entries‖set=primary&action_class=agonist&source=any&cluster=all⟧
entries) and SERM (⟦20‖S12_cov‖n_entries‖set=primary&action_class=SERM&source=any&cluster=all⟧) pass the coverage gate of ⟦5‖S12_sum‖only if each side has >= (\d+) entries⟧ entries per side, and the classes nearly coincide
with the clusters (agonist entries ⟦25‖S12_cov‖n_entries‖set=primary&action_class=agonist&source=any&cluster=cl7⟧ in cl7 and ⟦1‖S12_cov‖n_entries‖set=primary&action_class=agonist&source=any&cluster=cl4⟧ in cl4; SERM ⟦3‖S12_cov‖n_entries‖set=primary&action_class=SERM&source=any&cluster=cl7⟧ and ⟦19‖S12_cov‖n_entries‖set=primary&action_class=SERM&source=any&cluster=cl4⟧; SERD ⟦0‖S12_cov‖n_entries‖set=primary&action_class=SERD&source=any&cluster=cl7⟧ and ⟦4‖S12_cov‖n_entries‖set=primary&action_class=SERD&source=any&cluster=cl4⟧).
⟦1170‖S12_sum‖(\d+) activities on CHEMBL206⟧ ChEMBL activities on the ERalpha target were saved for SV15. The labels are provisional; when the
teammate's labels are ready they replace them and SV13-SV15 are re-run. Figure: none.

**SV13 - H12 by action class.** *Request: can strong H12 destabilisation by a ligand class be seen in the structures -
"Good to show the picture of H12 and the ligand and compare between several class types".* Every class comparison
fails the coverage gate (⟦0‖S13_n‖gate_rows_pass⟧ of ⟦22‖S13_n‖gate_rows_primary⟧ comparison rows testable): SERD vs SERM (SERD ⟦4‖S13_g‖entries_a‖set=primary&comparison=C1_SERD_vs_SERM&stratum=all⟧ entries, all in cl4), antagonist vs
SERM (⟦0‖S13_g‖entries_a‖set=primary&comparison=C1_antagonist_vs_SERM&stratum=all⟧ antagonist entries), partial agonist vs agonist (⟦0‖S13_g‖entries_a‖set=primary&comparison=C2a_partial_agonist_vs_agonist⟧ partial-agonist entries) and agonist affinity vs H12
within cl7 (⟦3‖S13_g‖entries_a‖set=primary&comparison=C2b_agonist_rho_dpX⟧ agonist entries with a tested SV15 dpX). No link between action class and H12 disorder or H12
mobility can therefore be stated from these data. Drawn, not tested: complete H12 ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ in agonist ⟦53‖S13_r‖n_complete_h12‖action_class=agonist⟧ of ⟦55‖S13_r‖n_chains‖action_class=agonist⟧
chains, SERM ⟦19‖S13_r‖n_complete_h12‖action_class=SERM⟧ of ⟦53‖S13_r‖n_chains‖action_class=SERM⟧ and SERD ⟦5‖S13_r‖n_complete_h12‖action_class=SERD⟧ of ⟦9‖S13_r‖n_chains‖action_class=SERD⟧ - but agonist entries are ⟦25‖S12_cov‖n_entries‖set=primary&action_class=agonist&source=any&cluster=cl7⟧ of ⟦27‖S13_r‖n_entries‖action_class=agonist⟧ in cl7 and SERM
entries ⟦19‖S12_cov‖n_entries‖set=primary&action_class=SERM&source=any&cluster=cl4⟧ of ⟦20‖S13_r‖n_entries‖action_class=SERM⟧ in cl4, so this restates the clusters. Main figures: the per-class H12 panels
`figures/sv/SV13/SV13_h12_class_panels.png` (one view for all panels; representatives `5WGD|B|1`, agonist, H12
⟦10‖S13_sum‖\| 5WGD\|B\|1 \|.*?\| (\d+) of 10 \|⟧ of ⟦10‖S13_sum‖\| 5WGD\|B\|1 \|.*?\| \d+ of (\d+) \|⟧ residues modelled; `2BJ4|A|1`, SERM, H12 not modelled; `7MSA|C|1`, SERD, ⟦7‖S13_sum‖\| 7MSA\|C\|1 \|.*?\| (\d+) of 10 \|⟧ of ⟦10‖S13_sum‖\| 7MSA\|C\|1 \|.*?\| \d+ of (\d+) \|⟧; drawn
with matplotlib because PyMOL is not available) and the forest plot `figures/sv/SV13/SV13_forest.png`, which shows the
gate reason in place of estimates. The redrawn MDS and H12 maps are backup.

**SV14 - ligand size vs pocket size, and pocket filling.** *Request: "Correlation between ligand size and pocket
size"; normalise pocket and ligand volume to ⟦100‖FB‖normalized to (\d+)%⟧ % and "calculate the percentage of filling".* A buriedness grid
(⟦0.4‖S14_ref‖grid/spacing_A⟧ A spacing; a point counts as buried when at least ⟦10‖S14_ref‖grid/buried_min_rays_primary⟧ of ⟦14‖S14_ref‖grid/rays⟧ rays hit the receptor) was computed for ⟦701‖S14_n‖primary_chains⟧ chains
with a complete ligand (⟦331‖S14_n‖primary_entries⟧ entries, ⟦287‖S14_n‖primary_ligands⟧ ligands). Ligand and pocket volume go together (Spearman rho ⟦0.44‖S14_c‖rho‖subset=all&level=ligand&kind=ligand vs pocket volume⟧
at ligand level, ⟦0.47‖S14_c‖rho‖subset=all&level=entry&kind=ligand vs pocket volume⟧ at entry level), within cl7 (⟦0.41‖S14_c‖rho‖subset=cl7&level=ligand&kind=ligand vs pocket volume⟧ / ⟦0.46‖S14_c‖rho‖subset=cl7&level=entry&kind=ligand vs pocket volume⟧) but not within cl4 (⟦0.11‖S14_c‖rho‖subset=cl4&level=ligand&kind=ligand vs pocket volume⟧ / ⟦0.14‖S14_c‖rho‖subset=cl4&level=entry&kind=ligand vs pocket volume⟧) nor
within agonist or SERM. Filling has a median of ⟦33.8‖S14_n‖filling_primary_median⟧ % and does not differ between cl7 and cl4 (entry medians
⟦33.0‖S14_t‖entry_median_cl7‖comparison=cluster&level=entry⟧ vs ⟦32.3‖S14_t‖entry_median_cl4‖comparison=cluster&level=entry⟧ %, q ⟦0.14‖S14_t‖q_bh‖comparison=cluster&level=entry⟧) or between agonist and SERM (q ⟦0.94‖S14_t‖q_bh‖comparison=action class&level=entry⟧): cl4 ligands are larger (median
⟦363‖S14_fg‖ligand_volume_A3_chain_median‖grouping=cluster&group=cl4⟧ vs ⟦250‖S14_fg‖ligand_volume_A3_chain_median‖grouping=cluster&group=cl7⟧ A^3), but so are their pockets (⟦929‖S14_fg‖pocket_volume_A3_chain_median‖grouping=cluster&group=cl4⟧ vs ⟦768‖S14_fg‖pocket_volume_A3_chain_median‖grouping=cluster&group=cl7⟧ A^3). Reference caveats (first line of the SUMMARY):
only ⟦3‖S14_ref‖apo_baseline/entries_passing⟧ apo entries pass, fewer than ⟦5‖S14_sum‖apo entries pass, < (\d+)\)⟧, so pocket ⟦100‖S14_sum‖pocket (100) % = wild-type⟧ % is the wild-type E2 pocket (⟦956‖S14_ref‖pocket_100pct/median_A3⟧ A^3) and ligand
⟦100‖S14_sum‖ligand (100) % = E2⟧ % the E2 volume (⟦224‖S14_ref‖ligand_100pct/median_A3⟧ A^3); the E2 pocket is bimodal (⟦505‖S14_n‖e2_ref_pocket_values‖‖min⟧-⟦579‖S14_n‖e2_ref_pocket_values/7⟧ vs ⟦956‖S14_n‖e2_ref_pocket_values/8⟧-⟦1027‖S14_n‖e2_ref_pocket_values‖‖max⟧ A^3)
because in ⟦9‖S14_sum‖in (\d+) of 17 chains it joins⟧ of ⟦17‖S14_ref‖ligand_100pct/n_chains⟧ reference chains it joins a buried region past E353/R394; and the pocket reaches the grid box
edge in ⟦183‖S14_n‖touches_box_b10_primary⟧ analysed chains, whose volumes are bounded by the box. Main figure:
`figures/sv/SV14/SV14_pocket_vs_ligand.png`; the filling figure is backup.

**SV15 - affinities across assay systems.** *Request: "Different affinities since of different assays ... you need to
think how to compare based on different systems".* Each ligand's binding pChEMBL is expressed relative to E2 measured
in the same assay (dpX = pX of the ligand minus pX of E2). ⟦14‖S15_n‖tested_ligands_primary⟧ ligands get a same-assay dpX (⟦31‖S15_n‖pairs⟧ ligand x assay pairs
in ⟦19‖S15_n‖paired_assays⟧ assays), ⟦48‖S15_n‖basis_counts_primary_set/absolute⟧ only an absolute value; ⟦9‖S15_n‖tested_single_assay_ligands⟧ of the ⟦14‖S15_n‖tested_ligands_primary⟧ rest on one assay. The E2 reference narrows the
between-assay spread for ⟦4‖S15_n‖spread_n_smaller⟧ of the ⟦5‖S15_n‖spread_ligands⟧ ligands with at least ⟦2‖S15_v‖with >= (\d+) paired assays⟧ paired assays (median SD ⟦0.86‖S15_n‖spread_median_sd_abs⟧ -> ⟦0.56‖S15_n‖spread_median_sd_dpx⟧ log
units). dpX does not differ significantly between cl7 and cl4 ligands (median ⟦-2.17‖S15_t‖median_1‖set=primary&comparison=dpX by SV10 ligand group&level=ligand⟧ vs ⟦-0.23‖S15_t‖median_2‖set=primary&comparison=dpX by SV10 ligand group&level=ligand⟧, q ⟦0.15‖S15_t‖q_bh‖set=primary&comparison=dpX by SV10 ligand group&level=ligand⟧; entries q
⟦0.057‖S15_t‖q_bh‖set=primary&comparison=dpX by SV10 ligand group&level=entry⟧) and does not follow pocket filling (rho ⟦-0.03‖S15_t‖statistic‖set=primary&comparison=dpX vs filling % (SV14)&level=ligand⟧); the class comparison is not testable (⟦3‖S15_n‖class_gate_primary/entries_in_sv15_tested_set/agonist⟧ agonist
entries with a dpX). The two agonists with a dpX are DES (⟦+0.63‖S15_a‖dpx_median‖ccd=DES⟧) and J3Z (⟦-0.78‖S15_a‖dpx_median‖ccd=J3Z⟧); E2 itself is the reference.
Caveat (first line of the SUMMARY): after the first run the pairing was corrected to one dpX per assay, which moved the
entry-level cluster test from q ⟦0.033‖S15_t1‖q_bh‖comparison=dpX by SV10 ligand group&level=entry⟧ to ⟦0.057‖S15_t‖q_bh‖set=primary&comparison=dpX by SV10 ligand group&level=entry⟧. These are in-vitro binding affinities, not cellular activity. Figure
(backup): `figures/sv/SV15/SV15_dpx_by_group.png`.

**SV08b - literature for the class questions.** ⟦19‖S08b_sum‖\*\*(\d+) distinct sources opened\*\*⟧ sources were opened for the three class questions (⟦20‖S08b_sum‖\*\*(\d+) claims\*\*⟧ claims; ⟦12‖S08b_sum‖\*\*(\d+) sources could not be opened\*\*⟧ more could not be opened). First,
in VDR crystal structures the position of H12 is the same for superagonists and agonists (three superagonist families
and the superpotent 2MD; C-alpha RMSD ⟦0.08‖S08b_lit‖Calpha RMSD (0\.08)-0\.38 A against⟧-⟦0.38‖S08b_lit‖Calpha RMSD 0\.08-(0\.38) A against⟧ A against the reference complex), but two sources claim
extra ligand-H12 contacts, so the supported statement is "same H12 position", not "no H12 difference". Second, in
ERalpha, H12 destabilisation by SERDs relative to SERMs is reported mainly in solution (HDX-MS, limited proteolysis),
SERD crystals still show H12 in the AF-2 cleft, and the step from H12 destabilisation to receptor degradation is
asserted rather than shown. Third, the published ERalpha pocket (⟦450‖S08b_lit‖volumes are (\d+) A3 for the pocket⟧ A^3) and E2 (about ⟦250‖S08b_lit‖~(\d+) A3 for E2⟧ A^3) volumes are quoted
second-hand without a method, so the SV14 references cannot be compared like for like; only the lower E2-pocket mode
(⟦505‖S14_n‖e2_ref_pocket_values‖‖min⟧-⟦579‖S14_n‖e2_ref_pocket_values/7⟧ A^3) is near ⟦450‖S08b_lit‖volumes are (\d+) A3 for the pocket⟧ A^3. Table: `docs/sv/literature_class.md`. Figure: none.

### Backup figures

`figures/sv/SV01a/SV01a_completeness_profile.png`, `figures/sv/SV01a/SV01a_forest_completeness.png` (first pass,
without lattice terms), `figures/sv/SV02/SV02_sse_heatmap.png`, the other plots in `figures/sv/SV02/rama/`,
`figures/sv/SV03/SV03_h12_map.png` (the map without the SV04 overlay), `figures/sv/SV03/SV03_descriptors_by_cluster.png`,
`figures/sv/SV05/SV05_core_v4.png`, `figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.png`,
`figures/sv/SV11/SV11_site_share.png`, `figures/sv/SV13/SV13_mds_map.png`, `figures/sv/SV13/SV13_h12_map.png`,
`figures/sv/SV14/SV14_filling.png`, `figures/sv/SV15/SV15_dpx_by_group.png`, `figures/sv/SV15/SV15_dpx_vs_filling.png`.
Every graph is also saved alone in `figures/panels/sv/<ID>/`. SV00, SV01c, SV08, SV08b, SV10 and SV12 made no figure.
