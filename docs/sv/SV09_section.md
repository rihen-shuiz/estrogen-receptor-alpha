## Supervisor requests 2026-10-28

Answers to the requests of the meeting on 2026-10-28 and to the notes received 2026-11-04,
produced by the unattended run `docs/sv/RUN_SV1028.md` (status `docs/sv/RUN_STATUS.md`). One paragraph and one figure
per step, written from each step's `results/sv/<ID>/SUMMARY.md` and verdict files; every number in this section is
matched to the file it comes from by `scripts/sv/SV09_number_check.py` (`results/sv/SV09/number_check.csv`). Fixed
choices: the H12 window stays 538-547 for all primary analyses, the superposition core stays v3, and cl7 / cl4
are the S08 clusters at 1.5 A. Crystal structures show H12 order and position, so results are worded as H12 disorder
(H12 not modelled) or H12 mobility (local B), never as receptor activity or degradation. SV00-SV11 describe ligands by
chemistry only; the action classes of SV12-SV15 are provisional ChEMBL-mechanism labels that the teammate's labels
will replace. New limitations are in `docs/limitations.md`.

### Main figures for the next meeting

| figure | step | what it shows |
|---|---|---|
| `figures/sv/SV01b/SV01b_forest_completeness.png`, `figures/sv/SV01b/SV01b_forest_mobility.png` | SV01b (the SV01 forest plot, redrawn with lattice terms) | odds ratios (complete H12) and coefficients (H12 mobility) per predictor, coloured by block Q / X / B, within cl7 and cl4; chain level coloured, entry level grey; each block's gain entered first and last |
| `figures/sv/SV02/SV02_helix_track.png` | SV02 | pydssp helix fraction per residue for all chains, cl7 and cl4, with the consensus helices and the historical H12 windows |
| `figures/sv/SV02/rama/SV02_rama_537.png` | SV02 | phi / psi of residue 537 by cluster (named a main figure in the run file) |
| `figures/sv/SV04/h12_map_partial.png` | SV03 with the SV04 overlay | the SV03 H12 map (H12 distance to the H11 axis vs H12-H3 angle) with the partial-H12 chains added, hollow, at their SV04 midpoint estimates |
| `figures/sv/SV06/SV06_outliers_map.png` | SV06 | consensus outliers on the shared MDS map (A) and their RMSD to their own vs the other complete-H12 medoid (B) |
| `figures/sv/SV07/SV07_spike_profile.png` | SV07 | the Y10 cl7-cl4 displacement profile with every spike labelled by its audit verdict, and the share of chains observing each position |
| `figures/sv/SV13/SV13_h12_class_panels.png` | SV13 | H12 and the ligand of one representative per provisional action class, with the H12 traces of all chains of the class, in one saved view |
| `figures/sv/SV13/SV13_forest.png` | SV13 | the planned class comparisons; every row shows its coverage-gate reason in place of an estimate |
| `figures/sv/SV14/SV14_pocket_vs_ligand.png` | SV14 | pocket vs ligand volume per chain as % of the E2 references, cluster as fill, action class as marker |

### Results per step

**SV00 - data restore and shared chain table.** All 526 mmCIF files of the frozen entry list were restored,
526 of them with the frozen sha256, and the reproduction gate passed (ARI 1.0 at all four cuts, max |dRMSD|
9.2e-15 A, 0 differing chains), so every later step works on the same 788 chains from 351 entries. The shared
chain table (`results/sv/SV00/chain_table.csv`, 53 columns) adds dimer partner, bound peptide and SIFTS
offsets per chain; crystallisation pH is missing for 285 chains, because only the mmCIF pH field and the RCSB Data API were read.
Figure: none.

**SV01a - H12 completeness, first pass.** *Request: "check distribution of completeness, and if completeness correlates
with any experimental conditions/parameters, or is it biology".* H12 538-547 is complete in 431 of 445 cl7 chains
(0.969) but in 120 of 329 cl4 chains (0.365); 207 of the 212 partial H12 lack only
C-terminal residues. In multi-chain entries the H12 label differs between chains in 59 of 317 entries, the cluster
in only 5. Logistic models within each cluster, with blocks Q (quality), X (experiment) and B (biology: peptide,
construct mutations, ligand chemistry), give the [DEFAULT] verdict **mixed**: entry level is artifact (Q+X vs B gain
cl7 0.613 vs 0.315, cl4 0.346 vs 0.117), but at chain level cl7 is a near tie (Q+X 0.224 vs B 0.225; cl4
0.208 vs 0.079). cl7 has only 14 incomplete chains for 16 predictors. Figure (backup):
`figures/sv/SV01a/SV01a_completeness_profile.png`.

**SV01b - with coordinates: lattice contacts and H12 mobility.** Two additions: heavy-atom lattice contacts of residues
528-550 with ASU partners and symmetry mates, added to block X, and a second outcome, local H12 mobility (C-alpha B
z-score of the observed residues of 528-547 against the same chain's core v3). The re-issued verdict is
**artifact** for both outcomes at chain and entry level: H12 complete, chain level, Q+X vs B gain cl7 0.414 vs 0.225 and
cl4 0.523 vs 0.079; H12 mobility cl7 0.184 vs 0.076 and cl4 0.210 vs 0.143; the
sensitivity runs (fraction observed, entries without TLS refinement) agree. The change from SV01a is the lattice
terms' share: per SD of symmetry-mate contact atoms the odds of a complete H12 rise 20.2-fold in cl4 (12.6-fold in cl7)
and H12 mobility falls by 0.268 z in cl4 (0.17 in cl7). Caveat (first line and Reading of the SUMMARY): for completeness
this share is partly mechanical, because contacts are counted on modelled atoms. In cl4, complete chains carry
24.3 symmetry-contact atoms on 528-550 on average and partial chains 9.1, but 12.3 of the former sit on 546-550, which
partial chains hardly model; on 538-545, modelled by both, the means are 8.22 vs 6.66. The data cannot tell whether a
contact holds the C-terminal end in place or exists because the end is modelled; the mobility outcome is not tied to
missing atoms in that way. The medoid clash test was not run (it runs only if the verdict is mixed). Main figure: the
SV01 forest plots `figures/sv/SV01b/SV01b_forest_completeness.png` and `figures/sv/SV01b/SV01b_forest_mobility.png`.

**SV01c - sensitivity to the H12 window.** On the data-derived window 538-544 (SV02 proposal), 204 chains that are partial
on 538-547 become complete, so an incomplete H12 becomes a rare event (cl7 14 -> 10 chains, cl4 209 -> 19),
mostly H12 not modelled at all. The completeness verdict moves from artifact to **mixed** at chain level only (cl7 Q+X
0.130 vs B 0.386, cl4 0.506 vs 0.216); entry level stays artifact but rests on 4 (cl7) and 6 (cl4)
incomplete entries, and H12 mobility stays artifact. Because the verdict was mixed, the medoid clash test ran: of 31
partial or not-modelled chains on 538-544, 3 have both medoid H12 positions blocked by lattice neighbours, 0 only the
cl7 and 0 only the cl4 position, 28 neither, so the lattice does not favour one H12 position. The
completeness verdict therefore depends on the H12 window, which the team has yet to sign off. Figure: none (tables
`results/sv/SV01c/side_by_side.csv` and `results/sv/SV01c/clash_summary.csv`).

**SV02 - consensus secondary structure.** *Request: "consensus alignment of secondary structures".* pydssp on all
788 chains (fallback for mkdssp, which is not installed; its 'H' covers every helix type) gives 8 global consensus
helices (at least 4 consecutive positions helical in at least 80 % of the observing chains). The data-derived H12
window is 538-544 (cl7 538-545, cl4 537-544), a proposal only: it matches none of the
historical windows (Jaccard 0.70 with each of 538-547, 537-546 and 536-545), and 538-547 stays
primary. The clusters also differ in where other helices end: H10/11 ends at 530 in cl7 and 525 in cl4, H9
starts at 473 in cl7 and 466 in cl4. Residue 537 alone separates the clusters: circular-mean psi
+156.6 deg in cl7 vs -46.8 deg in cl4 (separation 157.2 deg; overlap 0 at chain and 0 at entry
level), as recorded earlier (+157 / -47). Caveat: pydssp has no chain-break handling, so beta bridges
across breaks are missed (helices are not affected). Main figures: `figures/sv/SV02/SV02_helix_track.png` and
`figures/sv/SV02/rama/SV02_rama_537.png`; the SSE heatmap and the other Ramachandran plots are backup.

**SV03 - helix geometry.** *Request: "we can try geometrical approach since it's more sensitive for things like the
length of the helix".* 16 descriptors per chain (helix lengths, inter-helix angles, H12 distances to the H3 and
H11 axes, H11 kink, H11-H12 hinge) in the core v3 frame. cl7 vs cl4: 12 of 16 differ at q < 0.05 at chain and
entry level, all 7 H12 descriptors among them; the largest are the H12-H3 angle (median 79.05 vs 46.63 deg)
and the distance from the H12 centroid to the H11 axis (12.48 vs 22.26 A). Within clusters, 37 covariate
comparisons meet the [DEFAULT] rule (|r| >= 0.3 and q < 0.05 at both levels; 19 in cl7, 18 in cl4), so by that
rule the descriptors are more sensitive than RMSD for those subgroups. The effects are small (angles 0.5-5.4 deg,
lengths 1-4 residues) and 4 of the 37 rest on a group of at most 5 entries. Caveat: helix ends
come from the pydssp fallback, and cl4 H12 axes rest mostly on 7 C-alpha, where the axis fit on an ideal helix
is off by 4.5 deg, so H12 angle differences of a few degrees inside a cluster are within the method's range.
Figure (backup): `figures/sv/SV03/SV03_descriptors_by_cluster.png`; the SV03 H12 map is shown with the SV04 overlay.

**SV04 - helix midpoints for incomplete helices.** *Request (Ye): compare incomplete helices by "taking the middle
atoms".* Each helix is reduced to the midpoint of its fitted axis, validated by deleting residues from complete
helices; the H12 curve rests on the 75 chains whose whole 538-547 is helical. Clustering chains on midpoints
reproduces S08: ARI 0.954 at K = 2 (774 cl7 and cl4 chains) and 0.931 at K = 7 (788 chains);
entry level 0.977 / 0.959. This mostly restates the H12 split - the cl7 and cl4 H12 midpoints lie a median
12.58 A apart, every other helix at most 0.30 A - and is not independent evidence, since S08 used the
same observed H12 C-alpha. Of 212 partial-H12 chains, 199 have a usable H12 midpoint, and every one of the 191
in cl7 or cl4 lands on its own S08 side. Main figure (with SV03): `figures/sv/SV04/h12_map_partial.png`; panel (a)
places complete chains at SV03 positions and partial chains at SV04 estimates, panel (b) uses the SV04 estimator for all
chains, because midpoint and centroid differ by a cluster-dependent offset in x (complete cl4 +1.51 A, cl7 -0.68 A).

**SV05 - stable core v4 from complete helices.** *Request: "take helices without the problem of completeness, find a
stable core".* From 132 candidate positions (SV02 consensus helices outside 528-550, observed in at least 98 % of
chains), iterative superposition keeps 128 positions with spread at most 0.5 A as core v4: it adds H1 and H9,
drops 25 v3 positions at helix ends, and has Jaccard 0.64 with v3. Recomputed with v4, pairwise RMSDs barely
change (Spearman 0.9975 with the primary matrix; ARI with S08 0.990 at 1.5 A, 3 chains change cluster), and
98.5 % of the squared cl7-cl4 difference stays in 528-550 (v3 98.6 %). Whether to adopt v4 is a team decision;
SV03, SV04 and SV07 used v3. Figure (backup): `figures/sv/SV05/SV05_core_v4.png`.

**SV06 - outliers on the clustering map.** *Request: "color outliers in clustering (but no pattern is also
possible)".* Of the 22 consensus outliers (Y15), 2 are in cl7 (445 chains), 12 in cl4 (329) and 8 in the minor
clusters (14 chains). Between the main clusters this is a pattern (hypergeometric, BH over 8 tests): 2 in cl7
against 8.0 expected (q 0.0023), and at entry level 2 against 7.0 (q 0.0034). The minor-cluster
enrichment is expected by construction, because the detectors measure distance from the main clusters. None
(0 of 14) of the main-cluster outliers is closer to the other cluster's medoid than to its own (own/other
RMSD ratio 0.21-0.63), so they sit at the edge of their own cluster rather than between the two. Caveat:
panel B uses the complete-H12 medoids `5TM1|A|1` and `7UJ7|C|1` (3.68 A apart), because the S12 cl4 medoid has no
modelled H12. Main figure: `figures/sv/SV06/SV06_outliers_map.png`.

**SV07 - spike audit of the displacement profile.** *Request: "check spikes on figure in slide 5".* 19 positions
outside 528-550 exceed 1.0 A in the Y10 profile (cl7-cl4 C-alpha distance of the cluster means or of the medoid
pair). 14 are genuine loop differences (331-340, 464, 466-468): the distance stays above 1.0 A
after removing the 3 most influential entries, at chain and entry level. 5 are driven by at most 3 entries:
330, 397, 420 and 459 exceed only on the medoid pair, and 462 drops to 0.96 A without one entry. No
spike is a numbering error. Among the genuine spikes, ligand chemotype is associated (q < 0.05 at both levels) at
331-332 and 334-339, dimer-partner source at 334-340 and 467-468, and the bound peptide at
331 and 466-468; at 335 the cl4-cl7 coefficient falls from 2.74 to 0.51 A once chemotype is in
the model, but chemotype is unevenly spread over the clusters, so the adjusted value rests on few chains. Partly
BLOCKED: the 'sparse-sample effect' rule is degenerate, because Y10 already averages only the chains that observe a
position (restricted means equal Y10 within 2e-14 A). Main figure: `figures/sv/SV07/SV07_spike_profile.png`.

**SV08 - literature: is only H12 different?** *Request: "check literature if only H12 is significant? so far seems
so".* From 14 opened sources (25 claims; 8 sources could not be opened) the answer is no: they report differences
between ERalpha LBD structures in H11 (shifting with ligand chemistry), H3, H7, the H11-H12 loop and the
His524 / Leu525 rotamers, and one deposited structure leaves residues 331-340 and 462-463 unmodelled
(Li et al. 2011) - the two loop regions where the genuine SV07 spikes lie. Y537S, D538G and L536S are documented as
conformation-driving construct mutations, but no opened source gives the residue-level effect of L536S (carried by
243 of our chains), and none attributes a difference outside H12 to the coactivator peptide. Table:
`docs/sv/literature.md`. Figure: none.

**SV10 - ligand drug-likeness.** Of 21 drug-likeness parameters of the crystallised ligands, 10 add information
beyond ligand size (7 distinct quantities: H-bond donors and acceptors, rotatable bonds, TPSA, QED alerts, Brenk and
PAINS alerts) and 11 restate size (heavy atoms median cl7 24 vs cl4 34); table
`results/sv/SV10/druglikeness_by_group.md`.

**SV11 - pocket sites and steroid flag.** Part A: P2Rank, which saw only the receptor chain and its ligand (no peptide,
no partner chain), ranks the ligand pocket first in 788 of 788 chains, with no cl7-cl4 difference there. A
pocket in the AF-2 groove is found in 0.53 of cl7 chains (234 of 445) but 0.02 of cl4 chains (5 of 329;
entries 0.72 vs 0.01), and within cl7 only in chains with a bound peptide (234 of 414, 0 of
31 without); a dimer-interface pocket in 0.18 vs 0.31 of chains (q 0.0034 chain level, 0.022 entry level).
Both are low-probability P2Rank pockets (median probability 0.010 and 0.013 vs 0.924 for the ligand pocket),
and the groove rank contrasts rest on 5 cl4 chains. Part B: a steroid ring system is present in 23 of 306
ligands, more often in cl7 entries (41 of 217) than in cl4 entries (5 of 124; q 0.00011), and not
significantly at ligand level (q 0.071). Figure (backup): `figures/sv/SV11/SV11_site_share.png`.

#### Addition: notes received 2026-11-04

**SV12 - ligand identifiers and provisional action classes.** *Request: "6 different action class types and based
on their DOIs, show in the plot"; "get the SMILES, InChiKey from api and using chembl api".* SMILES and InChIKeys came
from RCSB for all 306 ligand codes and were matched to ChEMBL. The teammate's label file did not exist at run time, so
classes come from ChEMBL mechanism records only; the DOI route failed its validation gate (1 of 2 rows agree)
and is not used. 51 of 362 entry x ligand rows are labelled - agonist 27, SERM 20, SERD 4; partial agonist
0, antagonist 0, other 0 - and 311 stay unlabelled. Only agonist (27
entries) and SERM (20) pass the coverage gate of 5 entries per side, and the classes nearly coincide
with the clusters (agonist entries 25 in cl7 and 1 in cl4; SERM 3 and 19; SERD 0 and 4).
1170 ChEMBL activities on the ERalpha target were saved for SV15. The labels are provisional; when the
teammate's labels are ready they replace them and SV13-SV15 are re-run. Figure: none.

**SV13 - H12 by action class.** *Request: can strong H12 destabilisation by a ligand class be seen in the structures -
"Good to show the picture of H12 and the ligand and compare between several class types".* Every class comparison
fails the coverage gate (0 of 22 comparison rows testable): SERD vs SERM (SERD 4 entries, all in cl4), antagonist vs
SERM (0 antagonist entries), partial agonist vs agonist (0 partial-agonist entries) and agonist affinity vs H12
within cl7 (3 agonist entries with a tested SV15 dpX). No link between action class and H12 disorder or H12
mobility can therefore be stated from these data. Drawn, not tested: complete H12 538-547 in agonist 53 of 55
chains, SERM 19 of 53 and SERD 5 of 9 - but agonist entries are 25 of 27 in cl7 and SERM
entries 19 of 20 in cl4, so this restates the clusters. Main figures: the per-class H12 panels
`figures/sv/SV13/SV13_h12_class_panels.png` (one view for all panels; representatives `5WGD|B|1`, agonist, H12
10 of 10 residues modelled; `2BJ4|A|1`, SERM, H12 not modelled; `7MSA|C|1`, SERD, 7 of 10; drawn
with matplotlib because PyMOL is not available) and the forest plot `figures/sv/SV13/SV13_forest.png`, which shows the
gate reason in place of estimates. The redrawn MDS and H12 maps are backup.

**SV14 - ligand size vs pocket size, and pocket filling.** *Request: "Correlation between ligand size and pocket
size"; normalise pocket and ligand volume to 100 % and "calculate the percentage of filling".* A buriedness grid
(0.4 A spacing; a point counts as buried when at least 10 of 14 rays hit the receptor) was computed for 701 chains
with a complete ligand (331 entries, 287 ligands). Ligand and pocket volume go together (Spearman rho 0.44
at ligand level, 0.47 at entry level), within cl7 (0.41 / 0.46) but not within cl4 (0.11 / 0.14) nor
within agonist or SERM. Filling has a median of 33.8 % and does not differ between cl7 and cl4 (entry medians
33.0 vs 32.3 %, q 0.14) or between agonist and SERM (q 0.94): cl4 ligands are larger (median
363 vs 250 A^3), but so are their pockets (929 vs 768 A^3). Reference caveats (first line of the SUMMARY):
only 3 apo entries pass, fewer than 5, so pocket 100 % is the wild-type E2 pocket (956 A^3) and ligand
100 % the E2 volume (224 A^3); the E2 pocket is bimodal (505-579 vs 956-1027 A^3)
because in 9 of 17 reference chains it joins a buried region past E353/R394; and the pocket reaches the grid box
edge in 183 analysed chains, whose volumes are bounded by the box. Main figure:
`figures/sv/SV14/SV14_pocket_vs_ligand.png`; the filling figure is backup.

**SV15 - affinities across assay systems.** *Request: "Different affinities since of different assays ... you need to
think how to compare based on different systems".* Each ligand's binding pChEMBL is expressed relative to E2 measured
in the same assay (dpX = pX of the ligand minus pX of E2). 14 ligands get a same-assay dpX (31 ligand x assay pairs
in 19 assays), 48 only an absolute value; 9 of the 14 rest on one assay. The E2 reference narrows the
between-assay spread for 4 of the 5 ligands with at least 2 paired assays (median SD 0.86 -> 0.56 log
units). dpX does not differ significantly between cl7 and cl4 ligands (median -2.17 vs -0.23, q 0.15; entries q
0.057) and does not follow pocket filling (rho -0.03); the class comparison is not testable (3 agonist
entries with a dpX). The two agonists with a dpX are DES (+0.63) and J3Z (-0.78); E2 itself is the reference.
Caveat (first line of the SUMMARY): after the first run the pairing was corrected to one dpX per assay, which moved the
entry-level cluster test from q 0.033 to 0.057. These are in-vitro binding affinities, not cellular activity. Figure
(backup): `figures/sv/SV15/SV15_dpx_by_group.png`.

**SV08b - literature for the class questions.** 19 sources were opened for the three class questions (20 claims; 12 more could not be opened). First,
in VDR crystal structures the position of H12 is the same for superagonists and agonists (three superagonist families
and the superpotent 2MD; C-alpha RMSD 0.08-0.38 A against the reference complex), but two sources claim
extra ligand-H12 contacts, so the supported statement is "same H12 position", not "no H12 difference". Second, in
ERalpha, H12 destabilisation by SERDs relative to SERMs is reported mainly in solution (HDX-MS, limited proteolysis),
SERD crystals still show H12 in the AF-2 cleft, and the step from H12 destabilisation to receptor degradation is
asserted rather than shown. Third, the published ERalpha pocket (450 A^3) and E2 (about 250 A^3) volumes are quoted
second-hand without a method, so the SV14 references cannot be compared like for like; only the lower E2-pocket mode
(505-579 A^3) is near 450 A^3. Table: `docs/sv/literature_class.md`. Figure: none.

### Backup figures

`figures/sv/SV01a/SV01a_completeness_profile.png`, `figures/sv/SV01a/SV01a_forest_completeness.png` (first pass,
without lattice terms), `figures/sv/SV02/SV02_sse_heatmap.png`, the other plots in `figures/sv/SV02/rama/`,
`figures/sv/SV03/SV03_h12_map.png` (the map without the SV04 overlay), `figures/sv/SV03/SV03_descriptors_by_cluster.png`,
`figures/sv/SV05/SV05_core_v4.png`, `figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.png`,
`figures/sv/SV11/SV11_site_share.png`, `figures/sv/SV13/SV13_mds_map.png`, `figures/sv/SV13/SV13_h12_map.png`,
`figures/sv/SV14/SV14_filling.png`, `figures/sv/SV15/SV15_dpx_by_group.png`, `figures/sv/SV15/SV15_dpx_vs_filling.png`.
Every graph is also saved alone in `figures/panels/sv/<ID>/`. SV00, SV01c, SV08, SV08b, SV10 and SV12 made no figure.
