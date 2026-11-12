# ESR1 ligand-binding domain: how diverse are the crystal structures?

We collect the deposited X-ray structures of the oestrogen receptor α ligand-binding domain, group the chains by shape,
find what drives the groups and how much of the set repeats itself. "Lid" means helix 12 (H12). Numbers come from
`results/<stage>/numbers.json` and `evidence/INDEX.md` (claim IDs in brackets); `bash reproduce.sh` rebuilds them.

## 1. September baseline

- **Added.** 487 ESR1 and 39 ESR2 entries downloaded; 444 ESR1 entries hold a receptor chain (901 chain-models); 788
  chains from 351 entries pass the completeness screen (≥ 180 LBD, ≥ 60 core C-alpha); 247 of 351 entries pass
  resolution ≤ 2.5 Å, R-free ≤ 0.30 and mean B ≤ 50 Å².
- **Method.** UniProt numbering; one chain per entry superposed by Kabsch on an 84-residue core onto 7NFB chain B
  (1.33 Å, engineered variant); C-alpha RMSD over shared positions; average linkage: 9 groups at 1.5 Å. Tools: Python
  (gemmi, Biopython, NumPy, SciPy, scikit-learn), RDKit, PyMOL.
- **Professor's view.** Keep every eligible chain and treat quality as a label; try many clustering methods; report
  compactness; test whether outliers come from quality or conditions; colour maps by quality and solvent; explain why
  there are two groups.

## 2. October extension

- **Changed.** All 788 chains instead of 351 (quality tiers A 544, B 240, C 4) and a core of the five most stable helices
  by normalised B (123 positions); groups at 1.5 Å went from 9 to 7.
- **Two groups.** cl7 (445 chains; lid closed over the pocket) and cl4 (329; lid in the coactivator groove), plus 14
  chains in five minor groups. Median RMSD 0.67 Å within cl7, 0.87 Å within cl4, 3.49 Å between; residues 528-550 carry
  98.6 % of the squared difference, and without them the groups are 0.82 Å apart. Other clustering methods and five 2-D
  maps find the same split.
- **What goes with it.** Coactivator peptide in 414 of 445 cl7 vs 4 of 329 cl4 chains; Y537S in 405 cl7 chains, L536S in
  237 cl4 and no cl7 chains; cl4 ligands are larger (median 34 vs 24 heavy atoms) and carry a basic amine (74 % vs 1 %).
  P2Rank ranks the ligand pocket first in all 788 chains.
- **Outliers and redundancy.** 22 chains flagged by ≥ 3 of 5 detectors, 6 with a quality reason. Effective number of
  conformations 2.19 (Hill q1); 2 complete-H12 structures cover 90 % of chains within 1.0 Å (48 at 0.5 Å); copies of
  one crystal differ by a median 0.53 Å; no chain deposited since 2010 is > 1.5 Å from all earlier ones.
- **ESR2.** 86 chains split the same way: the main group (93 % with peptide) is 1.33 Å from ESR1 cl7, the second 2.31 Å
  from cl4, while the ESR1 groups are 3.49 Å apart.
- **Statistics and tools.** Entry-level tests beside chain-level ones, Benjamini-Hochberg, seed 42; UMAP, P2Rank, ChEMBL.
- **Professor's view (28 Oct).** Does H12 completeness follow conditions or biology; geometric descriptors; consensus
  secondary structure; a stable core; colour outliers; audit the profile spikes; literature on changes outside H12;
  helix midpoints for incomplete helices (Ye).

## 3. Requests of 28 Oct and 4 Nov

- **Lid order.** H12 538-547 is complete in 431 of 445 cl7 but 120 of 329 cl4 chains. Consensus secondary structure
  (pydssp) puts H12 at 538-544, where an incomplete lid becomes rare (cl4 209 → 19) [W1]. Per SD of symmetry-mate
  contacts the odds of a complete H12 rise 20.2-fold in cl4 (logistic models, entry-clustered errors) [C3].
- **Geometry.** Residue 537 switches: psi +156.6° in cl7 vs −46.8° in cl4, no overlap [W2]. The lid sits 12.48 vs
  22.26 Å from the H11 axis [L1]; H12-H3 angle 79.0° vs 46.6°. Helix midpoints move ≤ 0.91 Å with 4 end residues
  missing and reproduce the groups (ARI 0.954) [W3].
- **Outliers and spikes.** Outliers are rarer in cl7 than expected (2 vs 8.0 of 14, q 0.0023) [O1]; outside the lid only
  loops 331-340 and 464-468 move (≤ 3.6 Å), as the literature also reports (`docs/literature.md`).
- **Ligands (4 Nov).** Uncurated ChEMBL mechanism classes label 51 of 362 entry-ligand rows (agonist 27, SERM 20,
  SERD 4); agonists are 25 of 27 entries in cl7, SERMs 19 of 20 in cl4 [L7]. Pocket filling: median 33.8 % over 701
  chains, no group difference (q 0.14) [P2]. Referencing affinity to E2 in the same assay narrows the between-assay spread
  for 4 of 5 ligands (SD 0.86 → 0.56) [P3]. Size explains most drug-likeness differences; 7 quantities go beyond it [P1].
  In VDR, superagonists keep the agonist H12 position; SERD destabilisation of H12 is reported in solution.
- **Professor's view (11 Nov).** Look into cl4, the two cl7 chains and the grey outliers of the lid-to-H11 distance; try
  Sammon and self-organising maps and UMAP of PCA; hide H12; ask whether many crystallisation components make the
  receptor look stable; use phi/psi and centres of mass without superposition; separate ordinary variability from real
  effects (Dana).

## 4. Analyses of 11 Nov

- **Rule.** An effect counts only above ordinary variation: the larger of the within-group SD and the 95th percentile of
  differences between copies of one crystal (lid distance: 0.91 Å in cl7, 1.08 Å in cl4) [L2].
- **Lid outliers.** 15 chains exceed it [L3]: the two cl7 chains 5KCU B and 5KD9 B sit 1.5-1.6 Å further out; cl4 has no
  Tukey outlier, only three extremes; five grey chains hold the lid between the two positions [L4]. In 7 of 15 crystals
  the other copies are ordinary, so the unusual lid belongs to one copy [L5]. No ligand, construct, peptide or crystal
  property moves the lid inside a group beyond ordinary variation [L6].
- **Maps.** Sammon, a 12 × 12 self-organising map and UMAP on 6 PCA components show the same two groups; Sammon keeps
  distances best (Spearman 0.92) [M1, M2]; extra islands follow crystal form and stay within copy variation [M3].
- **Lid hidden.** Without 528-550 the chains form one ensemble (best split 784 / 4); the rest is loop 331-340, where
  copies of one crystal differ as much [N1, N2].
- **No superposition.** From phi/psi helix ranges and centres of mass in the deposited coordinates, only the 8 distances
  involving H12 change beyond ordinary variation (H7-H12 +5.45 Å); all others ≤ 0.20 Å [G1]. These features recover the
  groups (ARI 0.94) only with H12 [G2]; helix ends shift 3-5 residues, within copy variation [G3]; in cl4 a longer H12 is
  less mobile (ρ −0.37) [G4].
- **Crystal environment.** Additives, recipe and solvent content stay inside ordinary variation [C1]; within one crystal
  the cl4 copy with more lattice contacts at the lid has more of it modelled (38 vs 14 crystals) [C2]. Crystals cannot
  show dynamics in water; they show that the local packing, not the drop composition, goes with lid order.

## Open work

- Replace the uncurated classes with the teammate's hand labels and test class against lid position within each group.
- Sign off the H12 window and repeat the secondary structure with DSSP.
- Link structure to affinity or activity (asked on 4 Nov, not attempted) and compare with solution data (HDX, MD).
- Find the residue-level effect of L536S (243 chains) and an apo wild-type baseline, both absent from the literature read.
