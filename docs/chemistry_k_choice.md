# Choice of K for chemistry clusters

Status: prepared by Claude, pending human sign-off (task C08).

Input: ECFP4 bit vectors (radius 2, 2048 bits) for the 480 distinct ligands in `results/ligands.csv`. K-means with `n_init=20` and seed 42. Scores are in `results/kmeans_scores.csv`. Labels are in `results/ligand_clusters.csv`.

Silhouette by K (Euclidean, binary features):

| K | Silhouette |
|---|---|
| 3 | 0.103 |
| 4 | 0.112 |
| 5 | 0.120 |
| 6 | 0.126 |
| 7 | 0.122 |

Silhouette is highest at K=6, but every value is between 0.10 and 0.13. That indicates weak cluster structure at all K tested.

K=6 cluster sizes: 21, 105, 235, 53, 44, 22.

Chemical sense check, K=6:

- Cluster 3 (235 ligands, 49 %) contains crystallization additives and ions (MG, CL, GOL) alongside estradiol (EST). It is not a chemical class.
- Cluster 1 holds tetrahydronaphthalene and isoquinoline phenols with a basic aminoethoxy or azetidinyl side chain (SERM-like scaffolds).
- Cluster 2 mixes small molecules and additives (for example BME, SAH) with other compounds.
- Cluster 4 holds 7-oxabicyclo[2.2.1]heptene bis(4-hydroxyphenyl) sulfonates and sulfonamides.
- Cluster 5 holds acetylpiperidine-containing compounds.
- Cluster 6 holds carbamimidoylthiophene compounds.

Recommendation for the human: the silhouette favours K=6, but the margin is small and cluster 3 is driven by non-ligand species. Decide first whether crystallization additives (ions, glycerol, ethylene glycol and similar) belong in the chemistry set. If they are removed, rerun C07 and then choose K again. Until then, K=6 is the provisional value, not a chosen one.


## Auto-selection (final run, task Y08, 2026-10-25)

Status: **auto-selected**, still needs the human sign-off that C08 requires.

- Set: the 306 ligands of interest (C11), i.e. the ligand in the pocket of each matrix chain. Crystallisation
  additives are not in this set, so the additive question above no longer affects the choice.
- Zoo: 6 fingerprints (ECFP4, ECFP6, FCFP4, MACCS, atom pair, topological torsion) x 7 methods (Butina, PAM,
  spectral, average linkage, k-means on MDS, HDBSCAN, DBSCAN) and their settings = 276 runs
  (results/ext/chem_zoo_runs.csv). Each run: silhouette on its own Tanimoto distance, silhouette on ECFP4 distance,
  bootstrap stability (30 x 80%, ARI).
- Rule declared before running (docs/PLAN_EXTENDED.md): eligible = >= 2 clusters and noise <= 10%; score = mean of the
  percentile ranks of silhouette and stability.
- The rule's winner is `MACCS|hier_average|k=2`, which is degenerate: one cluster holds 99% of the ligands.
  A non-triviality constraint (largest cluster <= 60%) was added after seeing this, and is reported as post hoc.
- **Auto-selected: `MACCS|kmeans_mds|k=6`, K = 6** (silhouette 0.291, stability
  0.975). Labels: results/ext/chem_labels_best.csv, column chem_cluster_auto.
- Cluster sizes: 0: 66, 1: 42, 2: 68, 3: 37, 4: 55, 5: 38.
- Silhouettes are modest (best about 0.29), as in the earlier ECFP4 K-means (0.10-0.13): chemistry clusters are real
  but not sharply separated.
- Structural link (Y25 tests it formally): chemistry cluster (rows) x structural majority of the ligand's chains:

|   chem_cluster_auto |   cl4 |   cl7 |   minor |   mixed cl7/cl4 |   none |
|--------------------:|------:|------:|--------:|----------------:|-------:|
|                   0 |    20 |    43 |       0 |               1 |      2 |
|                   1 |     6 |    35 |       0 |               1 |      0 |
|                   2 |    63 |     3 |       2 |               0 |      0 |
|                   3 |     1 |    36 |       0 |               0 |      0 |
|                   4 |     2 |    52 |       0 |               1 |      0 |
|                   5 |    12 |    23 |       0 |               3 |      0 |

What the human must do: accept K = 6 (MACCS, k-means on MDS) or choose another run from
results/ext/chem_zoo_ranking.csv, then log C08 DONE --who human.
