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
