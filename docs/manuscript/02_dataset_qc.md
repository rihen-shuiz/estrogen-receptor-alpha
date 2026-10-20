# Dataset and quality control

## Methods

### Structure collection

The RCSB search API was queried for all X-ray entries that contain UniProt P03372 (ESR1) or Q92731 (ESR2) (data/pdb_ids_esr1.txt, data/pdb_ids_esr2.txt). Every entry was downloaded as mmCIF. Each file was validated by parsing and by checking the deposited entry identifier, and its SHA-256 hash was recorded (data/raw/download_manifest.csv, data/raw/sha256sums.txt).

### Receptor chains and the ligand-binding domain

Each polymer chain was classified as receptor, peptide, partner protein or DNA from its deposited accession and its sequence. Receptor chains were mapped residue by residue onto canonical UniProt numbering, and ESR2 chains also onto ESR1 numbering (results/lbd_extract/residue_mapping.csv). The LBD is ESR1 residues 305-550. A chain enters the structural analysis when it has at least 180 mapped LBD C-alpha atoms and at least 60 C-alpha atoms in the alignment core. This completeness screen is unchanged from the original project (results/lbd_chains.csv).

### Entry metadata and QC

For each entry we recorded resolution, R-free, mean B-factor, space group, deposition date, depositing authors and ligand codes (results/entries.csv). Three flags were set: resolution above 2.5 Å, R-free above 0.30 and mean B above 50 Å² fail. A missing value is unknown, never a pass. For the analysis, each entry also receives a graded label, qc_label, from three tiers (results/qc_tiers.csv; limits in esr_analysis/resources/config/qc_tiers.json):

- Tier A uses the flag cutoffs above.
- Tier B allows resolution up to 3.5 Å, R-free up to 0.40 and mean B up to 80 Å².
- Tier C is anything beyond tier B.
- An entry is labelled unknown when a metric is missing.

No entry or chain is removed for quality. Labels are carried into every downstream table, and the effect of restricting to better tiers is reported as a sensitivity analysis.

### Ligands

Every distinct non-water chemical component was tabulated with its SMILES from the PDB Chemical Component Dictionary (results/ligands.csv) and ranked by the number of entries that carry it (results/ligands_ranked.csv). Ions, buffers, solvents and cryoprotectants are marked is_additive and kept. For each eligible ESR1 chain, the ligand of interest is the deposited "subject of investigation" when one is declared. Otherwise it is any non-polymer within 4.5 Å of at least three pocket residues, and additives are never counted (results/chain_ligand_of_interest.csv).

### Redundancy and temporal bias

Entries that share the same ligand set, space group and depositing authors were flagged as near-duplicates (results/near_duplicates.csv). Deposits per year were counted from the deposition date (figures/deposits_per_year.png).

### QC questions

Two questions test whether quality drives the structural results.

1. Do structures in rare clusters (fewer than 10 chains at a cut) have worse resolution or R-free? Entries in rare clusters were compared with the rest by a two-sided Mann-Whitney U test, per entry, at each cut (results/qc_outliers_stats.csv).
2. Do high-B structures concentrate in particular clusters? Chains were split into tertiles of mean C-alpha B, and the tertiles were cross-tabulated against the 1.5 Å clusters. The test is chi-square, at the chain level and at the entry level using tier_b (results/qc_bfactor_clusters_tests.json).

## Results

### Collection

The query returned 487 ESR1 and 39 ESR2 X-ray entries. All 526 mmCIF files downloaded and validated (data/raw/download_manifest.csv). Deposition dates run from 1995 to 2026 and are uneven: 113 entries were deposited in 2022, 73 in 2015 and 70 in 2016 (results/entries.csv).

### Chains

The 526 entries hold 1,666 chain-models. Of these, 149 belong to the ESR2 entries, and 462 peptide, 146 partner and 8 DNA chains are not receptor chains (results/attrition_funnel.csv). That leaves 901 ESR1 receptor chains in 444 entries. 788 chains pass the completeness screen, and they come from 351 entries. Of the 113 that fail, 109 have no LBD residues, 2 are short constructs and 2 have an incomplete observed LBD (results/lbd_chains.csv). All 788 chains enter the RMSD matrix (results/chain_status.csv). The 39 ESR2 entries contribute 86 eligible ESR2 chains, which are analysed separately.

### QC flags and tiers

Over all 526 entries, 481 pass the resolution flag, 513 pass R-free (11 fail, 2 unknown) and 424 pass mean B (results/entries.csv). 397 pass all three, 128 fail at least one and 1 is unknown. In the 351 entries with eligible ESR1 chains, the counts are 315, 343 and 267, and 247 entries pass all three. By tier, the 526 entries are 397 A, 125 B, 2 C and 2 unknown. The 351 entries are 247 A, 103 B and 1 C (results/qc_tiers.csv). At the chain level, the 788 chains in the matrix are 544 tier A, 240 tier B and 4 tier C. The 244 tier B and C chains carry a written warning (results/chain_status.csv).

### Ligands

The ligand table has 480 distinct non-water codes (results/ligands.csv), of which 18 are marked as additives. Additives lead the raw ranking: magnesium is in 88 entries and chloride in 29. Estradiol (EST) is the most frequent real ligand, in 28 entries (results/ligands_ranked.csv). Among the 788 eligible chains, 765 have one ligand of interest and 38 rows come from chains with several. Seven chains have no pocket ligand (status apo_or_unclear) and are kept. In total, 306 distinct ligands of interest occur (results/chain_ligand_of_interest.csv).

### Redundancy

Twelve entries in 6 groups are near-duplicates (results/near_duplicates.csv). Most redundancy is not of this exact kind. It is structural: at 1.5 Å, the two largest clusters hold 98.2% of the 788 chains (results/representative_set.csv). The cluster-size entropy is 1.13 bits, against a maximum of 2.81 for 7 clusters (results/entropy.csv).

### QC question 1: rare clusters and resolution

At the 1.5 Å cut, 7 entries have a chain in a rare cluster, against 344 that do not. The median resolution is 2.40 Å against 2.10 Å (p = 0.21), and the median R-free is 0.236 against 0.238 (p = 0.75). No cut reaches p < 0.05. The lowest p is 0.061, for resolution at 2.5 Å, with only 3 rare entries (results/qc_outliers_stats.csv, docs/qc_outliers.md). Rare conformations in this set are not explained by worse resolution or R-free. The rare groups are small, so this is weak evidence either way.

### QC question 2: high-B structures

High-B structures concentrate rather than scatter. At 1.5 Å, 42% of the chains in the largest cluster (445 chains) are in the high-B tertile. In the second cluster (329 chains) the share is 21% (chain-level chi-square p = 5.6 × 10⁻¹²). The entry-level test agrees: 66 of 219 entries in the largest cluster are tier B for mean B, against 16 of 124 in the second (p = 7.4 × 10⁻⁴) (results/qc_bfactor_clusters.csv, results/qc_bfactor_clusters_tests.json). The two small clusters with 5 and 3 chains are mostly high-B (7 of 8 chains), but they are too small to test.

### Does QC change the structural conclusion?

Reclustering only tier A + B chains (784) gives the same partition as all chains at every cut (ARI 1.0). Reclustering only tier A chains (544 chains, 247 entries) keeps the large clusters, with ARI 0.98 at 1.5 Å. It merges some small ones (4 clusters instead of 7 at 1.5 Å), and one all-chain cluster has no tier-A member (results/qc_clustering_delta.csv). Restricting to R-free ≤ 0.30 changes the entropy by at most 0.06 bits (results/entropy_qc_compare.csv). The structural redundancy reported here does not depend on the poorer structures. Some rare conformations rest only on tier B structures.
