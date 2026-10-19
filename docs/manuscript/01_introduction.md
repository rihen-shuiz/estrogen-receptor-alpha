# Introduction

Estrogen receptor alpha (ESR1, UniProt P03372) and estrogen receptor beta (ESR2, UniProt Q92731) are nuclear receptors. Their ligand-binding domains (LBDs) bind estradiol and many synthetic compounds. Ligand binding changes the position of helix 12 at the C-terminal edge of the LBD. That position controls whether coactivator peptides can bind. Agonists stabilise an active helix 12 position. Antagonists, including 4-hydroxytamoxifen, place helix 12 over the coactivator groove. Two early crystal structures show these two states: an agonist complex and an antagonist complex of the ESR1 LBD (Brzozowski et al., 1997; Shiau et al., 1998).

Ligands for these receptors fall into several classes. Agonists activate the receptor. Antagonists block it. Selective estrogen receptor modulators (SERMs) act as agonists in some tissues and as antagonists in others. Medicinal chemistry has produced many series within each class. Whether a compound acts as an agonist, antagonist or SERM is not always clear from its chemical structure alone. Structural data help with that question, but the labels used in the literature are not always consistent.

The Protein Data Bank holds many X-ray structures of these LBDs. This project collected 487 ESR1 X-ray entries and 39 ESR2 X-ray entries (data/pdb_ids_esr1.txt and data/pdb_ids_esr2.txt). Of the 526 entries, 351 ESR1 entries have at least one chain that passes the historical coverage screen. Those 351 form the selected ESR1 set. Ligand identities are taken from the entries' CCD codes. The table contains 480 distinct non-water codes.

Redundancy matters for these data. Many entries are the same protein with the same ligand, deposited by the same group. Such entries add little independent information about conformational variation, but they weigh heavily in any count. Deposits are also uneven over time. In this collection, 113 entries were deposited in 2022, and 73 and 70 in 2015 and 2016 (results/entries.csv). Burst depositions like these can make a few groups dominate the apparent diversity. This project therefore flags near-duplicate entries (12 entries in 6 groups, results/near_duplicates.csv). Entries that fail QC are flagged, not removed (docs/qc_summary.md).

The analysis uses a single reference chain, 7NFB chain B, for coordinate frames. That structure is an engineered variant, so the reference is a coordinate frame and not a wild-type model (docs/reference.md). The aim of the work is to describe the conformational diversity of the deposited ESR1 LBD structures and to test whether that diversity matches the ligand classes assigned to them.

## Scope of the data

- 487 ESR1 X-ray entries and 39 ESR2 X-ray entries were queried from RCSB for UniProt P03372 and Q92731 (data/pdb_ids_esr1.txt, data/pdb_ids_esr2.txt).
- All 526 mmCIF files were downloaded and parsed without error (data/raw/download_manifest.csv).
- 351 ESR1 entries form the selected set. 247 of them pass the combined QC cutoffs (docs/qc_summary.md).
- The ligand table has 480 distinct non-water CCD codes across 756 entry-ligand pairs (results/ligands.csv).
- The reference chain is 7NFB chain B, at 1.33 A resolution (docs/reference.md).

## Questions

- Do deposited ESR1 LBD structures fall into distinct conformational groups?
- Do those groups agree with the ligand classes assigned to the entries?
- How much of the observed diversity remains after removing near-duplicate entries and entries that fail QC?
- Which canonical positions differ in conformation between ESR1 and ESR2 LBD structures?

## References

- Brzozowski AM, et al. Molecular basis of agonism and antagonism in the oestrogen receptor. *Nature*. 1997;389(6652):753-758. doi:10.1038/39645. PMID 9338790.
- Shiau AK, Barstad D, Loria PM, Cheng L, Kushner PJ, Agard DA, Greene GL. The structural basis of estrogen receptor/coactivator recognition and the antagonism of this interaction by tamoxifen. *Cell*. 1998;95(7):927-937. doi:10.1016/s0092-8674(00)81717-1. PMID 9875847.
