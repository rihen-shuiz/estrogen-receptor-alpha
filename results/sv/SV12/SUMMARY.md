Not done / substituted: teammate file data/sv/action_class_teammate.csv absent at run time -> labels come from route (b) ChEMBL mechanism only; route (c) FAILED its validation gate (1/2 rows agree = 0.50; needs >= 0.80 of >= 10 rows), so its labels are not used; class_agreement.csv not written (no teammate file). All labels are provisional (results/sv/SV12/ only; results/ligand_action_labels.csv and Q10/Q11 untouched; no literature searched; no class from structure, cluster, chemistry or a figure). Declared before labelling where the run file is open: (1) 'modulator' text excludes the phrase 'allosteric modulator' so the explicit NEGATIVE ALLOSTERIC MODULATOR -> antagonist mapping can apply; (2) ChEMBL document lookup uses doi__iexact (the plain doi filter is case-sensitive: the CHEMBL1145862 test DOI upper-cased gives 0 hits, doi__iexact 1); (3) connectivity matches with several parents take the no-stereo-layer parent if present, else the lowest ChEMBL id (all candidates listed); (4) ESR1 targets = all 9 ChEMBL targets with component P03372 (CHEMBL206 + family/complex/PPI targets) for routes (b) and (c); activities for SV15 are CHEMBL206 only; (5) conflict = different classes counting partial_agonist as agonist; (6) chain label = label of the chain's SV00 ligand_id_primary; an entry counts in a class x cluster cell if >= 1 of its chains there has the class. 23 of 351 entries have no primary-citation DOI in the mmCIF or the Data API. Supporting files not named in the run file: chain_class.csv, numbers.json, chembl_checks.json (results/sv/SV12/) and data/sv/SV12/*.

# SV12 - Ligand identifiers, provisional action classes and ChEMBL activities

**Result:** 362 entry x ligand rows (349 entries, 306 CCD codes). Final provisional classes (rows): agonist 27, partial_agonist 0, SERM 20, SERD 4, antagonist 0, other 0, unlabelled 311.
- Conflicts flagged: 0

## Identifiers (step 1)
- RCSB chemcomp: 306/306 records; stereo SMILES missing 0, InChIKey missing 0 (all 306 formal charge 0, standard InChIKeys).
- Primary-citation DOI: 328 entries from mmCIF _citation (id primary), 0 from the Data API, 23 none; 87 distinct DOIs.

## ChEMBL checks and match (step 2, 6)
- ChEMBL cache: 1237 responses; 1076 list pages of the pipeline (limit 1000), 0 with a non-null page_meta.next, so every list fitted one page (paging code in cached_http.chembl_get_all); 1 field probe (activity, limit 2) not used as data; 0 un-cached requests on a re-run.
- ChEMBL release: ChEMBL_37 (release date 2026-05-01; status resource cached data/sv/chembl_cache/status_6da02864d76f03ef892df65e2647831271838ebd.json).
- EST test: RCSB InChIKey VOXZDWNPVJITMN-ZBRFXRBCSA-N -> CHEMBL135 (ESTRADIOL), parent CHEMBL135; E2 reference for SV15 = CHEMBL135 (full match).
- CHEMBL206: SINGLE PROTEIN, Homo sapiens, components [['P03372', 'PROTEIN']] -> single-protein target for P03372 = True; the only single-protein target with P03372: ['CHEMBL206']; 9 targets contain P03372.
- Match per CCD: full 129, connectivity_only 31, no_chembl 146 (of 306); several candidate parents 3 (match types: connectivity_only); 157 distinct parents.
- Rows: full 172, connectivity_only 32 (excluded from the SV13/SV15 primary analyses, primary_set = False), no_chembl 158.

## Labels (steps 3-4)
- Route (b): 19 mechanism records for the parents, 16 on ESR1 targets, for 14 parents; parent labels: SERM 6, agonist 4, SERD 4; rows: unlabelled 311, agonist 27, SERM 20, SERD 4.
  - CHEMBL3581693 (7I0): SERD <- DEGRADER [Estrogen receptor degrader]
  - CHEMBL4297509 (I0V): SERD <- DEGRADER [Estrogen receptor alpha degrader]
  - CHEMBL4475463 (L5B): SERD <- DEGRADER [Estrogen receptor alpha degrader]
  - CHEMBL4650316 (ZNM): SERD <- DEGRADER [Estrogen receptor alpha degrader]
  - CHEMBL193482 (ESL): SERM <- MODULATOR [Estrogen receptor modulator]
  - CHEMBL328190 (C3D): SERM <- MODULATOR [Estrogen receptor modulator]
  - CHEMBL46740 (29S): SERM <- MODULATOR [Estrogen receptor modulator]
  - CHEMBL489 (OHT): SERM <- MODULATOR [Estrogen receptor modulator]
  - CHEMBL81 (RAL): SERM <- MODULATOR [Estrogen receptor modulator]
  - CHEMBL954 (53Q): SERM <- MODULATOR [Estrogen receptor alpha modulator]
  - CHEMBL1230314 (4OH): agonist <- AGONIST [Estrogen receptor agonist]
  - CHEMBL135 (EST): agonist <- AGONIST [Estrogen receptor alpha agonist]
  - CHEMBL1405 (J3Z): agonist <- AGONIST [Estrogen receptor alpha agonist]
  - CHEMBL411 (DES): agonist <- AGONIST [Estrogen receptor alpha agonist]
- Route (c): 31 of 87 DOIs found as ChEMBL documents (31 documents, 1766 ESR1-target activities); row outcomes: no_chembl_parent 158, doi_not_in_chembl 134, pattern_match 24, no_pattern_match 18, no_doi 16, ligand_has_no_esr1_activity_in_document 11, no_activity_with_pchembl_ge6 1; labels: unlabelled 338, SERD 19, antagonist 4, agonist 1.
- Validation gate: 2 rows labelled by both (2 ligands, 2 entries), agreement 1/2 = 0.50 -> FAIL ([DEFAULT] >= 0.80 of >= 10 rows).
  - disagreements (b vs c): 7MSA/ZNM SERD vs antagonist
- Final rows by source: chembl_mechanism SERD 4; chembl_mechanism SERM 20; chembl_mechanism agonist 27; none unlabelled 311.

## Coverage gate [DEFAULT] (SV13-SV15: a class comparison is tested only if each side has >= 5 entries)

| set | class | chains all | entries all | chains cl7 | entries cl7 | chains cl4 | entries cl4 | chains minor | entries minor | >= 5 entries (all) |
|---|---|---|---|---|---|---|---|---|---|---|
| primary | agonist | 55 | 27 | 52 | 25 | 1 | 1 | 2 | 1 | yes |
| primary | partial_agonist | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| primary | SERM | 53 | 20 | 4 | 3 | 48 | 19 | 1 | 1 | yes |
| primary | SERD | 9 | 4 | 0 | 0 | 9 | 4 | 0 | 0 | no |
| primary | antagonist | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| primary | other | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| primary | unlabelled | 602 | 271 | 350 | 174 | 242 | 94 | 10 | 4 |  |
| all | agonist | 55 | 27 | 52 | 25 | 1 | 1 | 2 | 1 | yes |
| all | partial_agonist | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| all | SERM | 53 | 20 | 4 | 3 | 48 | 19 | 1 | 1 | yes |
| all | SERD | 9 | 4 | 0 | 0 | 9 | 4 | 0 | 0 | no |
| all | antagonist | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| all | other | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| all | unlabelled | 671 | 301 | 389 | 193 | 271 | 105 | 11 | 5 |  |

- 788 SV00 chains; 7 without a ligand of interest are in 'unlabelled'; chains whose two ligand codes carry different classes: 0.

## Activities for SV15 (step 5)
- 1170 activities on CHEMBL206 for 128 of 157 matched parents (551 with pChEMBL; assay types B 982, F 174, A 14); E2 (CHEMBL135): 280 activities in 266 assays.

## Verdict
- The provisional label set covers 51 of 362 rows (14 of 306 ligands, 14 ChEMBL parents, ChEMBL max_phase 2-4); 311 rows stay unlabelled (no ESR1 mechanism record, no ChEMBL parent, or route (c) not used).
- Coverage gate (primary set, all clusters): >= 5 entries for agonist, SERM; below 5 entries: partial_agonist 0, SERD 4, antagonist 0, other 0. Classes and clusters overlap: agonist entries cl7 25 / cl4 1; SERM entries cl7 3 / cl4 19; SERD entries cl7 0 / cl4 4.
- Consequence for the class comparisons named in SV13 (by the entry counts above): SERD vs SERM, antagonist vs SERM and agonist vs partial_agonist each have a side below 5 entries; agonist and SERM are the only classes with >= 5 entries.
- Note: ESL (ESTRIOL, CHEMBL193482) is SERM because its only ESR1 mechanism record has action type MODULATOR ('Estrogen receptor modulator'); kept as the [DEFAULT] map gives it.

## Files
- results/sv/SV12/ligand_ids.csv (54.3 KB)
- results/sv/SV12/entry_doi.csv (63.3 KB)
- results/sv/SV12/action_class.csv (55.3 KB)
- results/sv/SV12/class_coverage.csv (3.8 KB)
- results/sv/SV12/route_c_validation.csv (34.1 KB)
- results/sv/SV12/chembl_activities.csv (521.7 KB)
- results/sv/SV12/chain_class.csv (59.6 KB)
- results/sv/SV12/numbers.json (2.4 KB)
- results/sv/SV12/chembl_checks.json (3.0 KB)
- results/sv/SV12/env.txt (0.9 KB)
- data/sv/SV12/ccd_ids.csv (145.3 KB)
- data/sv/SV12/chembl_info.json (0.5 KB)
- data/sv/SV12/chembl_stdout.txt (0.5 KB)
- data/sv/SV12/document_esr1_activities.csv (473.3 KB)
- data/sv/SV12/doi_documents.csv (5.3 KB)
- data/sv/SV12/esr1_targets.csv (1.0 KB)
- data/sv/SV12/ids_info.json (0.6 KB)
- data/sv/SV12/mechanism_records.csv (3.8 KB)
- data/sv/SV12/route_b_per_parent.csv (1.2 KB)
- data/sv/SV12/route_c_evidence.csv (13.9 KB)
- docs/EXPLORATION.md: appended 'Suggested, not done (SV12, 2026-11-08)' (rule 1.1)
- caches: data/sv/rcsb_chemcomp/, data/sv/chembl_cache/ (index.jsonl), data/sv/api_cache/sv12_*
- scripts: scripts/sv/SV12_chembl_test.py, SV12_ids.py, SV12_chembl.py, SV12_classes.py, SV12_env.py, SV12_summary.py; scripts/sv/common/cached_http.py
