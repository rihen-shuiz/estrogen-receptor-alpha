# ESR1 ligand-binding domain: structural diversity of the crystal structures

788 chains of the oestrogen receptor α ligand-binding domain (351 PDB entries) fall into two groups set by the position
of helix 12 (the lid): closed over the ligand pocket (cl7, 445 chains) or lying in the coactivator groove (cl4, 329).
This repository holds the data, the code that rebuilds every result from it, the figures and the checked evidence.
Summary of all stages: [docs/project_summary.md](docs/project_summary.md).

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
PYTHON=.venv/bin/python bash reproduce.sh
```

`reproduce.sh` runs the pipeline stages in order, redraws every figure and the slides, and re-checks every number of the evidence index
(about one minute). It ends with `all numbers match`.

## Layout (pipeline order)

| Folder | Content |
|---|---|
| `data/` | Inputs. `ids/` PDB IDs and download checksums; `structures/` entries, chain table, quality flags, additive counts; `coordinates/` deposited backbone atoms of the 788 chains, alignment core, reference; `lid/`, `outliers/`, `ligands/`, `maps/`, `esr2/` per-chain tables from earlier stages (helix geometry, H12 mobility, lattice contacts, outlier flags, ligand descriptors, action classes, pockets, affinities, ESR2 matrices) |
| `pipeline/01_structures.py` | Counts of entries, chains, quality tiers, constructs |
| `pipeline/02_lid_groups.py` | Superposition, C-alpha RMSD (with and without residues 528-550), average-linkage groups, where the groups differ |
| `pipeline/03_lid_outliers.py` | Outliers of the lid-to-H11 distance, one by one, against ordinary variation |
| `pipeline/04_outliers_consensus.py` | Consensus RMSD outliers by group |
| `pipeline/05_maps.py` | Sammon map, self-organising map, UMAP on PCA scores |
| `pipeline/06_lid_hidden.py` | Clustering with the lid region removed |
| `pipeline/07_geometry_free.py` | Helix ranges (phi/psi), lengths and centres of mass without superposition |
| `pipeline/08_crystal_environment.py` | Additives, solvent and lattice contacts against lid order and mobility |
| `pipeline/09_ligands_pockets.py` | Ligand size and chemistry, action classes, pockets, filling, affinities |
| `pipeline/10_redundancy.py` | Hill numbers, set cover, copies of one crystal, novelty over time |
| `pipeline/11_esr2.py` | ESR2 groups and their distance to the ESR1 groups |
| `pipeline/figures/` | One script per figure (`fig_<name>.py` → `figures/<name>.png` and `.svg`) |
| `results/<stage>/` | Tables and `numbers.json` written by each stage |
| `figures/` | All figures |
| `evidence/` | `INDEX.md` (claims, numbers, figures, sources), `data/` (tables behind the figures), `rebuild.py` |
| `docs/` | Project summary, meeting notes, literature checked; `recent_results.pdf`, slides of the results since the October report (built by `recent_results.py`) |

## Conventions

Groups: average linkage on the C-alpha RMSD at 1.5 Å after superposition on 123 core positions of the five most stable
helices onto chain 7NFB B. Ordinary variation: an effect counts only if larger than both the within-group SD and the
95th percentile of the difference between copies of one crystal. Seed 42. Distances in Å. Action classes are uncurated
ChEMBL mechanism labels (51 of 362 entry-ligand rows) until the hand-labelled set replaces them.

## Inputs not rebuilt here

The tables in `data/` come from the mmCIF files listed in `data/ids/` (checksums in `download_manifest.csv` and
`sha256sums.txt`) through extraction, pydssp secondary structure, lattice-contact counting, P2Rank pockets, a buriedness
grid, RDKit descriptors and ChEMBL queries. The pipeline starts from these tables and from the deposited backbone
coordinates; superposition, RMSD and groups are recomputed and checked against the frozen groups.
