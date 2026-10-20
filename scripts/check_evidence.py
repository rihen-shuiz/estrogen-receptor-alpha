"""Re-check the evidence in docs/EVIDENCE.md: one independent recount per task, the invariants, and tl.py check.

Read-only: reads result files with pandas/numpy/scipy, never imports or runs a producing script, writes nothing.
Every task in TASKS is a DONE task logged after CUTOFF. A recount is one expression (EXPR) evaluated with the helpers
below; its value is compared with EXPECTED, the claim taken from the task's log note. Floats match when they agree to
the number of decimals written in the claim. Prints PASS or FAIL per task and per invariant, then a total.
Exit code 0 only if everything passes.
Run from the repo root: .venv/Scripts/python scripts/check_evidence.py
"""
import ast
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu, spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score as ari

ROOT = Path(__file__).resolve().parents[1]
CUTOFF = "2026-10-18T23:56+05:00"
CUTS = ["cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5"]
PY = ".venv/Scripts/python"


# ---------- helpers available to every recount expression ----------
def csv(path, **kw):
    return pd.read_csv(ROOT / path, low_memory=False, **kw)


def npy(path):
    return np.load(ROOT / path)


def txt(path):
    return (ROOT / path).read_text(encoding="utf-8")


def up(m):
    """Upper triangle of a square matrix, diagonal excluded."""
    return m[np.triu_indices(len(m), 1)]


def nonwhite(path):
    a = np.asarray(Image.open(ROOT / path).convert("RGB"))
    return int((a.min(axis=2) < 245).sum())


def nclust(mask, cut):
    """Average-linkage clusters at a distance cut on the slice of results/rmsd_entry.npy given by a boolean mask."""
    m = npy("results/rmsd_entry.npy")[np.ix_(mask, mask)]
    return int(fcluster(linkage(squareform(m, checks=False), "average"), cut, "distance").max())


def medoid(cluster):
    """Chain with the smallest summed RMSD to the other chains of a cut_1.5 cluster."""
    c, m = csv("results/entry_clusters.csv"), npy("results/rmsd_entry.npy")
    idx = np.flatnonzero(c["cut_1.5"].to_numpy() == cluster)
    return c.chain_key[idx[np.argmin(m[np.ix_(idx, idx)].sum(axis=1))]]


def core_positions():
    ranges = json.loads(txt("esr_analysis/resources/config/alignment_core_v3.json"))["core_ranges"]
    return [p for a, b in ranges for p in range(a, b + 1)]


HELPERS = dict(pd=pd, np=np, json=json, csv=csv, npy=npy, txt=txt, up=up, nonwhite=nonwhite, nclust=nclust,
               medoid=medoid, core_positions=core_positions, CUTS=CUTS, fcluster=fcluster, linkage=linkage,
               squareform=squareform, spearmanr=spearmanr, kruskal=kruskal, mannwhitneyu=mannwhitneyu,
               chi2_contingency=chi2_contingency, ari=ari, PCA=PCA)


def T(name, outputs, inputs, code, cmd, claim, expr, expected, caveats=""):
    return dict(name=name, outputs=outputs, inputs=inputs, code=code, cmd=cmd, claim=claim, expr=expr,
                expected=expected, caveats=caveats)


def run(script):
    return f"PYTHONPATH=. {PY} {script}"


TASKS = {
    "S06": T("Note that distances use only C-alpha atoms", ["docs/limitations.md"],
             ["results/alignment_core_rmsd.csv", "results/rmsd_entry.npy"], [], "none (text edited by hand)",
             "788 chains / 351 entries, 123-residue core v3",
             "(npy('results/rmsd_entry.npy').shape[0], csv('results/alignment_core_rmsd.csv').entry_id.nunique(), json.loads(txt('esr_analysis/resources/config/alignment_core_v3.json'))['residues_total'], '788 x 788' in txt('docs/limitations.md'))",
             (788, 351, 123, True), "Text task: no producing script. The paragraph itself is unchanged except the first sentence."),
    "S08": T("Group all chains by structural similarity at four distance cuts", ["results/entry_clusters.csv", "results/linkage_entry.npy"],
             ["results/rmsd_entry.npy", "results/rmsd_entry_index.csv", "results/qc_tiers.csv"],
             [("scripts/cluster_selected.py", "main"), ("esr_analysis/cluster_structures.py", "structural_groups")],
             run("scripts/cluster_selected.py"), "clusters 16/7/4/3 at 1.0/1.5/2.0/2.5 A; qc_label A 544, B 240, C 4",
             "[csv('results/entry_clusters.csv')[c].nunique() for c in CUTS] + [int(fcluster(linkage(squareform(npy('results/rmsd_entry.npy'), checks=False), 'average'), 1.5, 'distance').max())] + [int((csv('results/entry_clusters.csv').qc_label == k).sum()) for k in 'ABC']",
             [16, 7, 4, 3, 7, 544, 240, 4], "The 5th value re-clusters the saved matrix independently with scipy (7 at 1.5 A). qc_label is the entry-level K05 label copied per chain."),
    "S09": T("Table of cluster counts and sizes per cut", ["results/cluster_stats.csv"], ["results/rmsd_entry.npy", "results/entry_clusters.csv"],
             [("scripts/cluster_stats.py", "main"), ("scripts/cluster_stats.py", "cut_stats")], run("scripts/cluster_stats.py"),
             "clusters 16/7/4/3; top3 chain share 0.954/0.989/0.997/1.000",
             "[csv('results/entry_clusters.csv')[c].nunique() for c in CUTS] + [round(csv('results/entry_clusters.csv')[c].value_counts().head(3).sum() / len(csv('results/entry_clusters.csv')), 3) for c in CUTS]",
             [16, 7, 4, 3, 0.954, 0.989, 0.997, 1.0], "Spec asks for the share of entries in the 3 largest clusters; the file gives the share of chains (rows are chains). Entry-level view is K07."),
    "S10": T("Tree drawing (dendrogram) of the clustering", ["figures/dendrogram.png"], ["results/linkage_entry.npy", "results/entry_clusters.csv"],
             [("scripts/figure_dendrogram.py", "main")], run("scripts/figure_dendrogram.py"), "788 leaves, cut 1.5 A = 7 clusters",
             "(npy('results/linkage_entry.npy').shape[0] + 1, int(fcluster(npy('results/linkage_entry.npy'), 1.5, 'distance').max()), nonwhite('figures/dendrogram.png') > 0)",
             (788, 7, True), "Static PNG; leaves unlabeled, table twin results/entry_clusters.csv."),
    "S12": T("Most central structure of each cluster", ["results/medoids.csv"], ["results/rmsd_entry.npy", "results/entry_clusters.csv"],
             [("scripts/medoids.py", "main")], run("scripts/medoids.py"), "sizes 5/3/3/329/1/2/445 chains; largest-cluster medoid 5TM1|A|1",
             "csv('results/medoids.csv').cluster_size.tolist() + [medoid(7)]", [5, 3, 3, 329, 1, 2, 445, "5TM1|A|1"],
             "Spec says 'the entry'; the file picks a chain (an entry can have chains in two clusters: 9B25 A and B are in clusters 6 and 3)."),
    "S13": T("Does keeping only good-quality structures change the clusters?", ["results/qc_clustering_delta.csv"],
             ["results/rmsd_entry.npy", "results/entry_clusters.csv"], [("scripts/qc_clustering_delta.py", "main")],
             run("scripts/qc_clustering_delta.py"), "tier A only clusters 14/4/3/2; tier A+B ARI 1.0 at every cut",
             "[nclust(csv('results/entry_clusters.csv').qc_label.eq('A').to_numpy(), t) for t in (1.0, 1.5, 2.0, 2.5)] + csv('results/qc_clustering_delta.csv').query(\"subset == 'tier_A+B'\").ari_vs_all.round(3).tolist()",
             [14, 4, 3, 2, 1.0, 1.0, 1.0, 1.0], "Tier-A cluster counts re-derived by re-clustering the matrix slice independently."),
    "S14": T("Plot of cluster count against distance cut", ["figures/structural_sensitivity.png"], ["results/linkage_entry.npy", "results/cluster_stats.csv"],
             [("scripts/figure_sensitivity.py", "main"), ("scripts/figure_sensitivity.py", "stats_at")], run("scripts/figure_sensitivity.py"),
             "reported cuts 16/7/4/3",
             "[int(fcluster(npy('results/linkage_entry.npy'), t, 'distance').max()) for t in (1.0, 1.5, 2.0, 2.5)] + [nonwhite('figures/structural_sensitivity.png') > 0]",
             [16, 7, 4, 3, True], "Static PNG; numbers in the figure checked against cluster_stats.csv by the producing script."),
    "S07": T("Compare C-alpha distances with all-atom distances on 50 entries", ["results/rmsd_allatom_subset.csv", "results/rmsd_allatom_subset_summary.json"],
             ["results/rmsd_entry.npy", "results/alignment/alignment_transforms.json", "results/lbd_extract/residue_mapping.csv"],
             [("scripts/allatom_rmsd_subset.py", "main"), ("scripts/allatom_rmsd_subset.py", "pair_rmsd")], run("scripts/allatom_rmsd_subset.py"),
             "1225 pairs over 50 entries; Spearman CA vs all-atom 0.985",
             "(lambda d: (len(d), round(spearmanr(d.ca_rmsd, d.allatom_rmsd)[0], 3), pd.concat([d.entry_a, d.entry_b]).nunique()))(csv('results/rmsd_allatom_subset.csv'))",
             (1225, 0.985, 50), "Uses the saved S03 superposition for both measures (no all-atom refit). Needs mmCIFs in data/raw (git-ignored). Subset = historical chain per entry, stratified by cluster."),
    "E05": T("Superimpose each matched ESR2 chain on its ESR1 partner", ["results/pair_superposition.csv"],
             ["results/matched_pairs.csv", "results/lbd_chains.csv", "results/lbd_extract/residue_mapping.csv"],
             [("scripts/pair_superposition.py", "main"), ("scripts/pair_superposition.py", "shared_rmsd")], run("scripts/pair_superposition.py"),
             "671 chain pairs, 135 primary; core median 0.524 (max 1.876), LBD median 1.434, H12 median 1.957 A",
             "(lambda d: (len(d), int(d.primary_pair.sum()), d.core_rmsd_A.median(), d.core_rmsd_A.max(), d.lbd_rmsd_A.median(), d.h12_rmsd_A.median()))(csv('results/pair_superposition.csv'))",
             (671, 135, 0.524, 1.876, 1.434, 1.957), "All chain combinations kept; primary_pair marks one pair per entry pair."),
    "R01": T("How evenly the chains spread over clusters (entropy)", ["results/entropy.csv"], ["results/entry_clusters.csv"],
             [("scripts/cluster_entropy.py", "main"), ("scripts/cluster_entropy.py", "entropy_bits")], run("scripts/cluster_entropy.py"),
             "entropy 1.60/1.13/1.06/1.03 bits",
             "(lambda d: [float(-(p * np.log2(p)).sum()) for p in (d[c].value_counts(normalize=True) for c in CUTS)])(csv('results/entry_clusters.csv'))",
             [1.60, 1.13, 1.06, 1.03]),
    "R02": T("How many clusters appear as structures are added at random", ["results/diversity_curve.csv"], ["results/entry_clusters.csv"],
             [("scripts/diversity_curve.py", "main"), ("scripts/provisional_common.py", "clusters_seen")], run("scripts/diversity_curve.py"),
             "3152 rows; n_90 = 520/488/315/147",
             "(lambda d: (len(d), [int(g.n[g.mean_clusters >= 0.9 * g.n_clusters_total].min()) for _, g in d.groupby('cut')]))(csv('results/diversity_curve.csv'))",
             (3152, [520, 488, 315, 147])),
    "E06": T("Per-residue shift between ESR1 and ESR2 for the first 5 pairs", ["results/first_five_displacements.csv"],
             ["results/pair_superposition.csv", "results/lbd_extract/residue_mapping.csv", "results/numbering_map_esr1_esr2.csv"],
             [("scripts/receptor_displacements.py", "main"), ("scripts/receptor_displacements.py", "displacements")],
             run("scripts/receptor_displacements.py --first 5"), "1125 rows, 5 pairs; per-pair median 0.53-0.73 A; >30 A only at 538-544",
             "(lambda d: (len(d), d.pair.nunique(), d.groupby('pair').displacement_A.median().min(), d.groupby('pair').displacement_A.median().max(), int(d.esr1_pos[d.displacement_A > 30].min()), int(d.esr1_pos[d.displacement_A > 30].max())))(csv('results/first_five_displacements.csv'))",
             (1125, 5, 0.53, 0.73, 538, 544), "All five first pairs use the same ESR1 chain (1A52|B|1, estradiol), so they are not a varied sample."),
    "E07": T("Per-residue shift for all matched pairs (position x pair table)", ["results/displacement_matrix.csv", "results/displacement_pairs.csv"],
             ["results/pair_superposition.csv", "results/lbd_extract/residue_mapping.csv"], [("scripts/receptor_displacements.py", "main")],
             run("scripts/receptor_displacements.py --all"), "246 positions x 135 pairs",
             "(lambda m: (m.shape[0], m.shape[1] - 1))(csv('results/displacement_matrix.csv'))", (246, 135),
             "Columns are the 135 primary pairs only, not all 671 chain pairs of E05."),
    "R04": T("Number of structures needed to see 90% of clusters", ["results/saturation.csv"], ["results/diversity_curve.csv"],
             [("scripts/saturation.py", "main")], run("scripts/saturation.py"), "n_90 = 520/488/315/147 (66/62/40/19% of chains)",
             "(lambda d: (d.n_90.tolist(), d.fraction_of_chains.tolist()))(csv('results/saturation.csv'))",
             ([520, 488, 315, 147], [0.66, 0.62, 0.40, 0.19])),
    "R03": T("Figure of the cumulative diversity curve with a 90% line", ["figures/diversity_curve.png"], ["results/diversity_curve.csv", "results/saturation.csv"],
             [("scripts/figure_diversity_curve.py", "main")], run("scripts/figure_diversity_curve.py"), "4 cuts drawn, non-blank figure",
             "(len(csv('results/saturation.csv')), nonwhite('figures/diversity_curve.png') > 0)", (4, True), "Static PNG."),
    "Q14": T("Are rare structural groups worse in resolution or R-free?", ["results/qc_outliers.csv", "docs/qc_outliers.md"],
             ["results/entry_clusters.csv", "results/qc_tiers.csv"], [("scripts/qc_outliers.py", "main")], run("scripts/qc_outliers.py"),
             "1.5 A: 7 rare vs 344 common entries; resolution 2.40 vs 2.10 A (p 0.21); R-free p 0.75",
             "(lambda e: (int(e.rare.sum()), int((~e.rare).sum()), e.resolution_A[e.rare].median(), e.resolution_A[~e.rare].median(), mannwhitneyu(e.resolution_A[e.rare], e.resolution_A[~e.rare]).pvalue, mannwhitneyu(e.r_free[e.rare], e.r_free[~e.rare]).pvalue))(csv('results/qc_outliers.csv').groupby('entry_id').agg(rare=('rare_1.5', 'any'), resolution_A=('resolution_A', 'first'), r_free=('r_free', 'first')))",
             (7, 344, 2.40, 2.10, 0.21, 0.75), "'Rare' = cluster with fewer than 10 chains: a threshold chosen in the script, not given by the spec."),
    "Q15": T("Do high B-factor structures gather in particular clusters?", ["results/qc_bfactor_clusters.csv", "results/qc_bfactor_clusters_tests.json"],
             ["results/entry_clusters.csv", "results/helix_bfactor_by_chain.csv", "results/qc_tiers.csv"], [("scripts/qc_bfactor_clusters.py", "main")],
             run("scripts/qc_bfactor_clusters.py"), "high-B share 42% (cluster 7) vs 21% (cluster 4), chain chi2 p 5.6e-12; clusters 1-2: 7 of 8 chains high-B",
             "(lambda d: ((d[d['cut_1.5'] == 7].b_group == 'high').mean(), (d[d['cut_1.5'] == 4].b_group == 'high').mean(), f\"{chi2_contingency(pd.crosstab(d[d['cut_1.5'].isin([4, 7])]['cut_1.5'], d.b_group))[1]:.1e}\", int((d[d['cut_1.5'].isin([1, 2])].b_group == 'high').sum()), int(d['cut_1.5'].isin([1, 2]).sum())))(csv('results/entry_clusters.csv').merge(csv('results/helix_bfactor_by_chain.csv')[['chain_key', 'b_group']], on='chain_key'))",
             (0.42, 0.21, "5.6e-12", 7, 8), "Chain-level counts repeat entries; the entry-level test (p 7.4e-4) is the safer reading. Only clusters with >= 20 chains are tested."),
    "K06": T("One status row for every ESR1 chain, none dropped", ["results/chain_status.csv"],
             ["results/chain_loss_reasons.csv", "results/qc_tiers.csv", "results/helix_bfactor_by_chain.csv", "results/h12_metric.csv"],
             [("scripts/chain_status.py", "main")], run("scripts/chain_status.py"), "901 chains / 444 entries; in matrix 788, out 113; qc A 628 / B 261 / C 4 / unknown 8",
             "(lambda d: (len(d), d.entry_id.nunique(), int(d.in_matrix.sum()), int((~d.in_matrix).sum()), [int((d.qc_label == k).sum()) for k in ['A', 'B', 'C', 'unknown']]))(csv('results/chain_status.csv'))",
             (901, 444, 788, 113, [628, 261, 4, 8])),
    "K09": T("Does the choice of alignment core change the clusters?", ["results/core_sensitivity.csv", "results/rmsd_old_core.npy"],
             ["results/alignment/aligned_coordinates_old_core.npz", "results/rmsd_entry.npy", "results/entry_clusters.csv"],
             [("scripts/core_sensitivity.py", "main")], run("scripts/core_sensitivity.py"),
             "old-core matrix 788x788; Spearman 0.988; median |diff| 0.044 A (max 2.22); old-core clusters 20/7/6/4",
             "(lambda a, b: (a.shape, spearmanr(up(a), up(b))[0], float(np.median(np.abs(up(a) - up(b)))), float(np.abs(up(a) - up(b)).max()), [int(fcluster(linkage(squareform(a, checks=False), 'average'), t, 'distance').max()) for t in (1.0, 1.5, 2.0, 2.5)]))(npy('results/rmsd_old_core.npy'), npy('results/rmsd_entry.npy'))",
             ((788, 788), 0.988, 0.044, 2.22, [20, 7, 6, 4]), "Spec allows this new matrix explicitly (the primary matrix is untouched)."),
    "E08": T("Heatmap of the ESR1/ESR2 per-residue shifts", ["figures/displacement_heatmap.png"], ["results/displacement_matrix.csv", "results/displacement_pairs.csv"],
             [("scripts/figure_displacement_heatmap.py", "main"), ("scripts/figure_displacement_heatmap.py", "draw")],
             run("scripts/figure_displacement_heatmap.py"), "246 positions drawn, non-blank figure",
             "(csv('results/displacement_matrix.csv').shape[0], nonwhite('figures/displacement_heatmap.png') > 0)", (246, True),
             "Static PNG; colour clipped at 5 A."),
    "E10": T("Positions that shift by more than 1.5 A in most pairs", ["results/consistent_positions.csv"],
             ["results/displacement_matrix.csv", "results/numbering_map_esr1_esr2.csv"], [("scripts/consistent_positions.py", "main")],
             run("scripts/consistent_positions.py"), "48 consistent positions; 0 in core, 7 in H12, 29 with a different residue",
             "(lambda m, nm: (lambda hit: (int(hit.sum()), int((hit & m.index.isin(core_positions())).sum()), int((hit & m.index.isin(range(538, 548))).sum()), int((hit & np.array([nm.get(p) == False for p in m.index])).sum())))(((m > 1.5).sum(axis=1) / m.notna().sum(axis=1) > 0.5) & (m.notna().sum(axis=1) >= 20)))(csv('results/displacement_matrix.csv', index_col='esr1_pos'), dict(zip(csv('results/numbering_map_esr1_esr2.csv').esr1_pos, csv('results/numbering_map_esr1_esr2.csv').same_residue)))",
             (48, 0, 7, 29), "'Most pairs' = more than 50% of pairs, with at least 20 observed pairs: thresholds chosen in the script."),
    "R06": T("Entropy with and without structures failing R-free", ["results/entropy_qc_compare.csv"], ["results/rmsd_entry.npy", "results/entry_clusters.csv", "results/qc_tiers.csv"],
             [("scripts/entropy_qc_compare.py", "main"), ("scripts/entropy_qc_compare.py", "n_90")], run("scripts/entropy_qc_compare.py"),
             "R-free set 767 chains / 343 entries, 21 out; clusters 17/6/4/3",
             "(lambda d, m: (int(m.sum()), d.entry_id[m].nunique(), int((~m).sum()), [nclust(m, t) for t in (1.0, 1.5, 2.0, 2.5)]))(csv('results/entry_clusters.csv'), csv('results/entry_clusters.csv').entry_id.map(csv('results/qc_tiers.csv').set_index('entry_id').r_free).le(0.30).to_numpy())",
             (767, 343, 21, [17, 6, 4, 3]), "A tier_A row was added beyond the spec. No tier A+B row (the S13 file has it). n_90 is unstable when clusters are few."),
    "K07": T("One-chain-per-entry view next to the all-chain view", ["results/entry_level_clusters.csv", "results/entry_vs_chain_clusters.csv"],
             ["results/rmsd_entry.npy", "results/entry_clusters.csv", "results/chain_status.csv", "results/helix_bfactor_by_chain.csv"],
             [("scripts/entry_level_clusters.py", "main")], run("scripts/entry_level_clusters.py"), "351 entries; clusters 18/6/3/3; ARI 0.84/0.97/0.98/0.99",
             "(lambda d: (len(d), bool(d.entry_id.is_unique), [d[c].nunique() for c in CUTS], [ari(d['chain_level_' + c], d[c]) for c in CUTS]))(csv('results/entry_level_clusters.csv'))",
             (351, True, [18, 6, 3, 3], [0.84, 0.97, 0.98, 0.99])),
    "S15": T("Smallest set of representative structures", ["results/representative_set.csv"], ["results/medoids.csv", "results/entry_clusters.csv", "results/rmsd_entry.npy"],
             [("scripts/representative_set.py", "main")], run("scripts/representative_set.py"), "7 representatives, 6 QC tier A; cluster 2 keeps 1R5K|B|1; top 2 clusters cover 98.2% of chains",
             "(lambda r, s: (len(r), int(r.qc_pass.sum()), r.chain_key[r.cluster == 2].item(), s.head(2).sum() / s.sum()))(csv('results/representative_set.csv'), csv('results/entry_clusters.csv')['cut_1.5'].value_counts())",
             (7, 6, "1R5K|B|1", 0.982), "Definition choice: '90% of diversity' read as 90% of clusters, so all 7 are needed; by chain share 2 would do. Needs a human decision."),
    "K08": T("Plain note for the professor on what data is kept", ["docs/data_retention.md"], ["results/chain_status.csv", "results/attrition_funnel.csv"],
             [], "none (text written by hand)", "244 chains with a warning, 113 outside the matrix; figures quoted in the note",
             "(lambda s, t: (int((s.in_matrix & (s.warning != 'none')).sum()), int((~s.in_matrix).sum()), all(x in t for x in ['244', '113', '1.255', '0.988'])))(csv('results/chain_status.csv'), txt('docs/data_retention.md'))",
             (244, 113, True), "Text task; only these figures are recounted here."),
    "X02": T("Distance matrix between all ESR2 structures", ["results/rmsd_esr2.npy", "results/rmsd_esr2_index.txt", "results/esr2_reference.json"],
             ["results/lbd_chains.csv", "results/esr2_entries.csv", "results/lbd_extract/residue_mapping.csv"], [("scripts/esr2_rmsd.py", "main")],
             run("scripts/esr2_rmsd.py"), "86 chains / 39 entries; reference 3OLL|A|1; pairwise median 0.947 A (max 6.64); core median 0.572 A",
             "(lambda m, ix, al: (m.shape[0], bool(np.allclose(m, m.T)), len({k.split('|')[0] for k in ix}), float(np.median(up(m))), float(up(m).max()), al.core_rmsd_A.median(), json.loads(txt('results/esr2_reference.json'))['chain_key']))(npy('results/rmsd_esr2.npy'), txt('results/rmsd_esr2_index.txt').split(), csv('results/esr2_alignment.csv'))",
             (86, True, 39, 0.947, 6.64, 0.572, "3OLL|A|1"), "Stated deviation: the reference core-fraction minimum (1.0) was not applied because no ESR2 chain reaches it (best 0.976)."),
    "X15": T("Plot of ChEMBL potency per structural cluster", ["figures/potency_by_cluster.png", "results/potency_by_cluster.csv"],
             ["results/chain_ligand_of_interest.csv", "results/chembl_pdb_matches.csv", "data/chembl_esr1_activities.csv", "results/entry_clusters.csv"],
             [("scripts/potency_by_cluster.py", "ligand_table"), ("scripts/potency_by_cluster.py", "figure")], run("scripts/potency_by_cluster.py"),
             "82 ligands; cluster 4 n=53, 7 n=27, 1 and 2 n=1",
             "(lambda d: (len(d), d.cluster.value_counts().to_dict(), nonwhite('figures/potency_by_cluster.png') > 0))(csv('results/potency_by_cluster.csv'))",
             (82, {4: 53, 7: 27, 1: 1, 2: 1}, True), "Unit = one PDB ligand placed in the cluster holding most of its chains; assays are mixed (IC50/Ki/Kd)."),
    "X16": T("Statistical tests: potency against cluster", ["results/potency_tests.csv"], ["results/potency_by_cluster.csv"],
             [("scripts/potency_by_cluster.py", "tests")], run("scripts/potency_by_cluster.py"),
             "Kruskal p 6.0e-5; median pChEMBL 8.31 vs 7.11; Spearman vs H12 RMSD 0.45 (n=79)",
             "(lambda d: (f'{kruskal(d.pchembl[d.cluster == 4], d.pchembl[d.cluster == 7]).pvalue:.1e}', d.pchembl[d.cluster == 4].median(), d.pchembl[d.cluster == 7].median(), spearmanr(d.dropna(subset=['median_h12_rmsd_A']).pchembl, d.dropna(subset=['median_h12_rmsd_A']).median_h12_rmsd_A)[0], int(d.median_h12_rmsd_A.notna().sum())))(csv('results/potency_by_cluster.csv'))",
             ("6.0e-05", 8.31, 7.11, 0.45, 79), "Descriptive: ligand class probably confounds the link (no ligand-class labels yet)."),
    "X03": T("Group the ESR2 structures the same way as ESR1", ["results/esr2_clusters.csv", "results/esr2_linkage.npy"], ["results/rmsd_esr2.npy", "results/esr2_entries.csv"],
             [("scripts/esr2_clusters.py", "main")], run("scripts/esr2_clusters.py"), "86 chains / 39 entries; clusters 7/5/4/3; sizes at 1.5 A 68/12/2/2/2",
             "(lambda d: (len(d), d.entry_id.nunique(), [d[c].nunique() for c in CUTS], d['cut_1.5'].value_counts().tolist(), int(fcluster(linkage(squareform(npy('results/rmsd_esr2.npy'), checks=False), 'average'), 1.5, 'distance').max())))(csv('results/esr2_clusters.csv'))",
             (86, 39, [7, 5, 4, 3], [68, 12, 2, 2, 2], 5)),
    "X04": T("Two-axis map (PCA) of the ESR2 structures", ["figures/esr2_pca_landscape.png", "results/esr2_pca_coords.csv"], ["results/rmsd_esr2.npy", "results/esr2_clusters.csv"],
             [("scripts/esr2_pca.py", "main")], run("scripts/esr2_pca.py"), "PC1 81.9%, PC2 13.0%",
             "tuple(PCA(n_components=2).fit(npy('results/rmsd_esr2.npy')).explained_variance_ratio_.tolist())", (0.819, 0.130),
             "Spec says 'same method as S11', which colours by ligand class; ESR2 has no ligand-class labels, so the figure colours by structural cluster."),
    "X05": T("ESR2 cluster sizes, with small clusters marked not interpretable", ["results/esr2_basin_counts.csv"], ["results/esr2_clusters.csv"],
             [("scripts/esr2_basin_counts.py", "main")], run("scripts/esr2_basin_counts.py"), "5 basins with 29/7/2/1/1 entries; 2 interpretable",
             "(lambda n: (sorted(n.tolist(), reverse=True), int((n >= 3).sum())))(csv('results/esr2_clusters.csv').groupby('cut_1.5').entry_id.nunique())",
             ([29, 7, 2, 1, 1], 2)),
    "W10": T("Manuscript section: dataset and quality control", ["docs/manuscript/02_dataset_qc.md"],
             ["results/attrition_funnel.csv", "results/entries.csv", "results/lbd_chains.csv"], [], "none (text written by hand)",
             "1,666 chain-models; 247 selected entries pass all three QC flags; >= 25 lines",
             "(lambda t, f, e, c: ('1,666' in t, int(f.n_before.iloc[0]), int(e[e.entry_id.isin(set(c[(c.identity == 'esr1') & c.eligible_historical].pdb_id))].qc_pass.eq('pass').sum()), '247' in t, len([l for l in t.splitlines() if l.strip()]) >= 25))(txt('docs/manuscript/02_dataset_qc.md'), csv('results/attrition_funnel.csv'), csv('results/entries.csv'), csv('results/lbd_chains.csv'))",
             (True, 1666, 247, True, True), "Text task; only these figures are recounted here."),
    "W09": T("Worked example notebook that rebuilds the main figure", ["notebooks/worked_example.ipynb"], ["results/rmsd_entry.npy", "results/entry_clusters.csv", "results/diversity_curve.csv", "results/saturation.csv"],
             [], "python -m jupyter nbconvert --to notebook --execute notebooks/worked_example.ipynb --output worked_example.ipynb",
             "executed without error; n_90 values from saturation.csv appear in its outputs",
             "(lambda nb: (all(c.get('execution_count') for c in nb['cells'] if c['cell_type'] == 'code'), any(o.get('output_type') == 'error' for c in nb['cells'] for o in c.get('outputs', [])), all(str(x) in json.dumps(nb) for x in csv('results/saturation.csv').n_90)))(json.loads(txt('notebooks/worked_example.ipynb')))",
             (True, False, True), "The script that generated the notebook was a scratch file, not committed. README was not updated, although the spec says 'runs from README'."),
    "W16": T("Script that collects the headline numbers for checking the manuscript", ["results/numbers.json", "docs/numbers.md"], ["results/*.csv", "results/*.json"],
             [("scripts/collect_numbers.py", "keep_all"), ("scripts/collect_numbers.py", "write")], f"{PY} scripts/collect_numbers.py",
             "279 numbers, 64 new keep-all numbers",
             "(lambda n: (len(n), sum(x['name'].startswith(('keep_', 'q14_', 'q15_', 's07_', 'x16_', 'x02_', 'x05_')) for x in n)))(json.loads(txt('results/numbers.json'))['numbers'])",
             (279, 64), "Needs Biopython (imports ligand_of_interest), so it runs only with the .venv Python."),
}


# ---------- comparison ----------
def decimals(x):
    s = repr(x)
    return len(s.split(".")[1]) if "." in s and "e" not in s else 0


def same(value, expected):
    if isinstance(expected, bool) or isinstance(value, (bool, np.bool_)):
        return bool(value) == expected and isinstance(expected, bool)
    if isinstance(expected, float):
        try:
            return abs(float(value) - expected) <= 0.5 * 10 ** -decimals(expected) + 1e-9
        except (TypeError, ValueError):
            return False
    if isinstance(expected, int):
        try:
            return float(value) == expected
        except (TypeError, ValueError):
            return False
    if isinstance(expected, (list, tuple)):
        value = list(value) if isinstance(value, (list, tuple, np.ndarray, pd.Series)) else value
        return isinstance(value, list) and len(value) == len(expected) and all(same(v, e) for v, e in zip(value, expected))
    if isinstance(expected, dict):
        return isinstance(value, dict) and set(map(int, value)) == set(expected) and all(same(value[k], expected[int(k)]) for k in value)
    return value == expected


def show(value):
    if isinstance(value, (np.floating, float)):
        return f"{float(value):.6g}"
    if isinstance(value, (np.integer,)):
        return str(int(value))
    if isinstance(value, (list, tuple)):
        inner = ", ".join(show(v) for v in value)
        return f"({inner})" if isinstance(value, tuple) else f"[{inner}]"
    if isinstance(value, dict):
        return "{" + ", ".join(f"{show(k)}: {show(v)}" for k, v in value.items()) + "}"
    return repr(value) if isinstance(value, str) else str(value)


def recount(tid):
    spec = TASKS[tid]
    try:
        value = eval(spec["expr"], dict(HELPERS))
    except Exception as exc:  # reported, never hidden
        return f"ERROR {type(exc).__name__}: {exc}", False
    return show(value), same(value, spec["expected"])


def tl_check(tid):
    out = subprocess.run([sys.executable, "tracker/tl.py", "check", tid], cwd=ROOT, capture_output=True, text=True).stdout
    lines = [l.strip() for l in out.splitlines() if l.strip().startswith(("PASS", "FAIL"))]
    return lines, bool(lines) and all(l.startswith("PASS") for l in lines)


# ---------- log and provenance ----------
def log_lines():
    rows = []
    for line in txt("tracker/esr1.log").splitlines():
        parts = [p.strip() for p in line.split(" | ")]
        if len(parts) >= 3 and re.match(r"\d{4}-\d\d-\d\dT", parts[0]):
            fields = dict(p.split("=", 1) for p in parts[3:] if "=" in p)
            rows.append(dict(time=parts[0], id=parts[1], state=parts[2], raw=line, **fields))
    return rows


def latest_states():
    last = {}
    for row in log_lines():
        if row["id"] != "-":
            last[row["id"]] = row
    return last


def sha16(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16]


def prov(tid):
    path = ROOT / "tracker/prov" / f"{tid}.json"
    if not path.exists():
        return "no prov file", False
    record = json.loads(path.read_text(encoding="utf-8"))
    drift = [item["path"] for item in record["ev"] if item["sha"] and (not (ROOT / item["path"]).is_file() or sha16(item["path"]) != item["sha"])]
    return (f"{len(record['ev'])} ev files, drift: {', '.join(drift)}" if drift else f"{len(record['ev'])} ev files, no drift"), not drift


# ---------- invariants ----------
ARCHIVE = ROOT / "results/archive_v2_strict"


def inv_archive_raw():
    m = pd.read_csv(ARCHIVE / "MANIFEST.csv")
    bad = [f for f, s in zip(m.file, m.sha256) if hashlib.sha256((ARCHIVE / f).read_bytes()).hexdigest() != s]
    return not bad, f"{len(m) - len(bad)}/{len(m)} byte-identical; differ: {', '.join(bad) or 'none'}"


def inv_archive_lf():
    m = pd.read_csv(ARCHIVE / "MANIFEST.csv")
    bad = [f for f, s in zip(m.file, m.sha256)
           if hashlib.sha256((ARCHIVE / f).read_bytes().replace(b"\r\n", b"\n")).hexdigest() != s
           and hashlib.sha256((ARCHIVE / f).read_bytes()).hexdigest() != s]
    return not bad, f"{len(m) - len(bad)}/{len(m)} match after CRLF->LF; differ: {', '.join(bad) or 'none'}"


def inv_matrix():
    m, ix = npy("results/rmsd_entry.npy"), csv("results/rmsd_entry_index.csv")
    ok = m.shape == (788, 788) and not np.isnan(m).any() and np.allclose(m, m.T) and len(ix) == 788 and ix.chain_key.is_unique
    return ok, f"shape {m.shape}, NaN {int(np.isnan(m).sum())}, symmetric {np.allclose(m, m.T)}, index rows {len(ix)}, unique {ix.chain_key.nunique()}"


CHAIN_TABLES = ["results/entry_clusters.csv", "results/rmsd_entry_index.csv", "results/alignment_core_rmsd.csv", "results/h12_metric.csv",
                "results/helix_bfactor_by_chain.csv", "results/chain_status.csv", "results/qc_outliers.csv",
                "results/alignment/global_alignment.csv", "results/attrition_funnel.csv:skip"]


def inv_chain_tables():
    counts = {t: len(csv(t)) for t in CHAIN_TABLES if not t.endswith(":skip")}
    loi = csv("results/chain_ligand_of_interest.csv").chain_key.nunique()
    counts["results/chain_ligand_of_interest.csv (unique chain_key)"] = loi
    return all(v >= 788 for v in counts.values()), "; ".join(f"{k} {v}" for k, v in counts.items())


def inv_clustering_changed():
    files = ["entry_clusters.csv", "cluster_stats.csv", "medoids.csv", "qc_clustering_delta.csv"]
    same_bytes = [f for f in files if (ROOT / "results" / f).read_bytes().replace(b"\r\n", b"\n") == (ARCHIVE / f).read_bytes().replace(b"\r\n", b"\n")]
    n = len(csv("results/entry_clusters.csv"))
    return not same_bytes and n == 788, f"identical to archive: {', '.join(same_bytes) or 'none'}; entry_clusters rows {n}"


def inv_subsets():
    s13 = set(csv("results/qc_clustering_delta.csv").subset)
    r06 = set(csv("results/entropy_qc_compare.csv").subset)
    main = len(csv("results/entry_clusters.csv")) == npy("results/rmsd_entry.npy").shape[0] == 788
    want = {"all", "tier_A", "tier_A+B"}
    return want <= s13 and want <= r06 and main, f"S13 subsets {sorted(s13)}; R06 subsets {sorted(r06)}; main file all 788 chains {main}"


def inv_no_row_loss():
    m = pd.read_csv(ARCHIVE / "MANIFEST.csv")
    rows, lost = [], []
    for f, src in zip(m.file, m.source):
        if f.endswith((".csv", ".txt", ".npy")):
            if f.endswith(".csv"):
                old, new = len(pd.read_csv(ARCHIVE / f)), len(pd.read_csv(ROOT / src))
            elif f.endswith(".txt"):
                old, new = len((ARCHIVE / f).read_text().split()), len((ROOT / src).read_text().split())
            else:
                old, new = np.load(ARCHIVE / f).shape[0], np.load(ROOT / src).shape[0]
            rows.append(f"{f} {old}->{new}")
            if new < old:
                lost.append(f)
    return not lost, "; ".join(rows) + f"; lost rows: {', '.join(lost) or 'none'}"


def inv_h12():
    h = csv("results/h12_metric.csv")
    counts = h.h12_label.value_counts().to_dict()
    ok = len(h) == 788 and sum(counts.get(k, 0) for k in ["complete", "partial", "not_modelled"]) == 788
    return ok, f"rows {len(h)}; {counts}"


def inv_esr2():
    entries = set(csv("results/esr2_entries.csv").entry_id)
    out = {"esr2_clusters": set(csv("results/esr2_clusters.csv").entry_id), "esr2_alignment": set(csv("results/esr2_alignment.csv").entry_id),
           "esr2_pca_coords": set(csv("results/esr2_pca_coords.csv").entry_id)}
    missing = {k: sorted(entries - v) for k, v in out.items() if entries - v}
    return not missing, f"{len(entries)} ESR2 entries; missing: {missing or 'none'}"


def inv_validate():
    out = subprocess.run([sys.executable, "tracker/tl.py", "validate"], cwd=ROOT, capture_output=True, text=True).stdout.strip().splitlines()[-1]
    done = [tid for tid, row in latest_states().items() if row["state"] == "DONE"]
    noprov = sorted(t for t in done if not (ROOT / "tracker/prov" / f"{t}.json").exists())
    return " 0 errors" in out, f"{out}; DONE tasks without prov file ({len(noprov)}): {', '.join(noprov) or 'none'}"


INVARIANTS = [
    ("archive hashes equal MANIFEST (raw bytes)", inv_archive_raw),
    ("archive hashes equal MANIFEST (after CRLF->LF)", inv_archive_lf),
    ("rmsd_entry.npy 788x788, no NaN, symmetric; index 788 unique", inv_matrix),
    ("chain-level tables have >= 788 rows", inv_chain_tables),
    ("clustering files differ from archive; entry_clusters 788 rows", inv_clustering_changed),
    ("S13 and R06 show all, tier A and tier A+B; main file uses all chains", inv_subsets),
    ("no archived file lost rows", inv_no_row_loss),
    ("H12 present: 788 rows, labels add up to 788", inv_h12),
    ("every ESR2 entry appears in the ESR2 outputs", inv_esr2),
    ("tl.py validate 0 errors; DONE tasks without prov", inv_validate),
]


def main():
    failures = 0
    for tid in TASKS:
        value, match = recount(tid)
        lines, chk = tl_check(tid)
        ok = match and chk
        failures += not ok
        print(f"{'PASS' if ok else 'FAIL'}  task {tid}: recount {'MATCH' if match else 'MISMATCH'} -> {value}; "
              f"expected {show(TASKS[tid]['expected'])}; tl.py check {'PASS' if chk else 'FAIL'} ({len(lines)} chk)")
    for name, fn in INVARIANTS:
        ok, detail = fn()
        failures += not ok
        print(f"{'PASS' if ok else 'FAIL'}  invariant: {name}: {detail}")
    total = len(TASKS) + len(INVARIANTS)
    print(f"TOTAL: {total - failures}/{total} PASS, {failures} FAIL")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
