"""Tracker Y40: independent recounts for the final-run tasks (Y..) and the tracker checks; writes docs/EVIDENCE_EXT.md.

Read-only on results: it never imports or runs a producing ext_* script. For every Y task logged DONE, one recount
expression (pandas / numpy / scipy / sklearn on the result files, several of them recomputed from the frozen inputs,
e.g. re-clustering the saved RMSD matrix) is compared with the claim written in the task's log note. Floats match
when they agree to the number of decimals written in the claim. Also runs `tl.py check` for the task.
Prints PASS / FAIL per task and writes docs/EVIDENCE_EXT.md. Exit 0 only if every recount and check passes.
Run from the repo root: python scripts/check_evidence_ext.py
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score as ari

ROOT = Path(__file__).resolve().parents[1]


def csv(p, **kw):
    return pd.read_csv(ROOT / p, low_memory=False, **kw)


def npy(p):
    return np.load(ROOT / p)


def png_ok(p):
    a = np.asarray(Image.open(ROOT / p).convert("L"))
    return bool(a.std() > 2)


def S08():
    return csv("results/entry_clusters.csv")["cut_1.5"].values


def peptide_by_cluster():
    cm = csv("results/ext/chain_metadata.csv").merge(csv("results/entry_clusters.csv")[["chain_key", "cut_1.5"]], on="chain_key")
    return int(cm[cm["cut_1.5"] == 7].peptide_present.sum()), int(cm[cm["cut_1.5"] == 4].peptide_present.sum())


def triangle_frac():
    d = npy("results/rmsd_entry.npy")
    rng = np.random.default_rng(42)
    i, j, k = rng.integers(0, len(d), (3, 200_000))
    ok = (i != j) & (j != k) & (i != k)
    i, j, k = i[ok], j[ok], k[ok]
    return float(((d[i, k] - d[i, j] - d[j, k]) > 1e-6).mean())


def hill_q1(lab):
    p = pd.Series(lab).value_counts(normalize=True).values
    return float(np.exp(-(p * np.log(p)).sum()))


def cover_at(r, frac):
    d = npy("results/rmsd_entry.npy")
    comp = csv("results/h12_metric.csv").set_index("chain_key").loc[csv("results/rmsd_entry_index.csv").chain_key, "h12_label"].values == "complete"
    within = d[:, comp] <= r
    covered = np.zeros(len(d), bool)
    n = 0
    while covered.mean() < frac:
        g = within[~covered].sum(0)
        k = int(np.argmax(g))
        if g[k] == 0:
            return None
        covered |= within[:, k]
        n += 1
    return n


H = dict(pd=pd, np=np, csv=csv, npy=npy, png_ok=png_ok, S08=S08, ari=ari, fcluster=fcluster, linkage=linkage,
         squareform=squareform, peptide_by_cluster=peptide_by_cluster, triangle_frac=triangle_frac, hill_q1=hill_q1,
         cover_at=cover_at)

# id: (plain name, claim from the log note, recount expression, expected)
TASKS = {
    "Y02": ("Experimental metadata from mmCIF", "351 entries; pH 307/351; 323 PEG; 17 space groups; peptide in 210 entries; cl7 414 / cl4 4 chains with peptide",
            "(lambda m: [len(m), int(m.ph.notna().sum()), int((m.additive_class == 'PEG').sum()), m.space_group.nunique(), int(m.peptide_present.sum())] + list(peptide_by_cluster()))(csv('results/ext/entry_metadata.csv'))",
            [351, 307, 323, 17, 210, 414, 4]),
    "Y03": ("Shared embeddings", "788 chains in every embedding; PCA 2-D 98.3% variance",
            "(lambda e: [len(e), int(e[['pca1','mds1','tsne1','umap1']].notna().all(axis=1).sum())])(csv('results/ext/embedding_rmsd.csv')) + [round(float(csv('results/ext/embedding_quality.csv').set_index('embedding').loc['pca','explained_2d']), 3)]",
            [788, 788, 0.983]),
    "Y04": ("Structural clustering zoo", "81 runs in 14 families; winner spectral sigma=2 K=2; ARI vs S08 0.94",
            "(lambda z, l: [len(z), z.method.nunique(), csv('results/ext/zoo_ranking.csv').run_id.iloc[0], round(ari(S08(), l['spectral|sigma=2.0,k=2']), 2)])(csv('results/ext/zoo_runs.csv'), csv('results/ext/zoo_labels.csv'))",
            [81, 14, "spectral|sigma=2.0,k=2", 0.94]),
    "Y05": ("Cluster compactness", "cl7 mean intra 0.70 A, cl4 0.91 A (recomputed from the matrix)",
            "(lambda d, c: [round(float(d[np.ix_(c == k, c == k)][np.triu_indices((c == k).sum(), 1)].mean()), 2) for k in (7, 4)])(npy('results/rmsd_entry.npy'), S08())",
            [0.70, 0.91]),
    "Y06": ("Method agreement", "selected spectral K=2 vs S08: 324/329 cl4 and 444/445 cl7 together",
            "(lambda l, c: sorted([int(pd.crosstab(c, l).loc[4].max()), int(pd.crosstab(c, l).loc[7].max())]))(csv('results/ext/zoo_labels.csv')['spectral|sigma=2.0,k=2'].values, S08())",
            [324, 444]),
    "Y07": ("Dendrogram and annotated heatmap", "788 leaves; both PNGs non-blank",
            "[len(csv('results/ext/leaf_order.csv')), png_ok('figures/ext/Y07_dendrogram_metadata.png'), png_ok('figures/ext/Y07_rmsd_heatmap_annotated.png')]",
            [788, True, True]),
    "Y08": ("Chemistry zoo", "276 runs; auto-selected MACCS k-means K=6; rule winner largest cluster 99%",
            "(lambda r, b: [len(r), int(b.chem_cluster_auto.nunique()), round(float(r.set_index('run_id').loc['MACCS|hier_average|k=2', 'largest_share']), 2)])(csv('results/ext/chem_zoo_runs.csv'), csv('results/ext/chem_labels_best.csv'))",
            [276, 6, 0.99]),
    "Y09": ("Ligand descriptors", "basic amine 1% of cl7-only vs 74% of cl4-only ligands; 169 scaffolds",
            "(lambda d: [round(float(d[d.structural_majority == 'cl7'].basic_amine.mean()) * 100), round(float(d[d.structural_majority == 'cl4'].basic_amine.mean()) * 100), d.murcko.nunique()])(csv('results/ext/ligand_descriptors.csv'))",
            [1, 74, 169]),
    "Y10": ("Displacement profile", "528-550 carry 99% of the squared difference; complete medoids 5TM1|A|1 / 7UJ7|C|1",
            "(lambda p: [round(float((p[p.position >= 528].d_mean_A ** 2).sum() / (p.d_mean_A ** 2).sum()), 2), p.medoid_cl7.iloc[0], p.medoid_cl4.iloc[0]])(csv('results/ext/cluster_displacement_profile.csv'))",
            [0.99, "5TM1|A|1", "7UJ7|C|1"]),
    "Y11": ("Structural features", "median K362-E542 19.3 vs 13.6 A; H12-H3 angle 82 vs 43 deg",
            "(lambda f: [round(float(f[f.cluster == 7].k362_e542_A.median()), 1), round(float(f[f.cluster == 4].k362_e542_A.median()), 1), round(float(f[f.cluster == 7].h12_h3_angle_deg.median())), round(float(f[f.cluster == 4].h12_h3_angle_deg.median()))])(csv('results/ext/structural_features.csv'))",
            [19.3, 13.6, 82, 43]),
    "Y12": ("phi/psi PCA", "190 positions; 81% of dPC1 loading at residue 537",
            "(lambda l: [len(l), round(float(l.set_index('position').load1.loc[537] / l.load1.sum()), 2)])(csv('results/ext/phipsi_loadings.csv'))",
            [190, 0.81]),
    "Y13": ("Supervised check", "RF balanced accuracy quality 0.80, construct 0.96, context 0.98",
            "(lambda s: [round(float(s.loc[(k, 'random_forest'), 'balanced_accuracy_mean']), 2) for k in ('quality', 'construct', 'context')])(csv('results/ext/supervised_scores.csv').set_index(['feature_set', 'model']))",
            [0.80, 0.96, 0.98]),
    "Y14": ("Hypothesis tests", "entry-level log OR peptide 6.7, Y537S 4.2, L536S -6.7",
            "(lambda s: [round(float(s.loc[(h, 'entry'), 'effect']), 1) for h in ('coactivator/other peptide present', 'Y537S construct', 'L536S construct')])(csv('results/ext/split_hypotheses.csv').set_index(['hypothesis', 'level']))",
            [6.7, 4.2, -6.7]),
    "Y15": ("Outlier detection", "kNN 13, medoid 13, LOF 42, iforest 94, Mahalanobis 339; 22 consensus",
            "(lambda o: [int(o[c].sum()) for c in ('knn_flag', 'medoid_flag', 'lof_flag', 'iforest_flag', 'mahal_flag')] + [int((o[['knn_flag', 'medoid_flag', 'lof_flag', 'iforest_flag', 'mahal_flag']].sum(axis=1) >= 3).sum())])(csv('results/ext/outlier_flags.csv'))",
            [13, 13, 42, 94, 339, 22]),
    "Y16": ("Outliers vs quality", "6 quality / 16 conformation verdicts",
            "(lambda v: [int((v.verdict == 'quality').sum()), int((v.verdict == 'conformation').sum())])(csv('results/ext/outlier_verdicts.csv'))",
            [6, 16]),
    "Y17": ("Metadata tests", "100 tests; entry-level space group V 0.73 vs cluster",
            "(lambda t: [len(t), round(float(t[(t.variable == 'space_group8') & (t.target == 'structural cluster') & (t.level == 'entry')].effect.iloc[0]), 2)])(csv('results/ext/metadata_tests.csv'))",
            [100, 0.73]),
    "Y18": ("Colour-by small multiples", "4 figures, 788 chains in the table",
            "[len(csv('results/ext/colour_by_table.csv'))] + [png_ok(f'figures/ext/Y18_colour_by_{k}.png') for k in ('quality', 'experiment', 'biology', 'pca')]",
            [788, True, True, True, True]),
    "Y19": ("Interactive embedding", "self-contained HTML (plotly inside, > 3 MB)",
            "[(ROOT / 'docs/interactive_embedding.html').stat().st_size > 3_000_000, 'plotly' in (ROOT / 'docs/interactive_embedding.html').read_text(encoding='utf8')[:200000].lower()]",
            [True, True]),
    "Y20": ("Coordinate PCA", "196 positions; mode 1 86.5%",
            "[len(csv('results/ext/coord_pca_loadings.csv')), round(float(csv('results/ext/coord_pca_variance.csv').explained_variance_ratio.iloc[0]) * 100, 1)]",
            [196, 86.5]),
    "Y21": ("B-factors and core", "per-chain normalised B core v3 0.928 vs old core 1.038",
            "(lambda c: [round(float(c.core_v3_norm_b.median()), 3), round(float(c.old_core_norm_b.median()), 3)])(csv('results/ext/core_b_v3.csv'))",
            [0.928, 1.038]),
    "Y22": ("ESR2 parallel", "86 chains; 39 entries; main basin 93% peptide, second 0%",
            "(lambda t: [len(t), t.entry_id.nunique(), round(float(t[t.x03_cluster == 2].peptide_present.mean()), 2), round(float(t[t.x03_cluster == 5].peptide_present.mean()), 2)])(csv('results/ext/esr2_chain_table.csv'))",
            [86, 39, 0.93, 0.0]),
    "Y23": ("Joint ESR1+ESR2", "874 x 874; ESR1 block equals S04; ESR1 cl7 - ESR2 main median 1.33 A",
            "(lambda d, i: [d.shape[0], bool(np.allclose(d[:788, :788], npy('results/rmsd_entry.npy'), atol=1e-9)), round(float(np.median(d[:788][S08() == 7][:, 788:][:, (i.cluster.values[788:] == 'ESR2 X03-2')])), 2)])(npy('results/ext/rmsd_joint.npy'), csv('results/ext/joint_index.csv'))",
            [874, True, 1.33]),
    "Y24": ("ESR1 vs ESR2 potency", "2683 ESR2 rows; 1615 compounds at both; 54% within 10-fold",
            "(lambda s: [len(csv('results/ext/chembl_esr2_activities.csv')), len(s), round(float((s.selectivity.abs() <= 1).mean()), 2)])(csv('results/ext/selectivity_ratios.csv'))",
            [2683, 1615, 0.54]),
    "Y25": ("Ligand-conformation link", "304 ligands; chi2 170; nearest-neighbour agreement 0.89",
            "(lambda r, n: [int(r.n.iloc[0]), round(float(r.statistic.iloc[0])), round(float(n.agree.mean()), 2)])(csv('results/ext/chem_struct_link.csv'), csv('results/ext/nn_agreement.csv'))",
            [304, 170, 0.89]),
    "Y26": ("Redundancy and diversity", "q1 2.19 at 1.5 A; 2 complete-H12 reps cover 90% at 1.0 A, 48 at 0.5 A (recomputed)",
            "[round(hill_q1(S08()), 2), cover_at(1.0, 0.9), cover_at(0.5, 0.9)]",
            [2.19, 2, 48]),
    "Y27": ("Time trends", "peak 2016 with 149 chains; 3 of 673 chains since 2010 novel > 1.5 A",
            "(lambda t, n: [int(t.loc[t.n_chains.idxmax(), 'release_year']), int(t.n_chains.max()), int((n[n.release_year >= 2010].novelty_A > 1.5).sum()), int((n.release_year >= 2010).sum())])(csv('results/ext/time_trends.csv'), csv('results/ext/novelty_by_chain.csv'))",
            [2016, 149, 3, 673]),
    "Y28": ("Robustness grid", "120 cells; agreement 0.98-1.00 on both real cores",
            "(lambda g: [len(g), round(float(g[g.core.isin(['core v3', 'old core'])].main_split_agreement.min()), 2)])(csv('results/ext/robustness_grid.csv'))",
            [120, 0.98]),
    "Y29": ("ESR2 core fraction", "3 policy candidates + ESR1 frame all ARI 1.00 vs X03",
            "(lambda a: [len(a), round(float(a.ari_vs_x03.min()), 2)])(csv('results/ext/esr2_reference_alternatives.csv'))",
            [4, 1.0]),
    "Y30": ("PyMOL renders", "4 raw PyMOL panels + 2 composites, non-blank",
            "[png_ok(f'figures/ext/pymol_panels/{p}.png') for p in ('medoids_front', 'medoids_side', 'bputty_cl7', 'bputty_cl4')] + [png_ok('figures/ext/Y30_pymol_medoids.png'), png_ok('figures/ext/Y30_pymol_bputty.png')]",
            [True] * 6),
    "Y31": ("Upstream rebuild", "11 identical, 7 equivalent, 3 differs",
            "(lambda r: [int((r.status == s).sum()) for s in ('identical', 'equivalent', 'differs')])(csv('results/rebuild/rebuild_report.csv'))",
            [11, 7, 3]),
    "Y42": ("Metric check", "triangle violations 5.4%; full recomputation equals S04; H12-only median cl7-cl4 14.0 A",
            "(lambda m: [round(triangle_frac() * 100, 1), float(m.loc['full', 'max_abs_diff_vs_s04']), round(float(m.loc['h12only', 'median_cl7_cl4_A']), 1)])(csv('results/ext/metric_check.csv').set_index('matrix'))",
            [5.4, 0.0, 14.0]),
    "Y43": ("NCS copies", "317 multi-chain entries; 5 split across clusters",
            "(lambda n: [len(n), int(n.mixed_clusters.sum())])(csv('results/ext/ncs_copies.csv'))",
            [317, 5]),
    "Y44": ("H12 crystal contacts", "median contacts cl7 1 vs cl4 10 atoms; 788 chains",
            "(lambda h: [len(h), float(h[h.cluster == 7].h12_contact_atoms.median()), float(h[h.cluster == 4].h12_contact_atoms.median())])(csv('results/ext/h12_crystal_contacts.csv'))",
            [788, 1.0, 10.0]),
    "Y45": ("cl4 substructure", "cl4-1 262 chains (L536S 84%), cl4-2 48 (L536S 0%)",
            "(lambda s: [int((s.sub_cluster == 1).sum()), round(float((s[s.sub_cluster == 1].mutation_status == 'L536S').mean()), 2), int((s.sub_cluster == 2).sum()), round(float((s[s.sub_cluster == 2].mutation_status == 'L536S').mean()), 2)])(csv('results/ext/cl4_subclusters.csv'))",
            [262, 0.84, 48, 0.0]),
}


def match(got, exp):
    if isinstance(exp, list):
        return isinstance(got, (list, tuple)) and len(got) == len(exp) and all(match(g, e) for g, e in zip(got, exp))
    if isinstance(exp, float):
        dec = len(repr(exp).split(".")[1]) if "." in repr(exp) else 0
        try:
            return round(float(got), dec) == round(exp, dec)
        except (TypeError, ValueError):
            return False
    return got == exp


def main():
    log = (ROOT / "tracker/esr1.log").read_text(encoding="utf8").splitlines()
    last = {}
    for line in log:
        p = line.split(" | ")
        if len(p) > 2 and p[1] in TASKS:
            last[p[1]] = p
    rows, ok_all = [], True
    for tid, (name, claim, expr, exp) in TASKS.items():
        st = last.get(tid)
        state = st[2] if st else "no log line"
        when = st[0] if st else ""
        ev = re.search(r"ev=([^|]+)", " | ".join(st)).group(1).strip() if st and "ev=" in " | ".join(st) else ""
        try:
            got = eval(expr, {**H, "ROOT": ROOT})
            rec = match(got, exp)
        except Exception as e:
            got, rec = f"error {e}", False
        chk = subprocess.run([sys.executable, str(ROOT / "tracker/tl.py"), "check", tid], capture_output=True, text=True, cwd=ROOT)
        chk_ok = chk.returncode == 0 and "FAIL" not in chk.stdout
        ok = rec and chk_ok and state == "DONE"
        ok_all &= ok
        rows.append((tid, name, state, when, ev, claim, expr, exp, got, "PASS" if chk_ok else "FAIL", "MATCH" if rec else "MISMATCH"))
        print(f"{tid} {state} check={'PASS' if chk_ok else 'FAIL'} recount={'MATCH' if rec else 'MISMATCH'} got={got}")
    md = ["# Evidence for the final-run tasks (Y01-Y45)", "",
          "Generated by `python scripts/check_evidence_ext.py` (read-only; it never runs a producing script). For each task: "
          "the claim from its log note, an independent recount from the result files (several recomputed from the frozen "
          "inputs, e.g. Y05, Y23, Y26, Y42), the `tl.py check` result. Floats match to the decimals written in the claim.", "",
          f"Tasks in scope: {len(rows)}. All pass: {'yes' if ok_all else 'NO'}.", "",
          "| id | plain name | state | log time | tl.py check | recount | claim | recount value |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[9]} | {r[10]} | {r[5]} | {str(r[8]).replace('|', '/')} |")
    md += ["", "## Recount expressions", ""]
    for r in rows:
        md += [f"- **{r[0]}** expected `{r[7]}`; evidence: {r[4] or '-'}", f"  `{r[6]}`"]
    (ROOT / "docs/EVIDENCE_EXT.md").write_text("\n".join(md) + "\n", encoding="utf8")
    print("ALL PASS" if ok_all else "SOME FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
