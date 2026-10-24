"""Tracker Y22: ESR2 parallel analysis - clustering zoo, outliers, metadata and maps on the 86-chain ESR2 matrix.

Inputs: results/rmsd_esr2.npy + results/rmsd_esr2_index.txt (X02; ESR2 chains fitted on ESR2 reference 3OLL|A|1,
ESR1 numbering), results/esr2_clusters.csv (X03, average linkage), results/esr2_entries.csv (X01 quality),
mmCIFs of the 39 ESR2 entries downloaded from RCSB into data/cache/esr2_cif/ (not git-tracked).
Metadata per entry (same fields as Y02; peptide = any polymer entity of <= peptide_max_length residues that is not
ESR2, read from _entity_poly in the mmCIF).
Zoo: the Y04 run list and the same declared rule (scripts/ext_zoo.py: fit, evaluate), 10-D classical MDS for
vector methods, bootstrap 50 x 80%. Small N (86 chains, 39 entries): results are descriptive.
Outliers: kNN (k=3) robust z > 3.5 and LOF (10 neighbours) > 1.5.
Outputs: results/ext/esr2_metadata.csv, results/ext/esr2_zoo_runs.csv, results/ext/esr2_chain_table.csv,
figures/ext/Y22_esr2_parallel.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_esr2.py
"""
import re
import sys
import warnings
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import LocalOutlierFactor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RED, RES, ROOT, SEED, SETTINGS, VIOLET,
                        mpl_setup, save_fig, write_numbers)
from ext_metadata import crystal_system, first, num, precipitant_flags
from ext_zoo import classical_mds, evaluate, run_list

warnings.filterwarnings("ignore")
CIF = ROOT / "data/cache/esr2_cif"


def meta_entry(path, rc):
    b = gemmi.cif.read(str(path)).sole_block()
    e = path.stem.upper()
    det = first(b, "_exptl_crystal_grow.pdbx_details") or ""
    row = dict(entry_id=e, exp_method=first(b, "_exptl.method") or "not_reported",
               cryst_method=(first(b, "_exptl_crystal_grow.method") or "not_reported").upper(),
               ph=num(first(b, "_exptl_crystal_grow.pH")), temp_K=num(first(b, "_exptl_crystal_grow.temp")),
               space_group=first(b, "_symmetry.space_group_name_H-M") or "not_reported")
    if np.isnan(row["ph"]):
        m = re.search(r"PH\s*=?\s*(\d{1,2}(?:\.\d+)?)", det.upper())
        row["ph"] = float(m.group(1)) if m and 2 <= float(m.group(1)) <= 12 else np.nan
    row["crystal_system"] = crystal_system(row["space_group"] if row["space_group"] != "not_reported" else None)
    row.update({k: v for k, v in precipitant_flags(det).items() if k in ("additive_class", "has_peg")})
    rev = list(b.find_values("_pdbx_audit_revision_history.revision_date"))
    row["release_year"] = int(min(rev)[:4]) if rev else np.nan
    ids = list(b.find_values("_entity_poly.entity_id"))
    seqs = list(b.find_values("_entity_poly.pdbx_seq_one_letter_code_can"))
    ent_ids = list(b.find_values("_entity.id"))
    desc = dict(zip(ent_ids, [gemmi.cif.as_string(x) for x in b.find_values("_entity.pdbx_description")]))
    rec_ent = set(rc[(rc.pdb_id == e) & (rc.identity == "esr2")].entity_id.astype(str))
    peps = [desc.get(i, "?") for i, s in zip(ids, seqs)
            if i not in rec_ent and len(gemmi.cif.as_string(s).replace("\n", "")) <= SETTINGS["peptide_max_length"]]
    row["peptide_present"] = bool(peps)
    row["peptide_names"] = "; ".join(sorted(set(peps))) if peps else "none"
    row["chains_asu"] = int(((rc.pdb_id == e) & (rc.identity == "esr2")).sum())
    row["r_work"] = num(first(b, "_refine.ls_R_factor_R_work"))
    return row


def main():
    d = np.load(ROOT / "results/rmsd_esr2.npy")
    keys = [l.strip() for l in (ROOT / "results/rmsd_esr2_index.txt").read_text().splitlines() if l.strip()]
    cl = pd.read_csv(ROOT / "results/esr2_clusters.csv").set_index("chain_key").loc[keys].reset_index()
    rc = pd.read_csv(ROOT / "results/lbd_extract/receptor_chains.csv", low_memory=False)
    meta = pd.DataFrame([meta_entry(p, rc) for p in sorted(CIF.glob("*.cif"))])
    q = pd.read_csv(ROOT / "results/esr2_entries.csv")[["entry_id", "resolution", "r_free", "mean_b", "qc_pass", "ligand_ids"]]
    meta = meta.merge(q, on="entry_id", how="left")
    meta.to_csv(RES / "esr2_metadata.csv", index=False)

    x = classical_mds(d, 10)
    runs = run_list()
    res = Parallel(n_jobs=-1)(delayed(evaluate)(i, f, p, d, x) for i, (f, p) in enumerate(runs))
    df = pd.DataFrame([r for r, _ in res])
    labs = {r["run_id"]: l for r, l in res}
    ref = cl["cut_1.5"].values
    df["ari_vs_x03"] = [adjusted_rand_score(ref, labs[r]) for r in df.run_id]
    elig = (df.n_clusters >= 2) & (df.noise_fraction <= 0.10) & df.silhouette.notna()
    e = df[elig].copy()
    e["score"] = (e.silhouette.rank(pct=True) + e.stability_ari.rank(pct=True)) / 2
    df = df.merge(e[["run_id", "score"]], on="run_id", how="left").sort_values("score", ascending=False)
    df.to_csv(RES / "esr2_zoo_runs.csv", index=False)
    best = df.iloc[0]

    dd = d.copy(); np.fill_diagonal(dd, np.inf)
    knn = np.sort(dd, 1)[:, :3].mean(1)
    z = (knn - np.median(knn)) / (1.4826 * np.median(np.abs(knn - np.median(knn))))
    lof = -LocalOutlierFactor(n_neighbors=10, metric="precomputed").fit(d).negative_outlier_factor_
    tab = cl.rename(columns={"cut_1.5": "x03_cluster"})[["chain_key", "entry_id", "x03_cluster"]].copy()
    tab["best_zoo_cluster"] = labs[best.run_id]
    tab["knn3_A"], tab["knn_flag"], tab["lof"], tab["lof_flag"] = knn, z > 3.5, lof, lof > 1.5
    tab = tab.merge(meta[["entry_id", "peptide_present", "space_group", "resolution", "r_free", "ph", "release_year"]],
                    on="entry_id", how="left")
    emb = MDS(2, dissimilarity="precomputed", random_state=SEED, n_init=4, normalized_stress="auto").fit_transform(d)
    tab["mds1"], tab["mds2"] = emb[:, 0], emb[:, 1]
    tab.to_csv(RES / "esr2_chain_table.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    sizes = pd.Series(ref).value_counts()
    main_c = sizes.index[:2]
    ccol = {main_c[0]: VIOLET, main_c[1]: CATEGORICAL[3]}
    ax = axs[0, 0]
    for c in sorted(set(ref)):
        m = ref == c
        ax.scatter(emb[m, 0], emb[m, 1], s=26 if c in ccol else 40, c=ccol.get(c, GREY),
                   marker="o" if c in ccol else "^", label=f"X03 cluster {c} (n={m.sum()})", lw=0.4, edgecolor="white")
    ax.legend(fontsize=7.5)
    ax.set(xlabel="MDS 1 (A)", ylabel="MDS 2 (A)", title="ESR2 (86 chains): X03 average-linkage clusters at 1.5 A")
    ax = axs[0, 1]
    ax.scatter(df.silhouette, df.stability_ari, s=16, c=GREY, lw=0, label=f"all runs (n={len(df)})")
    ax.scatter([best.silhouette], [best.stability_ari], s=110, facecolor="none", edgecolor=INK, lw=1.5,
               label=f"selected: {best.run_id}")
    ax.set(xlabel="silhouette (precomputed RMSD)", ylabel="bootstrap stability (ARI)", title="ESR2 clustering zoo")
    ax.legend(fontsize=7.5, loc="lower left")
    ax = axs[0, 2]
    pp = tab.peptide_present.fillna(False).astype(bool).values
    ax.scatter(emb[~pp, 0], emb[~pp, 1], s=22, c=GREY, label=f"no peptide (n={(~pp).sum()})", lw=0)
    ax.scatter(emb[pp, 0], emb[pp, 1], s=22, c=CATEGORICAL[6], label=f"peptide (n={pp.sum()})", lw=0)
    ax.legend(fontsize=8)
    ax.set(xlabel="MDS 1 (A)", ylabel="MDS 2 (A)", title="Coloured by co-crystallised peptide")
    ax = axs[1, 0]
    sc = ax.scatter(emb[:, 0], emb[:, 1], s=22, c=tab.resolution, cmap="esr_seq", lw=0)
    fig.colorbar(sc, ax=ax, label="resolution (A)")
    ax.set(xlabel="MDS 1 (A)", ylabel="MDS 2 (A)", title="Coloured by resolution")
    ax = axs[1, 1]
    fl = (tab.knn_flag | tab.lof_flag).values
    ax.scatter(emb[~fl, 0], emb[~fl, 1], s=18, c=GREY, lw=0, label="not flagged")
    ax.scatter(emb[fl, 0], emb[fl, 1], s=40, c=RED, lw=0.4, edgecolor="white", label=f"kNN or LOF outlier (n={fl.sum()})")
    for r in tab[fl].itertuples():
        ax.annotate(r.chain_key.split("|")[0], (r.mds1, r.mds2), fontsize=7, color=INK2, xytext=(3, 3), textcoords="offset points")
    ax.legend(fontsize=8)
    ax.set(xlabel="MDS 1 (A)", ylabel="MDS 2 (A)", title="Outliers (kNN k=3 robust z > 3.5, LOF > 1.5)")
    ax = axs[1, 2]
    vc = meta.space_group.value_counts().head(8)
    ax.barh(vc.index[::-1], vc.values[::-1], color=VIOLET)
    ax.set(xlabel="ESR2 entries", title=f"{meta.space_group.nunique()} space groups in {len(meta)} entries")
    ari = adjusted_rand_score(ref, labs[best.run_id])
    big = pd.Series(ref == main_c[0])
    pep_main = tab.peptide_present.fillna(False).astype(bool)[big.values].mean()
    pep_second = tab.peptide_present.fillna(False).astype(bool)[(ref == main_c[1])].mean()
    fig.suptitle(f"ESR2 in parallel: main basin (n={sizes.iloc[0]}, {pep_main:.0%} with peptide) vs second basin "
                 f"(n={sizes.iloc[1]}, {pep_second:.0%} with peptide), the same split as ESR1; rule winner {best.run_id} "
                 f"is degenerate (largest cluster {best.largest_fraction:.0%})",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y22_esr2_parallel", section="I. ESR2 and ESR1 vs ESR2", task="Y22",
             shows="ESR2 metric-MDS map coloured by cluster, peptide, resolution, outliers; ESR2 clustering zoo; space groups",
             takeaway=(f"ESR2 shows the ESR1 pattern: the main basin ({sizes.iloc[0]} chains) is peptide-bound, the second "
                       f"({sizes.iloc[1]}) is not; {fl.sum()} chains are kNN/LOF outliers; small N, descriptive only."),
             source="results/ext/esr2_zoo_runs.csv, results/ext/esr2_metadata.csv, results/ext/esr2_chain_table.csv",
             script="scripts/ext_esr2.py")
    write_numbers("Y22", dict(n_chains=len(keys), n_entries=len(meta), best_run=best.run_id, best_sil=float(best.silhouette),
                              best_stab=float(best.stability_ari), best_ari_x03=float(ari),
                              n_peptide_entries=int(meta.peptide_present.sum()), n_outliers=int(fl.sum()),
                              peptide_frac_main=float(pep_main), peptide_frac_second=float(pep_second),
                              basin_sizes=sizes.to_dict()))
    print(df.head(8)[["run_id", "n_clusters", "silhouette", "stability_ari", "ari_vs_x03", "score"]])
    print(meta[["entry_id", "peptide_present", "peptide_names", "space_group", "ph"]].to_string())


if __name__ == "__main__":
    main()
