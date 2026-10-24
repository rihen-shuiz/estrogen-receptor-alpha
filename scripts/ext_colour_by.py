"""Tracker Y18: one shared embedding drawn as small multiples, coloured by every variable of interest.

Map: metric MDS of the 788 x 788 RMSD matrix (results/ext/embedding_rmsd.csv, Y03; chosen by the declared rule), and
a PCA version for the key variables. Point = chain; every chain is drawn (missing value = grey).
Colours: structural clusters use the shared palette (blue cl7, orange cl4, aqua minor); other categorical variables
use palette slots 3-8 and grey (6 most frequent levels, the rest 'other'), so blue/orange always mean cl7/cl4;
continuous variables use the one-hue blue ramp with a colour bar.
Ligand action: only the 7 evidence-supported labels of block 07 (results/task_outputs/block_07/ligand_actions.csv,
supported_action_type) are shown; every other ligand is 'unlabelled' (Q10/Q11 labels are made by hand by a teammate).
Outputs: figures/ext/Y18_colour_by_quality.png, Y18_colour_by_experiment.png, Y18_colour_by_biology.png,
Y18_colour_by_pca.png; the merged per-chain table results/ext/colour_by_table.csv.
Run from the repo root: PYTHONPATH=. python scripts/ext_colour_by.py
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, CLUSTER_COLOR, CLUSTER_MARKER, GREY, INK2, ORANGE, RES, ROOT,
                        chain_table, mpl_setup, save_fig, scatter_clusters, write_numbers)

OTHER = [AQUA, CATEGORICAL[3], CATEGORICAL[4], CATEGORICAL[5], CATEGORICAL[6], CATEGORICAL[7]]


def table():
    t = chain_table()
    emb = pd.read_csv(RES / "embedding_rmsd.csv").drop(columns=["entry_id"])
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status", "peptide_class"]]
    em = pd.read_csv(RES / "entry_metadata.csv")[["entry_id", "ph", "additive_class", "crystal_system", "release_year",
                                                   "cryst_method", "chains_asu", "r_work", "peg_mw"]]
    o = pd.read_csv(RES / "outlier_flags.csv")[["chain_key", "n_methods"]]
    sf = pd.read_csv(RES / "structural_features.csv")[["chain_key", "h12_shift_A", "h12_h3_angle_deg", "k362_e542_A",
                                                        "pocket_sas_points", "ligand_heavy_atoms", "ligand_basic_amine"]]
    ch = pd.read_csv(RES / "chem_labels_best.csv")[["ligand_id", "chem_cluster_auto"]]
    pp = pd.read_csv(RES / "phipsi_pca.csv")[["chain_key", "dpc1"]]
    act = pd.read_csv(ROOT / "results/task_outputs/block_07/ligand_actions.csv")[["ligand_code", "supported_action_type"]]
    act = act.dropna().rename(columns={"ligand_code": "ligand_id", "supported_action_type": "action_supported"})
    t = (t.merge(emb, on="chain_key").merge(cm, on="chain_key").merge(em, on="entry_id", how="left")
         .merge(o, on="chain_key").merge(sf, on="chain_key").merge(ch, on="ligand_id", how="left")
         .merge(pp, on="chain_key").merge(act, on="ligand_id", how="left"))
    t["action_supported"] = t.action_supported.fillna("unlabelled (Q10 pending)")
    t["chem_cluster_auto"] = t.chem_cluster_auto.map(lambda v: "no ligand" if pd.isna(v) else f"chem {int(v)}")
    t["ligand_basic_amine"] = t.ligand_basic_amine.map({True: "basic amine", False: "no basic amine"}).fillna("no ligand")
    t["year"] = t.deposit_date.astype(str).str[:4].astype(float)
    return t


def panel(ax, fig, x, y, v, name, kind, order=None, fixed=None):
    v = pd.Series(v).reset_index(drop=True)
    if kind == "cluster":
        scatter_clusters(ax, x, y, v, s=6, legend=True)
    elif kind == "cat":
        s = v.fillna("not_reported").astype(str).replace({"nan": "not_reported"})
        levels = order or [l for l in s.value_counts().index if l not in ("not_reported", "other")]
        top = levels[:6]
        s = s.where(s.isin(top) | s.isin(["not_reported"]), "other")
        cmap = {l: (fixed or {}).get(l, OTHER[i % 6]) for i, l in enumerate(top)}
        cmap.update({"other": "#d4d3cd", "not_reported": GREY})
        for l in ["not_reported", "other"] + top[::-1]:
            m = (s == l).values
            if m.any():
                ax.scatter(np.asarray(x)[m], np.asarray(y)[m], s=6 if l in ("other", "not_reported") else 8,
                           c=cmap[l], lw=0, label=f"{l} ({m.sum()})")
        h, lab = ax.get_legend_handles_labels()
        ax.legend(h[::-1], lab[::-1], fontsize=6.5, markerscale=1.6, loc="best", handletextpad=0.2)
    else:
        vals = pd.to_numeric(v, errors="coerce").values
        m = np.isfinite(vals)
        ax.scatter(np.asarray(x)[~m], np.asarray(y)[~m], s=5, c=GREY, lw=0)
        lo, hi = np.nanpercentile(vals, [2, 98])
        sc = ax.scatter(np.asarray(x)[m], np.asarray(y)[m], s=7, c=vals[m], cmap="esr_seq", vmin=lo, vmax=hi, lw=0)
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.01)
        cb.ax.tick_params(labelsize=7)
    n_na = int(pd.isna(v).sum()) if kind == "num" else 0
    ax.set_title(name + (f" (grey = missing, n={n_na})" if n_na else ""), fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])


def grid(t, specs, title, xcol="mds1", ycol="mds2", ncol=4):
    plt = mpl_setup()
    nrow = int(np.ceil(len(specs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.1 * nrow))
    axs = np.atleast_1d(axs).ravel()
    for ax, sp in zip(axs, specs):
        panel(ax, fig, t[xcol], t[ycol], t[sp[0]], sp[1], sp[2], *(sp[3:] if len(sp) > 3 else []))
    for ax in axs[len(specs):]:
        ax.axis("off")
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold", fontsize=12)
    fig.tight_layout()
    return fig


def main():
    t = table()
    t.to_csv(RES / "colour_by_table.csv", index=False)
    yes_no = {"True": CATEGORICAL[6], "False": GREY}
    q = [("cluster", "S08 structural cluster (1.5 A)", "cluster"), ("r_free", "R-free", "num"),
         ("resolution_A", "resolution (A)", "num"), ("mean_b_A2", "entry mean B (A^2)", "num"),
         ("mean_b_chain", "chain mean C-alpha B (A^2)", "num"),
         ("b_group", "B-group (tertile of chain mean B)", "cat", ["low", "mid", "high"],
          {"low": CATEGORICAL[5], "mid": CATEGORICAL[3], "high": CATEGORICAL[7]}),
         ("tier", "QC tier", "cat", ["A", "B", "C"], {"A": CATEGORICAL[5], "B": CATEGORICAL[3], "C": CATEGORICAL[7]}),
         ("h12_label", "H12 label", "cat", ["complete", "partial", "not_modelled"],
          {"complete": CATEGORICAL[6], "partial": CATEGORICAL[3], "not_modelled": CATEGORICAL[7]}),
         ("n_methods", "outlier detectors flagging the chain (0-5)", "num"),
         ("r_work", "R-work", "num"), ("mapped_lbd_ca", "observed LBD C-alpha", "num"),
         ("core_rmsd", "core RMSD to reference (A)", "num")]
    fig = grid(t, q, "Colour by quality: R-free, resolution and B-factors are spread over both clusters; only H12 "
                     "completeness follows the split (shared metric-MDS map, 788 chains)")
    save_fig(fig, "Y18_colour_by_quality", section="F. Embeddings and colour-by", task="Y18",
             shows="metric-MDS map of the RMSD matrix coloured by cluster, R-free, resolution, B-factors, B-group, tier, H12 label, outlier count, R-work, completeness, core RMSD",
             takeaway="Quality variables are mixed across the map; cl4 even holds more high-resolution chains; H12 partial/not modelled concentrates in cl4.",
             source="results/ext/colour_by_table.csv", script="scripts/ext_colour_by.py")
    e = [("cluster", "S08 structural cluster", "cluster"), ("space_group", "space group", "cat"),
         ("crystal_system", "crystal system", "cat"), ("ph", "crystallisation pH", "num"),
         ("additive_class", "main precipitant class", "cat"), ("peptide_class", "co-crystallised peptide", "cat", None, {"none": "#d4d3cd"}),
         ("mutation_status", "construct mutation", "cat"), ("year", "deposition year", "num"),
         ("cryst_method", "crystallisation method", "cat"), ("chains_asu", "ESR1 chains per ASU", "num"),
         ("peg_mw", "PEG molecular weight (g/mol)", "num"), ("release_year", "release year", "num")]
    fig = grid(t, e, "Colour by experiment: peptide, construct mutation and space group paint the two clusters; "
                     "precipitant class does not")
    save_fig(fig, "Y18_colour_by_experiment", section="F. Embeddings and colour-by", task="Y18",
             shows="the same map coloured by space group, crystal system, pH, precipitant, peptide, mutation, year, method, ASU, PEG size",
             takeaway="P21 / NCOA2 peptide / Y537S cover cl7; P1, C2, P3x and L536S cover cl4; PEG is everywhere (no signal).",
             source="results/ext/colour_by_table.csv", script="scripts/ext_colour_by.py")
    b = [("cluster", "S08 structural cluster", "cluster"), ("chem_cluster_auto", "ligand chemistry cluster (Y08 auto)", "cat"),
         ("ligand_basic_amine", "ligand basic amine", "cat", ["basic amine", "no basic amine", "no ligand"],
          {"basic amine": CATEGORICAL[6], "no basic amine": CATEGORICAL[3], "no ligand": GREY}),
         ("ligand_heavy_atoms", "ligand heavy atoms", "num"),
         ("action_supported", "action label (7 evidence-supported only)", "cat", None,
          {"unlabelled (Q10 pending)": "#d4d3cd"}),
         ("pocket_sas_points", "P2Rank pocket size (SAS points)", "num"), ("h12_shift_A", "H12 shift from cl7 mean (A)", "num"),
         ("h12_h3_angle_deg", "H12-H3 axis angle (deg)", "num"), ("k362_e542_A", "K362-E542 C-alpha (A)", "num"),
         ("dpc1", "dihedral PC1 (psi of 537)", "num"), ("pocket_score", "P2Rank pocket score", "num"),
         ("h12_to_core_ratio", "H12 / core B ratio", "num")]
    fig = grid(t, b, "Colour by biology: ligand chemistry, H12 geometry and the 537 hinge dihedral all follow the "
                     "two clusters; action labels are pending (Q10)")
    save_fig(fig, "Y18_colour_by_biology", section="F. Embeddings and colour-by", task="Y18",
             shows="the same map coloured by ligand chemistry cluster, basic amine, ligand size, supported action label, pocket size, H12 geometry, 537 dihedral, H12/core B",
             takeaway="Basic-amine ligands, large ligands and large pockets sit in cl4; the H12 geometry variables change abruptly at the cluster border.",
             source="results/ext/colour_by_table.csv", script="scripts/ext_colour_by.py")
    k = [("cluster", "S08 structural cluster", "cluster"), ("peptide_class", "peptide", "cat", None, {"none": "#d4d3cd"}),
         ("mutation_status", "construct mutation", "cat"), ("ligand_basic_amine", "ligand basic amine", "cat"),
         ("r_free", "R-free", "num"), ("mean_b_chain", "chain mean B (A^2)", "num"), ("ph", "pH", "num"),
         ("space_group", "space group", "cat")]
    fig = grid(t, k, "Same story on the PCA of RMSD rows (PC1 90%): construct and ligand variables split along PC1, "
                     "quality variables do not", xcol="pca1", ycol="pca2")
    save_fig(fig, "Y18_colour_by_pca", section="F. Embeddings and colour-by", task="Y18",
             shows="PCA (RMSD rows) map coloured by cluster, peptide, mutation, basic amine, R-free, B, pH, space group",
             takeaway="The PCA map reproduces the MDS picture: PC1 is the H12 state, set by construct and ligand, not by quality.",
             source="results/ext/colour_by_table.csv, results/ext/embedding_rmsd.csv", script="scripts/ext_colour_by.py")
    write_numbers("Y18", dict(n_chains=len(t), n_supported_action=int((t.action_supported != "unlabelled (Q10 pending)").sum())))


if __name__ == "__main__":
    main()
