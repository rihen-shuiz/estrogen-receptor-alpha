"""Tracker Y13: supervised check - which variables predict the main structural cluster (cl7 vs cl4)?

Rows: the 774 chains in cl7 or cl4 (the 14 minor-cluster chains have no binary target; they are kept in every other
analysis). Target: 1 = cl7, 0 = cl4. Cross-validation: GroupKFold (5 folds) with entry_id as group, so chains of one
PDB entry are never split between train and test.
Feature sets (non-structural ones are the candidate explanations; the structural set is a positive control):
  quality    resolution, R-free, R-work, mean B, Wilson B, completeness, tier, B-group, observed LBD C-alpha
  crystal    space group (5 most common, else other), crystal system, ASU chains, pH, temperature, precipitant class,
             release year
  construct  peptide present, peptide class, construct mutation status (Y537S / L536S / D538G / WT / Cys / other)
  ligand     ligand heavy atoms, MW, cLogP, rotatable bonds, basic amine, phenol count (Y09)
  context    quality + crystal + construct + ligand
  structure  H12 shift, H12-H3 angle, K362-E542, L540-pocket, H12-H11 distance (Y11)
Models: decision tree (depth 3), random forest (300 trees), L2 logistic regression (standardised). Median imputation
(fitted inside each training fold). Metrics: balanced accuracy and ROC AUC (mean and SD over folds).
Importance: permutation importance of the random forest on each held-out fold (10 repeats, balanced accuracy),
averaged; mutual information of each single feature with the target (sklearn, seed 42).
Outputs: results/ext/supervised_scores.csv, results/ext/feature_importance.csv, figures/ext/Y13_feature_importance.png,
figures/ext/Y13_decision_tree.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_supervised.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, plot_tree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK2, ORANGE, RES, SEED, VIOLET, chain_table, mpl_setup, save_fig, write_numbers,
                        CATEGORICAL)


def build():
    t = chain_table()
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status", "peptide_present", "peptide_class"]]
    em = pd.read_csv(RES / "entry_metadata.csv")[["entry_id", "crystal_system", "chains_asu", "ph", "temp_K",
                                                   "additive_class", "release_year", "r_work", "wilson_b",
                                                   "completeness_pct"]]
    sf = pd.read_csv(RES / "structural_features.csv")[["chain_key", "h12_shift_A", "h12_h3_angle_deg", "k362_e542_A",
                                                        "l540_pocket_A", "h12_to_h11_A", "lbd_ca_observed"]]
    ld = pd.read_csv(RES / "ligand_descriptors.csv")[["ligand_id", "heavy_atoms", "mw", "clogp", "rotb", "basic_amine",
                                                       "n_phenol"]]
    t = t.merge(cm, on="chain_key").merge(em, on="entry_id", how="left").merge(sf, on="chain_key")
    t = t.merge(ld, on="ligand_id", how="left")
    t = t[t.cluster.isin((7, 4))].reset_index(drop=True)
    y = (t.cluster == 7).astype(int).values
    sg = t.space_group.where(t.space_group.isin(t.space_group.value_counts().index[:5]), "other")
    X = pd.DataFrame(index=t.index)
    sets = {}
    q = {"resolution_A": t.resolution_A, "r_free": t.r_free, "r_work": t.r_work, "mean_b": t.mean_b_A2,
         "wilson_b": t.wilson_b, "completeness_pct": t.completeness_pct, "tier_A": (t.tier == "A").astype(float),
         "b_group_high": (t.b_group == "high").astype(float), "b_group_low": (t.b_group == "low").astype(float),
         "lbd_ca_observed": t.lbd_ca_observed}
    c = {f"sg={s}": (sg == s).astype(float) for s in sg.unique()}
    c.update({f"system={s}": (t.crystal_system == s).astype(float) for s in t.crystal_system.unique()})
    c.update({"chains_asu": t.chains_asu, "ph": t.ph, "temp_K": t.temp_K, "release_year": t.release_year,
              "precip_PEG": (t.additive_class == "PEG").astype(float)})
    k = {"peptide_present": t.peptide_present.astype(float)}
    k.update({f"peptide={s}": (t.peptide_class == s).astype(float) for s in t.peptide_class.unique() if s != "none"})
    k.update({f"mut={s}": (t.mutation_status == s).astype(float) for s in t.mutation_status.unique()})
    l = {"lig_heavy_atoms": t.heavy_atoms, "lig_mw": t.mw, "lig_clogp": t.clogp, "lig_rotb": t.rotb,
         "lig_basic_amine": t.basic_amine.astype(float), "lig_n_phenol": t.n_phenol}
    s = {"h12_shift_A": t.h12_shift_A, "h12_h3_angle": t.h12_h3_angle_deg, "k362_e542_A": t.k362_e542_A,
         "l540_pocket_A": t.l540_pocket_A, "h12_to_h11_A": t.h12_to_h11_A}
    for name, dct in (("quality", q), ("crystal", c), ("construct", k), ("ligand", l), ("structure", s)):
        for kk, v in dct.items():
            X[kk] = pd.to_numeric(v, errors="coerce").values
        sets[name] = list(dct)
    sets["context"] = sets["quality"] + sets["crystal"] + sets["construct"] + sets["ligand"]
    return t, X, y, sets


def models():
    return {"decision_tree": lambda: make_pipeline(SimpleImputer(strategy="median"),
                                                   DecisionTreeClassifier(max_depth=3, random_state=SEED,
                                                                          class_weight="balanced")),
            "random_forest": lambda: make_pipeline(SimpleImputer(strategy="median"),
                                                   RandomForestClassifier(300, random_state=SEED, n_jobs=-1,
                                                                          class_weight="balanced", min_samples_leaf=2)),
            "logistic": lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                              LogisticRegression(max_iter=5000, class_weight="balanced"))}


def main():
    t, X, y, sets = build()
    g = t.entry_id.values
    cv = GroupKFold(n_splits=5)
    rows, perm = [], []
    for sname, cols in sets.items():
        Xs = X[cols].values
        for mname, mk in models().items():
            ba, auc = [], []
            for k, (tr, te) in enumerate(cv.split(Xs, y, g)):
                m = mk().fit(Xs[tr], y[tr])
                p = m.predict_proba(Xs[te])[:, 1]
                ba.append(balanced_accuracy_score(y[te], p > 0.5))
                auc.append(roc_auc_score(y[te], p))
                if sname == "context" and mname == "random_forest":
                    pi = permutation_importance(m, Xs[te], y[te], n_repeats=10, random_state=SEED,
                                                scoring="balanced_accuracy")
                    perm.append(pi.importances_mean)
            rows.append(dict(model=mname, feature_set=sname, n_features=len(cols), n_chains=len(y),
                             n_entries=len(set(g)), balanced_accuracy_mean=np.mean(ba), balanced_accuracy_sd=np.std(ba),
                             roc_auc_mean=np.mean(auc), roc_auc_sd=np.std(auc)))
    sc = pd.DataFrame(rows)
    sc.to_csv(RES / "supervised_scores.csv", index=False)
    cols = sets["context"] + sets["structure"]
    Xi = SimpleImputer(strategy="median").fit_transform(X[cols])
    mi = mutual_info_classif(Xi, y, random_state=SEED)
    fi = pd.DataFrame({"feature": cols, "mutual_info": mi})
    fi["perm_importance"] = np.nan
    fi.loc[fi.feature.isin(sets["context"]), "perm_importance"] = np.mean(perm, 0)
    fi["set"] = [next(s for s in ("quality", "crystal", "construct", "ligand", "structure") if c in sets[s]) for c in cols]
    fi = fi.sort_values("mutual_info", ascending=False)
    fi.to_csv(RES / "feature_importance.csv", index=False)

    plt = mpl_setup()
    setcol = {"quality": GREY, "crystal": CATEGORICAL[3], "construct": VIOLET, "ligand": CATEGORICAL[5],
              "structure": INK2}
    fig, axs = plt.subplots(1, 3, figsize=(18, 7), gridspec_kw={"width_ratios": [1.1, 1, 1]})
    ax = axs[0]
    piv = sc.pivot(index="feature_set", columns="model", values="balanced_accuracy_mean").loc[
        ["quality", "crystal", "ligand", "construct", "context", "structure"]]
    yy = np.arange(len(piv))
    for k, (m, col) in enumerate(zip(piv.columns, [CATEGORICAL[2], CATEGORICAL[4], CATEGORICAL[6]])):
        ax.barh(yy + (k - 1) * 0.26, piv[m], height=0.25, color=col, label=m)
        for y_, v in zip(yy, piv[m]):
            ax.text(v + 0.005, y_ + (k - 1) * 0.26, f"{v:.2f}", va="center", fontsize=7, color=INK2)
    ax.axvline(0.5, color="black", ls="--", lw=1)
    ax.set_yticks(yy, piv.index)
    ax.set_xlim(0.4, 1.08)
    ax.set(xlabel="balanced accuracy (5-fold GroupKFold by entry)", title="Which variable groups predict cl7 vs cl4?")
    ax.legend(fontsize=8, loc="lower right")
    ax = axs[1]
    top = fi[fi.perm_importance.notna()].sort_values("perm_importance", ascending=False).head(15)[::-1]
    ax.barh(top.feature, top.perm_importance, color=[setcol[s] for s in top.set])
    ax.set(xlabel="permutation importance (drop in balanced accuracy)",
           title="Random forest on all context variables")
    ax = axs[2]
    top = fi.head(18)[::-1]
    ax.barh(top.feature, top.mutual_info, color=[setcol[s] for s in top.set])
    ax.set(xlabel="mutual information with cl7/cl4 (nats)", title="Single-feature mutual information")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=v, label=k) for k, v in setcol.items()], fontsize=8, loc="lower right")
    s_ = sc.set_index(["feature_set", "model"]).balanced_accuracy_mean
    fig.suptitle(f"Construct (peptide + mutation) predicts the cluster with balanced accuracy "
                 f"{s_.loc[('construct', 'random_forest')]:.2f}, ligand chemistry {s_.loc[('ligand', 'random_forest')]:.2f}, "
                 f"quality only {s_.loc[('quality', 'random_forest')]:.2f}", x=0.01, ha="left", fontweight="bold",
                 fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y13_feature_importance", section="C. Why two clusters", task="Y13",
             shows="cross-validated balanced accuracy per feature group and model; permutation importance; mutual information",
             takeaway=(f"Quality variables predict the split less well ({s_.loc[('quality', 'random_forest')]:.2f}; quality is "
                       f"confounded with construct and crystal form, see Y14/Y16); "
                       f"construct ({s_.loc[('construct', 'random_forest')]:.2f}) and ligand chemistry "
                       f"({s_.loc[('ligand', 'random_forest')]:.2f}) predict it well; top single feature "
                       f"{fi.iloc[0].feature} (MI {fi.iloc[0].mutual_info:.2f})."),
             source="results/ext/supervised_scores.csv, results/ext/feature_importance.csv", script="scripts/ext_supervised.py")

    fig, ax = plt.subplots(figsize=(17, 7.5))
    cols_c = sets["context"]
    dt = make_pipeline(SimpleImputer(strategy="median"),
                       DecisionTreeClassifier(max_depth=3, random_state=SEED, class_weight="balanced")).fit(X[cols_c].values, y)
    plot_tree(dt[-1], feature_names=cols_c, class_names=["cl4", "cl7"], filled=False, impurity=False, proportion=False,
              rounded=True, fontsize=8, ax=ax)
    for txt in ax.texts:
        bb = txt.get_bbox_patch()
        if bb is not None:
            lab = txt.get_text()
            bb.set_facecolor("#dbe8f8" if "class = cl7" in lab else "#fbe1d6")
            bb.set_edgecolor(INK2)
    ba_dt = s_.loc[("context", "decision_tree")]
    ax.set_title(f"Depth-3 decision tree on the context variables (all 774 chains; CV balanced accuracy {ba_dt:.2f}); "
                 "blue boxes predict cl7, orange cl4", fontsize=11, loc="left")
    save_fig(fig, "Y13_decision_tree", section="C. Why two clusters", task="Y13",
             shows="depth-3 decision tree (class-balanced) on quality, crystal, construct and ligand variables",
             takeaway="The first splits use the peptide / mutation and the ligand's basic amine, not quality variables.",
             source="results/ext/supervised_scores.csv", script="scripts/ext_supervised.py")
    write_numbers("Y13", dict(scores=sc.round(4).to_dict("records"), top_mi=fi.head(10)[["feature", "mutual_info"]].round(4).to_dict("records")))
    print(sc.round(3).to_string())
    print(fi.head(12).round(3).to_string())


if __name__ == "__main__":
    main()
