"""Tracker Y24: ESR1 vs ESR2 potency and selectivity from ChEMBL.

ESR1: data/chembl_esr1_activities.csv (X13; CHEMBL206, IC50/Ki/Kd, relation '='). ESR2: the same query for CHEMBL242,
downloaded here over all pages from the ChEMBL REST API (cached in results/ext/chembl_esr2_activities.csv; the run
date is written to results/ext/chembl_esr2_meta.json). Per compound and receptor: median pChEMBL over its activity
rows (rows without pChEMBL dropped from the median only). Selectivity = pESR1 - pESR2 (positive = ESR1-selective).
Sensitivity: the same with Ki/Kd rows only (binding constants, less assay-format dependent).
PDB ligands: ChEMBL ids matched to PDB ligand codes by X14 (results/chembl_pdb_matches.csv), joined to the structural
cluster of the ligand's chains (Y08 structural_majority, >= 80% of its chains). Mixed assays: descriptive.
Outputs: results/ext/chembl_esr2_activities.csv, results/ext/selectivity_ratios.csv, results/ext/selectivity_tests.csv,
figures/ext/Y24_potency_esr1_vs_esr2.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_selectivity.py [--refresh]
"""
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import BLUE, GREY, INK, INK2, ORANGE, RES, ROOT, VIOLET, mpl_setup, save_fig, write_numbers

BASE = "https://www.ebi.ac.uk"
FIRST = ("/chembl/api/data/activity.json?target_chembl_id=CHEMBL242&standard_type__in=IC50,Ki,Kd"
         "&standard_relation__exact=%3D&limit=1000&offset=0")
COLS = ["activity_id", "assay_chembl_id", "molecule_chembl_id", "pchembl_value", "relation", "standard_relation",
        "standard_type", "standard_units", "standard_value", "target_chembl_id", "type", "units", "value",
        "canonical_smiles", "assay_type"]


def download():
    out = RES / "chembl_esr2_activities.csv"
    if out.exists() and "--refresh" not in sys.argv:
        return pd.read_csv(out)
    url, rows = BASE + FIRST, []
    while url:
        for attempt in range(5):
            try:
                r = requests.get(url, timeout=120)
                r.raise_for_status()
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
        page = r.json()
        rows += [{c: a.get(c) for c in COLS} for a in page["activities"]]
        nxt = page["page_meta"].get("next")
        url = BASE + nxt if nxt else None
    df = pd.DataFrame(rows)
    df = df[df.relation.astype(str).eq("=") | df.relation.isna()]
    df.to_csv(out, index=False)
    (RES / "chembl_esr2_meta.json").write_text(json.dumps({"downloaded": dt.datetime.now().isoformat(timespec="minutes"),
                                                           "rows": len(df), "query": BASE + FIRST}))
    return df


def per_molecule(a, types=None):
    a = a.assign(pchembl_value=pd.to_numeric(a.pchembl_value, errors="coerce")).dropna(subset=["pchembl_value"])
    if types:
        a = a[a.standard_type.isin(types)]
    return a.groupby("molecule_chembl_id").agg(p=("pchembl_value", "median"), n=("pchembl_value", "size"))


def main():
    e1 = pd.read_csv(ROOT / "data/chembl_esr1_activities.csv")
    e2 = download()
    rows = {}
    for name, types in (("all", None), ("KiKd", ["Ki", "Kd"])):
        a, b = per_molecule(e1, types), per_molecule(e2, types)
        j = a.join(b, lsuffix="_esr1", rsuffix="_esr2", how="inner")
        rows[name] = j
    j = rows["all"].rename(columns={"p_esr1": "p_esr1", "p_esr2": "p_esr2"})
    j["selectivity"] = j.p_esr1 - j.p_esr2
    kk = rows["KiKd"]
    j["selectivity_KiKd"] = (kk.p_esr1 - kk.p_esr2).reindex(j.index)
    m = pd.read_csv(ROOT / "results/chembl_pdb_matches.csv")
    lab = pd.read_csv(RES / "chem_labels_best.csv")[["ligand_id", "structural_majority"]]
    m = m.merge(lab, on="ligand_id", how="left")
    pdb = m.groupby("molecule_chembl_id").agg(ligand_id=("ligand_id", "first"), structural=("structural_majority", "first"))
    j = j.join(pdb, how="left").reset_index()
    j.to_csv(RES / "selectivity_ratios.csv", index=False)
    tests = []
    w = wilcoxon(j.selectivity)
    tests.append(dict(test="median selectivity != 0 (all compounds)", n=len(j), statistic=w.statistic, p_value=w.pvalue,
                      effect=j.selectivity.median()))
    rho = spearmanr(j.p_esr1, j.p_esr2)
    tests.append(dict(test="Spearman pESR1 vs pESR2", n=len(j), statistic=rho[0], p_value=rho[1], effect=rho[0]))
    a7 = j[j.structural == "cl7"].selectivity
    a4 = j[j.structural == "cl4"].selectivity
    if len(a7) >= 3 and len(a4) >= 3:
        u, p = mannwhitneyu(a7, a4)
        tests.append(dict(test="selectivity, PDB ligands cl7 vs cl4", n=len(a7) + len(a4), statistic=u, p_value=p,
                          effect=a7.median() - a4.median()))
    tests = pd.DataFrame(tests)
    tests.to_csv(RES / "selectivity_tests.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={"width_ratios": [1.2, 1, 0.8]})
    ax = axs[0]
    ax.scatter(j.p_esr1, j.p_esr2, s=6, c="#c9c8c2", lw=0, label=f"ChEMBL compounds with both (n={len(j)})")
    for key, col in (("cl7", BLUE), ("cl4", ORANGE)):
        s = j[j.structural == key]
        ax.scatter(s.p_esr1, s.p_esr2, s=34, c=col, lw=0.5, edgecolor="white", label=f"PDB ligands in {key} chains (n={len(s)})")
    lim = [3.5, 11.5]
    ax.plot(lim, lim, color=INK, lw=0.8)
    ax.plot(lim, [l - 1 for l in lim], color=INK, lw=0.6, ls="--")
    ax.plot(lim, [l + 1 for l in lim], color=INK, lw=0.6, ls="--")
    ax.set(xlim=lim, ylim=lim, xlabel="median pChEMBL at ESR1 (CHEMBL206)", ylabel="median pChEMBL at ESR2 (CHEMBL242)",
           title=f"Potency at the two receptors (Spearman {rho[0]:.2f}); dashed = 10-fold")
    ax.legend(fontsize=8, loc="upper left")
    ax = axs[1]
    ax.hist(j.selectivity, bins=np.arange(-4, 4.01, 0.25), color=VIOLET, edgecolor="white", label="all compounds")
    ax.hist(j.selectivity_KiKd.dropna(), bins=np.arange(-4, 4.01, 0.25), histtype="step", color=INK, lw=1.2,
            label=f"Ki/Kd only (n={j.selectivity_KiKd.notna().sum()})")
    ax.axvline(0, color=INK, lw=0.8)
    ax.set(xlabel="selectivity pESR1 - pESR2 (log units; + = ESR1-selective)", ylabel="compounds",
           title=f"Median selectivity {j.selectivity.median():+.2f} (Wilcoxon p={w.pvalue:.1g})")
    ax.legend(fontsize=8)
    ax = axs[2]
    rng = np.random.default_rng(42)
    for k, (vals, col, nm) in enumerate(((a7, BLUE, "cl7"), (a4, ORANGE, "cl4"))):
        ax.scatter(k + rng.uniform(-.15, .15, len(vals)), vals, s=20, c=col, lw=0)
        if len(vals):
            ax.hlines(vals.median(), k - .3, k + .3, color=INK, lw=1.5)
    ax.axhline(0, color=GREY, lw=0.8)
    ax.set_xticks([0, 1], [f"cl7 ligands\nn={len(a7)}", f"cl4 ligands\nn={len(a4)}"])
    pt = tests[tests.test.str.contains("cl7 vs cl4")]
    ttl = f"MWU p={pt.p_value.iloc[0]:.2g}" if len(pt) else "too few ligands"
    ax.set(ylabel="pESR1 - pESR2", title=f"PDB ligands by structural cluster ({ttl})")
    fig.suptitle(f"ESR1 vs ESR2 potency for {len(j)} ChEMBL compounds: potencies correlate (rho {rho[0]:.2f}); how many "
                 f"are within 10-fold? only {(j.selectivity.abs() <= 1).mean():.0%}; PDB ligands of cl4 chains lean ESR1-selective", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y24_potency_esr1_vs_esr2", section="I. ESR2 and ESR1 vs ESR2", task="Y24",
             shows="median pChEMBL at ESR1 vs ESR2 per compound, selectivity distribution (all assays and Ki/Kd only), selectivity of PDB ligands by structural cluster",
             takeaway=(f"{(j.selectivity.abs() <= 1).mean():.0%} of {len(j)} compounds are within 10-fold between receptors; "
                       f"{(j.selectivity > 1).sum()} are >10-fold ESR1-selective and {(j.selectivity < -1).sum()} >10-fold ESR2-selective."),
             source="results/ext/selectivity_ratios.csv, results/ext/selectivity_tests.csv, results/ext/chembl_esr2_activities.csv",
             script="scripts/ext_selectivity.py")
    write_numbers("Y24", dict(n_esr2_rows=len(e2), n_both=len(j), median_sel=float(j.selectivity.median()),
                              frac_within_10fold=float((j.selectivity.abs() <= 1).mean()),
                              n_esr1_sel=int((j.selectivity > 1).sum()), n_esr2_sel=int((j.selectivity < -1).sum()),
                              spearman=float(rho[0]), tests=tests.round(6).to_dict("records")))
    print(tests)
    print(len(e2), len(j))


if __name__ == "__main__":
    main()
