"""SV10 (RUN_SV1028): drug-likeness of the crystallised ligands of interest, cl7 vs cl4, with a heavy-atom size check.

Inputs (read only; paths confirmed through tracker/esr1.tasks, rule 1.10):
- results/ligands.csv (Q08 out): SMILES = canonical_isomeric_smiles, else smiles.
- results/ext/ligand_descriptors.csv (Y09 out): the ligand set (306 ligands of interest), structural_majority (cl7 / cl4
  when >= 80 % of the ligand's chains are in that S08 cluster; mixed cl7/cl4, minor or none otherwise) and heavy_atoms.
- results/potency_by_cluster.csv (X15/X16 evidence, tracker/prov/X16.json): pchembl, one value per ligand.
- results/sv/SV00/chain_table.csv (SV00 out): kept chains with entry_id, s08_cluster and ligand_id_primary (the first
  ligand-of-interest row of results/chain_ligand_of_interest.csv, Q18) for the entry level.

Step 1, descriptors per ligand (RDKit, definitions fixed by the run file):
  Lipinski MW (Descriptors.MolWt), cLogP (Crippen.MolLogP), HBD (Lipinski.NumHDonors), HBA (Lipinski.NumHAcceptors),
  violations = #(MW > 500, cLogP > 5, HBD > 5, HBA > 10); Veber rotatable bonds
  (rdMolDescriptors.CalcNumRotatableBonds), TPSA (rdMolDescriptors.CalcTPSA), violations = #(rotb > 10, TPSA > 140);
  QED.qed and the eight QED.properties; Brenk and PAINS A/B/C match counts with FilterCatalog plus the names of the
  matched alerts; LE = 1.37 * pChEMBL / heavy atoms and LipE = pChEMBL - cLogP for ligands with pChEMBL.
  A SMILES that does not parse gets an empty row with the reason.
Step 2, comparison cl7 vs cl4 at two levels:
  ligand level (primary): one row per ligand with structural_majority cl7 or cl4;
  entry level: one row per entry; the entry's group uses the same >= 80 % rule on the s08_cluster of its kept chains;
  its ligand of interest is the ligand_id_primary held by most of its kept chains (tie: of the tied ligands, the one on
  the entry's first chain in chain_table order); entries without a ligand of interest are dropped.
  Continuous parameters: median [IQR] per group, two-sided Mann-Whitney U. Counts of violations or alerts (Lipinski
  violations, Veber violations, QED ALERTS, Brenk matches, PAINS matches = A + B + C): share of units with >= 1,
  two-sided Fisher exact test. BH over all parameters within each level.
Step 3, [DEFAULT] size check: Spearman rho of each parameter with heavy atoms (ligand and entry level; the md table
  shows the ligand level); OLS rank(parameter) ~ group (cl4 = 1) + rank(heavy atoms) on the units of that level with
  the parameter present (average ranks within that subset), group p BH-corrected within each level.
  Label "adds beyond size" if the size-adjusted q < 0.05 at both levels, otherwise "restates size".
  Declared before running: a parameter with no variation among the compared units of a level is not tested at that
  level (excluded from both BH families) and its label is "BLOCKED: no variation".
Outputs: results/sv/SV10/ligand_druglikeness.csv, entry_druglikeness.csv, druglikeness_by_group.csv,
  druglikeness_by_group.md, verdict.md, numbers.json, env.txt. No figure.
Run: .venv/bin/python -B -I scripts/sv/SV10_druglikeness.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from rdkit import Chem, rdBase
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from scipy.stats import fisher_exact, mannwhitneyu, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV10"
GROUPS = ("cl7", "cl4")
SHARE_MIN = 0.8  # existing ligand analysis rule (scripts/ext_chem_zoo.py structural_majority)
QED_KEYS = ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

POTENCY_NOTE = ("pChEMBL from results/potency_by_cluster.csv mixes IC50, Ki and Kd from different assays "
                "(median per ChEMBL molecule, then per ligand)")
# parameter, kind, definition (kind: cont = Mann-Whitney on the value; count = Fisher on share with >= 1)
PARAMS = [
    ("mw", "cont", "Lipinski: molecular weight, Descriptors.MolWt (g/mol)"),
    ("clogp", "cont", "Lipinski: cLogP, Crippen.MolLogP"),
    ("hbd", "cont", "Lipinski: H-bond donors, Lipinski.NumHDonors"),
    ("hba", "cont", "Lipinski: H-bond acceptors, Lipinski.NumHAcceptors"),
    ("lipinski_violations", "count", "Lipinski violations: #(MW > 500, cLogP > 5, HBD > 5, HBA > 10)"),
    ("rotb", "cont", "Veber: rotatable bonds, rdMolDescriptors.CalcNumRotatableBonds"),
    ("tpsa", "cont", "Veber: TPSA, rdMolDescriptors.CalcTPSA (A^2)"),
    ("veber_violations", "count", "Veber violations: #(rotatable bonds > 10, TPSA > 140)"),
    ("qed", "cont", "QED.qed (weighted desirability, 0-1)"),
    ("qed_mw", "cont", "QED property MW"),
    ("qed_alogp", "cont", "QED property ALOGP"),
    ("qed_hba", "cont", "QED property HBA (QED acceptor SMARTS)"),
    ("qed_hbd", "cont", "QED property HBD"),
    ("qed_psa", "cont", "QED property PSA"),
    ("qed_rotb", "cont", "QED property ROTB"),
    ("qed_arom", "cont", "QED property AROM (aromatic rings)"),
    ("qed_alerts", "count", "QED property ALERTS (QED structural-alert matches)"),
    ("brenk_n", "count", "Brenk alerts: FilterCatalog BRENK match count"),
    ("pains_n", "count", "PAINS alerts: FilterCatalog PAINS_A + PAINS_B + PAINS_C match count"),
    ("le", "cont", f"Ligand efficiency LE = 1.37 * pChEMBL / heavy atoms; {POTENCY_NOTE}"),
    ("lipe", "cont", f"Lipophilic efficiency LipE = pChEMBL - cLogP; {POTENCY_NOTE}"),
]


def catalog(which):
    params = FilterCatalogParams()
    params.AddCatalog(which)
    return FilterCatalog(params)


CATALOGS = {"brenk": catalog(FilterCatalogParams.FilterCatalogs.BRENK),
            "pains_a": catalog(FilterCatalogParams.FilterCatalogs.PAINS_A),
            "pains_b": catalog(FilterCatalogParams.FilterCatalogs.PAINS_B),
            "pains_c": catalog(FilterCatalogParams.FilterCatalogs.PAINS_C)}


def describe(smiles):
    """Step 1 descriptors of one SMILES; (dict, reason) with reason empty when it parsed."""
    if not isinstance(smiles, str) or not smiles:
        return {}, "no SMILES"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}, "RDKit MolFromSmiles returned None"
    d = dict(mw=Descriptors.MolWt(mol), clogp=Crippen.MolLogP(mol), hbd=Lipinski.NumHDonors(mol),
             hba=Lipinski.NumHAcceptors(mol), rotb=rdMolDescriptors.CalcNumRotatableBonds(mol),
             tpsa=rdMolDescriptors.CalcTPSA(mol), rdkit_heavy_atoms=mol.GetNumHeavyAtoms())
    d["lipinski_violations"] = int(d["mw"] > 500) + int(d["clogp"] > 5) + int(d["hbd"] > 5) + int(d["hba"] > 10)
    d["veber_violations"] = int(d["rotb"] > 10) + int(d["tpsa"] > 140)
    d["qed"] = QED.qed(mol)
    props = QED.properties(mol)
    for k in QED_KEYS:
        d[f"qed_{k.lower()}"] = getattr(props, k)
    for name, cat in CATALOGS.items():
        matches = cat.GetMatches(mol)
        d[f"{name}_n"] = len(matches)
        d[f"{name}_names"] = "; ".join(m.GetDescription() for m in matches)
    d["pains_n"] = d["pains_a_n"] + d["pains_b_n"] + d["pains_c_n"]
    d["pains_names"] = "; ".join(x for x in (d["pains_a_names"], d["pains_b_names"], d["pains_c_names"]) if x)
    return d, ""


def group_of(clusters):
    """Existing ligand-analysis rule: cl7 / cl4 when >= 80 % of the chains are in that S08 cluster."""
    c = np.asarray(clusters)
    if len(c) == 0:
        return "none"
    if np.mean(c == 7) >= SHARE_MIN:
        return "cl7"
    if np.mean(c == 4) >= SHARE_MIN:
        return "cl4"
    return "mixed cl7/cl4" if np.isin(c, (7, 4)).any() else "minor"


def ligand_table():
    lig = pd.read_csv(ROOT / "results/ligands.csv", keep_default_na=False, na_values=[""], dtype={"ligand_id": str})
    desc = pd.read_csv(ROOT / "results/ext/ligand_descriptors.csv", keep_default_na=False, na_values=[""],
                       dtype={"ligand_id": str})
    pot = pd.read_csv(ROOT / "results/potency_by_cluster.csv", keep_default_na=False, na_values=[""],
                      dtype={"ligand_id": str})
    assert desc.ligand_id.is_unique and pot.ligand_id.is_unique and lig.ligand_id.is_unique
    base = desc[["ligand_id", "n_chains", "structural_majority", "heavy_atoms"]].rename(
        columns={"structural_majority": "group"})
    base = base.merge(lig[["ligand_id", "ligand_name", "canonical_isomeric_smiles", "smiles"]], on="ligand_id",
                      how="left", indicator=True)
    base["smiles_column"] = np.where(base.canonical_isomeric_smiles.notna(), "canonical_isomeric_smiles",
                                     np.where(base.smiles.notna(), "smiles", ""))
    base["smiles_used"] = base.canonical_isomeric_smiles.fillna(base.smiles)
    missing = set(base.ligand_id[base["_merge"] == "left_only"])
    rows = []
    for r in base.itertuples():
        d, reason = describe(r.smiles_used)
        if r.ligand_id in missing:
            d, reason = {}, "ligand_id not in results/ligands.csv"
        rows.append(dict(ligand_id=r.ligand_id, parsed=reason == "", parse_reason=reason, **d))
    out = base.drop(columns=["canonical_isomeric_smiles", "smiles", "_merge"]).merge(pd.DataFrame(rows), on="ligand_id")
    out = out.merge(pot[["ligand_id", "pchembl"]], on="ligand_id", how="left")
    out["le"] = 1.37 * out.pchembl / out.heavy_atoms
    out["lipe"] = out.pchembl - out.clogp
    out["compared"] = out.group.isin(GROUPS) & out.parsed
    lead = ["ligand_id", "ligand_name", "group", "compared", "n_chains", "heavy_atoms", "smiles_column", "smiles_used",
            "parsed", "parse_reason", "mw", "clogp", "hbd", "hba", "lipinski_violations", "rotb", "tpsa",
            "veber_violations", "qed"] + [f"qed_{k.lower()}" for k in QED_KEYS] + [
            "brenk_n", "brenk_names", "pains_a_n", "pains_b_n", "pains_c_n", "pains_n", "pains_names",
            "pains_a_names", "pains_b_names", "pains_c_names", "pchembl", "le", "lipe", "rdkit_heavy_atoms"]
    return out[lead], pot


def entry_table(lt):
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""],
                     dtype={"ligand_id_primary": str, "entry_id": str})
    rows = []
    for entry, g in ct.groupby("entry_id", sort=False):
        ligs = g.ligand_id_primary.dropna()
        counts = ligs.value_counts()
        if len(counts):
            top = counts[counts == counts.max()].index
            ligand = next(x for x in ligs if x in top)  # tie: the ligand of the entry's first chain
        else:
            ligand = None
        rows.append(dict(entry_id=entry, n_chains=len(g), n_cl7_chains=int((g.s08_cluster == 7).sum()),
                         n_cl4_chains=int((g.s08_cluster == 4).sum()), group=group_of(g.s08_cluster.values),
                         n_ligands_of_interest=int(len(counts)), ligand_id=ligand,
                         ligand_tie=bool(len(counts) and (counts == counts.max()).sum() > 1),
                         ligands_all=";".join(counts.index)))
    et = pd.DataFrame(rows)
    keep = ["ligand_id", "group", "parsed", "heavy_atoms"] + [p for p, _, _ in PARAMS]
    et = et.merge(lt[keep].rename(columns={"group": "ligand_group"}), on="ligand_id", how="left")
    et["compared"] = et.group.isin(GROUPS) & et.ligand_id.notna() & et.parsed.fillna(False).astype(bool)
    # check: the ligand groups of Y09 recomputed from the SV00 chain table
    regroup = ct.dropna(subset=["ligand_id_primary"]).groupby("ligand_id_primary").s08_cluster.agg(
        lambda s: group_of(s.values))
    chk = lt.set_index("ligand_id").group
    same = int((regroup.reindex(chk.index).fillna("none") == chk).sum())
    return et, ct, dict(ligand_groups_recomputed_equal=same, ligand_groups_total=int(len(chk)))


def fmt(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.3g}" if abs(x) < 1000 else f"{x:.0f}"


def compare(df, level):
    """Rows of the group table for one level (df = compared units of that level)."""
    rows = []
    for name, kind, _ in PARAMS:
        sub = df.dropna(subset=[name])
        a, b = sub[sub.group == "cl7"], sub[sub.group == "cl4"]
        r = dict(parameter=name, level=level, n_cl7=len(a), n_cl4=len(b))
        values = sub[name].astype(float)
        varies = values.nunique() > 1 and len(a) > 0 and len(b) > 0
        if kind == "cont":
            for g, x in (("cl7", a), ("cl4", b)):
                q1, med, q3 = np.percentile(x[name].astype(float), [25, 50, 75]) if len(x) else (np.nan,) * 3
                r.update({f"median_{g}": med, f"q1_{g}": q1, f"q3_{g}": q3,
                          f"summary_{g}": f"{fmt(med)} [{fmt(q1)}-{fmt(q3)}]"})
            if varies:
                res = mannwhitneyu(a[name].astype(float), b[name].astype(float), alternative="two-sided")
                r.update(test="Mann-Whitney U", statistic=res.statistic, p=res.pvalue)
        else:
            for g, x in (("cl7", a), ("cl4", b)):
                k = int((x[name] >= 1).sum())
                r.update({f"n_ge1_{g}": k, f"share_{g}": k / len(x) if len(x) else np.nan,
                          f"summary_{g}": f"{k}/{len(x)} ({100 * k / len(x):.0f} %)" if len(x) else "NA"})
            if varies:
                tab = [[r["n_ge1_cl7"], len(a) - r["n_ge1_cl7"]], [r["n_ge1_cl4"], len(b) - r["n_ge1_cl4"]]]
                odds, p = fisher_exact(tab, alternative="two-sided")
                r.update(test="Fisher exact (>= 1)", statistic=odds, p=p)
        r["varies"] = bool(varies)
        if varies:
            ha = sub.heavy_atoms.astype(float)
            r["rho_heavy_atoms"] = spearmanr(values, ha).statistic
            X = pd.DataFrame({"const": 1.0, "group_cl4": (sub.group == "cl4").astype(float).values,
                              "rank_heavy_atoms": rankdata(ha)})
            fit = sm.OLS(rankdata(values), X).fit()
            r.update(adj_coef_group_cl4=fit.params["group_cl4"], adj_p=fit.pvalues["group_cl4"],
                     adj_coef_rank_ha=fit.params["rank_heavy_atoms"], adj_n=int(fit.nobs))
        rows.append(r)
    t = pd.DataFrame(rows)
    for col, qcol in (("p", "q"), ("adj_p", "adj_q")):
        t[qcol] = np.nan
        ok = t[col].notna()
        if ok.any():
            t.loc[ok, qcol] = multipletests(t.loc[ok, col], method="fdr_bh")[1]
    return t


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lt, pot = ligand_table()
    et, ct, chk = entry_table(lt)
    lt.to_csv(OUT / "ligand_druglikeness.csv", index=False)
    et.to_csv(OUT / "entry_druglikeness.csv", index=False)

    lig = compare(lt[lt.compared], "ligand")
    ent = compare(et[et.compared], "entry")
    defs = {p: d for p, _, d in PARAMS}
    kinds = {p: k for p, k, _ in PARAMS}
    wide = []
    for name, _, _ in PARAMS:
        L = lig.set_index("parameter").loc[name]
        E = ent.set_index("parameter").loc[name]
        if not (L.varies and E.varies):
            label = "BLOCKED: no variation"
        else:
            label = "adds beyond size" if (L.adj_q < 0.05 and E.adj_q < 0.05) else "restates size"
        row = dict(parameter=name, kind=kinds[name], definition=defs[name])
        for pre, X in (("lig", L), ("ent", E)):
            row.update({f"{pre}_n_cl7": X.n_cl7, f"{pre}_cl7": X.summary_cl7, f"{pre}_n_cl4": X.n_cl4,
                        f"{pre}_cl4": X.summary_cl4, f"{pre}_test": X.get("test"), f"{pre}_statistic": X.get("statistic"),
                        f"{pre}_p": X.p, f"{pre}_q": X.q, f"{pre}_rho_heavy_atoms": X.get("rho_heavy_atoms"),
                        f"{pre}_adj_coef_group_cl4": X.get("adj_coef_group_cl4"), f"{pre}_adj_p": X.adj_p,
                        f"{pre}_adj_q": X.adj_q})
            for g in GROUPS:
                for k in ("median", "q1", "q3", "share", "n_ge1"):
                    if f"{k}_{g}" in X.index:
                        row[f"{pre}_{k}_{g}"] = X[f"{k}_{g}"]
        row["label"] = label
        signs = {np.sign(L.get("adj_coef_group_cl4", np.nan)), np.sign(E.get("adj_coef_group_cl4", np.nan))}
        row["adj_direction"] = ("cl4 higher at equal heavy atoms" if signs == {1.0} else
                                "cl4 lower at equal heavy atoms" if signs == {-1.0} else "differs between levels")
        wide.append(row)
    w = pd.DataFrame(wide)
    # QED properties that are numerically the same as a Lipinski / Veber parameter for every ligand
    pairs = [("qed_mw", "mw"), ("qed_alogp", "clogp"), ("qed_hba", "hba"), ("qed_hbd", "hbd"), ("qed_psa", "tpsa"),
             ("qed_rotb", "rotb")]
    parsed = lt[lt.parsed]
    n_equal = {q: int(np.isclose(parsed[q].astype(float), parsed[b].astype(float), rtol=0, atol=1e-9).sum())
               for q, b in pairs}
    dup_of = {q: b for q, b in pairs if n_equal[q] == len(parsed)}
    w["identical_to"] = w.parameter.map(dup_of).fillna("")
    w.to_csv(OUT / "druglikeness_by_group.csv", index=False)

    def pq(p):
        return "NA" if pd.isna(p) else f"{p:.2g}"

    md = ["# SV10 drug-likeness by structural group (cl7 vs cl4)", "",
          "Unit: ligand level = one row per ligand of interest whose chains are >= 80 % in one group; entry level = one "
          "row per entry (its majority group by the same rule, the ligand held by most of its chains). Continuous: "
          "median [IQR], Mann-Whitney; counts: units with >= 1 (share), Fisher exact. q = BH over the parameters "
          "within each level. rho = Spearman with heavy atoms (ligand level). Size-adjusted q: group term of OLS "
          "rank(parameter) ~ group + rank(heavy atoms), BH within each level (ligand / entry). Label 'adds beyond "
          "size' when the size-adjusted q < 0.05 at both levels.", "",
          "These parameters describe oral drug-likeness, not binding, conformation or activity; the set includes the "
          "natural hormone, fragments and tool compounds. LE and LipE use pChEMBL, which mixes assay types "
          "(IC50, Ki, Kd from different assays).", "",
          "| parameter | definition | cl7 n | cl7 median [IQR] or share | cl4 n | cl4 median [IQR] or share | "
          "p ligand | q ligand | p entry | q entry | rho heavy atoms | size-adj q ligand / entry | label |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in w.itertuples():
        md.append(f"| {r.parameter} | {r.definition} | {r.lig_n_cl7} | {r.lig_cl7} | {r.lig_n_cl4} | {r.lig_cl4} | "
                  f"{pq(r.lig_p)} | {pq(r.lig_q)} | {pq(r.ent_p)} | {pq(r.ent_q)} | "
                  f"{'NA' if pd.isna(r.lig_rho_heavy_atoms) else f'{r.lig_rho_heavy_atoms:+.2f}'} | "
                  f"{pq(r.lig_adj_q)} / {pq(r.ent_adj_q)} | {r.label} |")
    md += ["", f"QED properties equal to a Lipinski / Veber row for all {len(parsed)} parsed ligands (same RDKit "
           f"function): " + ", ".join(f"{q} = {b}" for q, b in dup_of.items()) + "; kept as separate rows as the run "
           "file lists them. " + "; ".join(f"{q} equals {b} for {n_equal[q]} of {len(parsed)}"
                                            for q, b in pairs if q not in dup_of) + ".",
           "Direction of the size-adjusted group term (cl4 vs cl7 at equal heavy atoms): " +
           "; ".join(f"{r.parameter} {r.adj_direction.replace(' at equal heavy atoms', '')}"
                     for r in w.itertuples() if r.label == "adds beyond size") + ".",
           "", "Entry level: n and median [IQR] or share per group are in druglikeness_by_group.csv (ent_* columns).",
           "Source: scripts/sv/SV10_druglikeness.py; inputs results/ligands.csv, results/ext/ligand_descriptors.csv, "
           "results/potency_by_cluster.csv, results/sv/SV00/chain_table.csv."]
    (OUT / "druglikeness_by_group.md").write_text("\n".join(md) + "\n")

    # numbers for the verdict and SUMMARY
    cmp_l, cmp_e = lt[lt.compared], et[et.compared]

    def ha_summary(df):
        return {g: dict(n=int((df.group == g).sum()),
                        median=float(df[df.group == g].heavy_atoms.median()),
                        q1=float(df[df.group == g].heavy_atoms.quantile(0.25)),
                        q3=float(df[df.group == g].heavy_atoms.quantile(0.75))) for g in GROUPS}

    ha_p = {lvl: float(mannwhitneyu(d[d.group == "cl7"].heavy_atoms, d[d.group == "cl4"].heavy_atoms).pvalue)
            for lvl, d in (("ligand", cmp_l), ("entry", cmp_e))}
    beyond = w[w.label == "adds beyond size"]
    restates = w[w.label == "restates size"]
    blocked = w[w.label.str.startswith("BLOCKED")]
    nums = dict(
        n_ligands=int(len(lt)), n_parsed=int(lt.parsed.sum()), n_unparsed=int((~lt.parsed).sum()),
        unparsed=lt.loc[~lt.parsed, ["ligand_id", "parse_reason"]].to_dict("records"),
        smiles_column_counts=lt.smiles_column.value_counts().to_dict(),
        heavy_atoms_rdkit_equal_y09=int((lt.rdkit_heavy_atoms == lt.heavy_atoms).sum()),
        ligand_group_counts=lt.group.value_counts().to_dict(),
        n_ligands_compared={g: int((cmp_l.group == g).sum()) for g in GROUPS},
        n_ligands_with_pchembl=int(lt.pchembl.notna().sum()),
        n_ligands_with_pchembl_compared={g: int(cmp_l[cmp_l.group == g].pchembl.notna().sum()) for g in GROUPS},
        n_potency_rows=int(len(pot)),
        n_entries=int(len(et)), entry_group_counts=et.group.value_counts().to_dict(),
        n_entries_without_ligand=int(et.ligand_id.isna().sum()),
        n_entries_multiple_ligands=int((et.n_ligands_of_interest > 1).sum()),
        n_entries_ligand_tie=int(et.ligand_tie.sum()), entries_ligand_tie=et.loc[et.ligand_tie, "entry_id"].tolist(),
        n_entries_group_differs_from_ligand_group=int((cmp_e.group != cmp_e.ligand_group).sum()),
        n_entries_compared={g: int((cmp_e.group == g).sum()) for g in GROUPS},
        n_distinct_ligands_entry_level=int(cmp_e.ligand_id.nunique()),
        heavy_atoms_ligand=ha_summary(cmp_l), heavy_atoms_entry=ha_summary(cmp_e), heavy_atoms_mw_p=ha_p,
        n_parameters=len(PARAMS),
        n_q_lt_005={"ligand": int((w.lig_q < 0.05).sum()), "entry": int((w.ent_q < 0.05).sum())},
        n_adj_q_lt_005={"ligand": int((w.lig_adj_q < 0.05).sum()), "entry": int((w.ent_adj_q < 0.05).sum())},
        adds_beyond_size=beyond.parameter.tolist(), restates_size=restates.parameter.tolist(),
        blocked=blocked.parameter.tolist(),
        qed_property_n_equal=n_equal, qed_property_identical_to=dup_of,
        n_distinct_beyond_size=int((beyond.identical_to == "").sum()),
        n_ligands_brenk_ge1=int((lt.brenk_n >= 1).sum()), n_ligands_pains_ge1=int((lt.pains_n >= 1).sum()),
        brenk_names_compared=pd.Series([x for s in cmp_l.brenk_names.dropna() for x in s.split("; ") if x])
        .value_counts().head(10).to_dict(),
        pains_names_all=pd.Series([x for s in lt.pains_names.dropna() for x in s.split("; ") if x])
        .value_counts().to_dict(),
        **chk)
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1, default=float) + "\n")

    # verdict (<= 5 lines)
    twins = {b: q for q, b in dup_of.items()}

    def item(r):
        name = r.parameter + (f" = {twins[r.parameter]}" if r.parameter in twins and
                              twins[r.parameter] in set(beyond.parameter) else "")
        return (f"{name}: {r.adj_direction.replace(' at equal heavy atoms', '')} (cl7 {r.lig_cl7} vs cl4 {r.lig_cl4}; "
                f"size-adj q {pq(r.lig_adj_q)} / {pq(r.ent_adj_q)})")

    shown = beyond[beyond.identical_to == ""]
    v = [f"# SV10 verdict ([DEFAULT] size check): {len(beyond)} of {len(w)} drug-likeness parameters add beyond size "
         f"(size-adjusted q < 0.05 at ligand and entry level; {len(shown)} distinct quantities, since "
         f"{len(beyond) - len(shown)} are QED copies); {len(restates)} restate size"
         + (f"; {len(blocked)} not testable ({', '.join(blocked.parameter)}: no variation)" if len(blocked) else "") + ".",
         "- Adds beyond size, direction at equal heavy atoms (ligand-level values; q ligand / entry): "
         + ("; ".join(item(r) for r in shown.itertuples()) if len(shown) else "none") + ".",
         "- Restates size: " + (", ".join(restates.parameter) if len(restates) else "none") + ".",
         f"- Size context: heavy atoms median cl7 {nums['heavy_atoms_ligand']['cl7']['median']:.0f} vs cl4 "
         f"{nums['heavy_atoms_ligand']['cl4']['median']:.0f} (ligands, Mann-Whitney p {ha_p['ligand']:.2g}).",
         "- These parameters describe oral drug-likeness of the crystallised ligands (hormone, fragments, tool "
         "compounds included), not binding, conformation or activity; LE/LipE use pChEMBL that mixes assay types."]
    (OUT / "verdict.md").write_text("\n".join(v) + "\n")

    tz = datetime.now().astimezone()
    env = [f"# SV10 environment, {tz.isoformat(timespec='minutes')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}",
           f"rdkit {rdBase.rdkitVersion} (metadata {metadata.version('rdkit')})"]
    for name in ("numpy", "pandas", "scipy", "statsmodels"):
        env.append(f"{name} {metadata.version(name)}")
    env += ["nothing installed, upgraded or downgraded by SV10", "network: none"]
    (OUT / "env.txt").write_text("\n".join(env) + "\n")
    print(w[["parameter", "lig_cl7", "lig_cl4", "lig_q", "ent_q", "lig_rho_heavy_atoms", "lig_adj_q", "ent_adj_q",
             "label"]].to_string())
    print(json.dumps({k: nums[k] for k in ("n_ligands", "n_parsed", "n_ligands_compared", "n_entries_compared",
                                           "adds_beyond_size", "blocked", "ligand_groups_recomputed_equal")},
                     default=float))


if __name__ == "__main__":
    main()
