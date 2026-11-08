"""SV11 Part B (RUN_SV1028): steroid ring-system flag and cl7 vs cl4 share of steroid ligands.

Inputs (rule 1.10):
- Q08 out results/ligands.csv: SMILES as in SV10 = canonical_isomeric_smiles, else smiles.
- SV10 out results/sv/SV10/ligand_druglikeness.csv: the ligand set (306 ligands of interest), ligand group (cl7 / cl4
  when >= 80 % of the ligand's chains are in that S08 cluster) and `compared`; results/sv/SV10/entry_druglikeness.csv:
  entry group (same rule on the entry's kept chains), the entry's ligand (ligand_id) and `compared`.
- SV00 out results/sv/SV00/chain_table.csv: ligand_chemotype per ligand_id_primary, copied unchanged for reference.
Step 6: steroid ring system = RDKit substructure match of the run-file SMARTS (carbon atoms, any bond order).
  Check before use: EST must match; DES, GEN, OHT and RAL must not; otherwise Part B is BLOCKED (no flag, no test).
Step 7: Fisher exact test (two-sided) of the steroid share, cl7 vs cl4, at ligand level (SV10 compared ligands) and
  entry level (SV10 compared entries; the entry's flag = the flag of its SV10 ligand). BH over these two tests.
Step 8: the flag is its own column (steroid_ring_system); ligand_chemotype is not changed. The flag names the ring
  system only (rule 1.9: no action classes).
Outputs: results/sv/SV11/steroid_smarts_check.csv, ligand_steroid.csv, steroid_tests.csv.
Run: .venv/bin/python -B -I scripts/sv/SV11_steroid.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV11"
SMARTS = ("[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2")
CONTROLS = {"EST": True, "DES": False, "GEN": False, "OHT": False, "RAL": False}
GROUPS = ("cl7", "cl4")


def smiles_of(lig):
    col = np.where(lig.canonical_isomeric_smiles.notna(), "canonical_isomeric_smiles",
                   np.where(lig.smiles.notna(), "smiles", ""))
    return pd.Series(col, index=lig.index), lig.canonical_isomeric_smiles.fillna(lig.smiles)


def flag(patt, smi):
    if not isinstance(smi, str) or not smi:
        return np.nan, "no SMILES"
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return np.nan, "SMILES did not parse"
    return bool(m.HasSubstructMatch(patt)), ""


def fisher_row(level, df):
    a, b = df[df.group == "cl7"], df[df.group == "cl4"]
    k7, k4 = int(a.steroid_ring_system.sum()), int(b.steroid_ring_system.sum())
    odds, p = fisher_exact([[k7, len(a) - k7], [k4, len(b) - k4]])
    return dict(level=level, n_cl7=len(a), steroid_cl7=k7, share_cl7=k7 / len(a), n_cl4=len(b), steroid_cl4=k4,
                share_cl4=k4 / len(b), odds_ratio_cl7_vs_cl4=odds, p=p, test="Fisher exact, two-sided",
                steroid_ligands_cl7=";".join(sorted(set(a.ligand_id[a.steroid_ring_system]))),
                steroid_ligands_cl4=";".join(sorted(set(b.ligand_id[b.steroid_ring_system]))))


def main():
    RDLogger.DisableLog("rdApp.*")
    OUT.mkdir(parents=True, exist_ok=True)
    patt = Chem.MolFromSmarts(SMARTS)
    lig = pd.read_csv(ROOT / "results/ligands.csv", keep_default_na=False, na_values=[""], dtype={"ligand_id": str})
    lig["smiles_column"], lig["smiles_used"] = smiles_of(lig)
    lig = lig.set_index("ligand_id")
    chk = []
    for code, expect in CONTROLS.items():
        got, reason = flag(patt, lig.smiles_used.get(code)) if code in lig.index else (np.nan, "not in ligands.csv")
        chk.append(dict(ligand_id=code, expected_match=expect, observed_match=got, reason=reason,
                        passed=bool(got == expect) if reason == "" else False))
    chk = pd.DataFrame(chk)
    chk.to_csv(OUT / "steroid_smarts_check.csv", index=False)
    if patt is None or not chk.passed.all():
        print("SMARTS check failed: Part B BLOCKED")
        print(chk.to_string())
        return

    ld = pd.read_csv(ROOT / "results/sv/SV10/ligand_druglikeness.csv", keep_default_na=False, na_values=[""],
                     dtype={"ligand_id": str})
    ed = pd.read_csv(ROOT / "results/sv/SV10/entry_druglikeness.csv", keep_default_na=False, na_values=[""],
                     dtype={"ligand_id": str})
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""],
                     dtype={"ligand_id_primary": str})
    chemo = ct.dropna(subset=["ligand_id_primary"]).groupby("ligand_id_primary").ligand_chemotype.agg(
        lambda s: ";".join(sorted(set(s))))
    rows = []
    for r in ld.itertuples():
        smi = lig.smiles_used.get(r.ligand_id)
        f, reason = flag(patt, smi)
        rows.append(dict(ligand_id=r.ligand_id, ligand_name=r.ligand_name, group=r.group, compared=r.compared,
                         n_chains=r.n_chains, heavy_atoms=r.heavy_atoms,
                         smiles_column=lig.smiles_column.get(r.ligand_id, ""), smiles_used=smi,
                         steroid_ring_system=f, flag_reason=reason,
                         ligand_chemotype_sv00=chemo.get(r.ligand_id, "")))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "ligand_steroid.csv", index=False)

    lc = out[out.compared.astype(bool) & out.steroid_ring_system.notna()].copy()
    lc["steroid_ring_system"] = lc.steroid_ring_system.astype(bool)
    ec = ed[ed.compared.astype(bool)].merge(out[["ligand_id", "steroid_ring_system"]], on="ligand_id", how="left")
    ec = ec[ec.steroid_ring_system.notna()].copy()
    ec["steroid_ring_system"] = ec.steroid_ring_system.astype(bool)
    t = pd.DataFrame([fisher_row("ligand", lc), fisher_row("entry", ec)])
    t["q"] = multipletests(t.p, method="fdr_bh")[1]
    t["n_ligands_flagged_all"] = int(out.steroid_ring_system.fillna(False).astype(bool).sum())
    t["n_ligands_all"] = len(out)
    t.to_csv(OUT / "steroid_tests.csv", index=False)
    print(chk.to_string())
    print("steroid ligands:", sorted(out.ligand_id[out.steroid_ring_system.fillna(False).astype(bool)]))
    print(t.drop(columns=["steroid_ligands_cl7", "steroid_ligands_cl4"]).to_string())


if __name__ == "__main__":
    main()
