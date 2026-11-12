"""Stage 09 - ligands, action classes, pockets and affinities against the lid groups.

Inputs (data/ligands/): ligand_descriptors.csv (RDKit descriptors per ligand code and the group most of its chains
fall in, >= 80 %), chain_class.csv and action_class.csv (uncurated ChEMBL-mechanism classes), p2rank_pockets.csv
(P2Rank on each chain with its ligand), pocket_sites.csv (where every P2Rank pocket sits: ligand pocket, coactivator
groove, dimer interface), pocket_filling_per_chain.csv (buriedness-grid pocket and ligand volumes),
affinity_per_ligand.csv (ChEMBL binding affinity relative to E2 in the same assay).
Outputs (results/09_ligands_pockets/): ligand_size_by_group.csv, class_by_group.csv, pocket_sites_by_group.csv,
numbers.json.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("09_ligands_pockets")
t = c.chains()

# ligand chemistry by the group most of a ligand's chains fall in
ld = c.read("ligand_descriptors")
ld = ld[ld.structural_majority.isin(["cl7", "cl4"])]
size = ld.groupby("structural_majority").agg(ligands=("ligand_id", "size"), heavy_atoms_median=("heavy_atoms", "median"),
                                             mw_median=("mw", "median"),
                                             basic_amine_share=("basic_amine", lambda s: s.astype(str).eq("True").mean()))
size.to_csv(OUT / "ligand_size_by_group.csv")
p_size = mannwhitneyu(ld[ld.structural_majority == "cl7"].heavy_atoms, ld[ld.structural_majority == "cl4"].heavy_atoms).pvalue

# classes (entries with >= 1 chain of the class in the group)
cc = c.read("chain_class")
cc["action_class"] = cc.action_class.fillna("unlabelled")
cls = cc.groupby(["action_class", "cluster"]).entry_id.nunique().unstack(fill_value=0)
cls.to_csv(OUT / "class_by_group.csv")
ac = c.read("action_class")

# P2Rank pockets: is the top pocket the ligand pocket; groove and dimer-interface pockets by group
ps = c.read("pocket_sites")
top = ps[ps["rank"] == 1]
sites = ps.groupby(["group", "chain_key"]).label.apply(set).reset_index()
site_share = sites.groupby("group").label.apply(
    lambda s: pd.Series({"chains": len(s), "groove": float(np.mean(["AF-2 groove" in x or "groove" in x for x in s])),
                         "dimer": float(np.mean(["dimer interface" in x or "dimer" in x for x in s]))})).unstack()
site_share.to_csv(OUT / "pocket_sites_by_group.csv")

# pocket filling
pf = c.read("pocket_filling")
pf = pf[pf.ligand_incomplete.astype(str).eq("False") & pf.filling_pct.notna()]   # chains with a complete ligand
ft = c.read("pocket_filling_tests")
q_fill = float(ft.query("comparison == 'cluster' and level == 'entry'").q_bh.iloc[0])
rho_vol = spearmanr(pf.groupby("ligand_id_primary").ligand_volume_A3.median(),
                    pf.groupby("ligand_id_primary").pocket_volume_A3.median())[0]

af = c.read("affinity")
c.write_numbers(OUT, {
    "ligands_compared": size.ligands.to_dict(), "heavy_atoms_median": size.heavy_atoms_median.to_dict(),
    "basic_amine_share": size.basic_amine_share.to_dict(), "heavy_atoms_p": p_size,
    "class_rows": ac.action_class.value_counts().to_dict(), "class_rows_total": len(ac),
    "class_entries_by_group": cls.to_dict(),
    "top_pocket_is_ligand_pocket": int(top.label.astype(str).str.contains("ligand").sum()), "chains_with_pockets": int(top.chain_key.nunique()),
    "pocket_sites_share": site_share.to_dict(),
    "filling_chains": len(pf), "filling_median_pct": float(pf.filling_pct.median()), "filling_q_cl7_cl4_entry": q_fill,
    "ligand_vs_pocket_volume_rho_ligand_level": float(rho_vol),
    "affinity_ligands_relative_to_e2": int(af.affinity_basis.eq("relative").sum()),
    "affinity_ligands_absolute_only": int(af.affinity_basis.eq("absolute").sum()),
})
print(size, cls, site_share, sep="\n")
