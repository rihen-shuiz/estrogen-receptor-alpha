"""Stage 01 - the structure set: entries, receptor chains, quality tiers, constructs and peptides.

Inputs: data/structures/attrition_funnel.csv (chain counts at each selection step), entries.csv (526 entries with
resolution / R-free / mean-B flags), chain_table.csv (788 analysed chains), september_cluster_stats.csv (one chain per
entry, old 84-residue core) and october_cluster_stats.csv (all 788 chains, core v3).
Output: results/01_structures/numbers.json and counts.csv.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("01_structures")

att = c.read("attrition").set_index("step")
ent = c.read("entries")
ct = c.read("chain_table")
sep = c.read("september_clusters").set_index("cut")
octo = c.read("october_clusters").set_index("cut")
esr1_ids = (c.DATA / "ids/pdb_ids_esr1.txt").read_text().split()
esr2_ids = (c.DATA / "ids/pdb_ids_esr2.txt").read_text().split()

sel = ent[ent.entry_id.isin(ct.entry_id.unique())]
pep_entries = ct.groupby("entry_id").peptide_present_entry.max()
n = {
    "esr1_entries_queried": len(esr1_ids),
    "esr2_entries_queried": len(esr2_ids),
    "esr1_entries_with_receptor_chain": int(att.loc["ESR1 entries without an ESR1 receptor chain", "n_after"]),
    "esr1_chain_models": int(att.loc["DNA chains in ESR1 entries", "n_after"]),
    "chains_analysed": len(ct),
    "entries_analysed": int(ct.entry_id.nunique()),
    "selected_entries_qc_pass_all_three": int(sel.qc_pass.astype(str).eq("pass").sum()),
    "chains_tier": ct.qc_label.value_counts().to_dict(),
    "entries_with_peptide": int(pep_entries.astype(bool).sum()),
    "chains_with_peptide_bound": int(ct.peptide_bound.astype(bool).sum()),
    "chains_by_construct": ct.mutation_status.value_counts().to_dict(),
    "september_matrix_chains": 351,
    "september_groups_at_1p5": int(sep.loc[1.5, "n_clusters"]),
    "september_groups_by_cut": sep.n_clusters.to_dict(),
    "october_groups_by_cut": octo.n_clusters.to_dict(),
}
c.write_numbers(OUT, n)
pd.DataFrame([{"quantity": k, "value": v} for k, v in n.items() if not isinstance(v, dict)]).to_csv(
    OUT / "counts.csv", index=False)
print(n)
