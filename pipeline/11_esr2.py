"""Stage 11 - ESR2: does the second oestrogen receptor show the same two lid positions?

Inputs (data/esr2/): rmsd_esr2.npy + rmsd_esr2_index.txt (86 ESR2 chains from 39 entries, superposed and measured
the same way as ESR1, ESR1 numbering), esr2_chain_table.csv (peptide per chain), rmsd_joint.npy + joint_index.csv
(788 ESR1 + 86 ESR2 chains in one frame on the same reference).
Method: average linkage at 1.5 A, as for ESR1; the main and second ESR2 groups; median joint RMSD between the ESR1
groups and the ESR2 groups.
Outputs (results/11_esr2/): esr2_groups.csv, joint_block_medians.csv, numbers.json.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("11_esr2")
d2 = np.load(c.p("rmsd_esr2"))
keys = c.p("rmsd_esr2_index").read_text().split()
lab = fcluster(linkage(squareform(d2, checks=False), "average"), t=1.5, criterion="distance")
e2 = pd.DataFrame({"chain_key": keys, "entry_id": [k.split("|")[0] for k in keys], "group_1p5": lab})
tab = c.read("esr2_chains")[["chain_key", "group_1p5_saved", "peptide_present"]]
e2 = e2.merge(tab, on="chain_key", how="left", validate="one_to_one")
assert (e2.group_1p5 == e2.group_1p5_saved).all(), "ESR2 groups differ from the saved October groups"
e2.to_csv(OUT / "esr2_groups.csv", index=False)
sizes = e2.group_1p5.value_counts()
main, second = sizes.index[0], sizes.index[1]
pep = e2.groupby("group_1p5").peptide_present.apply(lambda s: s.astype(str).eq("True").mean())

dj = np.load(c.p("rmsd_joint"))
ji = c.read("joint_index")
assert len(ji) == len(dj)
lab_j = ji.cluster.to_numpy()
blocks = {
    ("ESR1 cl7", "ESR2 main"): ("ESR1 cl7", f"ESR2 group {main}"), ("ESR1 cl4", "ESR2 main"): ("ESR1 cl4", f"ESR2 group {main}"),
    ("ESR1 cl7", "ESR2 second"): ("ESR1 cl7", f"ESR2 group {second}"),
    ("ESR1 cl4", "ESR2 second"): ("ESR1 cl4", f"ESR2 group {second}"), ("ESR1 cl7", "ESR1 cl4"): ("ESR1 cl7", "ESR1 cl4"),
}
rows = []
for (a, b), (la, lb) in blocks.items():
    rows.append(dict(group_a=a, group_b=b, median_rmsd_A=float(np.median(dj[np.ix_(lab_j == la, lab_j == lb)]))))
bm = pd.DataFrame(rows)
bm.to_csv(OUT / "joint_block_medians.csv", index=False)
c.write_numbers(OUT, {
    "esr2_chains": len(e2), "esr2_entries": int(e2.entry_id.nunique()), "group_sizes": sizes.to_dict(),
    "peptide_share_main": float(pep[main]), "peptide_share_second": float(pep[second]),
    "block_medians": {f"{r.group_a} - {r.group_b}": r.median_rmsd_A for r in bm.itertuples()},
})
print(sizes.to_dict(), pep.to_dict(), bm, sep="\n")
