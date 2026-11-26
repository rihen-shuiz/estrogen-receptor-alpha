"""Stage 11 - ESR2: does the second oestrogen receptor show the same two lid positions?

Inputs (data/esr2/): rmsd_esr2.npy + rmsd_esr2_index.txt (86 ESR2 chains from 39 entries, superposed and measured
the same way as ESR1, ESR1 numbering), esr2_chain_table.csv (peptide per chain), rmsd_joint.npy + joint_index.csv
(788 ESR1 + 86 ESR2 chains in one frame on the same reference), esr1_esr2_displacement_profile.csv (per ESR1 position:
median C-alpha distance over 135 matched ESR1/ESR2 chain pairs with the same ligand after a core fit, and the residue
pair) with esr1_esr2_matched_pairs.csv.
Method: average linkage at 1.5 A, as for ESR1; the main and second ESR2 groups; median joint RMSD between the ESR1
groups and the ESR2 groups; classical MDS of the joint matrix (the ESR1+ESR2 landscape); the pocket selectivity
positions of Pike et al. 1999 (ESR1 L384 / ESR2 M336, ESR1 M421 / ESR2 I373) on the displacement profile.
Outputs (results/11_esr2/): esr2_groups.csv, joint_block_medians.csv, joint_map.csv, numbers.json.
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
# classical MDS of the joint matrix: the shared ESR1 + ESR2 landscape
n_j = len(dj)
jc = np.eye(n_j) - 1.0 / n_j
bmat = -0.5 * jc @ (dj ** 2) @ jc
w, v = np.linalg.eigh(bmat)
top = np.argsort(w)[::-1][:2]
xy = v[:, top] * np.sqrt(w[top])
xy *= np.sign(xy[np.abs(xy).argmax(axis=0), [0, 1]])  # fix the sign of each axis
jm = ji.assign(mds1=xy[:, 0], mds2=xy[:, 1])
jm.to_csv(OUT / "joint_map.csv", index=False)
prof = c.read("esr2_displacement")
full = prof[prof.n_pairs == prof.n_pairs.max()]
sel = prof.set_index("esr1_position").loc[[384, 421]]

c.write_numbers(OUT, {
    "esr2_chains": len(e2), "esr2_entries": int(e2.entry_id.nunique()), "group_sizes": sizes.to_dict(),
    "peptide_share_main": float(pep[main]), "peptide_share_second": float(pep[second]),
    "block_medians": {f"{r.group_a} - {r.group_b}": r.median_rmsd_A for r in bm.itertuples()},
    "joint_mds_variance_share": float(w[top].sum() / w[w > 0].sum()),
    "displacement_pairs": int(prof.n_pairs.max()), "displacement_positions": int(len(prof)),
    "displacement_median_all_positions_A": float(full.median_displacement_A.median()),
    "displacement_positions_over_1p5_A": int((full.median_displacement_A > 1.5).sum()),
    "differing_residues_shared_positions": int((~full.same_residue.astype(bool)).sum()),
    "selectivity_positions_A": {f"{r.esr1_aa}{i}/{r.esr2_aa}{int(r.esr2_position)}": float(r.median_displacement_A)
                                for i, r in sel.iterrows()},
})
print(sizes.to_dict(), pep.to_dict(), bm, sep="\n")
