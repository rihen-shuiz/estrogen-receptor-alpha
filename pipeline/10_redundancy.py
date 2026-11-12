"""Stage 10 - redundancy: how many distinct conformations the 788 chains hold.

From the stage-02 RMSD matrix and groups:
  * Hill numbers of the group partition at each cut: q0 = number of groups, q1 = exp(Shannon), q2 = 1 / sum p^2.
  * Greedy set cover: a chain is covered at radius r if a representative lies within r A; representatives are chains
    with a complete H12 (chains without H12 look close to everyone). Representatives needed for 90 % of all chains.
  * Copies of one crystal: RMSD between chains of the same entry and group.
  * Novelty over time: for chains deposited from 2010 on, the RMSD to the nearest chain deposited in an earlier year.
Outputs (results/10_redundancy/): hill_numbers.csv, set_cover.csv, novelty.csv, numbers.json.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("10_redundancy")
d, idx = c.rmsd()
t = c.chains()
gr = c.read("groups")


def hill(lab):
    p = pd.Series(lab).value_counts(normalize=True).to_numpy()
    return len(p), float(np.exp(-(p * np.log(p)).sum())), float(1 / (p ** 2).sum())


hn = pd.DataFrame([dict(cut=cut, q0=q[0], q1=q[1], q2=q[2]) for cut in ("1.0", "1.5", "2.0", "2.5")
                   for q in [hill(gr[f"cut_{cut}"])]])
hn.to_csv(OUT / "hill_numbers.csv", index=False)


def greedy_cover(r, cand):
    within = d[:, cand] <= r
    cand_idx = np.where(cand)[0]
    covered = np.zeros(len(d), bool)
    cov = []
    while covered.mean() < 1.0:
        gain = within[~covered].sum(0)
        k = int(np.argmax(gain))
        if gain[k] == 0:
            break
        covered |= within[:, k]
        cov.append(covered.mean())
    return cov


cand = (t.h12_label == "complete").to_numpy()
sc = []
for r in (0.5, 0.75, 1.0, 1.5, 2.0):
    cov = greedy_cover(r, cand)
    for target in (0.90, 0.95):
        hit = [i + 1 for i, x in enumerate(cov) if x >= target - 1e-12]
        sc.append(dict(radius_A=r, coverage_target=target, representatives=hit[0] if hit else np.nan,
                       max_coverage=cov[-1] if cov else 0.0))
sc = pd.DataFrame(sc)
sc.to_csv(OUT / "set_cover.csv", index=False)
# coverage curve at 1.0 A for the figure
pd.DataFrame({"representatives": np.arange(1, len(greedy_cover(1.0, cand)) + 1),
              "coverage": greedy_cover(1.0, cand)}).to_csv(OUT / "cover_curve_1A.csv", index=False)

grp = t.group.to_numpy()
same = (t.entry_id.to_numpy()[:, None] == t.entry_id.to_numpy()[None, :]) & (grp[:, None] == grp[None, :])
copies = d[np.triu(same, 1)]

ent = c.read("entries")[["entry_id", "deposit_date"]]
year = t.merge(ent, on="entry_id", how="left").deposit_date.str[:4].astype(int).to_numpy()
nov = []
for i in np.where(year >= 2010)[0]:
    earlier = year < year[i]
    nov.append(dict(chain_key=t.chain_key[i], year=year[i], nearest_earlier_A=float(d[i, earlier].min())))
nov = pd.DataFrame(nov)
nov.to_csv(OUT / "novelty.csv", index=False)
first = {g: int(year[grp == g].min()) for g in ("cl7", "cl4")}

s = sc.set_index(["radius_A", "coverage_target"]).representatives
c.write_numbers(OUT, {
    "hill_1p5": hn[hn.cut == "1.5"].iloc[0].to_dict(), "hill_by_cut": hn.set_index("cut").q1.to_dict(),
    "reps_90pct_at_1A": s.loc[(1.0, 0.90)], "reps_90pct_at_0p5A": s.loc[(0.5, 0.90)],
    "copy_rmsd_median": float(np.median(copies)), "copy_pairs": int(len(copies)),
    "novel_share_since_2010_gt_1p5A": float((nov.nearest_earlier_A > 1.5).mean()),
    "novelty_median_since_2010": float(nov.nearest_earlier_A.median()), "first_year": first,
})
print(hn, sc, sep="\n")
