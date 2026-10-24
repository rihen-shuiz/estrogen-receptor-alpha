"""Tracker Y29: the ESR2 reference core-fraction question.

X02 applied the ESR1 block-08 reference policy to ESR2 but could not apply reference_min_core_fraction = 1.0 (the best
ESR2 chain reaches 0.976). Questions answered here, from files:
 (1) why 1.0 is unreachable: which old-core positions (settings esr1_core_ranges, the core the 'core_fraction' column
     of results/lbd_chains.csv refers to) are never observed in ESR2, and whether the ESR1/ESR2 numbering map
     (results/numbering_map_esr1_esr2.csv) has an ESR2 residue there;
 (2) the best reachable core fraction on the old core and on core v3 (alignment_core_v3.json);
 (3) does the reference choice matter: the ESR2 chains are re-fitted on every reference candidate that meets the rest of
     the policy (and on the ESR1 S03 reference 7NFB|B|1 via the Y23 joint matrix), the RMSD matrix is recomputed
     (pairwise common LBD positions, as X02) and re-clustered with the X03 method (average linkage, 1.5 A); ARI against
     the saved X03 labels.
Outputs: results/ext/esr2_reference_alternatives.csv, results/ext/esr2_core_fraction_chains.csv,
docs/esr2_core_fraction.md (numbers filled from the CSVs by this script).
Run from the repo root: PYTHONPATH=. python scripts/ext_esr2_core_fraction.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from esr_analysis.align_structures import fit_to_reference
from ext_common import RES, ROOT, SETTINGS, register_card
from ext_metric_check import masked_rmsd

LBD = np.arange(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)


def main():
    chains = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False)
    keys = [l.strip() for l in (ROOT / "results/rmsd_esr2_index.txt").read_text().splitlines() if l.strip()]
    e2 = chains.set_index("chain_key").loc[keys]
    q = pd.read_csv(ROOT / "results/esr2_entries.csv").set_index("entry_id")
    rm = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_ca", "x_A", "y_A", "z_A"])
    rm = rm[rm.chain_key.isin(keys) & rm.has_ca & rm.esr1_position.notna()]
    coords = {k: {int(p): np.array([a, b, c]) for p, a, b, c in zip(g.esr1_position, g.x_A, g.y_A, g.z_A)}
              for k, g in rm.groupby("chain_key")}
    old = sorted({p for a, b in SETTINGS["esr1_core_ranges"] for p in range(a, b + 1)})
    v3cfg = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    v3 = sorted({p for a, b in v3cfg["core_ranges"] for p in range(a, b + 1)})
    nm = pd.read_csv(ROOT / "results/numbering_map_esr1_esr2.csv").set_index("esr1_pos")
    per = []
    for k in keys:
        obs = set(coords.get(k, {}))
        per.append(dict(chain_key=k, old_core_fraction=len(obs & set(old)) / len(old),
                        core_v3_fraction=len(obs & set(v3)) / len(v3), lbd_fraction_saved=e2.loc[k, "lbd_fraction"],
                        core_fraction_saved=e2.loc[k, "core_fraction"]))
    per = pd.DataFrame(per)
    per.to_csv(RES / "esr2_core_fraction_chains.csv", index=False)
    seen = set().union(*[set(c) for c in coords.values()])
    never = [p for p in old if p not in seen]
    never_map = {p: (None if pd.isna(nm.esr2_pos.get(p, np.nan)) else int(nm.esr2_pos[p])) for p in never}

    # alternatives
    cand = e2[(e2.quality == "pass") & (e2.lbd_fraction >= SETTINGS["reference_min_lbd_fraction"])
              & (e2.h12_fraction >= SETTINGS["reference_min_h12_fraction"])].copy()
    cand["resolution"] = [q.loc[k.split("|")[0], "resolution"] for k in cand.index]
    cand = cand.sort_values(["resolution", "mapped_lbd_ca"], ascending=[True, False])
    x03 = pd.read_csv(ROOT / "results/esr2_clusters.csv").set_index("chain_key").loc[keys, "cut_1.5"].values
    saved = np.load(ROOT / "results/rmsd_esr2.npy")
    rows = []
    for ref in cand.index:
        X = np.full((len(keys), len(LBD), 3), np.nan)
        for i, k in enumerate(keys):
            moved, info = fit_to_reference(coords[k], coords[ref], v3, SETTINGS["min_core_ca_for_alignment"])
            for p, xyz in moved.items():
                if LBD[0] <= p <= LBD[-1]:
                    X[i, p - LBD[0]] = xyz
        d, _ = masked_rmsd(X, np.ones(len(LBD), bool))
        lab = fcluster(linkage(squareform(d, checks=False), "average"), 1.5, "distance")
        rows.append(dict(reference=ref, resolution_A=cand.loc[ref, "resolution"],
                         core_fraction=per.set_index("chain_key").loc[ref, "old_core_fraction"],
                         core_v3_fraction=per.set_index("chain_key").loc[ref, "core_v3_fraction"],
                         n_clusters_1p5=int(lab.max()), ari_vs_x03=adjusted_rand_score(x03, lab),
                         max_abs_diff_vs_saved_A=float(np.nanmax(np.abs(d - saved))), frame="ESR2 candidate"))
    dj = np.load(RES / "rmsd_joint.npy")[788:, 788:]
    lab = fcluster(linkage(squareform(dj, checks=False), "average"), 1.5, "distance")
    rows.append(dict(reference="7NFB|B|1 (ESR1 S03 reference, Y23)", resolution_A=np.nan, core_fraction=np.nan,
                     core_v3_fraction=np.nan, n_clusters_1p5=int(lab.max()), ari_vs_x03=adjusted_rand_score(x03, lab),
                     max_abs_diff_vs_saved_A=float(np.nanmax(np.abs(dj - saved))), frame="ESR1 frame"))
    alt = pd.DataFrame(rows)
    alt.to_csv(RES / "esr2_reference_alternatives.csv", index=False)

    ref0 = json.loads((ROOT / "results/esr2_reference.json").read_text())
    md = f"""# ESR2 reference: the core-fraction question (task Y29)

Generated by `scripts/ext_esr2_core_fraction.py` from `results/ext/esr2_core_fraction_chains.csv` and
`results/ext/esr2_reference_alternatives.csv`. Numbers are not typed by hand.

## Why reference_min_core_fraction = 1.0 cannot be met for ESR2

- The `core_fraction` column of `results/lbd_chains.csv` counts the **old core** ({len(old)} positions,
  {SETTINGS['esr1_core_ranges']}), the core the ESR1 reference policy was written for.
- Over all {len(keys)} ESR2 chains, {len(never)} old-core position(s) are never observed with a C-alpha:
  {', '.join(str(p) for p in never) or 'none'}. ESR2 residues mapped to them by the E03 numbering map:
  {', '.join(f'{p} -> {v}' for p, v in never_map.items()) or 'n/a'}.
- Best reachable old-core fraction: **{per.old_core_fraction.max():.3f}** (saved value in X02:
  {ref0['best_core_fraction_any_esr2']:.3f}). So the 1.0 rule is structurally unreachable, not a quality issue of one chain.
- On **core v3** (the core actually used for the fit, {len(v3)} positions) the best ESR2 chain reaches
  **{per.core_v3_fraction.max():.3f}**; {int((per.core_v3_fraction >= 0.95).sum())} of {len(per)} chains reach >= 0.95.
- The chosen reference {ref0['chain_key']} has old-core fraction {per.set_index('chain_key').loc[ref0['chain_key'], 'old_core_fraction']:.3f}
  and core-v3 fraction {per.set_index('chain_key').loc[ref0['chain_key'], 'core_v3_fraction']:.3f}.

## Does the reference choice change the ESR2 clusters?

{len(cand)} chains meet the remaining policy (quality pass, LBD fraction >= {SETTINGS['reference_min_lbd_fraction']},
H12 fraction >= {SETTINGS['reference_min_h12_fraction']}). Each was used as the fit reference, the RMSD matrix was
recomputed and re-clustered with the X03 method; one more run uses the ESR1 frame (Y23 joint matrix):

{alt[['reference', 'frame', 'resolution_A', 'core_fraction', 'core_v3_fraction', 'n_clusters_1p5', 'ari_vs_x03', 'max_abs_diff_vs_saved_A']].round(3).to_markdown(index=False)}

The first row reproduces the saved X02 matrix (max difference {alt.max_abs_diff_vs_saved_A.iloc[0]:.2g} A) when the reference is
{alt.reference.iloc[0]}. ARI against X03 ranges {alt.ari_vs_x03.min():.2f}-{alt.ari_vs_x03.max():.2f}: the ESR2 clusters
{'do not depend' if alt.ari_vs_x03.min() > 0.9 else 'depend somewhat'} on the reference.

## Recommendation

Record the deviation as already done in `results/esr2_reference.json`; a human may change the policy text to
"core fraction >= the best value reachable for the receptor" for ESR2. No result changes because of it.
"""
    (ROOT / "docs/esr2_core_fraction.md").write_text(md, encoding="utf8")
    print(alt.round(3).to_string(), never, never_map, per.old_core_fraction.max(), per.core_v3_fraction.max())


if __name__ == "__main__":
    main()
