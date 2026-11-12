"""Stage 08 - crystallisation additives, solvent content and lattice contacts against lid (H12) order and lid mobility.

Question from the meeting: with many components added for crystallisation the receptor should look stable, and it
should be more dynamic in water. Crystal structures cannot show solution dynamics; they can show whether lid order
and lid mobility (local B) follow what is in or around the crystal.
Declared before running:
  * Predictors per chain: (1) bound additive copies per ESR1 chain (mmCIF non-polymer molecules whose code is on
    the project's additive list: ions, solvents, buffers, cryoprotectants); (2) other bound heterogens per chain
    (non-water, not the chain's primary ligand, not on the additive list); (3) recipe component classes in the
    crystallisation text (PEG, ammonium sulfate, MPD, DMSO, glycerol, other salt; 0-6); (4) solvent content of
    the crystal (mmCIF _exptl_crystal.density_percent_sol); (5) symmetry-mate contact atoms on residues 528-550;
    (6) residues of the whole chain touching a symmetry mate.
  * Outcomes: lid order = observed C-alpha of 538-547 (0-10); lid mobility = H12 C-alpha B z-score against the
    same chain's core (28 Oct lid mobility, defined chains only).
  * Between crystals: within cl7 and within cl4, entry means, Spearman; BH over all tests of this step. Effect =
    difference of entry medians of the outcome between the top and bottom quartile of the predictor; real only if
    q < 0.05 and |effect| exceeds both the within-group SD and the copy-to-copy 95th percentile of the outcome.
  * Within crystals (copies share additives, solvent and recipe, but not lattice contacts): for every entry with
    two or more chains in the same group that differ in contacts, does the copy with more contacts have the more
    ordered / less mobile lid? Sign test (binomial) over entries.
Outputs: results/08_crystal_environment/additives_per_entry.csv, additives_lattice_tests.csv, lattice_within_crystal.csv,
numbers_additives_lattice.json.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("08_crystal_environment")

RECIPE = ["has_peg", "has_ammonium_sulfate", "has_mpd", "has_dmso", "has_glycerol", "has_other_salt"]


def main():
    t = c.chains()
    # per-entry counts parsed from the deposited mmCIF files (_pdbx_nonpoly_scheme, _exptl_crystal) and the
    # recipe flags of the crystallisation text; additive = code on the project's additive list
    ent = c.read("additives")
    ent["recipe_classes"] = ent[RECIPE].astype(str).eq("True").sum(1)
    ent.to_csv(OUT / "additives_per_entry.csv", index=False)

    x = t.merge(ent[["entry_id", "additives_per_chain", "other_het_per_chain", "recipe_classes", "solvent_pct"]],
                on="entry_id", how="left")
    x["lid_order"] = x.n_h12_ca
    preds = {"additives_per_chain": "bound additive copies per chain",
             "other_het_per_chain": "other bound heterogens per chain",
             "recipe_classes": "recipe component classes (0-6)",
             "solvent_pct": "solvent content of the crystal (%)",
             "sym_atoms": "symmetry-mate contact atoms on 528-550",
             "n_residues_sym_contact": "residues touching a symmetry mate"}
    outs = {"lid_order": "observed C-alpha of 538-547", "mobility_h12_z": "lid mobility (B z-score)"}
    res = []
    for g in ("cl7", "cl4"):
        d = x[x.group == g]
        e = d.groupby("entry_id")[list(preds) + list(outs)].mean()
        for o in outs:
            ov = c.ordinary_variation(x, o, g)
            for pr in preds:
                ee = e[[pr, o]].dropna()
                if ee[pr].nunique() < 3:
                    res.append(dict(group=g, outcome=o, predictor=pr, n_entries=len(ee), rho=np.nan, p=np.nan,
                                    effect=np.nan, within_sd=ov["within_sd"], copy_p95=ov["copy_p95"]))
                    continue
                rho, pv = spearmanr(ee[pr], ee[o])
                q1, q3 = ee[pr].quantile([0.25, 0.75])
                if q3 > q1:   # top vs bottom quartile of the predictor
                    eff = float(ee[ee[pr] >= q3][o].median() - ee[ee[pr] <= q1][o].median())
                else:         # coarse predictor: above vs at-or-below its median, if both sides hold >= 5 entries
                    hi, lo = ee[ee[pr] > ee[pr].median()], ee[ee[pr] <= ee[pr].median()]
                    eff = float(hi[o].median() - lo[o].median()) if len(hi) >= 5 and len(lo) >= 5 else np.nan
                res.append(dict(group=g, outcome=o, predictor=pr, n_entries=len(ee), rho=rho, p=pv, effect=eff,
                                within_sd=ov["within_sd"], copy_p95=ov["copy_p95"]))
    rt = pd.DataFrame(res)
    rt["q"] = c.bh(rt.p)
    rt["threshold"] = rt[["within_sd", "copy_p95"]].astype(float).max(axis=1)
    rt["real"] = (rt.q < 0.05) & (rt.effect.abs() > rt.threshold)
    rt["predictor_label"] = rt.predictor.map(preds)
    rt["outcome_label"] = rt.outcome.map(outs)
    rt.to_csv(OUT / "additives_lattice_tests.csv", index=False)

    # ---------------- within one crystal: the copy with more contacts
    wc = []
    for contact in ("sym_atoms", "n_residues_sym_contact"):
        for o, better in (("lid_order", 1), ("mobility_h12_z", -1)):
            for g in ("cl7", "cl4"):
                agree = disagree = ties = 0
                for _, d in x[(x.group == g)].dropna(subset=[contact, o]).groupby("entry_id"):
                    if len(d) < 2 or d[contact].nunique() < 2:
                        continue
                    hi = d.loc[d[contact].idxmax()]
                    lo = d.loc[d[contact].idxmin()]
                    diff = (hi[o] - lo[o]) * better
                    if diff > 0:
                        agree += 1
                    elif diff < 0:
                        disagree += 1
                    else:
                        ties += 1
                n = agree + disagree
                wc.append(dict(contact=contact, outcome=o, group=g, entries=n + ties, more_contacts_better=agree,
                               more_contacts_worse=disagree, ties=ties,
                               p_sign=binomtest(agree, n).pvalue if n else np.nan))
    wc = pd.DataFrame(wc)
    wc["q"] = c.bh(wc.p_sign)
    wc.to_csv(OUT / "lattice_within_crystal.csv", index=False)

    c.write_numbers(OUT, {
        "n_entries": len(ent), "entries_with_additive": int((ent.n_additive_copies > 0).sum()),
        "median_additives_per_chain": float(ent.additives_per_chain.median()),
        "max_additives_per_chain": float(ent.additives_per_chain.max()),
        "solvent_pct_median": float(ent.solvent_pct.median()), "solvent_pct_n": int(ent.solvent_pct.notna().sum()),
        "n_tests": int(rt.p.notna().sum()), "n_q_lt_005": int((rt.q < 0.05).sum()), "n_real": int(rt.real.sum()),
        "significant": rt[rt.q < 0.05][["group", "outcome", "predictor", "rho", "effect", "threshold", "q", "real"]]
        .to_dict("records"),
        "within_crystal": wc.to_dict("records"),
    })
    pd.set_option("display.width", 250)
    print(ent.describe().round(2).to_string())
    print(rt[["group", "outcome", "predictor", "n_entries", "rho", "effect", "threshold", "q", "real"]].round(3)
          .to_string())
    print(wc.round(4).to_string())


if __name__ == "__main__":
    main()
