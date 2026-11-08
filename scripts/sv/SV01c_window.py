"""SV01c (RUN_SV1028) step 9: SV01b re-run on the SV02 data-derived H12 window, verdict side by side with SV01b.

SV02 (results/sv/SV02/helix_boundaries.csv, row H12, global consensus) proposes the H12 window 538-544, which
differs from the primary window 538-547, so step 9 runs. This script re-runs SV01b on the proposed window by
importing the SV01b model code (scripts/sv/SV01b_models.py) and the SV01a code it builds on
(scripts/sv/SV01a_completeness.py) and setting their window parameters; neither file is edited and no SV01a/SV01b
output is written. The lattice contacts come from the saved SV01b table (scripts/sv/common/lattice_contacts.py output,
results/sv/SV01b/lattice_contacts_per_chain.csv); they are not recomputed.

How the window change propagates (fixed before any SV01c model was fitted):
- Outcome 1, H12 complete: all 7 C-alpha of 538-544 observed in results/sv/SV00/observed_mask.npy
  (SV01a_completeness.H12 = proposed window; SV01b: all 10 of 538-547). H12 status complete / partial /
  not modelled is re-derived on 538-544. Sensitivity outcome: fraction of 538-544 observed (fractional logit).
- Outcome 2, local mobility: mean C-alpha B z-score (against the same chain's core v3 C-alpha B) over every observed
  residue of 528-544 = loop 528-537 + the observed H12 residues of the proposed window 538-544
  (SV01b_models.MOB = (528, window end); SV01b: 528-547). Residues 545-547 leave the mobility mean; the loop 528-537,
  the core v3 reference (528-531 in both core and window), C-alpha B, altloc policy and ddof 0 are unchanged.
  TLS sensitivity unchanged (chains of TLS-refined entries removed).
- Lattice terms (block X) keep residues 528-550: step 6 defines them on 528-550 as a fixed segment, not as the H12
  window, and the step asks for the same lattice terms; the SV01b per-chain counts are reused unchanged
  (lattice_asu, lattice_sym = heavy atoms of 528-550 within 4.0 A of an ASU / symmetry copy of another chain).
- Everything else is SV01b unchanged: chains (788 kept, models within cl7 and within cl4), predictors and blocks
  Q / X (+ lattice terms) / B, median imputation + 'not reported' indicators, z-scoring within the cluster, chain level
  with entry-clustered errors and entry level on the first matrix chain of each entry, drop rule (constant,
  separation, collinearity, joint separation), McFadden R2 gains (logit) / OLS R2 gains (mobility), block tests of
  the block entered last, Wald 95 % CI.
- BH: one family over every p-value of SV01c (as SV01b: one family per step). The reproduction check below is not
  part of the family.
- [DEFAULT] verdict rule per outcome and level as SV01b; step verdict "artifact" / "biology" only if the condition
  holds for both outcomes (complete, mobility) at both levels, otherwise "mixed". Fraction observed and mobility
  without TLS are sensitivity verdicts. The SV01c verdict is reported next to the primary SV01b verdict
  (results/sv/SV01b/numbers.json, verdict_table.csv).
- Step 8 (medoid H12 clash test) only if the SV01c step verdict is "mixed": then for each chain partial or not
  modelled on 538-544, medoid residues 538-544 of 5TM1|A|1 and 7UJ7|C|1 are superposed on the chain's core and
  clashes counted with the SV01b_clash.py rules (clash 3.0 A, blocked = >= 1 backbone atom in clash).
- Check: the same code run with the primary window 538-547 must reproduce results/sv/SV01b/model_blocks.csv,
  odds_ratios.csv and mobility_per_chain.csv (reproduce_SV01b.json), so that only the window differs.

Outputs (results/sv/SV01c/): model_blocks.csv, odds_ratios.csv, verdict_table.csv, side_by_side.csv, verdict.md,
chain_outcomes.csv, status_transitions.csv, mobility_per_chain.csv, lattice_contacts_by_group.csv,
predictor_scaling.csv, model_data.csv, reproduce_SV01b.json, numbers.json, env.txt.
Run: .venv/bin/python -B scripts/sv/SV01c_window.py
"""
from __future__ import annotations

import json
import platform
import sys
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV01b_models as B  # noqa: E402  (imports SV01a_completeness as B.A)

A = B.A
OUT = ROOT / "results/sv/SV01c"
OUT.mkdir(parents=True, exist_ok=True)
SV01B = ROOT / "results/sv/SV01b"
PRIMARY = (538, 547)
LOOP_START = 528


def proposed_window():
    hb = pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv").set_index("helix_name")
    return int(hb.loc["H12", "global_start"]), int(hb.loc["H12", "global_end"])


def prepare(window):
    """Chains, outcomes and coded predictors for one H12 window (SV01b main(), parameterised)."""
    A.H12 = window
    B.MOB = (LOOP_START, window[1])
    t, mask, pos = A.load()
    mob, _, b_check, core_ranges = B.mobility(t, mask, pos)
    mob = mob.rename(columns={"n_obs_538_547": "n_obs_h12_window"})  # SV01b column name = H12 part of the window
    lat = pd.read_csv(SV01B / "lattice_contacts_per_chain.csv")
    assert (lat.chain_key.values == t.chain_key.values).all()
    t["mobility_z"] = mob.mobility_z.values
    t["mobility_defined"] = mob.mobility_defined.values
    t["mobility_defined_no_tls"] = (mob.mobility_defined & ~mob.tls_refined).values
    dcoded = A.code_predictors(t)
    dcoded["lattice_asu"] = lat.asu_atoms.astype(float).values
    dcoded["lattice_sym"] = lat.sym_atoms.astype(float).values
    return t, mob, lat, dcoded, b_check, core_ranges


def fit(window):
    t, mob, lat, dcoded, b_check, core_ranges = prepare(window)
    mb, orr, scale, mdata, n_tests = B.run_models(t, dcoded, B.X_SV01B, B.SPECS)
    return dict(t=t, mob=mob, lat=lat, mb=mb, orr=orr, scale=scale, mdata=mdata, n_tests=n_tests,
                b_check=b_check, core_ranges=core_ranges)


def reproduce_sv01b(r):
    """The parameterised code on 538-547 against the saved SV01b tables."""
    b_mb, b_or = pd.read_csv(SV01B / "model_blocks.csv"), pd.read_csv(SV01B / "odds_ratios.csv")
    b_mob = pd.read_csv(SV01B / "mobility_per_chain.csv")
    k = ["outcome", "cluster", "level", "block", "order"]
    j = r["mb"].merge(b_mb, on=k, suffixes=("", "_b"))
    ko = ["outcome", "cluster", "level", "predictor"]
    jo = r["orr"].merge(b_or, on=ko, suffixes=("", "_b"))
    f = jo[jo.status == "fitted"]
    m = r["mob"].merge(b_mob, on="chain_key", suffixes=("", "_b"))
    return dict(window=list(PRIMARY), model_block_rows=len(j), model_block_rows_sv01b=len(b_mb),
                max_abs_diff_gain=float(np.nanmax(np.abs(j.gain - j.gain_b))),
                max_abs_diff_block_p=float(np.nanmax(np.abs(j.p_value - j.p_value_b))),
                max_abs_diff_block_q=float(np.nanmax(np.abs(j.q_value - j.q_value_b))),
                odds_rows=len(jo), odds_rows_sv01b=len(b_or), status_equal=int((jo.status == jo.status_b).sum()),
                max_abs_diff_coef=float(np.nanmax(np.abs(f.coef - f.coef_b))),
                max_abs_diff_se=float(np.nanmax(np.abs(f.se - f.se_b))),
                max_abs_diff_q=float(np.nanmax(np.abs(f.q_value - f.q_value_b))),
                n_tests_bh=int(r["n_tests"]), n_tests_bh_sv01b=json.loads((SV01B / "numbers.json").read_text())[
                    "n_tests_bh"],
                mobility_rows=len(m), mobility_rows_sv01b=len(b_mob),
                max_abs_diff_mobility_z=float(np.nanmax(np.abs(m.mobility_z - m.mobility_z_b))),
                mobility_defined_equal=int((m.mobility_defined == m.mobility_defined_b).sum()))


def main():
    win = proposed_window()
    wtag, ptag = f"{win[0]}_{win[1]}", f"{PRIMARY[0]}_{PRIMARY[1]}"
    assert win != PRIMARY, "windows equal: SV01c is SKIPPED ('same window')"

    # check: primary window reproduces SV01b
    rp = fit(PRIMARY)
    rep = reproduce_sv01b(rp)
    (OUT / "reproduce_SV01b.json").write_text(json.dumps(rep, indent=1))
    print("reproduce SV01b:", rep)

    # SV01c: proposed window
    r = fit(win)
    mb, orr = r["mb"], r["orr"]
    mb.to_csv(OUT / "model_blocks.csv", index=False)
    orr.to_csv(OUT / "odds_ratios.csv", index=False)
    r["scale"].to_csv(OUT / "predictor_scaling.csv", index=False)
    r["mdata"].to_csv(OUT / "model_data.csv", index=False)
    mob = r["mob"].rename(columns={"mobility_h12_z": f"mobility_h12_{wtag}_z"})
    mob.to_csv(OUT / "mobility_per_chain.csv", index=False)
    t = r["t"]
    B.lattice_by_group(t, r["lat"]).to_csv(OUT / "lattice_contacts_by_group.csv", index=False)
    outcomes = [s["outcome"] for s in B.SPECS]
    vt, per_outcome, step, lv = B.verdicts(mb, outcomes)
    vt.to_csv(OUT / "verdict_table.csv", index=False)

    # per-chain outcomes on both windows
    tp = rp["t"]
    co = pd.DataFrame({"chain_key": t.chain_key, "entry_id": t.entry_id, "s08_cluster": t.s08_cluster,
                       "cluster_group": t.cluster_group, "tls_refined": rp["mob"].tls_refined.values,
                       f"h12_status_{ptag}": tp.h12_status, f"n_obs_{ptag}": tp.n_obs_window,
                       f"complete_{ptag}": tp.complete,
                       f"frac_obs_{ptag}": tp.frac_obs, f"mobility_z_528_{PRIMARY[1]}": tp.mobility_z,
                       f"h12_status_{wtag}": t.h12_status, f"n_obs_{wtag}": t.n_obs_window,
                       f"complete_{wtag}": t.complete,
                       f"frac_obs_{wtag}": t.frac_obs, f"mobility_z_528_{win[1]}": t.mobility_z})
    co.to_csv(OUT / "chain_outcomes.csv", index=False)
    tr = []
    for g, d in list(co.groupby("cluster_group")) + [("all", co)]:
        for lvl, dd in (("chain", d), ("entry", d.drop_duplicates("entry_id"))):
            ct = dd.groupby([f"h12_status_{ptag}", f"h12_status_{wtag}"]).size()
            for (a, b), n in ct.items():
                tr.append({"cluster_group": g, "level": lvl, f"h12_status_{ptag}": a, f"h12_status_{wtag}": b,
                           "n": int(n)})
    pd.DataFrame(tr).to_csv(OUT / "status_transitions.csv", index=False)

    # side by side with the primary SV01b verdict
    b_vt = pd.read_csv(SV01B / "verdict_table.csv")
    b_num = json.loads((SV01B / "numbers.json").read_text())
    k = ["outcome", "level", "cluster", "order"]
    sbs = b_vt[k + ["gain_QX", "gain_B", "level_verdict", "outcome_verdict", "step_verdict"]].merge(
        vt[k + ["gain_QX", "gain_B", "level_verdict", "outcome_verdict", "step_verdict"]], on=k,
        suffixes=(f"_{ptag}", f"_{wtag}"))
    sbs.to_csv(OUT / "side_by_side.csv", index=False)

    def g(df, outcome, level, c, order="first", suf=""):
        m = df[(df.outcome == outcome) & (df.level == level) & (df.cluster == c) & (df.order == order)].iloc[0]
        return m[f"gain_QX{suf}"], m[f"gain_B{suf}"]

    def cell(outcome, level):
        s = []
        for c in A.MAIN:
            qp, bp = g(sbs, outcome, level, c, suf=f"_{ptag}")
            qw, bw = g(sbs, outcome, level, c, suf=f"_{wtag}")
            s.append(f"cl{c} Q+X {qp:.3f} vs B {bp:.3f} -> {qw:.3f} vs {bw:.3f}")
        return "; ".join(s)

    bl = b_num["level_verdicts"]
    lines = [
        f"# SV01c verdict (H12 window {win[0]}-{win[1]}, SV02 proposal): **{step}**; primary SV01b ({PRIMARY[0]}-"
        f"{PRIMARY[1]}): **{b_num['step_verdict']}** ([DEFAULT] rule, lattice terms in X, chain and entry level)",
        f"- Complete H12: {b_num['outcome_verdicts']['complete']} -> **{per_outcome['complete']}** (chain "
        f"{bl['complete|chain']} -> {lv[('complete', 'chain')]}, entry {bl['complete|entry']} -> "
        f"{lv[('complete', 'entry')]}).",
        f"- Mobility (528-{PRIMARY[1]} -> 528-{win[1]}): {b_num['outcome_verdicts']['mobility']} -> "
        f"**{per_outcome['mobility']}** (chain {bl['mobility|chain']} -> {lv[('mobility', 'chain')]}, entry "
        f"{bl['mobility|entry']} -> {lv[('mobility', 'entry')]}).",
        f"- Sensitivity: fraction observed {b_num['outcome_verdicts']['fraction_observed']} -> "
        f"{per_outcome['fraction_observed']}; mobility without TLS {b_num['outcome_verdicts']['mobility_no_tls']} -> "
        f"{per_outcome['mobility_no_tls']}.",
        f"- Complete, gain entered first {ptag.replace('_', '-')} -> {wtag.replace('_', '-')}, chain: "
        + cell("complete", "chain") + ".",
        "- Complete, entry: " + cell("complete", "entry") + ".",
        "- Mobility, chain: " + cell("mobility", "chain") + "; entry: " + cell("mobility", "entry") + ".",
        "- Entered last gives the same Q+X vs B ordering as entered first (same rows in every model); both orders in "
        "side_by_side.csv.",
        f"- Lattice terms keep 528-550 (SV01b counts reused); only the outcomes moved. Step 8 "
        + ("runs (verdict mixed)." if step == "mixed" else "not run (verdict not mixed)."),
    ]
    (OUT / "verdict.md").write_text("\n".join(lines) + "\n")

    def n_inc(tt):
        return {g_: int((tt[tt.cluster_group == g_].complete == 0).sum()) for g_ in ("cl7", "cl4")}
    nums = dict(window_proposed=list(win), window_primary=list(PRIMARY), mobility_range=[LOOP_START, win[1]],
                lattice_range=[528, 550], n_chains=len(t), n_entries=int(t.entry_id.nunique()),
                n_tests_bh=int(r["n_tests"]), step_verdict=step, outcome_verdicts=per_outcome,
                level_verdicts={f"{a}|{b}": v for (a, b), v in lv.items()},
                sv01b_step_verdict=b_num["step_verdict"], sv01b_outcome_verdicts=b_num["outcome_verdicts"],
                sv01b_level_verdicts=b_num["level_verdicts"], run_step8=step == "mixed",
                n_incomplete_chains={ptag: n_inc(tp), wtag: n_inc(t)},
                n_incomplete_entries={ptag: {g_: int((tp[(tp.cluster_group == g_)].drop_duplicates("entry_id")
                                                      .complete == 0).sum()) for g_ in ("cl7", "cl4")},
                                      wtag: {g_: int((t[(t.cluster_group == g_)].drop_duplicates("entry_id")
                                                     .complete == 0).sum()) for g_ in ("cl7", "cl4")}},
                mobility_defined={g_: int(mob[mob.cluster_group == g_].mobility_defined.sum()) for g_ in ("cl7", "cl4")},
                mobility_defined_no_tls={g_: int((mob[mob.cluster_group == g_].mobility_defined
                                                 & ~mob[mob.cluster_group == g_].tls_refined).sum())
                                         for g_ in ("cl7", "cl4")},
                mobility_undefined_chains=mob.chain_key[~mob.mobility_defined & mob.s08_cluster.isin(A.MAIN)].tolist(),
                b_crosscheck=r["b_check"], core_ranges=r["core_ranges"])
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1))

    tz = timezone(timedelta(hours=5))
    import gemmi
    import matplotlib
    env = [f"# SV01c environment, {datetime.now(tz).strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}",
           f"numpy {np.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}",
           f"statsmodels {statsmodels.__version__}", f"matplotlib {matplotlib.__version__} (no figure drawn)",
           f"gemmi {gemmi.__version__} (step 8 clash test, SV01c_clash.py; lattice contacts reused from results/sv/SV01b)",
           "nothing installed, upgraded or downgraded by SV01c; no network used"]
    (OUT / "env.txt").write_text("\n".join(env) + "\n")
    print(json.dumps({k_: v for k_, v in nums.items() if k_ != "mobility_undefined_chains"}, indent=1))
    print(sbs.to_string())


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        main()
