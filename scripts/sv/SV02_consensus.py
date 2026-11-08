"""SV02 steps 3-6 (RUN_SV1028): per-position SSE fractions, consensus helices, per-chain helix table, H12 window.

Inputs: results/sv/SV02/sse_codes_per_chain.csv and data/sv/SV02/dssp_residues.csv.gz (SV02_dssp.py, pydssp),
results/sv/SV00/chain_table.csv (s08_cluster: 7 = cl7, 4 = cl4), results/helix_table.csv (K02 card, existing helix
table) and the project's names for its rows (scripts/ext_common.HELIX_NAMES: H1, H3, H4/5, H6, H7, H8, H9, H10/11,
H12; the table itself names only H12).

Step 3: for each position 305-550 and group (global = 788 kept chains, cl7, cl4), among the chains that observe the
position (all four backbone atoms, i.e. a pydssp state): fraction H (pydssp 'H' = any helix), 3_10 and pi = not
available (pydssp fallback), other = 1 - H.
Step 4 [DEFAULT], implemented as fixed before looking at the output:
  - helical position = H in >= 80 % of observing chains;
  - runs of consecutive helical positions; low-support positions (observed in < 50 % of the group's chains) are
    trimmed from both ends of a run, so they never form a start or end (they may lie inside a helix);
  - a consensus helix = >= 4 consecutive positions after trimming;
  - name = the helix-table row with the largest overlap (ties: earlier row), named with ext_common.HELIX_NAMES, or
    '<row> (unnamed)' for a row without a project name, 'unmatched' without overlap; two helices with the same name
    in one group get suffixes a, b, ... in sequence order.
Step 5: per chain x global consensus helix. Chain helix runs = maximal runs of consecutive UniProt positions with
pydssp 'H' inside one pydssp segment (full chain, so flanks outside 305-550 count). Runs overlapping the consensus
range: start = first position of the first run, end = last position of the last run; censored start/end = one of the
two flanking residues (start - 1, start - 2 / end + 1, end + 2) is not observed with all backbone atoms or lies across
a chain break. Two residues, not one: under DSSP's H-bond rule the residue right next to an H run is always observed
(it closes the last i -> i+4 H-bond), so a one-residue test is False by construction (it flagged 0 of 6277 chain
helices on the first run); the two-residue test equals SV03's 'real end = next 2 residues observed and not helical'.
Step 6: data-derived H12 window = global consensus helix named H12; compared with 538-547, 537-546 and 536-545.

Outputs (results/sv/SV02/): sse_per_position.csv, consensus_{global,cl7,cl4}.csv, consensus_runs_all.csv,
helix_boundaries.csv, sse_per_chain.csv, h12_window.md, numbers_consensus.json.
Run: .venv/bin/python -B -I scripts/sv/SV02_consensus.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from ext_common import HELIX_NAMES  # noqa: E402  (project names for the helix-table rows)

OUT = ROOT / "results/sv/SV02"
POS = np.arange(305, 551)
H_MIN, SUPPORT_MIN, MIN_LEN = 0.80, 0.50, 4
WINDOWS = {"538-547 (settings h12_range, primary)": (538, 547), "537-546 (helix table AC4, 7NFB struct_conf)": (537, 546),
           "536-545": (536, 545)}


def groups(ct):
    cl = ct.set_index("chain_key").s08_cluster
    return {"global": cl.index.tolist(), "cl7": cl[cl == 7].index.tolist(), "cl4": cl[cl == 4].index.tolist()}


def per_position(mat, ht, grp):
    rows = []
    for g, keys in grp.items():
        m = mat.loc[keys]
        n = len(keys)
        for p in POS:
            col = m[str(p)]
            nobs = int(col.notna().sum())
            nh = int((col == "H").sum())
            fh = nh / nobs if nobs else np.nan
            rows.append({"position": p, "group": g, "n_chains": n, "n_observing": nobs, "frac_observing": nobs / n,
                         "low_support": nobs / n < SUPPORT_MIN, "n_H": nh, "frac_alpha_H": fh,
                         "frac_3_10_G": "not available", "frac_pi_I": "not available",
                         "frac_other": (1 - fh) if nobs else np.nan,
                         "helical_ge_80": bool(nobs and fh >= H_MIN)})
    d = pd.DataFrame(rows)
    row_of = {}
    for r in ht.itertuples():
        for p in range(r.start, r.end + 1):
            row_of[p] = r.helix
    d.insert(1, "helix_table_row", d.position.map(row_of))
    d.insert(2, "helix_table_name", d.helix_table_row.map(lambda x: HELIX_NAMES.get(x, f"{x} (unnamed)")
                                                          if isinstance(x, str) else ""))
    return d


def runs(flags, positions):
    out, start = [], None
    for p, f in zip(positions, flags):
        if f and start is None:
            start = p
        if not f and start is not None:
            out.append((start, prev))
            start = None
        prev = p
    if start is not None:
        out.append((start, positions[-1]))
    return out


def name_helix(s, e, ht):
    ov = [(min(e, r.end) - max(s, r.start) + 1, -r.start, r.helix) for r in ht.itertuples()]
    ov = [o for o in ov if o[0] > 0]
    if not ov:
        return "unmatched", "", ""
    best = max(ov)
    row = best[2]
    return (HELIX_NAMES.get(row, f"{row} (unnamed)"), row,
            ";".join(f"{h}:{n}" for n, _, h in sorted(ov, key=lambda x: -x[1])))


def consensus(pp, g, ht):
    d = pp[pp.group == g].set_index("position")
    allruns, kept = [], []
    for s, e in runs(d.helical_ge_80.to_numpy(), POS):
        seg = d.loc[s:e]
        ok = seg.index[~seg.low_support.to_numpy()]
        rec = {"group": g, "run_start": s, "run_end": e, "run_length": e - s + 1}
        if len(ok) == 0:
            rec.update(status="dropped: all positions low support")
            allruns.append(rec)
            continue
        ts, te = int(ok.min()), int(ok.max())
        L = te - ts + 1
        rec.update(start=ts, end=te, length=L, trimmed_start=ts - s, trimmed_end=e - te,
                   status="kept" if L >= MIN_LEN else f"dropped: < {MIN_LEN} positions after trimming")
        allruns.append(rec)
        if L >= MIN_LEN:
            inner = d.loc[ts:te]
            name, row, ovs = name_helix(ts, te, ht)
            kept.append({"group": g, "helix_name": name, "helix_table_row": row, "table_overlaps": ovs,
                         "start": ts, "end": te, "length": L, "run_start_before_trim": s, "run_end_before_trim": e,
                         "n_low_support_inside": int(inner.low_support.sum()),
                         "min_frac_alpha_H": round(float(inner.frac_alpha_H.min()), 4),
                         "mean_frac_alpha_H": round(float(inner.frac_alpha_H.mean()), 4),
                         "min_frac_observing": round(float(inner.frac_observing.min()), 4)})
    k = pd.DataFrame(kept)
    if len(k):
        dup = k.helix_name.duplicated(keep=False)
        for nm in k.helix_name[dup].unique():
            idx = k.index[k.helix_name == nm]
            for j, i in enumerate(idx):
                k.loc[i, "helix_name"] = f"{nm}{'abcdefgh'[j]}"
    return k, pd.DataFrame(allruns)


def chain_runs(res):
    """H runs per chain in UniProt numbering inside one pydssp segment; observed positions with their segment."""
    out, obs = {}, {}
    for key, g in res.groupby("chain_key"):
        g = g[g.esr1_position.notna()].sort_values("esr1_position")
        pos = g.esr1_position.astype(int).to_numpy()
        seg = g.segment.to_numpy()
        h = (g.sse == "H").to_numpy()
        obs[key] = dict(zip(pos, seg))
        rr, start = [], None
        for i in range(len(pos)):
            cont = i > 0 and pos[i] == pos[i - 1] + 1 and seg[i] == seg[i - 1]
            if start is not None and (not h[i] or not cont):
                rr.append((start, pos[i - 1]))
                start = None
            if h[i] and start is None:
                start = pos[i]
        if start is not None:
            rr.append((start, pos[-1]))
        out[key] = rr
    return out, obs


def per_chain(cons, res, ct, mat):
    rr, obs = chain_runs(res)
    cl = ct.set_index("chain_key").s08_cluster
    rows = []
    for key in mat.index:
        o = obs.get(key, {})
        for h in cons.itertuples():
            ov = [r for r in rr.get(key, []) if r[1] >= h.start and r[0] <= h.end]
            n_obs_c = sum(1 for p in range(h.start, h.end + 1) if p in o)
            n_h_c = int((mat.loc[key, [str(p) for p in range(h.start, h.end + 1)]] == "H").sum())
            rec = {"chain_key": key, "entry_id": key.split("|")[0], "s08_cluster": cl[key], "helix_name": h.helix_name,
                   "consensus_start": h.start, "consensus_end": h.end, "n_observed_in_consensus": n_obs_c,
                   "n_helical_in_consensus": n_h_c, "n_runs": len(ov)}
            if ov:
                s, e = ov[0][0], ov[-1][1]

                def censored(p, nbs):
                    return any(not (nb in o and o[nb] == o[p]) for nb in nbs)
                rec.update(status="helix", start=s, end=e, n_helical_observed=int(sum(b - a + 1 for a, b in ov)),
                           start_censored=censored(s, (s - 1, s - 2)), end_censored=censored(e, (e + 1, e + 2)))
            else:
                rec.update(status="not observed" if n_obs_c == 0 else "not helical", start=np.nan, end=np.nan,
                           n_helical_observed=0, start_censored=np.nan, end_censored=np.nan)
            rows.append(rec)
    return pd.DataFrame(rows)


def window_md(cg, c7, c4, pp, ht):
    h12 = cg[cg.helix_name == "H12"]
    lines = ["# SV02 - data-derived H12 window (proposal only)", "",
             "Status: **proposed, not adopted.** The H12 window for all primary analyses stays 538-547 (run file section 2);",
             "a human signs off any change. SV01c tests the proposal as a sensitivity analysis.", ""]
    src = ("pydssp 0.9.1 (fallback for mkdssp, not installed); 'H' = any helix (3-, 4- and 5-turn helices together); "
           "[DEFAULT] consensus helix = >= 4 consecutive positions with H in >= 80 % of the observing chains, "
           "low-support positions (observed in < 50 % of chains) never a boundary.")
    lines += [f"Method: {src}", ""]
    if len(h12) != 1:
        lines += [f"No single global consensus helix named H12 (found {len(h12)}); no window proposed."]
        return "\n".join(lines) + "\n", None
    h = h12.iloc[0]
    prop = (int(h.start), int(h.end))
    lines += [f"**Data-derived H12 window (global consensus, 788 chains): {prop[0]}-{prop[1]}** "
              f"({prop[1] - prop[0] + 1} positions; run before trimming {int(h.run_start_before_trim)}-"
              f"{int(h.run_end_before_trim)}; low-support positions inside: {int(h.n_low_support_inside)}; "
              f"min H fraction {h.min_frac_alpha_H:.3f}, min fraction observing {h.min_frac_observing:.3f}).", ""]
    for nm, c in (("cl7", c7), ("cl4", c4)):
        x = c[c.helix_name == "H12"]
        lines.append(f"- {nm} consensus H12: " + ("; ".join(f"{int(r.start)}-{int(r.end)}" for r in x.itertuples())
                                                  if len(x) else "none (no H12 run meets the rule)"))
    lines += ["", "## Comparison with the historical windows", "",
              "| window | positions | shared with proposal | only in window | only in proposal | Jaccard | start shift | end shift |",
              "|---|---|---|---|---|---|---|---|"]
    P = set(range(prop[0], prop[1] + 1))
    comp = []
    for lab, (s, e) in WINDOWS.items():
        W = set(range(s, e + 1))
        j = len(P & W) / len(P | W)
        comp.append({"window": lab, "start": s, "end": e, "shared": len(P & W), "only_window": sorted(W - P),
                     "only_proposal": sorted(P - W), "jaccard": j, "start_shift": prop[0] - s, "end_shift": prop[1] - e,
                     "equal": P == W})
        lines.append(f"| {lab} | {s}-{e} | {len(P & W)} | {sorted(W - P) or '-'} | {sorted(P - W) or '-'} | {j:.2f} | "
                     f"{prop[0] - s:+d} | {prop[1] - e:+d} |")
    eq = [c["window"] for c in comp if c["equal"]]
    lines += ["", ("The proposal equals " + "; ".join(eq) + ".") if eq else
              "The proposal equals none of the three historical windows.",
              ("It differs from the primary window 538-547, so SV01c re-runs SV01b on it (run file SV01c step 9)."
               if not comp[0]["equal"] else "It equals the primary window 538-547, so SV01c is SKIPPED ('same window')."),
              "", "## Per-position support, 530-550", "",
              "| position | helix table | observing global | H global | H cl7 | H cl4 | obs cl7 | obs cl4 |", "|---|---|---|---|---|---|---|---|"]
    w = pp.pivot(index="position", columns="group", values=["frac_alpha_H", "frac_observing"])
    for p in range(530, 551):
        r = pp[(pp.group == "global") & (pp.position == p)].iloc[0]
        lab = r.helix_table_name or "loop"

        def f(v):
            return "NA" if pd.isna(v) else f"{v:.3f}"
        lines.append(f"| {p} | {lab} | {f(w.loc[p, ('frac_observing', 'global')])} | {f(w.loc[p, ('frac_alpha_H', 'global')])} | "
                     f"{f(w.loc[p, ('frac_alpha_H', 'cl7')])} | {f(w.loc[p, ('frac_alpha_H', 'cl4')])} | "
                     f"{f(w.loc[p, ('frac_observing', 'cl7')])} | {f(w.loc[p, ('frac_observing', 'cl4')])} |")
    lines += ["", "Values from results/sv/SV02/sse_per_position.csv and consensus_{global,cl7,cl4}.csv."]
    return "\n".join(lines) + "\n", {"proposal": prop, "comparisons": comp}


def main():
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv")
    mat = pd.read_csv(OUT / "sse_codes_per_chain.csv", index_col=0, keep_default_na=False, na_values=[""])
    res = pd.read_csv(ROOT / "data/sv/SV02/dssp_residues.csv.gz")
    ht = pd.read_csv(ROOT / "results/helix_table.csv")
    grp = groups(ct)
    pp = per_position(mat, ht, grp)
    pp.to_csv(OUT / "sse_per_position.csv", index=False)

    cons, allr = {}, []
    for g in grp:
        k, a = consensus(pp, g, ht)
        cons[g] = k
        allr.append(a)
        k.to_csv(OUT / f"consensus_{g}.csv", index=False)
    pd.concat(allr).to_csv(OUT / "consensus_runs_all.csv", index=False)

    # helix boundaries side by side (outer join on name), with the helix-table range
    hb = None
    for g in grp:
        x = cons[g][["helix_name", "helix_table_row", "start", "end", "length"]].rename(
            columns={"start": f"{g}_start", "end": f"{g}_end", "length": f"{g}_length",
                     "helix_table_row": f"{g}_table_row"})
        hb = x if hb is None else hb.merge(x, on="helix_name", how="outer")
    hb["helix_table_row"] = hb[[f"{g}_table_row" for g in grp]].bfill(axis=1).iloc[:, 0]
    hb = hb.drop(columns=[f"{g}_table_row" for g in grp])
    hts = ht.set_index("helix")
    hb["table_start"] = hb.helix_table_row.map(hts.start)
    hb["table_end"] = hb.helix_table_row.map(hts.end)
    hb["global_minus_table_start"] = hb.global_start - hb.table_start
    hb["global_minus_table_end"] = hb.global_end - hb.table_end
    hb = hb.sort_values(["global_start", "cl7_start", "cl4_start"]).reset_index(drop=True)
    cols = ["helix_name", "helix_table_row", "table_start", "table_end"] + \
           [f"{g}_{c}" for g in grp for c in ("start", "end", "length")] + \
           ["global_minus_table_start", "global_minus_table_end"]
    hb = hb[cols]
    hb.to_csv(OUT / "helix_boundaries.csv", index=False)

    pc = per_chain(cons["global"], res, ct, mat)
    pc.to_csv(OUT / "sse_per_chain.csv", index=False)

    md, win = window_md(cons["global"], cons["cl7"], cons["cl4"], pp, ht)
    (OUT / "h12_window.md").write_text(md)

    num = {"groups": {g: len(k) for g, k in grp.items()},
           "n_consensus": {g: int(len(cons[g])) for g in grp},
           "consensus": {g: [f"{r.helix_name} {r.start}-{r.end}" for r in cons[g].itertuples()] for g in grp},
           "runs_dropped": {g: int((a.status != "kept").sum()) for g, a in zip(grp, allr)},
           "low_support_positions": {g: [int(p) for p in pp[(pp.group == g) & pp.low_support].position] for g in grp},
           "h12_proposal": list(win["proposal"]) if win else None,
           "h12_comparisons": win["comparisons"] if win else None}
    st = pc.groupby("helix_name").agg(n_helix=("status", lambda s: int((s == "helix").sum())),
                                     n_not_obs=("status", lambda s: int((s == "not observed").sum())),
                                     n_not_helical=("status", lambda s: int((s == "not helical").sum())),
                                     start_cens=("start_censored", lambda s: int((s == True).sum())),  # noqa: E712
                                     end_cens=("end_censored", lambda s: int((s == True).sum())))  # noqa: E712
    num["per_chain_status"] = st.reset_index().to_dict(orient="records")
    (OUT / "numbers_consensus.json").write_text(json.dumps(num, indent=1, default=int))
    print(json.dumps({k: num[k] for k in ("n_consensus", "consensus", "runs_dropped", "h12_proposal")}, indent=1))
    print(hb.to_string())
    print(st.to_string())


if __name__ == "__main__":
    main()
