"""SV07 (RUN_SV1028): spike audit of the Y10 cl7-vs-cl4 displacement profile (slide 5).

Step text: docs/sv/briefs/SV07.md. Superposition core stays v3 (run file section 2). Diagnostic only: nothing in
results/ext or figures/ext is rewritten, the main analysis keeps all chains.

Inputs (rule 1.10):
- Y10 card out results/ext/cluster_displacement_profile.csv (spike rule is applied to its d_mean_A / d_medoid_A).
- Core v3 C-alpha: SV00 gate fit data/sv/gate/aligned_coordinates.npz (scripts/sv/common/core_frame.load_ca_core),
  checked against the S03 file results/alignment/aligned_coordinates.npz that Y10 read.
- SV00 results/sv/SV00/observed_mask.npy (+ _rows/_cols), chain_table.csv (s08_cluster, entry, peptide_bound,
  ligand_chemotype, dimer_partner_source, uniprot_author_offset, offset_constant, offset_source), chain_offsets.csv.
- Y10 medoid inputs: results/rmsd_entry.npy + rmsd_entry_index.csv (S04/S05), results/chain_status.csv h12_label.
- SV02 results/sv/SV02/sse_per_position.csv, sse_codes_per_chain.csv, helix_boundaries.csv.
- SV01b B z definition (scripts/sv/SV01b_models.mobility): C-alpha B from data/sv/backbone/backbone.npz (SV00),
  z = (B - mean core v3 B) / SD core v3 B of the same chain (valid B finite and >= 0, SD ddof 0, core needs >= 2 valid
  B and SD > 0). SV01b saved z only for 528-547 (results/sv/SV01b/mobility_per_residue.csv); the same formula is
  applied here to the spike positions and checked against that file on 528-547.
- SV01b module scripts/sv/common/lattice_contacts.py via scripts/sv/SV07_lattice.py
  (results/sv/SV07/lattice_contacts_spikes_per_residue.csv).
- Residue identity per chain and position: data/sv/lbd_extract/residue_mapping.csv (SV00 regenerated Q04 mapping;
  residue_status match / deposited_difference / unexplained_difference against UniProt P03372).

Definitions fixed before any diagnostic was computed:
1. [DEFAULT] spike = position outside 528-550 with Y10 d_mean_A > 1.0 A or d_medoid_A > 1.0 A (NaN never exceeds).
   region = maximal run of consecutive spike positions (labels only; every position gets its own verdict).
2. Per spike (cl7 = S08 cluster 7, cl4 = cluster 4, as Y10):
   a. observing chains / entries per cluster = SV00 mask.
   b. restricted recomputation = the Y10 formulas on the chains that observe the position, in both clusters:
      d_mean_restricted_A = distance between the cl7 and cl4 means of the observing chains;
      d_medoid_restricted_A = distance between the complete-H12 medoids picked (Y10 rule: smallest summed S04 RMSD
      within the set) among the complete-H12 chains of each cluster that observe the position.
      Entry level (rule 1.4): d_mean_entry_A = distance between the cl7 and cl4 means of entry means (unit = entry x
      cluster, mean C-alpha of the unit's observing chains).
   c. entries driving the mean: leave-one-entry-out influence = d_mean(all) - d_mean(without the entry's chains);
      top 5 listed with their numbering offsets. Greedy removal: up to 3 times remove the entry whose removal gives
      the smallest d_mean; d after 1, 2, 3 removals (chain and entry level).
   d. SSE context from SV02 (H fraction per group, table helix, consensus helix per group); B z per cluster.
   e. lattice and partner contacts (SV01b module): share of observing chains with >= 1 heavy atom of the residue
      within 4.0 A of a symmetry copy (lattice) or of another chain's ASU copy (partner), by neighbour kind, and of
      the chain's own bound peptide copy.
3. Numbering-error flag of a chain at a position: the residue observed there differs from UniProt with no deposited
   explanation (residue_status unexplained_difference), or the chain's UniProt-author offset is not one constant over
   the LBD (SV00 offset_constant False) or does not come from SIFTS (SV00 offset_source).
4. [DEFAULT] verdict, first rule that holds (E = measures above 1.0 A in Y10: mean and/or medoid; 'gone' = a finite
   value <= 1.0 A):
   numbering error      >= 1 flagged chain observes the position (or is a medoid) and every measure in E is gone
                         when the flagged chains are removed (mean recomputed, medoids re-picked);
   sparse-sample effect every measure in E is gone in the restricted recomputation (2b);
   driven by <= 3 entries  mean in E: d_mean <= 1.0 A after removing <= 3 entries (greedy, 2c); medoid-only spike
                         (d_mean <= 1.0 A with all chains): the excess rests on the two medoid entries;
   genuine loop difference  otherwise (d_mean > 1.0 A persists after removing the 3 most influential entries).
   Primary verdict = chain level (the Y10 unit); the same rules on entry-level means give the entry-level verdict.
5. Genuine spikes: per chain the signed C-alpha coordinate along the cl7->cl4 mean displacement,
   s = (x - mean_cl7) . u, u = (mean_cl4 - mean_cl7) / |...| (all observing chains). For each covariate (peptide_bound,
   ligand_chemotype, dimer_partner_source; SV00 chain table):
   chain level: OLS s ~ C(cluster) + C(covariate), entry-clustered (sandwich) errors, Wald chi2 of the covariate;
   entry level: unit = entry x cluster, s = unit mean, covariate = the unit's shared value (first matrix chain of the
   unit if its chains differ), OLS s ~ C(cluster) + C(covariate), nested F-test of the covariate.
   The cl4-vs-cl7 coefficient with and without the covariate is reported. One BH family over all covariate tests.
   'associated' = q < 0.05 at chain and entry level.
Outputs: results/sv/SV07/{spikes.csv, verdict.csv, top_entries.csv, restricted_profile.csv, covariate_tests.csv,
numbers.json} (+ SV07_lattice.py outputs). Run after SV07_lattice.py:
  .venv/bin/python -B scripts/sv/SV07_spikes.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))

OUT = ROOT / "results/sv/SV07"
THRESH_A = 1.0          # [DEFAULT] step 1
EXCLUDE = (528, 550)    # step 1
MAX_ENTRIES = 3         # [DEFAULT] verdict 'driven by <= 3 entries'
TOP_N = 5               # entries listed per spike
Q_ALPHA = 0.05
COVARIATES = ["peptide_bound", "ligand_chemotype", "dimer_partner_source"]
FLAG_STATUS = "unexplained_difference"


def spike_positions(y10: pd.DataFrame) -> pd.DataFrame:
    """Step-1 [DEFAULT] rule on the Y10 table."""
    out = y10[(y10.position < EXCLUDE[0]) | (y10.position > EXCLUDE[1])]
    s = out[(out.d_mean_A > THRESH_A) | (out.d_medoid_A > THRESH_A)].copy()
    s["exceeds_mean"] = s.d_mean_A > THRESH_A
    s["exceeds_medoid"] = s.d_medoid_A > THRESH_A
    return s.reset_index(drop=True)


def regions(positions):
    runs = []
    for p in sorted(positions):
        if runs and p == runs[-1][-1] + 1:
            runs[-1].append(p)
        else:
            runs.append([p])
    return {p: (f"{r[0]}-{r[-1]}" if len(r) > 1 else f"{r[0]}") for r in runs for p in r}


# ------------------------------------------------------------------------------------------------ displacement
def d_chain(xs, cs, keep):
    a, b = keep & (cs == 7), keep & (cs == 4)
    if not a.any() or not b.any():
        return np.nan
    return float(np.linalg.norm(xs[a].mean(0) - xs[b].mean(0)))


def d_entry(xs, cs, es, keep):
    """Entry level: unit = entry x cluster (mean of the unit's chains), cluster mean of unit means."""
    mus = {}
    for c in (7, 4):
        sel = keep & (cs == c)
        if not sel.any():
            return np.nan
        u, inv = np.unique(es[sel], return_inverse=True)
        sums = np.zeros((len(u), 3))
        np.add.at(sums, inv, xs[sel])
        cnt = np.bincount(inv, minlength=len(u))[:, None]
        mus[c] = (sums / cnt).mean(0)
    return float(np.linalg.norm(mus[7] - mus[4]))


def medoid(D, members):
    m = np.where(members)[0]
    if len(m) == 0:
        return None
    return int(m[np.argmin(D[np.ix_(m, m)].sum(1))])


def greedy(dfun, entries, k=MAX_ENTRIES):
    """Remove up to k entries, each time the one whose removal gives the smallest d. Returns [(entry, d), ...]."""
    removed, path = [], []
    for _ in range(k):
        best = None
        for e in entries:
            if e in removed:
                continue
            d = dfun(removed + [e])
            if np.isfinite(d) and (best is None or d < best[1]):
                best = (e, d)
        if best is None:
            break
        removed.append(best[0])
        path.append(best)
    return path


def gone(v):
    return bool(np.isfinite(v) and v <= THRESH_A)


# ------------------------------------------------------------------------------------------------ B z (SV01b)
def b_z(keys, pos, mask):
    bb = np.load(ROOT / "data/sv/backbone/backbone.npz")
    assert (bb["chain_keys"].astype(str) == keys).all() and (bb["positions"] == pos).all()
    ca = list(bb["atoms"]).index("CA")
    B = np.where(mask, bb["b_factor"][:, :, ca], np.nan)
    cfg = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    core = np.zeros(len(pos), bool)
    for lo, hi in cfg["core_ranges"]:
        core |= (pos >= lo) & (pos <= hi)
    Z = np.full(B.shape, np.nan)
    for i in range(len(keys)):
        cb = B[i, core]
        cb = cb[np.isfinite(cb) & (cb >= 0)]
        if len(cb) >= 2 and cb.std(ddof=0) > 0:
            z = (B[i] - cb.mean()) / cb.std(ddof=0)
            Z[i] = np.where(np.isfinite(B[i]) & (B[i] >= 0), z, np.nan)
    ref = pd.read_csv(ROOT / "results/sv/SV01b/mobility_per_residue.csv")
    ri = pd.Series(np.arange(len(keys)), index=keys)
    mine = Z[ri[ref.chain_key].values, ref.esr1_position.values - pos[0]]
    check = dict(n_compared=int(len(ref)), max_abs_diff=float(np.nanmax(np.abs(mine - ref.b_z.values))),
                 n_nan_mismatch=int((np.isnan(mine) != np.isnan(ref.b_z.values)).sum()))
    return Z, check


# ------------------------------------------------------------------------------------------------ covariate tests
def covariate_tests(df, pos_label):
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    rows = []
    df = df.copy()
    df["cluster"] = df.cluster.map({7: "cl7", 4: "cl4"})
    df["ent_code"] = pd.factorize(df.entry_id)[0]
    for cov in COVARIATES:
        df["cov"] = df[cov].astype(str)
        desc = "; ".join(f"{c}: " + ", ".join(f"{k} {v}" for k, v in g["cov"].value_counts().sort_index().items())
                         for c, g in df.groupby("cluster", sort=False))
        mean_s = "; ".join(f"{c}/{k}: {v:.2f}" for (c, k), v in df.groupby(["cluster", "cov"]).s.mean().items())
        # entry level units
        units = []
        for (e, c), g in df.groupby(["entry_id", "cluster"], sort=False):
            vals = g["cov"].unique()
            units.append(dict(entry_id=e, cluster=c, s=g.s.mean(), cov=vals[0] if len(vals) == 1 else
                              g.sort_values("row")["cov"].iloc[0], mixed=len(vals) > 1))
        u = pd.DataFrame(units)
        for level, d in (("chain", df), ("entry", u)):
            row = dict(position=pos_label, covariate=cov, level=level, n_units=len(d),
                       n_entries=int(d.entry_id.nunique()), n_levels=int(d["cov"].nunique()),
                       chain_counts_by_cluster=desc, chain_mean_s_by_cluster_level_A=mean_s,
                       n_entry_units_with_mixed_covariate=int(u.mixed.sum()) if level == "entry" else np.nan)
            base = smf.ols("s ~ C(cluster, Treatment('cl7'))", data=d)
            if level == "chain":
                rb = base.fit(cov_type="cluster", cov_kwds={"groups": d.ent_code})
            else:
                rb = base.fit()
            bname = [n for n in rb.params.index if n.startswith("C(cluster")][0]
            row["b_cl4_without_cov_A"] = float(rb.params[bname])
            if d["cov"].nunique() < 2:
                row.update(test="not run: covariate constant", stat=np.nan, df=np.nan, p=np.nan,
                           b_cl4_with_cov_A=np.nan)
                rows.append(row)
                continue
            full = smf.ols("s ~ C(cluster, Treatment('cl7')) + C(cov)", data=d)
            if level == "chain":
                rf = full.fit(cov_type="cluster", cov_kwds={"groups": d.ent_code})
                names = [n for n in rf.params.index if n.startswith("C(cov)")]
                R = np.zeros((len(names), len(rf.params)))
                for i, n in enumerate(names):
                    R[i, list(rf.params.index).index(n)] = 1
                w = rf.wald_test(R, use_f=False, scalar=True)
                row.update(test="Wald chi2, entry-clustered SE", stat=float(w.statistic), df=len(names),
                           p=float(w.pvalue))
            else:
                rf = full.fit()
                a = anova_lm(rb, rf)
                row.update(test="nested F", stat=float(a.F.iloc[1]), df=float(a.df_diff.iloc[1]),
                           p=float(a["Pr(>F)"].iloc[1]))
            row["b_cl4_with_cov_A"] = float(rf.params[bname])
            rows.append(row)
    return rows


def bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        from statsmodels.stats.multitest import multipletests
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


# ------------------------------------------------------------------------------------------------ main
def main():
    from common.core_frame import load_ca_core
    from SV05_core_v4 import y10_profile

    OUT.mkdir(parents=True, exist_ok=True)
    numbers = {}
    y10 = pd.read_csv(ROOT / "results/ext/cluster_displacement_profile.csv")
    keys, pos, X = load_ca_core()
    X = X.astype(float)
    mask = np.load(ROOT / "results/sv/SV00/observed_mask.npy").astype(bool)
    rows = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_rows.csv")
    mpos = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_cols.csv").esr1_position.to_numpy()
    assert list(keys) == rows.chain_key.tolist() and np.array_equal(pos, mpos)
    assert np.array_equal(np.isfinite(X[..., 0]), mask)
    s03 = np.load(ROOT / "results/alignment/aligned_coordinates.npz", allow_pickle=True)
    assert [str(k) for k in s03["chain_keys"]] == list(keys)
    numbers["gate_vs_s03_coords"] = dict(max_abs_diff_A=float(np.nanmax(np.abs(X - s03["coords"]))),
                                         same_gaps=bool(np.array_equal(np.isnan(X), np.isnan(s03["coords"]))))
    t = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    assert t.chain_key.tolist() == list(keys)
    t["row"] = np.arange(len(t))
    cl = t.s08_cluster.to_numpy().astype(int)
    ent = t.entry_id.to_numpy().astype(str)
    h12 = pd.read_csv(ROOT / "results/chain_status.csv").set_index("chain_key").loc[list(keys)].h12_label.to_numpy()
    D = np.load(ROOT / "results/rmsd_entry.npy")
    assert pd.read_csv(ROOT / "results/rmsd_entry_index.csv").chain_key.tolist() == list(keys)
    karr = np.asarray(keys)

    # Y10 reproduction on the core v3 gate coordinates (SV05 re-implementation of the Y10 formulas)
    p3 = y10_profile(X, pos, D, cl, h12, karr)
    numcols = ["d_mean_A", "d_medoid_A", "sd_cl7_A", "sd_cl4_A", "obs_cl7", "obs_cl4", "share_sq"]
    numbers["y10_reproduction_max_abs_diff"] = float(np.nanmax(np.abs(p3[numcols].to_numpy()
                                                                      - y10[numcols].to_numpy())))
    numbers["y10_reproduction_same_nan"] = bool(np.array_equal(p3[numcols].isna().to_numpy(),
                                                               y10[numcols].isna().to_numpy()))
    numbers["y10_medoids"] = [str(y10.medoid_cl7.iloc[0]), str(y10.medoid_cl4.iloc[0])]
    m7 = int(np.where(karr == y10.medoid_cl7.iloc[0])[0][0])
    m4 = int(np.where(karr == y10.medoid_cl4.iloc[0])[0][0])

    # numbering flags
    off_flag = ((t.offset_constant.astype(str) != "True") | (t.offset_source != "sifts_pdbe_api")).to_numpy()
    mp = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "expected_aa", "observed_aa", "residue_status",
                              "difference_annotation", "has_ca"])
    mp = mp[mp.chain_key.isin(set(keys)) & mp.has_ca.astype(bool) & mp.esr1_position.between(pos[0], pos[-1])]
    mp = mp.drop_duplicates(["chain_key", "esr1_position"])
    ri = pd.Series(np.arange(len(keys)), index=list(keys))
    status = np.full(mask.shape, "", dtype=object)
    obs_aa = np.full(mask.shape, "", dtype=object)
    exp_aa = np.full(mask.shape, "", dtype=object)
    ii, jj = ri[mp.chain_key].values, (mp.esr1_position.values - pos[0]).astype(int)
    status[ii, jj] = mp.residue_status.values
    obs_aa[ii, jj] = mp.observed_aa.fillna("").values
    exp_aa[ii, jj] = mp.expected_aa.fillna("").values
    numbers["residue_mapping_rows_vs_mask"] = dict(rows=int(len(mp)), mask_observed=int(mask.sum()),
                                                   mapped_and_masked=int((status != "")[mask].sum()))
    res_flag = status == FLAG_STATUS
    flag = (res_flag | off_flag[:, None]) & mask
    numbers["n_chains_offset_flag"] = int(off_flag.sum())
    numbers["chains_offset_flag"] = karr[off_flag].tolist()
    numbers["n_unexplained_residues_305_550"] = int(res_flag[mask].sum())

    # B z (SV01b definition)
    Z, zcheck = b_z(np.asarray(keys), pos, mask)
    numbers["b_z_check_vs_SV01b"] = zcheck

    # SSE (SV02)
    sse = pd.read_csv(ROOT / "results/sv/SV02/sse_per_position.csv")
    codes = pd.read_csv(ROOT / "results/sv/SV02/sse_codes_per_chain.csv", dtype=str, keep_default_na=False)
    codes = codes.set_index("chain_key").loc[list(keys)]
    hb = pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv")

    def helix_at(p, grp):
        lo, hi = ("table_start", "table_end") if grp == "table" else (f"{grp}_start", f"{grp}_end")
        h = hb[(hb[lo] <= p) & (hb[hi] >= p)].helix_name.tolist()
        return h[0] if h else "loop"

    # lattice (SV07_lattice.py)
    lat = pd.read_csv(OUT / "lattice_contacts_spikes_per_residue.csv")

    # -------------------------------------------------------------- full restricted profile (all 246 positions)
    main = np.isin(cl, (7, 4))
    complete = h12 == "complete"
    prof = []
    for j, p in enumerate(pos):
        obs = mask[:, j] & main
        xs, cs, es = X[obs, j], cl[obs], ent[obs]
        keep = np.ones(len(xs), bool)
        a = medoid(D, complete & (cl == 7) & mask[:, j])
        b = medoid(D, complete & (cl == 4) & mask[:, j])
        prof.append(dict(position=int(p), d_mean_y10_A=y10.d_mean_A.iloc[j], d_medoid_y10_A=y10.d_medoid_A.iloc[j],
                         d_mean_restricted_A=d_chain(xs, cs, keep),
                         d_medoid_restricted_A=float(np.linalg.norm(X[a, j] - X[b, j])) if a is not None and
                         b is not None else np.nan,
                         medoid_cl7_restricted=karr[a] if a is not None else "",
                         medoid_cl4_restricted=karr[b] if b is not None else "",
                         d_mean_entry_A=d_entry(xs, cs, es, keep),
                         n_obs_cl7=int((cs == 7).sum()), n_obs_cl4=int((cs == 4).sum()),
                         obs_cl7=float(mask[cl == 7, j].mean()), obs_cl4=float(mask[cl == 4, j].mean()),
                         n_entries_obs_cl7=int(len(set(es[cs == 7]))), n_entries_obs_cl4=int(len(set(es[cs == 4])))))
    prof = pd.DataFrame(prof)
    prof["excluded_528_550"] = prof.position.between(*EXCLUDE)
    prof.to_csv(OUT / "restricted_profile.csv", index=False)
    diff_mean = (prof.d_mean_restricted_A - prof.d_mean_y10_A).abs()
    diff_med = (prof.d_medoid_restricted_A - prof.d_medoid_y10_A).abs()
    numbers["restricted_vs_y10"] = dict(
        max_abs_diff_d_mean_A=float(diff_mean.max()),
        n_positions_d_medoid_changed=int((diff_med > 1e-9).sum()),
        n_positions_d_medoid_defined_only_restricted=int((prof.d_medoid_y10_A.isna()
                                                          & prof.d_medoid_restricted_A.notna()).sum()),
        positions_restricted_medoid_gt_1_outside_528_550_where_y10_nan=prof.position[
            prof.d_medoid_y10_A.isna() & (prof.d_medoid_restricted_A > THRESH_A) & ~prof.excluded_528_550].tolist(),
        max_d_mean_entry_outside_528_550_A=float(prof.d_mean_entry_A[~prof.excluded_528_550].max()),
        positions_d_mean_entry_gt_1_outside_528_550=prof.position[(prof.d_mean_entry_A > THRESH_A)
                                                                  & ~prof.excluded_528_550].tolist())

    # the other reading of 'in both clusters' (entries holding chains of both clusters): how many such entries
    both = t[np.isin(cl, (7, 4))].groupby("entry_id").s08_cluster.nunique()
    numbers["n_entries_with_cl7_and_cl4_chains"] = int((both == 2).sum())
    numbers["n_entries_cl7_or_cl4"] = int(len(both))

    # -------------------------------------------------------------- spikes
    sp = spike_positions(y10)
    reg = regions(sp.position.astype(int).tolist())
    numbers["n_spikes"] = int(len(sp))
    numbers["spike_positions"] = sp.position.astype(int).tolist()
    spikes, tops, verdicts, cov_rows = [], [], [], []
    for r in sp.itertuples():
        p = int(r.position)
        j = p - pos[0]
        obs = mask[:, j] & main
        idx = np.where(obs)[0]
        xs, cs, es = X[idx, j], cl[idx], ent[idx]
        allk = np.ones(len(idx), bool)
        E = [m for m, ex in (("mean", r.exceeds_mean), ("medoid", r.exceeds_medoid)) if ex]
        d_all, d_ent = d_chain(xs, cs, allk), d_entry(xs, cs, es, allk)
        pr = prof.iloc[j]
        row = dict(position=p, region=reg[p], exceeds=" + ".join(E), d_mean_y10_A=r.d_mean_A,
                   d_medoid_y10_A=r.d_medoid_A, sd_cl7_y10_A=r.sd_cl7_A, sd_cl4_y10_A=r.sd_cl4_A,
                   n_chains_cl7=int((cl == 7).sum()), n_chains_cl4=int((cl == 4).sum()),
                   n_obs_cl7=int(pr.n_obs_cl7), n_obs_cl4=int(pr.n_obs_cl4), obs_cl7=pr.obs_cl7, obs_cl4=pr.obs_cl4,
                   n_entries_obs_cl7=int(pr.n_entries_obs_cl7), n_entries_obs_cl4=int(pr.n_entries_obs_cl4),
                   d_mean_restricted_A=d_all, d_medoid_restricted_A=pr.d_medoid_restricted_A,
                   medoid_cl7_restricted=pr.medoid_cl7_restricted, medoid_cl4_restricted=pr.medoid_cl4_restricted,
                   d_mean_entry_A=d_ent)
        # leave-one-entry-out influence and greedy removal
        ents = sorted(set(es))
        loo = {e: d_chain(xs, cs, es != e) for e in ents}
        loo_e = {e: d_entry(xs, cs, es, es != e) for e in ents}
        path_c = greedy(lambda rem: d_chain(xs, cs, ~np.isin(es, rem)), ents)
        path_e = greedy(lambda rem: d_entry(xs, cs, es, ~np.isin(es, rem)), ents)
        for k in range(MAX_ENTRIES):
            row[f"d_mean_minus_{k + 1}_entries_A"] = path_c[k][1] if k < len(path_c) else np.nan
            row[f"d_mean_entry_minus_{k + 1}_entries_A"] = path_e[k][1] if k < len(path_e) else np.nan
        row["entries_removed_greedy"] = ";".join(e for e, _ in path_c)
        row["entries_removed_greedy_entry_level"] = ";".join(e for e, _ in path_e)
        # numbering flags
        fl = flag[idx, j]
        row["n_flagged_obs"] = int(fl.sum())
        row["flagged_chains"] = ";".join(karr[idx[fl]])
        row["n_unexplained_residue_obs"] = int(res_flag[idx, j].sum())
        row["d_mean_without_flagged_A"] = d_chain(xs, cs, ~fl) if fl.any() else np.nan
        row["d_mean_entry_without_flagged_A"] = d_entry(xs, cs, es, ~fl) if fl.any() else np.nan
        a = medoid(D, complete & (cl == 7) & ~(off_flag | flag[:, j]))
        b = medoid(D, complete & (cl == 4) & ~(off_flag | flag[:, j]))
        row["medoid_flagged"] = bool(off_flag[m7] or off_flag[m4] or flag[m7, j] or flag[m4, j])
        row["d_medoid_without_flagged_A"] = float(np.linalg.norm(X[a, j] - X[b, j]))
        # SSE context
        for g in ("global", "cl7", "cl4"):
            s = sse[(sse.position == p) & (sse.group == g)].iloc[0]
            row[f"frac_H_{g}"] = s.frac_alpha_H
            row[f"consensus_helix_{g}"] = helix_at(p, g)
        row["helix_table"] = helix_at(p, "table")
        cc = codes[str(p)].to_numpy()
        for c in (7, 4):
            sel = obs & (cl == c)
            row[f"frac_H_obs_cl{c}"] = float((cc[sel] == "H").mean()) if sel.any() else np.nan
        # B z
        for c in (7, 4):
            z = Z[obs & (cl == c), j]
            z = z[np.isfinite(z)]
            row[f"b_z_median_cl{c}"] = float(np.median(z)) if len(z) else np.nan
            row[f"b_z_mean_cl{c}"] = float(z.mean()) if len(z) else np.nan
            row[f"n_b_z_cl{c}"] = int(len(z))
        # lattice and partner contacts (observing chains only)
        lp = lat[lat.esr1_position == p].drop_duplicates("chain_key").set_index("chain_key")
        for c in (7, 4):
            ck = karr[obs & (cl == c)]
            q = lp.reindex(ck)
            row[f"n_lattice_rows_cl{c}"] = int(q.n_heavy_atoms.notna().sum())
            for col, name in (("sym_atoms", "lattice_sym"), ("asu_atoms", "partner_asu"),
                              ("sym_receptor_atoms", "sym_receptor"), ("sym_peptide_atoms", "sym_peptide"),
                              ("sym_other_atoms", "sym_other"), ("asu_receptor_atoms", "asu_receptor"),
                              ("asu_peptide_atoms", "asu_peptide"), ("asu_other_atoms", "asu_other"),
                              ("own_peptide_atoms", "own_peptide")):
                row[f"frac_{name}_cl{c}"] = float((q[col] > 0).mean()) if len(q) else np.nan
            row[f"mean_sym_atoms_cl{c}"] = float(q.sym_atoms.mean()) if len(q) else np.nan
            row[f"mean_asu_atoms_cl{c}"] = float(q.asu_atoms.mean()) if len(q) else np.nan

        # ---------------- verdict (chain level primary, entry level with the same rules)
        def verdict(level):
            dm = d_all if level == "chain" else d_ent
            path = path_c if level == "chain" else path_e
            dnf = row["d_mean_without_flagged_A"] if level == "chain" else row["d_mean_entry_without_flagged_A"]
            mean_ex = dm > THRESH_A
            med_ex = "medoid" in E
            ex = [m for m, v in (("mean", mean_ex), ("medoid", med_ex)) if v]
            if not ex:
                return "below threshold at this level", "d_mean <= 1.0 A at this level and no medoid excess"
            n_fl = row["n_flagged_obs"] + int(row["medoid_flagged"])
            if n_fl > 0:
                ok = all(gone(dnf) if m == "mean" else gone(row["d_medoid_without_flagged_A"]) for m in ex)
                if ok:
                    return "numbering error", f"{n_fl} flagged; without them " + ", ".join(
                        f"{m} {dnf if m == 'mean' else row['d_medoid_without_flagged_A']:.2f} A" for m in ex)
            if all(gone(dm if m == "mean" else row["d_medoid_restricted_A"]) for m in ex):
                return "sparse-sample effect", "gone on chains observing the position in both clusters"
            if mean_ex:
                for k, (e, d) in enumerate(path):
                    if gone(d):
                        return "driven by <= 3 entries", (f"d_mean {dm:.2f} -> {d:.2f} A after removing "
                                                          f"{k + 1} entries ({';'.join(x for x, _ in path[:k + 1])})")
                last = path[-1][1] if path else np.nan
                return "genuine loop difference", (f"d_mean {dm:.2f} A, still {last:.2f} A after removing "
                                                   f"{len(path)} entries")
            return "driven by <= 3 entries", (f"d_mean {dm:.2f} A <= 1.0; medoid excess "
                                              f"{r.d_medoid_A:.2f} A rests on the 2 medoid entries")

        vc, why_c = verdict("chain")
        ve, why_e = verdict("entry")
        verdicts.append(dict(position=p, region=reg[p], exceeds=row["exceeds"], verdict=vc, reason=why_c,
                             verdict_entry_level=ve, reason_entry_level=why_e, levels_agree=vc == ve))
        spikes.append(row)

        # top entries
        infl = sorted(((e, d_all - loo[e]) for e in ents), key=lambda kv: -kv[1])[:TOP_N]
        for rank, (e, inf) in enumerate(infl, 1):
            ce = idx[es == e]
            tops.append(dict(position=p, rank=rank, entry_id=e,
                             clusters=";".join(sorted({f"cl{c}" for c in cl[ce]})),
                             n_chains_obs=len(ce), chains=";".join(karr[ce]),
                             d_mean_without_A=loo[e], influence_A=inf,
                             d_mean_entry_without_A=loo_e[e], influence_entry_A=d_ent - loo_e[e],
                             in_greedy_removal=e in [x for x, _ in path_c],
                             dist_to_own_cluster_mean_A=";".join(
                                 f"{np.linalg.norm(X[i, j] - X[obs & (cl == cl[i]), j].mean(0)):.2f}" for i in ce),
                             uniprot_author_offset=";".join(str(t.uniprot_author_offset.iloc[i]) for i in ce),
                             offset_constant=";".join(str(t.offset_constant.iloc[i]) for i in ce),
                             offset_source=";".join(str(t.offset_source.iloc[i]) for i in ce),
                             residue_observed=";".join(obs_aa[i, j] for i in ce),
                             residue_expected=";".join(exp_aa[i, j] for i in ce),
                             residue_status=";".join(status[i, j] for i in ce),
                             numbering_flag=";".join(str(bool(flag[i, j])) for i in ce),
                             sse_code=";".join(cc[i] or "." for i in ce),
                             b_z=";".join(f"{Z[i, j]:.2f}" for i in ce),
                             peptide_bound=";".join(str(t.peptide_bound.iloc[i]) for i in ce),
                             ligand=";".join(str(t.ligand_id_primary.iloc[i]) for i in ce),
                             ligand_chemotype=";".join(str(t.ligand_chemotype.iloc[i]) for i in ce),
                             dimer_partner_source=";".join(str(t.dimer_partner_source.iloc[i]) for i in ce)))

        # covariate tests for genuine spikes (chain-level verdict)
        if vc == "genuine loop difference":
            mu7, mu4 = xs[cs == 7].mean(0), xs[cs == 4].mean(0)
            u = (mu4 - mu7) / np.linalg.norm(mu4 - mu7)
            df = t.iloc[idx][["row", "chain_key", "entry_id"] + COVARIATES].copy()
            df["cluster"] = cs
            df["s"] = (xs - mu7) @ u
            cov_rows += covariate_tests(df, p)

    spikes = pd.DataFrame(spikes)
    verd = pd.DataFrame(verdicts)
    cov = pd.DataFrame(cov_rows)
    if len(cov):
        cov["q_bh"] = bh(cov.p)
        numbers["n_covariate_tests"] = int(cov.p.notna().sum())
        assoc = []
        for (p, c), g in cov.groupby(["position", "covariate"]):
            qs = g.set_index("level").q_bh
            if (qs.reindex(["chain", "entry"]) < Q_ALPHA).all():
                assoc.append((int(p), c, float(qs["chain"]), float(qs["entry"])))
        numbers["covariates_associated"] = assoc
        summ = []
        for p in verd.position:
            g = cov[cov.position == p]
            if not len(g):
                summ.append(("", ""))
                continue
            a = [f"{c} (q chain {qc:.2g}, entry {qe:.2g})" for pp, c, qc, qe in assoc if pp == p]
            ch = g[g.level == "chain"]
            att = "; ".join(f"{x.covariate}: cl4 coef {x.b_cl4_without_cov_A:.2f} -> {x.b_cl4_with_cov_A:.2f} A"
                            for x in ch.itertuples() if np.isfinite(x.b_cl4_with_cov_A))
            summ.append((", ".join(a) if a else "none (q >= 0.05 at chain or entry level)", att))
        verd["covariates_associated"] = [s[0] for s in summ]
        verd["cluster_coef_with_covariate"] = [s[1] for s in summ]
        cov.to_csv(OUT / "covariate_tests.csv", index=False)
    else:
        numbers["n_covariate_tests"] = 0
    spikes.to_csv(OUT / "spikes.csv", index=False)
    verd.to_csv(OUT / "verdict.csv", index=False)
    pd.DataFrame(tops).to_csv(OUT / "top_entries.csv", index=False)
    numbers["verdict_counts_chain"] = verd.verdict.value_counts().to_dict()
    numbers["verdict_counts_entry"] = verd.verdict_entry_level.value_counts().to_dict()
    numbers["n_levels_disagree"] = int((~verd.levels_agree).sum())
    (OUT / "numbers.json").write_text(json.dumps(numbers, indent=1, default=str))
    print(json.dumps(numbers, indent=1, default=str))
    print(verd.to_string())


if __name__ == "__main__":
    main()
