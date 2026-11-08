"""SV11 Part A (RUN_SV1028): which site does each P2Rank pocket occupy, and do cl7 and cl4 differ per site?

Inputs (paths confirmed through the `out` / evidence of the tracker cards, rule 1.10):
- X08 out results/p2rank/<entry>_<label_asym>.pdb_predictions.csv (every pocket, all ranks; residue_ids = author chain
  + author number) and ..._residues.csv; X08 evidence scripts/p2rank_batch.py; X07 out docs/p2rank_check.md.
- S17 out data/clean/*.pdb (scripts/clean_chain_pdbs.py, S17 evidence), the files P2Rank read; data/clean/resname_map.csv.
- Q18 out results/chain_ligand_of_interest.csv, column pocket_residues, written by Q18 out scripts/ligand_of_interest.py
  (CONTACT_CUTOFF_A = 4.5 A, heavy atoms, over its 19 canonical POCKET_POSITIONS only); used as they are.
- SV00 out results/sv/SV00/chain_table.csv: entry_id, auth_asym_id, s08_cluster, uniprot_author_offset
  (= results/sv/SV00/chain_offsets.csv); data/sv/lbd_extract/residue_mapping.csv (SV00) for the wild-type residue letter
  and an offset cross-check.
- results/sv/SV11/chain_site_contacts.csv (scripts/sv/SV11_contacts.py, SV00 structures): groove and dimer contacts.
- SV10 out results/sv/SV10/entry_druglikeness.csv: entry groups (the entry groups below must equal them).

Step 1, P2Rank input check (files): one author chain id per data/clean file; HETATM records after TER are only the
  chain's Q18 ligand(s) (resname_map names); every residue in P2Rank residues.csv and in every pocket's residue_ids is
  a polymer residue (before TER) of that file and carries the file's chain id; sha256 (first 16 hex) of the three S17
  evidence files equals tracker/prov/S17.json. The script reading is reported in the SUMMARY.
Step 2, site residue sets [DEFAULT]: a residue belongs to a site if it is in contact in >= 50 % of the relevant chains:
  ligand pocket  = Q18 pocket_residues (union over a chain's ligand rows), among kept chains with a ligand of interest;
  AF-2 groove    = chain_site_contacts.csv, among kept chains with a bound peptide (<= 30 residues, 4.0 A);
  dimer interface= chain_site_contacts.csv, among kept chains with an SV00 dimer partner (4.0 A).
  site_residues.csv keeps every residue with support > 0 and flags in_site; sets are fixed before any labelling.
Step 3, pocket label (all ranks, UniProt numbering = author + offset): share(site) = pocket residues in the site set /
  pocket residues; label = the site with the largest share if that share >= 0.5 [DEFAULT], else "other surface";
  equal largest shares -> both sites, joined by " + " (the pocket then counts for each of them).
Step 4, per chain and label: any pocket with that label, n such pockets, best (lowest) rank, best probability, SAS
  points of the best-ranked such pocket (P2Rank probability is monotone in score, so it is the same pocket).
Step 5, cl7 vs cl4 per site (ligand pocket, AF-2 groove, dimer interface); declared before running:
  outcomes = share with a pocket at the site; best rank and best probability among the units that have one.
  Chain level (s08_cluster 7 vs 4): GEE, Gaussian identity link, independence working correlation, clusters = entries,
    cluster-robust Wald test of the cl4 term. The binary outcome is a linear-probability model (coefficient = difference
    in shares, defined under separation); rank and probability enter as average ranks within the compared chains.
  Entry level: entry group by the SV10 rule (>= 80 % of the entry's kept chains in cl7 or cl4); the entry has the site
    if >= 50 % of its kept chains have a pocket there (Fisher exact, two-sided); entry best rank / probability = median
    over its chains that have one (Mann-Whitney U, two-sided).
  No variation among the compared units, or < 2 units in a group -> not tested (excluded from BH).
  BH over all tested rows of Part A (both levels). A site "differs" when q < 0.05 at chain AND entry level.
  'other surface' is a label, not a site: described (rows with tested = False), not tested.
Outputs (results/sv/SV11/): p2rank_input_check.json, site_residues.csv, pocket_labels.csv, site_by_chain.csv,
  site_tests.csv, numbers_partA.json. verdict.md and SUMMARY.md: scripts/sv/SV11_summary.py.
Run: .venv/bin/python -B -I scripts/sv/SV11_sites.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV11"
SUPPORT_MIN = 0.5      # [DEFAULT] site residue: contact in >= 50 % of relevant chains
LABEL_MIN = 0.5        # [DEFAULT] pocket label: largest share >= 0.5
ENTRY_SHARE_MIN = 0.8  # SV10 entry-group rule
ENTRY_HAS_MIN = 0.5    # entry has the site if >= 50 % of its chains do (declared above)
LIG, GROOVE, DIMER, OTHER = "ligand pocket", "AF-2 groove", "dimer interface", "other surface"
SITES = (LIG, GROOVE, DIMER)
GROUPS = ("cl7", "cl4")


def read_pdb(path):
    """Polymer residues before the first TER (chain, resnum+icode, resname) and HETATM residue names after TER."""
    chains, poly, lig = set(), set(), set()
    after_ter = False
    for line in path.read_text().splitlines():
        if line.startswith("TER"):
            after_ter = True
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chains.add(line[21])
        if after_ter:
            lig.add(line[17:20].strip())
        else:
            poly.add((line[21], line[22:27].strip()))
    return chains, poly, lig


def read_predictions(path):
    d = pd.read_csv(path, skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    d["name"] = d["name"].str.strip()
    return d


def p2rank_check(ct, q18):
    resmap = pd.read_csv(ROOT / "data/clean/resname_map.csv", dtype=str, keep_default_na=False)
    resmap = dict(zip(resmap.ligand_id, resmap.pdb_resname))
    lig_by_chain = q18[q18.ligand_id != ""].groupby("chain_key").ligand_id.apply(lambda s: {resmap[x] for x in s})
    n = dict(files=0, one_chain_id=0, chain_id_is_author=0, het_after_ter_only_q18=0, apo_files_no_het=0,
             p2rank_residues_in_polymer=0, pocket_residues_in_polymer=0, n_p2rank_residues=0,
             n_p2rank_residues_not_in_polymer=0, n_pocket_tokens=0, n_pocket_tokens_not_in_polymer=0)
    bad = []
    for c in ct.itertuples():
        ent, label = c.chain_key.split("|")[:2]
        pdb = ROOT / "data/clean" / f"{ent}_{label}.pdb"
        chains, poly, lig = read_pdb(pdb)
        n["files"] += 1
        n["one_chain_id"] += len(chains) == 1
        n["chain_id_is_author"] += chains == {c.auth_asym_id}
        expect = lig_by_chain.get(c.chain_key, set())
        n["het_after_ter_only_q18"] += lig == expect
        n["apo_files_no_het"] += (not expect) and (not lig)
        rs = pd.read_csv(ROOT / "results/p2rank" / f"{ent}_{label}.pdb_residues.csv", skipinitialspace=True,
                         dtype=str)
        rs.columns = [x.strip() for x in rs.columns]
        keys = set(zip(rs["chain"].str.strip(), rs["residue_label"].str.strip()))
        miss = keys - poly
        n["n_p2rank_residues"] += len(keys)
        n["n_p2rank_residues_not_in_polymer"] += len(miss)
        n["p2rank_residues_in_polymer"] += not miss
        pr = read_predictions(ROOT / "results/p2rank" / f"{ent}_{label}.pdb_predictions.csv")
        toks = {tuple(t.split("_", 1)) for r in pr.residue_ids for t in r.split()}
        pmiss = toks - poly
        n["n_pocket_tokens"] += len(toks)
        n["n_pocket_tokens_not_in_polymer"] += len(pmiss)
        n["pocket_residues_in_polymer"] += not pmiss
        if len(chains) != 1 or lig != expect or miss or pmiss:
            bad.append(dict(chain_key=c.chain_key, chains=sorted(chains), het=sorted(lig), expected=sorted(expect),
                            p2rank_not_in_polymer=sorted(miss)[:5], pocket_not_in_polymer=sorted(pmiss)[:5]))
    prov = json.loads((ROOT / "tracker/prov/S17.json").read_text())
    sha = []
    for e in prov["ev"]:
        if e["path"].startswith("data/clean/") and e["path"].endswith(".pdb"):
            h = hashlib.sha256((ROOT / e["path"]).read_bytes()).hexdigest()[:16]
            sha.append(dict(path=e["path"], prov_sha=e["sha"], now_sha=h, equal=h == e["sha"]))
    info = json.loads((ROOT / "results/p2rank/run_info.json").read_text())
    return dict(counts=n, problems=bad, s17_sha_check=sha, p2rank_run_info=info,
                script_reading={
                    "scripts/p2rank_batch.py": "P2Rank 'predict' on a dataset listing every data/clean/*.pdb "
                                               "(default model, -visualizations 0); no other input file.",
                    "scripts/clean_chain_pdbs.py": "per chain: polymer atoms of its own label_asym_id (model of the "
                                                   "chain) + its Q18 ligand instance(s) as HETATM; waters, ions, "
                                                   "additives, partner chains and peptides are not written."})


def site_sets(ct, q18, contacts, wt):
    kept = set(ct.chain_key)
    rows = []
    # ligand pocket: Q18 pocket_residues as they are (4.5 A, the 19 canonical pocket positions)
    lq = q18[(q18.ligand_id != "") & q18.chain_key.isin(kept)]
    per_chain = lq.groupby("chain_key").pocket_residues.apply(
        lambda s: {int(t[1:]) for x in s for t in x.split(";") if t})
    sets = {LIG: per_chain}
    for site in (GROOVE, DIMER):
        s = contacts[(contacts.site == site) & contacts.has_partner & contacts.chain_key.isin(kept)]
        sets[site] = s.set_index("chain_key").contact_residues_uniprot.apply(
            lambda x: {int(t) for t in str(x).split(";") if t and t != "nan"})
    members = {}
    for site, pc in sets.items():
        n_rel = len(pc)
        cnt = pd.Series([p for s in pc for p in s]).value_counts()
        for pos, k in cnt.sort_index().items():
            rows.append(dict(site=site, uniprot_position=int(pos), wt_aa=wt.get(int(pos), ""), n_chains_contact=int(k),
                             n_chains_relevant=n_rel, support=k / n_rel, in_site=k / n_rel >= SUPPORT_MIN))
        members[site] = {r["uniprot_position"] for r in rows if r["site"] == site and r["in_site"]}
    sr = pd.DataFrame(rows)
    rel = {site: set(pc.index) for site, pc in sets.items()}
    return sr, members, rel


def label_pockets(ct, members):
    rows, chain_mismatch = [], 0
    for c in ct.itertuples():
        ent, label = c.chain_key.split("|")[:2]
        pr = read_predictions(ROOT / "results/p2rank" / f"{ent}_{label}.pdb_predictions.csv")
        off = float(c.uniprot_author_offset)
        for p in pr.itertuples():
            toks = [t.split("_", 1) for t in p.residue_ids.split()]
            chain_mismatch += sum(ch != c.auth_asym_id for ch, _ in toks)
            res = sorted({int(int(num) + off) for _, num in toks})
            k = {s: len(set(res) & members[s]) for s in SITES}
            best = max(k.values())
            share = best / len(res)
            if share >= LABEL_MIN:
                lab = " + ".join(s for s in SITES if k[s] == best)
            else:
                lab = OTHER
            rows.append(dict(chain_key=c.chain_key, entry_id=c.entry_id, s08_cluster=c.s08_cluster,
                             group=c.group, pocket=p.name, rank=int(p.rank), label=lab, share=share,
                             probability=float(p.probability), score=float(p.score), sas_points=int(p.sas_points),
                             n_residues=len(res), n_ligand_pocket=k[LIG], n_groove=k[GROOVE], n_dimer=k[DIMER],
                             share_ligand_pocket=k[LIG] / len(res), share_groove=k[GROOVE] / len(res),
                             share_dimer=k[DIMER] / len(res), tie=" + " in lab,
                             residues_uniprot=";".join(map(str, res))))
    return pd.DataFrame(rows), chain_mismatch


def by_chain(ct, pl, rel):
    rows = []
    for c in ct.itertuples():
        p = pl[pl.chain_key == c.chain_key]
        for site in SITES + (OTHER,):
            hit = p[p.label.str.split(" + ", regex=False).apply(lambda x: site in x)]
            best = hit.sort_values("rank").iloc[0] if len(hit) else None
            rows.append(dict(chain_key=c.chain_key, entry_id=c.entry_id, s08_cluster=c.s08_cluster, group=c.group,
                             site=site,
                             has_site_partner=(c.chain_key in rel[site]) if site in rel else np.nan,
                             has_pocket=len(hit) > 0, n_pockets=len(hit),
                             best_rank=best["rank"] if best is not None else np.nan,
                             best_probability=hit.probability.max() if len(hit) else np.nan,
                             sas_points=best["sas_points"] if best is not None else np.nan,
                             best_pocket=best["pocket"] if best is not None else "",
                             best_share=best["share"] if best is not None else np.nan))
    return pd.DataFrame(rows)


def group_of(clusters):
    c = np.asarray(clusters)
    if np.mean(c == 7) >= ENTRY_SHARE_MIN:
        return "cl7"
    if np.mean(c == 4) >= ENTRY_SHARE_MIN:
        return "cl4"
    return "mixed cl7/cl4" if np.isin(c, (7, 4)).any() else "minor"


def iqr_text(x):
    x = pd.Series(x).dropna()
    if not len(x):
        return "NA"
    return f"{x.median():.3g} [{x.quantile(0.25):.3g}-{x.quantile(0.75):.3g}]"


def gee(y, cl4, groups):
    X = sm.add_constant(np.asarray(cl4, float))
    res = sm.GEE(np.asarray(y, float), X, groups=np.asarray(groups), family=sm.families.Gaussian(),
                 cov_struct=sm.cov_struct.Independence()).fit()
    lo, hi = res.conf_int()[1]
    return float(res.params[1]), float(lo), float(hi), float(res.params[1] / res.bse[1]), float(res.pvalues[1])


def tests(sb, entry_groups):
    rows = []
    chains = sb[sb.group.isin(GROUPS)]
    for site in SITES + (OTHER,):
        s = chains[chains.site == site]
        tested_site = site in SITES
        # chain level
        a, b = s[s.group == "cl7"], s[s.group == "cl4"]
        base = dict(site=site, level="chain", n_cl7=len(a), n_cl4=len(b))
        r = dict(base, outcome="share with pocket", k_cl7=int(a.has_pocket.sum()), k_cl4=int(b.has_pocket.sum()),
                 value_cl7=a.has_pocket.mean(), value_cl4=b.has_pocket.mean(),
                 test="GEE linear probability, entry clusters (robust Wald)")
        r["effect_cl4_minus_cl7"] = r["value_cl4"] - r["value_cl7"]
        if not tested_site:
            r.update(tested=False, note="label, not a site: descriptive")
        elif s.has_pocket.nunique() < 2:
            r.update(tested=False, note="no variation")
        else:
            coef, lo, hi, z, p = gee(s.has_pocket.astype(float), s.group == "cl4", s.entry_id)
            r.update(tested=True, effect_cl4_minus_cl7=coef, ci_low=lo, ci_high=hi, statistic=z, p=p)
        rows.append(r)
        for col, name in (("best_rank", "best rank"), ("best_probability", "best probability")):
            h = s[s.has_pocket]
            ha, hb = h[h.group == "cl7"], h[h.group == "cl4"]
            r = dict(base, outcome=name, n_cl7=len(ha), n_cl4=len(hb), value_cl7=ha[col].median(),
                     value_cl4=hb[col].median(), iqr_cl7=iqr_text(ha[col]), iqr_cl4=iqr_text(hb[col]),
                     mean_cl7=ha[col].mean(), mean_cl4=hb[col].mean(),
                     test="GEE on average ranks, entry clusters (robust Wald)")
            if not tested_site:
                r.update(tested=False, note="label, not a site: descriptive")
            elif min(len(ha), len(hb)) < 2:
                r.update(tested=False, note="fewer than 2 chains in a group")
            elif h[col].nunique() < 2:
                r.update(tested=False, note="no variation")
            else:
                coef, lo, hi, z, p = gee(rankdata(h[col]), h.group == "cl4", h.entry_id)
                r.update(tested=True, effect_cl4_minus_cl7=coef, ci_low=lo, ci_high=hi, statistic=z, p=p,
                         effect_unit="average rank")
            rows.append(r)
        # entry level
        e = sb[sb.site == site].merge(entry_groups, on="entry_id")
        e = e[e.entry_group.isin(GROUPS)]
        agg = e.groupby(["entry_id", "entry_group"]).agg(
            share=("has_pocket", "mean"), best_rank=("best_rank", "median"),
            best_probability=("best_probability", "median")).reset_index()
        agg["has"] = agg.share >= ENTRY_HAS_MIN
        a, b = agg[agg.entry_group == "cl7"], agg[agg.entry_group == "cl4"]
        base = dict(site=site, level="entry", n_cl7=len(a), n_cl4=len(b))
        r = dict(base, outcome="share with pocket", k_cl7=int(a.has.sum()), k_cl4=int(b.has.sum()),
                 value_cl7=a.has.mean(), value_cl4=b.has.mean(), test="Fisher exact (entry has site: >= 50 % chains)")
        r["effect_cl4_minus_cl7"] = r["value_cl4"] - r["value_cl7"]
        if not tested_site:
            r.update(tested=False, note="label, not a site: descriptive")
        elif agg.has.nunique() < 2:
            r.update(tested=False, note="no variation")
        else:
            odds, p = fisher_exact([[r["k_cl7"], len(a) - r["k_cl7"]], [r["k_cl4"], len(b) - r["k_cl4"]]])
            r.update(tested=True, statistic=odds, p=p, effect_unit="share; statistic = odds ratio cl7/cl4")
        rows.append(r)
        for col, name in (("best_rank", "best rank"), ("best_probability", "best probability")):
            ha, hb = a[col].dropna(), b[col].dropna()
            r = dict(base, outcome=name, n_cl7=len(ha), n_cl4=len(hb), value_cl7=ha.median(), value_cl4=hb.median(),
                     iqr_cl7=iqr_text(ha), iqr_cl4=iqr_text(hb), mean_cl7=ha.mean(), mean_cl4=hb.mean(),
                     test="Mann-Whitney U on entry medians")
            if not tested_site:
                r.update(tested=False, note="label, not a site: descriptive")
            elif min(len(ha), len(hb)) < 2:
                r.update(tested=False, note="fewer than 2 entries in a group")
            elif pd.concat([ha, hb]).nunique() < 2:
                r.update(tested=False, note="no variation")
            else:
                u = mannwhitneyu(ha, hb, alternative="two-sided")
                r.update(tested=True, statistic=float(u.statistic), p=float(u.pvalue),
                         effect_cl4_minus_cl7=hb.median() - ha.median(), effect_unit="median")
            rows.append(r)
    t = pd.DataFrame(rows)
    t["tested"] = t.tested.fillna(False).astype(bool)
    t["q"] = np.nan
    m = t.tested & t.p.notna()
    t.loc[m, "q"] = multipletests(t.loc[m, "p"], method="fdr_bh")[1]
    cols = ["site", "outcome", "level", "n_cl7", "n_cl4", "k_cl7", "k_cl4", "value_cl7", "value_cl4", "iqr_cl7",
            "iqr_cl4", "mean_cl7", "mean_cl4", "effect_cl4_minus_cl7", "effect_unit", "ci_low", "ci_high",
            "statistic", "p", "q", "tested", "test", "note"]
    return t.reindex(columns=cols)


def main():
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    ct["group"] = ct.s08_cluster.map({7: "cl7", 4: "cl4"}).fillna("minor")
    q18 = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv", dtype=str, keep_default_na=False)
    contacts = pd.read_csv(OUT / "chain_site_contacts.csv", keep_default_na=False, na_values=[""],
                           dtype={"contact_residues_uniprot": str})
    rm = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "expected_aa", "has_residue", "auth_seq_id"])
    wt = rm.dropna(subset=["esr1_position"]).groupby("esr1_position").expected_aa.agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) else "")
    wt = {int(k): v for k, v in wt.items()}
    # offset cross-check against the per-residue mapping
    ob = rm[rm.has_residue & rm.chain_key.isin(ct.chain_key)].dropna(subset=["auth_seq_id"])
    ob = ob.merge(ct[["chain_key", "uniprot_author_offset"]], on="chain_key")
    off_ok = int(((ob.auth_seq_id.astype(float) + ob.uniprot_author_offset) == ob.esr1_position).sum())

    check = p2rank_check(ct, q18)
    check["offset_vs_residue_mapping"] = dict(residues=int(len(ob)), equal=off_ok)
    (OUT / "p2rank_input_check.json").write_text(json.dumps(check, indent=1, default=str) + "\n")

    sr, members, rel = site_sets(ct, q18, contacts, wt)
    sr.to_csv(OUT / "site_residues.csv", index=False)

    pl, mismatch = label_pockets(ct, members)
    pl.to_csv(OUT / "pocket_labels.csv", index=False)
    sb = by_chain(ct, pl, rel)
    sb.to_csv(OUT / "site_by_chain.csv", index=False)

    eg = pd.DataFrame([dict(entry_id=e, entry_group=group_of(g.s08_cluster.values))
                       for e, g in ct.groupby("entry_id", sort=False)])
    sv10 = pd.read_csv(ROOT / "results/sv/SV10/entry_druglikeness.csv", keep_default_na=False, na_values=[""])
    same = int((eg.set_index("entry_id").entry_group == sv10.set_index("entry_id").group.reindex(eg.entry_id)).sum())
    t = tests(sb, eg)
    t.to_csv(OUT / "site_tests.csv", index=False)

    lab_counts = pl.label.value_counts().to_dict()
    nums = dict(
        n_chains=int(len(ct)), n_entries=int(ct.entry_id.nunique()), n_pockets=int(len(pl)),
        chain_groups=ct.group.value_counts().to_dict(), entry_groups=eg.entry_group.value_counts().to_dict(),
        entry_groups_equal_sv10=same, residue_chain_mismatch=int(mismatch),
        relevant_chains={s: len(rel[s]) for s in SITES},
        site_size={s: len(members[s]) for s in SITES},
        site_members={s: sorted(members[s]) for s in SITES},
        overlap={f"{a} & {b}": sorted(members[a] & members[b]) for i, a in enumerate(SITES) for b in SITES[i + 1:]},
        label_counts=lab_counts, n_tie_pockets=int(pl.tie.sum()),
        label_counts_rank1=pl[pl["rank"] == 1].label.value_counts().to_dict(),
        n_pockets_per_chain=dict(min=int(pl.groupby("chain_key").size().min()),
                                 median=float(pl.groupby("chain_key").size().median()),
                                 max=int(pl.groupby("chain_key").size().max())),
        q18_cutoff_A=4.5, q18_positions=19, contact_cutoff_A=4.0,
        offset_vs_residue_mapping=check["offset_vs_residue_mapping"])
    (OUT / "numbers_partA.json").write_text(json.dumps(nums, indent=1, default=int) + "\n")
    print(json.dumps({k: nums[k] for k in ("n_pockets", "site_size", "label_counts", "label_counts_rank1",
                                            "entry_groups_equal_sv10", "residue_chain_mismatch")}, default=int))
    print(check["counts"], check["s17_sha_check"])
    print(t[["site", "outcome", "level", "n_cl7", "n_cl4", "value_cl7", "value_cl4", "p", "q", "note"]].to_string())


if __name__ == "__main__":
    main()
